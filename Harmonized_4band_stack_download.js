// ======================================================
// FULL SCRIPT (LIGHTWEIGHT / NO FREEZE VERSION)
// UPDATED MASKING:
//   - Landsat: mask fill + snow + (cloud/cirrus/shadow dilated), keep sat==0, SR_B1>0
//   - S2: use SCL+QA60 on 20 m grid, dilate cloud/shadow/cirrus, mask SCL no-data/sat/snow,
//         require B11>0 (valid20)
//
//
// Export: 4-band stacks per scene with band names/order:
//   red, nir, swir1, swir2
// Resolution: Landsat 30 m, S2 20 m
// Cloud metadata: 30%
// Chunking: export only a slice ic.toList(count, startIndex) (no full list / no size getInfo)
// ======================================================

// =======================
// USER SETTINGS
// =======================
var aoi = ee.FeatureCollection('projects/s2cloudless-kach/assets/Jesenik_aoi');
var aoiGeom = aoi.geometry();

var START_YEAR = 1984;
var END_YEAR   = 2024;

var START_MONTH = 5;
var END_MONTH   = 9;

var CLOUD_MAX_LANDSAT = 30;
var CLOUD_MAX_S2      = 30;

// Strict edge removal (DILATE BAD PIXELS)
var DILATE_PIXELS_LANDSAT = 1; // 30 m pixels
var DILATE_PIXELS_S2_20M  = 1; // 20 m pixels

// Export year window
var EXPORT_START_YEAR = 2012;
var EXPORT_END_YEAR   = 2012;

// Chunking
var MAX_IMAGES_PER_GROUP = 300;

// OPTIONAL: set small chunking (e.g. 20) and run multiple times by changing START_INDEX: 0, 20, 40, 60, ...
// var START_INDEX = 0;

var DRIVE_FOLDER = 'GEE_Exports';
var OUT_CRS = 'EPSG:32633';

// Fixed: use B8A (HLS coefficients are defined for B8A)
var nirMode = 'B8A';

// Optional visualization toggle (OFF by default to keep it fast)
var DO_PREVIEW = false;

// =======================
// DATE FILTERS
// =======================
var startDate = ee.Date.fromYMD(START_YEAR, 1, 1);
var endDate   = ee.Date.fromYMD(END_YEAR, 12, 31);
var monthFilter = ee.Filter.calendarRange(START_MONTH, END_MONTH, 'month');

// =======================
// HELPERS
// =======================

// Ensure band order + names exactly: red, nir, swir1, swir2
function forceBandOrder(img) {
  return img.select(['red','nir','swir1','swir2'], ['red','nir','swir1','swir2']);
}

// Landsat C2 L2 scaling (SR scale/offset)
function scaleLandsatC2L2(img) {
  var sr = img.select('SR_B.*').multiply(2.75e-05).add(-0.2);
  return img.addBands(sr, null, true);
}

// -----------------------
// Landsat mask
// -----------------------
function maskLandsatC2L2(img) {
  var qa = img.select('QA_PIXEL');

  var fill         = qa.bitwiseAnd(1 << 0).neq(0);
  var dilatedCloud = qa.bitwiseAnd(1 << 1).neq(0);
  var cirrus       = qa.bitwiseAnd(1 << 2).neq(0);
  var cloud        = qa.bitwiseAnd(1 << 3).neq(0);
  var shadow       = qa.bitwiseAnd(1 << 4).neq(0);
  var snow         = qa.bitwiseAnd(1 << 5).neq(0);

  // Dilate cloud/cirrus/shadow edges (do NOT "erode everything")
  var badEdge = dilatedCloud.or(cirrus).or(cloud).or(shadow).unmask(0);
  var badEdgeDilated = badEdge.focal_max({radius: DILATE_PIXELS_LANDSAT, units: 'pixels'});

  var bad = fill.or(snow).or(badEdgeDilated);

  // Saturation and basic validity
  var sat = img.select('QA_RADSAT').eq(0);
  var valid = img.select('SR_B1').gt(0);

  return img.updateMask(bad.not()).updateMask(sat).updateMask(valid);
}

// Landsat band renames
function selectRenameL57(img) {
  // L5/L7: Red=B3, NIR=B4, SWIR1=B5, SWIR2=B7
  return img.select(['SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'],
                    ['red',   'nir',   'swir1','swir2']);
}
function selectRenameL89(img) {
  // L8/L9: Red=B4, NIR=B5, SWIR1=B6, SWIR2=B7
  return img.select(['SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'],
                    ['red',   'nir',   'swir1','swir2']);
}

// Roy et al. (2016) ETM+ -> OLI adjustment (applied to L5/L7)
function harmonizeL57_to_OLI(img) {
  var red   = img.select('red').multiply(0.8474).add(0.0029);
  var nir   = img.select('nir').multiply(0.8349).add(0.0227);
  var swir1 = img.select('swir1').multiply(0.9372).add(0.0027);
  var swir2 = img.select('swir2').multiply(0.8339).add(0.0009);

  return ee.Image.cat([red, nir, swir1, swir2])
    .rename(['red','nir','swir1','swir2'])
    .copyProperties(img, img.propertyNames())
    .set('HARMONIZED_TO', 'OLI');
}

// -----------------------
// Sentinel-2 mask
// Works on 20 m grid (B11 projection), dilates cloud/shadow/cirrus,
// masks SCL no-data/sat/defective/snow, and requires B11 > 0.
// -----------------------
function s2MasksAt20m(img) {
  var proj20 = img.select('B11').projection(); // 20 m grid
  var scl = img.select('SCL');
  var qa  = img.select('QA60');

  var edgeSCL = scl.eq(3)  // cloud shadow
    .or(scl.eq(8))         // cloud medium
    .or(scl.eq(9))         // cloud high
    .or(scl.eq(10));       // cirrus

  var edgeQA = qa.bitwiseAnd(1 << 10).neq(0)   // opaque clouds
    .or(qa.bitwiseAnd(1 << 11).neq(0));        // cirrus

  var edgeBad = edgeSCL.or(edgeQA).unmask(0).reproject(proj20);
  var edgeBadDilated = edgeBad.focal_max({radius: DILATE_PIXELS_S2_20M, units: 'pixels'});

  var baseBad = scl.eq(0)   // no-data
    .or(scl.eq(1))          // sat/defective
    .or(scl.eq(11))         // snow/ice
    .unmask(0).reproject(proj20);

  var valid20 = img.select('B11').gt(0).unmask(0).reproject(proj20);

  var bad = baseBad.or(edgeBadDilated).or(valid20.not());
  return img.updateMask(bad.not());
}

// OFFICIAL NASA HLS COEFFICIENTS (MSI -> OLI): OLI = slope * MSI + offset
var HLS = {
  'Sentinel-2A': {
    red:   {slope: 0.9765, offset:  0.0009},
    nir:   {slope: 0.9983, offset: -0.0001}, // B8A
    swir1: {slope: 0.9987, offset: -0.0011},
    swir2: {slope: 1.0030, offset: -0.0012}
  },
  'Sentinel-2B': {
    red:   {slope: 0.9761, offset:  0.0010},
    nir:   {slope: 0.9966, offset:  0.0000}, // B8A
    swir1: {slope: 1.0000, offset: -0.0003},
    swir2: {slope: 0.9867, offset:  0.0004}
  }
};

function harmonizeS2_to_OLI(img) {
  var sc = ee.String(img.get('SPACECRAFT_NAME'));
  var key = ee.String(ee.Algorithms.If(sc.compareTo('Sentinel-2A').eq(0), 'Sentinel-2A', 'Sentinel-2B'));

  var coeffs = ee.Dictionary(HLS).get(key);
  coeffs = ee.Dictionary(coeffs);

  var redC = ee.Dictionary(coeffs.get('red'));
  var nirC = ee.Dictionary(coeffs.get('nir'));
  var sw1C = ee.Dictionary(coeffs.get('swir1'));
  var sw2C = ee.Dictionary(coeffs.get('swir2'));

  // IMPORTANT: Use scaled reflectance (divide by 10000).
  // We do it inline here so we don't add extra bands.
  var b4  = img.select('B4').divide(10000);
  var b8a = img.select('B8A').divide(10000);
  var b11 = img.select('B11').divide(10000);
  var b12 = img.select('B12').divide(10000);

  var red   = b4.multiply(redC.getNumber('slope')).add(redC.getNumber('offset'));
  var nir   = b8a.multiply(nirC.getNumber('slope')).add(nirC.getNumber('offset'));
  var swir1 = b11.multiply(sw1C.getNumber('slope')).add(sw1C.getNumber('offset'));
  var swir2 = b12.multiply(sw2C.getNumber('slope')).add(sw2C.getNumber('offset'));

  return ee.Image.cat([
      red.rename('red'),
      nir.rename('nir'),
      swir1.rename('swir1'),
      swir2.rename('swir2')
    ])
    .copyProperties(img, img.propertyNames())
    .set('HARMONIZED_TO', 'OLI')
    .set('NIR_MODE', nirMode);
}

// Sensor tag helper
function addSensorTag(img, tag) {
  return img.set('SENSOR', tag);
}

// =======================
// BUILD COLLECTIONS
// =======================
function prepL57(colId, tag) {
  return ee.ImageCollection(colId)
    .filterBounds(aoiGeom)
    .filterDate(startDate, endDate)
    .filter(monthFilter)
    .filter(ee.Filter.lte('CLOUD_COVER', CLOUD_MAX_LANDSAT))
    .map(maskLandsatC2L2)
    .map(scaleLandsatC2L2)
    .map(selectRenameL57)
    .map(harmonizeL57_to_OLI)
    .map(forceBandOrder)
    .map(function(img){ return addSensorTag(img, tag); });
}

function prepL89(colId, tag) {
  return ee.ImageCollection(colId)
    .filterBounds(aoiGeom)
    .filterDate(startDate, endDate)
    .filter(monthFilter)
    .filter(ee.Filter.lte('CLOUD_COVER', CLOUD_MAX_LANDSAT))
    .map(maskLandsatC2L2)
    .map(scaleLandsatC2L2)
    .map(selectRenameL89)
    .map(forceBandOrder)
    .map(function(img){
      return addSensorTag(img.set('HARMONIZED_TO', 'OLI'), tag);
    });
}

var l5 = prepL57('LANDSAT/LT05/C02/T1_L2', 'L5');
var l7 = prepL57('LANDSAT/LE07/C02/T1_L2', 'L7');
var l8 = prepL89('LANDSAT/LC08/C02/T1_L2', 'L8');
var l9 = prepL89('LANDSAT/LC09/C02/T1_L2', 'L9');

var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoiGeom)
  .filterDate(startDate, endDate)
  .filter(monthFilter)
  .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', CLOUD_MAX_S2))
  .map(s2MasksAt20m)
  .map(harmonizeS2_to_OLI)
  .map(forceBandOrder)
  .map(function(img){ return addSensorTag(img, 'S2'); });

var harmonized = l5.merge(l7).merge(l8).merge(l9).merge(s2)
  .sort('system:time_start');

print('AOI features:', aoi.size());
print('Harmonized collection (filtered):', harmonized);

// =======================
// OPTIONAL PREVIEW (OFF by default)
// =======================
if (DO_PREVIEW) {
  Map.centerObject(aoi, 11);
  var example = ee.Image(harmonized.filterDate('2024-05-01','2024-09-30').first()).clip(aoiGeom);
  Map.addLayer(example, {bands:['swir2','swir1','red'], min:0.0, max:0.4}, 'Harmonized example');
}

// =======================
// EXPORT (CHUNKED) — LIGHTWEIGHT
// - No ic.size().getInfo()
// - No toList(n) of whole year
// - Only toList(count, startIndex) for a small chunk
// =======================
function exportChunk(yearStart, yearEnd, startIndex, count) {
  startIndex = startIndex || 0;
  count = count || 20;

  var ic = harmonized
    .filter(ee.Filter.calendarRange(yearStart, yearEnd, 'year'))
    .sort('system:time_start');

  // Only pull a small slice
  var chunkList = ic.toList(count, startIndex);

  print('Exporting chunk: yearStart=', yearStart,
        'yearEnd=', yearEnd,
        'startIndex=', startIndex,
        'count=', count);

  for (var j = 0; j < count; j++) {
    var img = ee.Image(chunkList.get(j));

    // If we went past the end, system:time_start becomes null; stop.
    var t = img.get('system:time_start');
    if (t === null) break;

    img = forceBandOrder(img);

    // Small getInfo calls per image (OK for 20 images)
    var sensor = ee.String(img.get('SENSOR')).getInfo(); // L5/L7/L8/L9/S2
    var dateStr = ee.Date(t).format('YYYYMMdd').getInfo();

    var scale = (sensor === 'S2') ? 20 : 30;

    Export.image.toDrive({
      image: img.clip(aoiGeom),
      description: 'HARM_SR_' + sensor + '_' + dateStr + '_' + (startIndex + j),
      folder: DRIVE_FOLDER,
      fileNamePrefix: 'HARM_SR_' + sensor + '_' + dateStr,
      region: aoiGeom,
      crs: OUT_CRS,
      scale: scale,
      maxPixels: 1e13
    });
  }
}

// =======================
// RUN EXPORT (one chunk)
// - Run again with START_INDEX = 20, 40, 60, ...
// =======================
exportChunk(EXPORT_START_YEAR, EXPORT_END_YEAR, START_INDEX, MAX_IMAGES_PER_GROUP);