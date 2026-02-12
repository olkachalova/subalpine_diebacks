# subalpine_diebacks
Coding scripts and supplementary materials to the research paper Kachalova et al., 2026 "Recurrent Climate-Driven Dieback of Subalpine Grasslands in Central Europe: Detecting from Multi-Decadal Landsat and Sentinel-2 Time Series"
Coding scripts are free to use without citations.
While using data, please cite this work as: Kachalova, Olha; Řezník, Tomáš; Houška, Jakub; Řehoř, Jan; Trnka, Miroslav; Balek, Jan; and Hédl, Radim. Recurrent Climate-Driven Dieback of Subalpine Grasslands in Central Europe: Detecting from Multi-Decadal Landsat and Sentinel-2 Time Series. Remote Sensing, 2026 (submitted).
Preprint of the article is available from: https://www.preprints.org/manuscript/202602.0992

Scripts:

GEE_harmonized_4band_stack_download.js (Building the harmonized 40-year time series - GEE)
Download harmonized (by spectral reflectance) 4-band (Red, NIR, SWIR1, SWIR2) stacks based on Landsat-5,7,8,9 and Sentinel-2 collections 1984 - 2024 in Google Earth Engine

Field_valid_stats.ipynb (to obtain Figure 3, Table 2)
Calculating ANOVA with Tukey's post-hoc, plotting differences in NPV pixel values between 4 land cover types. Nine NPV indices were tested. At the and, Cohen's d calculated to explain differentiation power between classes

ArcGIS_Pro_thresholding_vectorization.ipynb (to obltain dieback maps: Figure 4, Figure 5; Fisure S1)
Python workflow for ArcGIS Pro (uses arcpy) to NBR index images automatic classification according to threshold values, vectorization and additional vector masking. Output - vector files (*.shp) of diebacks/regeneration extent for each date.

Seasonal_curves_plot.ipynb (to oblain Figure 8)
Python script to visualize seasonal trajectories of NBR in diebacks years (2000, 2003, 2012, 2019) vs baseline (1984-1999). Values vere extracted from pixels where total diebacks had been detected.

SPEI_gpaph_vs_diebacks.ipynb (to obtain Figure 10)
Processing SPEI extracted data, plotting SPEI dataseries together with estimated dieback areas

