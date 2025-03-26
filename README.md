# How to download and use PalaeoIce in ArcGIS or ArcGIS Pro
The github site includes an old PalaeoIce1.0 folder, an updated PalaeoIce2.0 folder, and a testdata folder, including the test datasets for this toolbox from the Daxi Valley, eastern Tian Shan, China. The PalaeoIce2.0 folder includes an ArcGIS toolbox (tbx) file and a python folder, including all python source codes associated with these tools. The user can click "Code" (green color) on the right side of the github page and choose Download Zip.

![image](https://user-images.githubusercontent.com/24683137/191109319-50965523-61dc-42cb-b977-7d463f9bcc27.png)

A zip file of the whole github folder will be downloaded to the local computer. Unzip this file will create a PalaeoIce-main folder with all folders and files. The reason to include two tbx files is because some tools created in ArcGIS cannot be fully transferred into ArcGIS Pro. The user can use the toolbox corresponding to ArcGIS or ArcGIS Pro, check the source codes, and continue improving this toolbox. Note that the source code file has not been imported to each tool in the current version, so that the toolbox cannot be run only with the tbx file. The test datasets in the file geodatabase include a SRTM DEM (UTM projection), a feature class of 23 extant glaciers, a feature class of 13 glacier outlines during the Little Ice Age, a feature class of a terminal moraine during MIS 2, a feature class of trimlines during MIS 2, and the measured ice thickness from Glacier #1 in this area in 2006.   

The toolboxes and tools have been tested successfully in ArcGIS 10.7, 10.8 and ArcGIS Pro 2.8, 2.9 and 3.0. Errors may occur if using old versions of ArcGIS or ArcGIS Pro. 

# Updated ELA calculation tool (3/26/2026)
Updated the ELA calculation tool: 1) provide the output of the ELA lines for each glacier; 2) Implemented Numba Jit to speed up the ELA calcualtion, achieving about 4 times faster than the old tool.

# ELA calculation tool is added in PalaeoIce 2.0 (1/3/2024)
A new toolbox, ELA.tbx, is added to the PalaeoIce 2.0 folder. This toolbox is revised from the original ArcGIS toolbox of Pellitero et al. (2015) (Pellitero, R., Rea, B.R., Spagnolo, M., Bakke, J., Hughes, P., Ivy-Ochs, S., Lukas, S., Ribolini, A., 2015. A GIS tool for automatic calculation of glacier equilibrium-line altitudes. Comput. Geosci. 82, 55–62). The revised toolbox provides one tool interface to derive the four ELAs (AAR, MGE, AA, and AABR) for each glacier based on the ice surface topography and glacier outline. The user can also specify the contour interval and the AAR and AABR ratios for the calculation. The ELA values are saved in the attrbute table of the glacier outlines. The following is the interface of the ELA Calculation tool. 

![image](https://github.com/yingkui2003/PalaeoIce/assets/24683137/25a6d982-215f-4036-bb7d-c57f30a88ed6)

# PalaeoIce 2.0 is available (6/26/2023)
PalaeoIce 2.0 is the updated version of the PalaeoIce toolbox. The major changes of PalaeoIce 2.0 include:
1.	Added a set of new tools for extant glacier centerlines, extant glacier thickness, flowline quality check, and bare-earth topography preparation.
2.	Reorganized the toolbox and removed the two fully automated tools because most reconstructions need to check separated steps and there are more tools for each step of the reconstruction.
3.	PalaeoIce used one shear stress to derive the ice thickness in all flowlines of a glacier. This may not be suitable for a glacier with many tributaries because different tributaries may have different shear stress values. PalaeoIce 2.0 revised the shear stress optimization method and each flowline will have its corresponding optimized shear stress based on the best fit with the target elevations along the flowline.
4.	Improved the reconstruction of pediment glaciers with much wider termini. The watershed-based approach may not cover the whole terminal area. The revised PalaeoIce toolbox extends the glacier coverage to the user-specified terminal area beyond the DEM-derived watershed. In this way, the whole terminal part can be reconstructed. The only requirement is to provide a detailed terminal part of the glacier outlines as the input.
5.	Revised the two PalaeoIce reconstruction tools to provide six surface interpretation methods. PalaeoIce used the TopoToRaster method to interpolate the ice thickness raster first and then add the ice thickness raster to the bed topography to derive the ice surface topography. This created a smoothed ice thickness raster, but the ice surface elevation raster has high roughness corresponding to the variations in bed topography. PalaeoIce 2.0 interpolates the ice surface elevation raster first with the selection of one of the six methods. This will create a relatively smoother ice surface elevation raster. Then, the ice thickness raster is generated based on the difference between the ice surface elevation raster and bed topography. This revision is more consistent with the field observation that glaciers usually have a smooth ice surface topography.

Because Esri will retire and no longer support ArcGIS 10, the updated PalaeoIce 2.0 toolbox is only for ArcGIS Pro 2.8 and newer.

The PalaeoIce 2.0 toolbox includes 13 tools that are organized in three toolsets: 1) ‘Flowline Creation’; 2) ‘Bare Earth Topography Preparation’; and 3) ‘Palaeoice Reconstruction’. Additional sub-toolsets are created within each toolset for the organizing purpose. The following is the screenshot of this updated toolbox. 

![image](https://github.com/yingkui2003/PalaeoIce/assets/24683137/e48c9eba-fef7-4cc4-a2c9-b7db8bc48754)


In the "Flowline Creation" toolset, two new tools are added in addition to the three tools in PalaeoIce 1.0. 

## Connect OGGM centerlines
This tool connects the modern glacier centerlines derived by the OGGM model.  Maussion et al., 2019. used the OGGM model to derived the centerlines of all RGI 6.0 glaciers (https://docs.oggm.org/en/stable/assets.html). Therefore, these centerlines can be used in PalaeoIce. However, the OGGM centerlines from different tributaries within a glacier are not connected. This tool connect the centerlines within a glacier to generate a connected centerline system for PalaeoIce Reconstruction. 

![image](https://github.com/yingkui2003/PalaeoIce/assets/24683137/e170f8f1-1603-4f39-9cb2-905df913b5b5)

## Remove Overlapping Flowlines
This tool removes the potential overlapping flowlines (one flowline is the identical or a small section of another flowline). The overlapping flowlines may cause the errors in the palaeoice reconstruction.

![image](https://github.com/yingkui2003/PalaeoIce/assets/24683137/d6699596-f96c-4ea5-a88b-d7d13a819276)

PalaeoIce 2.0 added more tools to prepare the bare earth topography, including extracting modern glaicer thickness data from Farinotti et al. (2019), interpreting and removing lake topography, and smoothing the topography of some flowline sections.

## Extract Glacier Thickness Data from Farinotti et al. (2019) Global Dataset
This tool extracts and merges the individual glacier thickness tif files generated by Farinotti et al. (2019): A consensus estimate for the ice thickness distribution of all glaciers on Earth. Nature Geoscience 12, 168–173 and its associated dataset link: https://www.research-collection.ethz.ch/handle/20.500.11850/315707. The individual glacier thickness tif file is named based on the RGI ID of each glacier. This tool creates a merged raster file for all input RGI glaciers, so that the data can be used to remove glacier topography from DEM or compared with the simulated ice thickness data by Volta, PalaeoIce, and other glacier reconstruction models.

![image](https://github.com/yingkui2003/PalaeoIce/assets/24683137/c1f5c124-b358-473c-b067-adb48e9cebae)

## Lake Extraction from DEM
This tool extracts lakes from a DEM based on a minimum area threshold. The lake in a DEM should have zero slopes within the lake area. Therefore, slope analysis can be used to extract lakes. A minimum area is used to remove the noises. 

![image](https://github.com/yingkui2003/PalaeoIce/assets/24683137/112ec9a7-ea81-4270-a8fe-881b659af5a7)

## Remove Lake Topography
This tool is designed to estimate lake bathymetric topography based on lake outlines and specified depth points (required), contours (optional) or simply a maximum depth in the lake polygon attribute table. Then, the estimated lake thickness is removed from the DEM to generate a bare-earth DEM.

![image](https://github.com/yingkui2003/PalaeoIce/assets/24683137/3c1f6bf6-2869-40e2-801e-6048d956cd1c)

## Smooth Elevations Along Certain Flowline Sections
This tool smooths the elevations along the specified flowline section that may be affected by the down-cut of rivers under the ice or after glacier retreat. This tool first uses a buffer of the flowline section to extract the elevations within the buffer boundary lines. Then, interpret the elevation surface within the buffer zone based on the points from the buffer boundary lines. Finally, the DEM is adjusted by replacing the interpreted elevations within the buffer zone.

![image](https://github.com/yingkui2003/PalaeoIce/assets/24683137/9559a928-6ba4-4902-9f38-852bf6b1a963)

PalaeoIce 2.0 also revised the two palaeoice reconstruction tools to provide six surface interpretation methods to interpolate palaeoice surface topography. The original PalaeoIce tools used the TopoToRaster method to interpolate ice thickness raster first and then add the ice thickness raster to the bed topography to derive the ice surface topography. This created a smoothed ice thickness raster, but the ice surface elevation raster has high roughness corresponding to the variations in bed topography. PalaeoIce 2.0 interpolates the ice surface elevation raster first with the selection of one of the six methods. This will create a relatively smoother ice surface elevation raster. Then, ice thickness raster is generated based on the difference between the ice surface elevation raster and bed topography. This revision is more consistent with the field observation that glaciers usually has a smooth ice surface topography.

## Palaeo Ice Reconstruction With Palaeoice Outline(s)
This tool calculates the palaeoglacier thickness (points) along flowlines and produces ice thickness and surface elevation rasters based on the Excel flowline model introduced by Benn and Houlton (2010) and the ArcGIS model, GLaRe, developed by Pellitero et al.(2016). This tool automatically adjusts shear stress and shape factors based on the DEM and the palaeoglacier outline(s). Because a large glacier may include many tributaries, this updated tool derives the shear stress value for each flowline of the reconstructed glacier. The optimization of the shear stress is implemented by an iteration process: 1) a default shear stress value is first used to derive the ice thickness points along the flowline(s); 2) the ice thickness points are then compared with the elevations extracted from the input palaeoglacier outlines; and 3) the shear stress is then optimized based on best-fit between these two sets of elevations. 

This tool automatically adjusts the shape factors based on the ice cross sections along the flowlines.  Two options are available for the shape factor calculation: one is based on the cross-section area, ice-contact perimeter, and ice thickness; the other is based on the fit of the polynomial function introduced by Li et al (2012). The updated tool also provides the six methods to interpret the palaeoice surface elevations: TopoToRaster, Kriging, IDW, Trend, NaturalNeighbor, and Spine. The default method is TopoToRaster, which generates a hydrological corrected ice surface raster. The ice thickness raster is then derived by subtracting the ice surface raster and the bare-ground DEM.

![image](https://github.com/yingkui2003/PalaeoIce/assets/24683137/6e8380b2-2ffe-4423-978c-4ae41eeb1a68)

## Palaeo Ice Reconstruction Without Palaeoice Outline(s)

This tool calculates the palaeoglacier thickness (points) along flowlines and produces the reconstructed palaeoglacier outlines (polygons), and the ice thickness and surface elevation rasters based on the Excel flowline model introduced by Benn and Houlton (2010) and the ArcGIS model, GLaRe, developed by Pellitero et al.(2016). This tool automatically adjusts shear stress and shape factors based on the DEM and the optional target geomorphic features (usually linear features, but can also be point and polygon features), such as trimlines and vegetation breaks, to constrain the palaeoice boundaries. Because a large glacier may include many tributaries, this updated tool derives the shear stress value for each flowline of the reconstructed glacier. The optimization of the shear stress is implemented by an iteration process: 1) a default shear stress value is first used to derive the ice thickness points along the flowline(s) and the ice thickness points is used to interpolate the palaeoice surface elevation; 2) then, the shear stress is updated based on the derived ice elevation distribution and compared to the previous shear stress value. The above iteration is stopped until the shear stress reaches a stable value. If target geomorphic features are provided, the shear stress is adjusted based on the best fit of the calculated ice thickness elevations and target elevations and the minimum distance offset between the reconstructed ice boundary and the target geomorphic features. 

This tool automatically adjusts the shape factors based on the ice cross sections along the flowlines. Two options are available for the shape factor calculation: one is based on the cross-section area, ice-contact perimeter, and ice thickness; the other is based on the fit of the polynomial function introduced by Li et al (2012). The updated tool also provides the six methods to interpret the palaeoice surface elevations: TopoToRaster, Kriging, IDW, Trend, NaturalNeighbor, and Spine. The default method is TopoToRaster, which generates a hydrological corrected ice surface raster. The ice thickness raster is then derived by subtracting the ice surface raster and the bare-ground DEM.

![image](https://github.com/yingkui2003/PalaeoIce/assets/24683137/1b3eb822-2ad5-4742-a352-154bd5fcee58)


# PalaeoIce 1.0
PalaeoIce is an automated method to reconstruct palaeoglaciers based on a digital elevation model (DEM) and geomorphic constraints. Coded in python, PalaeoIce provides a set of tools with user-friendly interfaces to generate glacial flowlines, optimize shear stress, derive shape factors, calculate ice thickness values along flowlines, and interpret palaeo ice surfaces based on geomorphic constraints. The GIS datasets used for these tools, such as the DEM, moraines/cross sections, and extant glacier outline(s), are recommended to use a UTM projection to ensure the correct calculations of the lengths, angles, and areas related to palaeoglacier reconstruction. The toolbox and tools have been tested successfully for ArcGIS 10.7 and 10.8 and ArcGIS Pro 2.8, 2.9 and 3.0. The PalaeoIce toolbox and the related python source codes for each tool are available on https://github.com/yingkui2003/PalaeoIce. 

The PalaeoIce toolbox includes nine GIS tools: Seven tools are organized in three toolsets: 1) ‘Flowline Creation’; 2) ‘Extant Glacier Topography Removal’; and 3) ‘PalaeoIce Reconstruction’. These tools provide the individual steps, allowing for the users to check and adjust the output(s) of each step for palaeoglacier reconstruction. The PalaeoIce toolbox also provides two fully automated palaeoglacier reconstruction tools with and without the constraints of palaeoglacier boundaries to automatically process individual steps without user’s interventions. 

![image](https://user-images.githubusercontent.com/24683137/191107923-6cf0cfe1-fed2-4021-87d4-72bca1d9b511.png)

The ‘Flowline Creation’ toolset includes three GIS tools. ‘Extant Glacier Centerlines (Revised from VOLTA)’ is developed to delineate glacier centerline(s) based on glacier outline(s) and a DEM. This tool is modified from the VOLTA centerline delineation tool (James and Carrivick, 2016). The inputs include extant glacier outline(s), a DEM, the minimum tributary to main valley ratio and the minimum area ratio of a tributary to the whole glacier. The latter two parameters define the two thresholds used for delineating the centerlines from tributaries. The output is the generated extant glacier centerline(s).

![image](https://user-images.githubusercontent.com/24683137/191108164-882a7e80-afe1-4b54-9bdb-01c88a37098f.png)

‘Flowlines from Stream Network’ is to delineate palaeoglacier flowlines based on the GIS stream network analysis. This tool includes four inputs: a DEM, moraines or cross sections (linear features), the minimum source area to start a stream, the minimum total contributing area of a stream before joining another stream, and the minimum tributary to main valley ratio. The outputs include the generated flowlines and an optional upstream watershed boundary for the moraines or cross sections.

![image](https://user-images.githubusercontent.com/24683137/191108235-c26d88da-e3cb-4c5b-9657-6c514b9c4739.png)

‘Combine Flowlines With Centerlines’ is to combine the flowlines derived using the stream network analysis with glacier centerlines derived from extant glacier outlines. This tool requires the inputs of the flowlines, centerlines from extant glaciers, extant glacier outlines, the watershed(s), a DEM, and a search distance (in meters) that is used to determine the linkage between the centerlines and flowlines. The output is the combined flowlines.

![image](https://user-images.githubusercontent.com/24683137/191108310-faf9aa5f-8af6-407e-81d7-efbf2a1e3f0a.png)

Two GIS tools are included in the ‘Extant Glacier Topography Removal’ toolset to determine the bed topography. ‘Extant Glacier Thickness (Revised from VOLTA)’ is a revised tool from VOLTA (James and Carrivick, 2016) to estimate ice thickness based on glacier centerlines, ice surface topography and glacier outlines. This tool incudes a set of inputs: centerlines from extant glaciers, a DEM, extant glacier outlines, ice density, effective slope limit, minimum slope limit, user-specified or default point resolution (cell size of DEM), user-specified or automatically derived shear stress. The outputs include the derived ice thickness points along the centerlines and the ice thickness raster for extant glaciers. Detailed description of this tool can be found in James and Carrivick (2016).

![image](https://user-images.githubusercontent.com/24683137/191108427-522d507b-d41d-4293-b41d-846b6da1be1d.png)

‘Subtract Ice Thickness From DEM’ is to derive bed ground topography by subtracting ice thickness from the DEM. This tool includes two inputs of the original DEM and the derived ice thickness raster. The output is the derived DEM representing the bare ground topography.

![image](https://user-images.githubusercontent.com/24683137/191108571-873a587a-9511-44d2-ad46-33e92312aa0b.png)

Two GIS tools are developed in the ‘PalaeoIce Reconstruction’ toolset for palaeoglacier reconstruction of two scenarios. If the whole palaeoglacier outline(s) can be defined with confidence, only the distributions of palaeo ice thickness and surface elevation are needed to be reconstructed. However, for most cases, the whole outlines of palaeoglaciers are unavailable. In these cases, the outlines, ice thickness, and ice surface of palaeoglaciers can be reconstructed based on terminal moraines or cross sections to define the low limits of the glaciers and optional target geomorphic features (i.e., trimlines) from some sections. 

The ‘Palaeo Ice Reconstruction With PalaeoIce Outline(s)’ includes a set of inputs: a bare ground DEM, flowlines, palaeoice outlines, ice thickness point resolution along flowlines, default shear stress, the minimum and maximum range for the shear stress, the method to derive shape (F) factors. The outputs include palaeoice thickness points along the flowlines and the reconstructed palaeo ice thickness and surface elevation rasters.

![image](https://user-images.githubusercontent.com/24683137/191108709-e31de288-fffb-4c66-b3f8-0c14ba14e646.png)

The ‘Palaeo Ice Reconstruction Without Palaeoice Outline(s)’ tool includes a set of inputs: a bare ground DEM, flowlines, ice thickness point resolution along flowlines, watershed boundaries (optional), target features constraining ice outlines (optional), default shear stress, the minimum and maximum range for the shear stress, the method to derive F factors. The outputs include palaeoice thickness points along the flowlines, the reconstructed palaeoglacier outlines, and the reconstructed palaeo ice thickness and surface elevation rasters. 

![image](https://user-images.githubusercontent.com/24683137/191108796-718fd3b8-966e-49ef-8608-7db665c092b1.png)

In addition to the above tools developed for each major steps of the reconstruction, two GIS tools are also developed to automatically facilitate the whole processes for the two scenarios: with and without palaeoice outline(s). Note that although these two fully automated tools are capable to reconstruct palaeoglaciers without user’s interventions for individual steps, it is strongly recommended that the users run the induvial steps first to obtain the suitable parameters for the reconstruction in a specific region before using the fully automated tools.

The ‘Whole PalaeoIce Model With PalaeoIce Outline(s)’ tool includes a set of inputs: a DEM, moraines or cross sections to define the low limits of palaeoglaciers, palaeoglacier outlines, extant glacier outlines (optional), the minimum source area to start a stream, the minimum total contributing area of a stream before joining another stream, the minimum tributary to main valley ratio, ice thickness point resolution along flowlines, several inputs to delineate extant glacier centerlines and estimate ice thickness, default shear stress, and the method to derive F factors. The outputs include palaeoice thickness points along the flowlines and the reconstructed palaeo ice thickness and surface elevation rasters.

![image](https://user-images.githubusercontent.com/24683137/191108964-596a92cc-888a-4450-8a67-c1d1956dd88a.png)

The ‘Whole PalaeoIce Model Without PalaeoIce Outline(s)’ tool includes a set of inputs of a DEM, moraines or cross sections to define the low limits of palaeoglaciers, target features (optional), extant glacier outlines (optional), the minimum source area to start a stream, the minimum total contributing area of a stream before joining another stream, the minimum tributary to main valley ratio, ice thickness point resolution along flowlines, several inputs to delineate extant glacier centerlines and estimate ice thickness, default shear stress, and the method to derive F factors. The outputs include palaeoice thickness points along the flowlines, the reconstructed palaeoglacier outlines, and the reconstructed palaeo ice thickness and surface elevation rasters. 

![image](https://user-images.githubusercontent.com/24683137/191109055-3ff98418-1c79-4c42-8a43-0595e5cea96f.png)



# How to avoid some potential errors
(1) The code file has not been imported to each tool (for continuous development purpose), so that the toolbox cannot be run only with the tbx file. The code file can be imported to each tool when no further improvement is needed in the future.  


(2) Make sure that the default setting for the scratch workspace in ArcGIS does not include space in the path or folder names. The space in the path or folder names may cause unexpected errors in some raster functions. 


(3) If using this toolbox or other ArcGIS functions (especially raster functions) many times in ArcGIS, ArcGIS may have memory issues, causing unexpected errors. These errors can be solved by restarting the ArcGIS program or the computer. It seems that ArcGIS Pro has a better memory management with few memory-related issues.    

(4) It is recommended to divide a large DEM to small DEMs (less than 1500 X 1500 cells) to avoid the memory issues in the palaeoglacier reconstruction. 

Please report any errors or questions to Yingkui Li (yli32@utk.edu).

# Cite this work
Li Y., 2023. PalaeoIce: an automated method to reconstruct palaeoglacier using geomorphic evidence and digital elevation models. Geomorphology 421, 108523. https://doi.org/10.1016/j.geomorph.2022.108523.

# Most recent updates
01/10/2023: Fixed some bugs in the program.


06/26/2023: PalaeoIce 2.0 is available. The previous toolbox was saved as PalaeoIce 1.0 on Github for the history backup purpose. This updated toolbox is only available for ArcGIS Pro 2.8 or newer.


07/06/2023: Removed the requirements of the UTM projection for some tools, so that the toolbox can be used for any projected coordinate systems;
            Fixed the errors of Volta centerline and Connect OGGM centerlines when using float DEM and when the resolution of DEM is a float number;
            Fixed the errors when the path and filename of the ArcGIS Pro project include empty space(s)
10/23/2023: Fixed an error in volta centerline (the output Features are same as overlay Clip Features).
1/3/2024:   Updated some codes related to OGGM centerline connection and volta ice thickness calculation.

# Contact info
Yingkui Li

Department of Geography & Sustainability

University of Tennessee

Knoxville, TN 37996

Email: yli32@utk.edu

Website: https://geography.utk.edu/about-us/faculty/dr-yingkui-li/

Google Scholar: https://scholar.google.com/citations?user=JoNuyCMAAAAJ&hl=en&oi=ao
