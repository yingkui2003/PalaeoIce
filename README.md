# Introduction
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

# How to download and use this toolbox in ArcGIS or ArcGIS Pro
The github site includes two ArcGIS toolbox (tbx) files (one is for ArcGIS 10.7 and 10.8, the other for ArcGIS Pro 2.8 or newer), a python folder, including all python source codes associated with these tools, and a testdata folder, including the test datasets for this toolbox from the Daxi Valley, eastern Tian Shan, China. The user can click "Code" (green color) on the right side of the github page and choose Download Zip.

![image](https://user-images.githubusercontent.com/24683137/191109319-50965523-61dc-42cb-b977-7d463f9bcc27.png)

A zip file of the whole github folder will be downloaded to the local computer. Unzip this file will create a PalaeoIce-main folder with the two tbx files, the python folder with all python code files, and the testdata folder with a geodatabase of test datasets. The reason to include two tbx files is because some tools created in ArcGIS cannot be fully transferred into ArcGIS Pro. The user can use the toolbox corresponding to ArcGIS or ArcGIS Pro, check the source codes, and continue improving this toolbox. Note that the source code file has not been imported to each tool in the current version, so that the toolbox cannot be run only with the tbx file. The test datasets in the file geodatabase include a SRTM DEM (UTM projection), a feature class of 23 extant glaciers, a feature class of 13 glacier outlines during the Little Ice Age, a feature class of a terminal moraine during MIS 2, a feature class of trimlines during MIS 2, and the measured ice thickness from Glacier #1 in this area in 2006.   

The toolboxes and tools have been tested successfully in ArcGIS 10.7, 10.8 and ArcGIS Pro 2.8, 2.9 and 3.0. Errors may occur if using old versions of ArcGIS or ArcGIS Pro. 

# How to avoid some potential errors
(1) The code file has not been imported to each tool (for continuous development purpose), so that the toolbox cannot be run only with the tbx file. The code file can be imported to each tool when no further improvement is needed in the future.  


(2) Make sure that the default setting for the scratch workspace in ArcGIS does not include space in the path or folder names. The space in the path or folder names may cause unexpected errors in some raster functions. 


(3) If using this toolbox or other ArcGIS functions (especially raster functions) many times in ArcGIS, ArcGIS may have memory issues, causing unexpected errors. These errors can be solved by restarting the ArcGIS program or the computer. It seems that ArcGIS Pro has a better memory management with few memory-related issues.    

Please report any errors or questions to Yingkui Li (yli32@utk.edu).

# Cite this work
Li Y., 2023. PalaeoIce: an automated method to reconstruct palaeoglacier using geomorphic evidence and digital elevation models. Geomorphology 421, 108523. https://doi.org/10.1016/j.geomorph.2022.108523.

# Most recent updates
01/10/2023: Fixed some bugs in the program.

# Contact info
Yingkui Li

Department of Geography & Sustainability

University of Tennessee

Knoxville, TN 37996

Email: yli32@utk.edu

Website: https://geography.utk.edu/about-us/faculty/dr-yingkui-li/

Google Scholar: https://scholar.google.com/citations?user=JoNuyCMAAAAJ&hl=en&oi=ao
