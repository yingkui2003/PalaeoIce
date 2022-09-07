# Introduction
PalaeoIce is an automated method to reconstruct palaeoglaciers based on a digitial elevation model (DEM) and geomorphic constraints. Coded in python, PalaeoIce provides a set of tools with user-friendly interfaces to generate glacial flowlines, optimize shear stress, derive shape factors, calculate ice thickness values along flowlines, and interpret palaeo ice surfaces based on geomorphic constraints. The GIS datasets used for these tools, such as the DEM, moraines/cross sections, and extant glacier outline(s), are recommended to use a UTM projection to ensure the correct calculations of the lengths, angles, and areas related to palaeoglacier reconstruction. The toolbox and tools have been tested successfully for ArcGIS 10.7 and 10.8 and ArcGIS Pro 2.8 and 2.9. The PalaeoIce toolbox and the related python source codes for each tool are available on https://github.com/yingkui2003/PalaeoIce. 

The PalaeoIce toolbox includes nine GIS tools: Seven tools are organized in three toolsets: 1) ‘Flowline Creation’; 2) ‘DEM Adjustment’; and 3) ‘PalaeoIce Reconstruction’. These tools provide the individual steps, allowing for the users to check and adjust the output(s) of each step for palaeoglacier reconstruction. The PalaeoIce toolbox also provides two fully automated palaeoglacier reconstruction tools with and without the constraints of palaeoglacier boundaries to automatically process individual steps without user’s interventions. 

![image](https://user-images.githubusercontent.com/24683137/175660173-3a09a7a6-6e08-4a24-986b-0632af2ce230.png)

The ‘Flowline Creation’ toolset includes three GIS tools. ‘Glacier Centerlines (Revised from VOLTA)’ is developed to delineate glacier centerline(s) based on glacier outline(s) and a DEM. This tool is modified from the VOLTA centerline delineation tool (James and Carrivick, 2016). The inputs include glacier outline(s), a DEM, the minimum tributary to main valley ratio and the minimum area ratio of a tributary to the whole glacier. The latter two parameters define the two thresholds used for delineating the centerlines from tributaries. The output is the generated glacier centerline(s).

![image](https://user-images.githubusercontent.com/24683137/175660330-ed16dd6b-da48-47c9-b681-6fe0274add8c.png)

‘Flowlines from Stream Network’ is to delineate palaeoglacier flowlines based on the GIS stream network analysis. This tool includes four inputs: a DEM, moraines or cross sections (linear features), the minimum source area to start a stream, the minimum total contributing area of a stream before joining another stream, and the minimum tributary to main valley ratio. The outputs include the generated flowlines and an optional upstream watershed boundary for the moraines or cross sections.

![image](https://user-images.githubusercontent.com/24683137/175660463-a99c883b-fec7-40b8-b1e4-745bc141a5f0.png)

‘Combine Flowlines With Centerlines’ is to combine the flowlines derived using the stream network analysis with glacier centerlines derived from extant glacier outlines. This tool requires the inputs of the flowlines, centerlines, glacier outlines, the watershed boundary, a DEM, and a search distance that is used to determine the linkage between the centerlines and flowlines. The output is the combined flowlines.

![image](https://user-images.githubusercontent.com/24683137/175660590-5934ee1f-69ba-484c-81a6-c336ae1dc45f.png)

Two GIS tools are included in the ‘DEM Adjustment’ toolset to determine the bed topography. ‘Glacier Thickness (Revised from VOLTA)’ is a revised tool from VOLTA (James and Carrivick, 2016) to estimate ice thickness based on glacier centerlines, ice surface topography and glacier outlines. This tool incudes a set of inputs: glacier centerlines, a DEM, glacier outlines, ice density, effective slope limit, minimum slope limit, user-specified or default point resolution (cell size of DEM), user-specified or automatically derived shear stress. The outputs include the derived ice thickness points along the centerlines and the ice thickness raster for extant glaciers. Detailed description of this tool can be found in James and Carrivick (2016).

![image](https://user-images.githubusercontent.com/24683137/175660729-02d2a0de-0abc-46c4-8e81-cfb6dfda1059.png)

‘Adjust DEM by subtracting ice thickness’ is to derive bed topography by subtracting ice thickness from the DEM. This tool includes two inputs of the original DEM and the derived ice thickness raster. The output is the adjusted DEM representing the bare ground topography.

![image](https://user-images.githubusercontent.com/24683137/175660815-70aeec3f-7c28-4236-8dc6-4bb3ae4665a7.png)
 
Two GIS tools are developed in the ‘PalaeoIce Reconstruction’ toolset for palaeoglacier reconstruction of two scenarios. If the whole palaeoglacier boundaries can be defined with confidence, only the distributions of palaeo ice thickness and surface elevation are needed to be reconstructed. However, for most cases, the whole boundaries of palaeoglaciers are unavailable. In these cases, the outlines, ice thickness, and ice surface of palaeoglaciers can be reconstructed based on terminal moraines or cross sections to define the low limits of the glaciers and optional target geomorphic features (i.e., trimlines) from some sections. 

The ‘Palaeo Ice Reconstruction With PalaeoIce Boundary’ includes a set of inputs: a DEM, flowlines, palaeoice boundaries, interval distance, default shear stress, the minimum and maximum range for the shear stress, the method to derive f factors. The outputs include ice thickness points along the flowlines and the reconstructed palaeo ice thickness and surface elevation rasters.

![image](https://user-images.githubusercontent.com/24683137/175660988-db73a388-b03f-481c-8052-7bca27dd9881.png)

The ‘Palaeo Ice Reconstruction Without Palaeoice Boundary’ tool includes a set of inputs: a DEM, flowlines, interval distance, watershed boundaries (optional), target features constraining ice boundaries (optional), default shear stress, the minimum and maximum range for the shear stress, the method to derive f factors. The outputs include ice thickness points along the flowlines, the reconstructed palaeoglacier outlines, and the reconstructed palaeo ice thickness and surface elevation rasters. 

![image](https://user-images.githubusercontent.com/24683137/175661073-e69fc83a-6475-48d3-a0f9-804194bffbfc.png)

In addition to the above tools developed for each major steps of the reconstruction, two GIS tools are also developed to automatically facilitate the whole processes for the two scenarios: with and without palaeoice boundaries. Note that although these two fully automated tools are capable to reconstruct palaeoglaciers without user’s interventions for individual steps, it is strongly recommended that the users run the induvial steps first to obtain the suitable parameters for the reconstruction in a specific region before using the fully automated tools.

The ‘Whole PalaeoIce Model With PalaeoIce Boundary’ tool includes a set of inputs: a DEM, moraines or cross sections to define the low limits of palaeoglaciers, palaeoglacier boundaries, extant glacier outlines (optional), the minimum source area to start a stream, the minimum total contributing area of a stream before joining another stream, the minimum tributary to main valley ratio, interval distance, the inputs to delineate extant glacier centerlines and estimate ice thickness, default shear stress, and the method to derive f factors. The outputs include ice thickness points along the flowlines and the reconstructed palaeo ice thickness and surface elevation rasters.

![image](https://user-images.githubusercontent.com/24683137/175661206-7168114c-4e00-49aa-a8e1-6c7666e213ad.png)
 
The ‘Whole PalaeoIce Model Without PalaeoIce Boundary’ tool includes a set of inputs of a DEM, moraines or cross sections to define the low limits of palaeoglaciers, target features (optional), extant glacier outlines (optional), the minimum source area to start a stream, the minimum total contributing area of a stream before joining another stream, the minimum tributary to main valley ratio, interval distance, the inputs to delineate extant glacier centerlines and estimate ice thickness, default shear stress, and the method to derive f factors. The outputs include ice thickness points along the flowlines, the reconstructed palaeoglacier outlines, and the reconstructed palaeo ice thickness and surface elevation rasters. 

![image](https://user-images.githubusercontent.com/24683137/175661270-9b1f1c6b-67a0-4236-9a4c-395a620fe86b.png)

# How to download and use this toolbox in ArcGIS or ArcGIS Pro
The github site includes two ArcGIS toolbox (tbx) file and a python folder, including all python source codes associated with these tools. The user can click "Code" (green color) on the right side of the github page and choose Download Zip.

![image](https://user-images.githubusercontent.com/24683137/175660065-06f763f9-580d-4f6a-8020-d19fd7620ac2.png)

A zip file of the while github folder will be downloaded to the local computer. Unzip this file will create a PalaeoIce-main folder with both the tbx file and the python folder and three code files. The user can use this toolbox, check the codes, and continue improving this toolbox. Note that the codes for each tool are not imported, so that the toolbox cannot be run just with the tbx file.

The toolboxes and tools have been tested in ArcGIS 10.7, 10.8 and ArcGIS Pro 2.8 and 2.9. Errors may occur if using old versions of ArcGIS or ArcGIS Pro. 

Please report any errors or questions to Yingkui Li (yli32@utk.edu).

# Cite this work
Li Y., in review. PalaeoIce: an automated method to reconstruct palaeoglacier using geomorphic evidence and digital elevation models.

# Contact info
Yingkui Li

Department of Geography & Sustainability

University of Tennessee

Knoxville, TN 37996

Email: yli32@utk.edu

Website: https://geography.utk.edu/about-us/faculty/dr-yingkui-li/

Google Scholar: https://scholar.google.com/citations?user=JoNuyCMAAAAJ&hl=en&oi=ao
