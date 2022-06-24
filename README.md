# PalaeoIce
PalaeoIce is an automated method to reconstruct palaeoglaciers based on digitial elevation models (DEM) and geomorphic constraints. Coded in python, PalaeoIce provides a set of tools with user-friendly interfaces to generate glacial flowlines, optimize shear stress, derive shape factors, calculate ice thickness values along flowlines, and interpret palaeo ice surfaces based on geomorphic constraints. The GIS datasets used for these tools, such as the DEM, moraines/cross sections, and extant glacier outline(s), are recommended to use a UTM projection to ensure the correct calculations of the lengths, angles, and areas related to palaeoglacier reconstruction. The toolbox and tools have been tested successfully for ArcGIS 10.7 and 10.8 and ArcGIS Pro 2.8 and 2.9. The PalaeoIce toolbox and the related python source codes for each tool are available on https://github.com/yingkui2003/PalaeoIce. 

The PalaeoIce toolbox includes nine GIS tools: Seven tools are organized in three toolsets: 1) ‘Flowline Creation’; 2) ‘DEM Adjustment’; and 3) ‘PalaeoIce Reconstruction’. These tools provide the individual steps, allowing for the users to check and adjust the output(s) of each step for palaeoglacier reconstruction. The PalaeoIce toolbox also provides two fully automated palaeoglacier reconstruction tools with and without the constraints of palaeoglacier boundaries to automatically process individual steps without user’s interventions. 


The ‘Flowline Creation’ toolset includes three GIS tools. ‘Glacier Centerlines (Revised from VOLTA)’ is developed to delineate glacier centerline(s) based on glacier outline(s) and a DEM (Fig. 5b). This tool is modified from the VOLTA centerline delineation tool (James and Carrivick, 2016, Please check 3.1.1 for detailed modifications). The input parameters include glacier outline(s), a DEM, the minimum tributary to main valley ratio and the minimum area ratio of a tributary to the whole glacier. The latter two parameters define the two thresholds used for delineating the centerlines from tributaries (check 3.1.1). The output is the generated glacier centerline(s).
 
Fig. 5 The PalaeoIce toolbox and three GIS tools in the ‘Flowline Creation’ toolset. (A) The PalaeoIce toolbox in ArcGIS. (B) The interface for the ‘Glacier Centerlines (Revised from VOLTA)’ tool. (C) The interface for the ‘Flowlines from Stream Network’ tool. (D) The interface for the ‘Combine Flowlines With Centerlines’ tool.

‘Flowlines from Stream Network’ is to delineate palaeoglacier flowlines based on the GIS stream network analysis (Section 3.1.2, Fig. 5c). This tool includes four input parameters: a DEM, moraines or cross sections (linear features), the minimum source area to start a stream, the minimum total contributing area of a stream before joining another stream, and the minimum tributary to main valley ratio. The output parameters include the generated flowlines and an optional upstream watershed boundary for the moraines or cross sections.

‘Combine Flowlines With Centerlines’ is to combine the flowlines derived using the stream network analysis with glacier centerlines derived from extant glacier outlines (Section 3.1.3, Fig. 5d). This tool requires the input parameters of the flowlines, centerlines, glacier outlines, the watershed boundary, a DEM, and a search distance that is used to determine the linkage between the centerlines and flowlines. The output parameter is the combined flowlines.

Two GIS tools are included in the ‘DEM Adjustment’ toolset to determine the bed topography. ‘Glacier Thickness (Revised from VOLTA)’ is a revised tool from VOLTA (James and Carrivick, 2016) to estimate ice thickness based on glacier centerlines, ice surface topography and glacier outlines (Section 3.2.1, Fig. 6a). This tool incudes a set of input parameters: glacier centerlines, a DEM, glacier outlines, ice density, effective slope limit, minimum slope limit, user-specified or default point resolution (cell size of DEM), user-specified or automatically derived shear stress. The output parameters include the derived ice thickness points along the centerlines and the ice thickness raster for extant glaciers (Fig. 6a). Detailed description of this tool can be found in James and Carrivick (2016).

‘Adjust DEM by subtracting ice thickness’ is to derive bed topography by subtracting ice thickness from the DEM (Section 3.2.1, Fig. 6b). This tool includes two input parameters of the original DEM and the derived ice thickness raster. The output is the adjusted DEM representing the bare ground topography.
 
Fig. 6 The two GIS tools in the ‘DEM Adjustment’ toolset. (a) The interface for the ‘Glacier Thickness (Revised from VOLTA)’ tool. (b) The interface for the ‘Adjust DEM by subtracting ice thickness’ tool.

Two GIS tools are developed in the ‘PalaeoIce Reconstruction’ toolset for palaeoglacier reconstruction of two scenarios. If the whole palaeoglacier boundaries can be defined with confidence, only the distributions of palaeo ice thickness and surface elevation are needed to be reconstructed. However, for most cases, the whole boundaries of palaeoglaciers are unavailable. In these cases, the outlines, ice thickness, and ice surface of palaeoglaciers can be reconstructed based on terminal moraines or cross sections to define the low limits of the glaciers and optional target geomorphic features (i.e., trimlines) from some sections (Section 3.3.3 and 3.3.4). 

Fig. 7a shows the interface of the ‘Palaeo Ice Reconstruction With PalaeoIce Boundary’ tool. This tool includes a set of input parameters: a DEM, flowlines, palaeoice boundaries, interval distance, default shear stress, the minimum and maximum range for the shear stress, the method to derive f factors. The output parameters include ice thickness points along the flowlines and the reconstructed palaeo ice thickness and surface elevation rasters. 
 
Fig. 7 The two GIS tools in the ‘PalaeoIce Reconstruction’ toolset. (a) The interface for the ‘Palaeo Ice Reconstruction With PalaeoIce Boundary’ tool. (b) The interface for the ‘Palaeo Ice Reconstruction Without PalaeoIce Boundary’ tool.

Fig. 7b illustrates the interface of the ‘Palaeo Ice Reconstruction Without Palaeoice Boundary’. This tool includes a set of input parameters: a DEM, flowlines, interval distance, watershed boundaries (optional), target features constraining ice boundaries (optional), default shear stress, the minimum and maximum range for the shear stress, the method to derive f factors. The output parameters include ice thickness points along the flowlines, the reconstructed palaeoglacier outlines, and the reconstructed palaeo ice thickness and surface elevation rasters. 

In addition to the above tools developed for each major steps of the reconstruction, two GIS tools are also developed to automatically facilitate the whole processes for the two scenarios: with and without palaeoice boundaries (Fig. 8). Note that although these two fully automated tools are capable to reconstruct palaeoglaciers without user’s interventions for individual steps, it is strongly recommended that the users run the induvial steps first to obtain the suitable parameters for the reconstruction in a specific region before using the fully automated tools.

The ‘Whole PalaeoIce Model With PalaeoIce Boundary’ tool (Fig. 8a) includes a set of input parameters: a DEM, moraines or cross sections to define the low limits of palaeoglaciers, palaeoglacier boundaries, extant glacier outlines (optional), the minimum source area to start a stream, the minimum total contributing area of a stream before joining another stream, the minimum tributary to main valley ratio, interval distance, the parameters to delineate extant glacier centerlines and estimate ice thickness, default shear stress, and the method to derive f factors. The output parameters include ice thickness points along the flowlines and the reconstructed palaeo ice thickness and surface elevation rasters.

 
Fig. 8 The two fully automated GIS tools for palaeoglacier reconstruction. (a) The interface for the ‘Whole PalaeoIce Model With PalaeoIce Boundary’ tool. (b) The interface for the ‘Whole PalaeoIce Model Without PalaeoIce Boundary’ tool.

The ‘Whole PalaeoIce Model Without PalaeoIce Boundary’ tool (Fig. 8b) includes a set of input parameters of a DEM, moraines or cross sections to define the low limits of palaeoglaciers, target features (optional), extant glacier outlines (optional), the minimum source area to start a stream, the minimum total contributing area of a stream before joining another stream, the minimum tributary to main valley ratio, interval distance, the parameters to delineate extant glacier centerlines and estimate ice thickness, default shear stress, and the method to derive f factors. The output parameters include ice thickness points along the flowlines, the reconstructed palaeoglacier outlines, and the reconstructed palaeo ice thickness and surface elevation rasters. 

# Contact info
Yingkui Li

Department of Geography

University of Tennessee

Knoxville, TN 37996

Email: yli32@utk.edu

Website: https://geography.utk.edu/about-us/faculty/dr-yingkui-li/

Google Scholar: https://scholar.google.com/citations?user=JoNuyCMAAAAJ&hl=en&oi=ao
