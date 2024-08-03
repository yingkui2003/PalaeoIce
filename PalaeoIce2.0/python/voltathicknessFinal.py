#-------------------------------------------------------------------------------
# Name: VOLTA Ice Thickness Reconstruction
# Purpose: This tool derives ice thickness for modern glaciers based on the revised ice thickness function proposed by Li et al (2012). This method
# was implemented in ArcGIS by James and Carrivick (2016) and called VOLTA. 
# 
# The codes here are modified based on the original codes by James and Carrivick(2016). The imporvements include:
# 1)refining the logic;
# 2) fixed the errors bugs;
# 3) optimize the code to speed up the process
# 4) revised the shear stress code based on the elevation distrbution of ice surface to represent closely
#    the original 1986 paper and include the last part of the highest contourline and the highest points in the shear stress calculation. 
# 
# Author:      Yingkui Li
# Created:     08/30/2020 to 12/16/2020
# Copyright:   (c) Yingkui Li 2020
# Department of Geography, University of Tennessee
# Knoxville, TN 37996, USA
#-------------------------------------------------------------------------------

#import SharedFunctions  
from SharedFunctions import *
import time
temp_workspace = "in_memory" 
#if ArcGISPro:
#    temp_workspace = "memory" 
#------------------------------------------------------------------------------------------------------------
# This function derives the shear stress of a glacier based on the primary flowline and ice surface elevations.
# To make sure to extract the valid elevations for the start and end points of the flowline, this tool applied a
# small buffer around these two points (2 cellsize) and extract the minimum and maximum elevations of the buffer.
#
# This function is revised based on the codes by James and Carrivick (2016) in VOLTA.
# Specifically, the contour step is changed to 300 m (1000ft) as proposed by the original paper (1986)
# and counted the missing part from the highest contour to the highest elevation of the flowline.
#------------------------------------------------------------------------------------------------------------
def shear_stress_calculation_bak(mainflowline, outline, icedem, min_ss, max_ss): 
    #arcpy.AddMessage("Start shear stress calculation...")
    cellsize = arcpy.GetRasterProperties_management(icedem,"CELLSIZEX")
    cellsize_value = float(cellsize.getOutput(0)) ### should change to float ?????

    primary_flowline = arcpy.CopyFeatures_management(mainflowline, "in_memory\\primary_flowline")
    ##Need to just keep the longest flowline
    lengthArr = arcpy.da.FeatureClassToNumPyArray(primary_flowline, ('SHAPE@LENGTH'))
    lengths = np.array([item[0] for item in lengthArr])
    max_length = max(lengths)
    ##Get the flowline length
    lengths = []
    with arcpy.da.UpdateCursor(primary_flowline,["SHAPE@LENGTH"]) as cursor:
        for row in cursor:
            if row[0] < max_length:
                #arcpy.AddMessage("Delete one line")
                cursor.deleteRow()
    del row, cursor
    
    flowlinelength = max_length

    ##4/21/2023: Can also used the extract value to  with the interpretation option
    arcpy.FeatureVerticesToPoints_management(primary_flowline, "in_memory\\startendpnt", "BOTH_ENDS")
    ExtractValuesToPoints("in_memory\\startendpnt", icedem, "in_memory\\startendpnt_ele", "INTERPOLATE", "VALUE_ONLY")
    elev = []
    with arcpy.da.SearchCursor("in_memory\\startendpnt_ele",["RASTERVALU"]) as cursor:
        for row in cursor:
            elev.append(row[0])
    del row, cursor

    #arcpy.AddMessage(elev)
    startz = min(elev)
    endz = max(elev)
    Zdiff = float(endz - startz)

    extractDEM = ExtractByMask(icedem, outline)    

    ##Revise the codes using raster calculations and clip analysis
    elev_bin_polygon = "in_memory\\elev_bin_polygon"
    subset_flowline = "in_memory\\subset_flowlines"
    #if Zdiff >= 300.0: ### the orginal paper was 300; so change to 300 m
    contour_interval = 300.0
    zone_num = int(Zdiff / contour_interval)
    sum_a_cos = 0
    z_start = startz
    sum_a_cos_sin = 0
    while (z_start < endz):
        #arcpy.AddMessage("Z_start is: " + str(z_start))
        z_end = min(z_start + contour_interval, endz)
        #arcpy.AddMessage("z_end is: " + str(z_end))
        ##get flowline length within the bin
        outCon = Con(((extractDEM > z_start) & (extractDEM < z_end)), 1)
        arcpy.RasterToPolygon_conversion(outCon, elev_bin_polygon, "SIMPLIFY", "VALUE")
        polyArr = arcpy.da.FeatureClassToNumPyArray(elev_bin_polygon, ('SHAPE@AREA'))
        areas = np.array([item[0] for item in polyArr])
        #arcpy.AddMessage(areas)
        if len(areas) > 0:
            area = max(areas)
            #arcpy.AddMessage("Area is:" + str(area))
            arcpy.Clip_analysis(primary_flowline, elev_bin_polygon, subset_flowline) ### Use clip analysis to replace the select by location
            arcpy.MultipartToSinglepart_management(subset_flowline, "in_memory\\singlepart_subset_flowline")
            lengthArr = arcpy.da.FeatureClassToNumPyArray("in_memory\\singlepart_subset_flowline", ('SHAPE@LENGTH'))
            lengths = np.array([item[0] for item in lengthArr])
            #arcpy.AddMessage("flowline Lengths is:")
            #arcpy.AddMessage(lengths)
            if len(lengths) > 0:
                max_length = max(lengths)
                if max_length > 0:
                    angle_rad = math.atan((z_end - z_start)/max_length)
                    area_cos = area/(math.cos(angle_rad))
                    sum_a_cos += area_cos
                    sum_a_cos_sin += area_cos * (math.sin(angle_rad))

        z_start = z_end

    shear_stress = 27000* sum_a_cos**0.106  ### this equation was from Glacier volume estimation on Cascade volcanoes 1986 paper in Annual of Glaciology
    #shear_stress = 21605* sum_a_cos**0.106  ### this equation was from Glacier volume estimation on Cascade volcanoes 1986 paper in Annual of Glaciology
    if shear_stress < min_ss:
        shear_stress = min_ss
        arcpy.AddMessage("WARNING: The calculated shear stress is smaller than the specified minimum value, using the specified minimum value instead")
    if shear_stress > max_ss: 
        shear_stress = max_ss
        arcpy.AddMessage("WARNING: The calculated shear stress is larger than the specified maximum value, using the specified maximum value instead")
            
    #else:
    #    arcpy.AddMessage("WARNING: glacier altitudinal extent < 300 m. using the default value instead")
    #    shear_stress = min_ss
    
    return shear_stress

def shear_stress_calculation(mainflowline, outline, icedem, min_ss, max_ss): 
    ##New codes
    points = temp_workspace + "\\points"
    primary_flowline = temp_workspace +"\\primary_flowline"
    flowline3D = temp_workspace +"\\flowline3D"
    elev_bin_polygon = temp_workspace + "\\elev_bin_polygon"
    subset_flowline = temp_workspace + "\\subset_flowlines"
    singlepart_subset_flowline = temp_workspace + "\\singlepart_subset_flowline"

    primary_flowline = arcpy.CopyFeatures_management(mainflowline, primary_flowline)
    ##Need to just keep the longest flowline
    lengthArr = arcpy.da.FeatureClassToNumPyArray(primary_flowline, ('SHAPE@LENGTH'))
    lengths = np.array([item[0] for item in lengthArr])
    max_length = max(lengths)
    ##Get the flowline length
    lengths = []
    with arcpy.da.UpdateCursor(primary_flowline,["SHAPE@LENGTH"]) as cursor:
        for row in cursor:
            if row[0] < max_length:
                cursor.deleteRow()
    del row, cursor
    
    flowlinelength = max_length
    
    area = 0
    with arcpy.da.SearchCursor(outline,["SHAPE@AREA"]) as cursor:
        for row in cursor:
            area = row[0]
    del row, cursor
    cellsizeX = icedem.meanCellWidth
    
    arcpy.InterpolateShape_3d(icedem, primary_flowline, flowline3D, cellsizeX)
    arcpy.FeatureVerticesToPoints_management(flowline3D, points, "ALL")
    pnt3DArray = arcpy.da.FeatureClassToNumPyArray(points,["SHAPE@Z"])
    Zs = np.array([item[0] for item in pnt3DArray])

    extractDEM = ExtractByMask(icedem, outline)

    max_elev = np.max(Zs)
    min_elev = np.min(Zs)
    contour_arr = np.arange(min_elev, max_elev, 300)
    sum_a_cos = 0
    if len(contour_arr) == 1:
        elev_diff = max_elev - min_elev
        #segment_length = dis.sum()
        segment_length = flowlinelength
        angle_rad = math.atan(elev_diff/segment_length)
        sum_a_cos = area/(math.cos(angle_rad))
    else: ## if the elevation difference is > 300 m
        if (max_elev - np.max(contour_arr)) > 1: ## if the max elevation is 10 m higher than the contour arr max
            np.append(contour_arr, max_elev)

        for i in range(1, len(contour_arr)):
            z_start = contour_arr[i-1]
            z_end = contour_arr[i]

            outCon = Con(((extractDEM > z_start) & (extractDEM < z_end)), 1)
            arcpy.RasterToPolygon_conversion(outCon, elev_bin_polygon, "SIMPLIFY", "VALUE")
            polyArr = arcpy.da.FeatureClassToNumPyArray(elev_bin_polygon, ('SHAPE@AREA'))
            areas = np.array([item[0] for item in polyArr])
            if len(areas) > 0:
                area = max(areas)
                arcpy.Clip_analysis(primary_flowline, elev_bin_polygon, subset_flowline) ### Use clip analysis to replace the select by location
                arcpy.MultipartToSinglepart_management(subset_flowline, singlepart_subset_flowline)
                lengthArr = arcpy.da.FeatureClassToNumPyArray(singlepart_subset_flowline, ('SHAPE@LENGTH'))
                lengths = np.array([item[0] for item in lengthArr])
                if len(lengths) > 0:
                    max_length = max(lengths)
                    if max_length > 0:
                        angle_rad = math.atan((z_end - z_start)/max_length)
                        area_cos = area/(math.cos(angle_rad))
                        sum_a_cos += area_cos
            z_start = z_end

    shear_stress = 27000* sum_a_cos**0.106  ### this equation was from Glacier volume estimation on Cascade volcanoes 1986 paper in Annual of Glaciology
    if shear_stress < min_ss:
        shear_stress = min_ss
    if shear_stress > max_ss: 
        shear_stress = max_ss
    return shear_stress

'''
#------------------------------------------------------------------------------------------------------------
# This function check each line in the line feature and make sure the line is from low elevation to high
# elevation (Glacier flowline needs from low to high elevation in order to reconstruct the paleo ice thickness).
# It is revised from the codes by Pellitero et al.(2016) in GlaRe.
#------------------------------------------------------------------------------------------------------------
def Check_If_Flip_Line_Direction(line, dem):
    
    #cellsize = dem.meanCellWidth
    #cellsize_int = int(cellsize)
    cellsize = arcpy.GetRasterProperties_management(dem,"CELLSIZEX")
    cellsize_int = int(float(cellsize.getOutput(0)))
    #arcpy.AddMessage("cellsize_int: " + str(cellsize_int))

    line3d = arcpy.env.scratchGDB + "\\line3d"
    arcpy.AddField_management(line, "Flip", "Long", "", "", "", "", "", "", "")

    arcpy.InterpolateShape_3d(dem, line, line3d, cellsize_int*3) 

    flip_list = []
    i = 0
    with arcpy.da.SearchCursor(line3d,["Shape@"]) as cursor:
        for row in cursor:
            startZ = row[0].firstPoint.Z
            #arcpy.AddMessage("startZ: " + str(startZ))
            endZ = row[0].lastPoint.Z
            #arcpy.AddMessage("endZ: " + str(endZ))

            if startZ >= endZ:  ##Flip = True use equal in case the start and end point are the same
                flip_list.append(1)
            else:  ##Flip = False
                flip_list.append(0)
            i += 1 


    del cursor
    if i>0:
        del row

    #arcpy.AddMessage(flip_list)
    #arcpy.AddMessage(str(sum(flip_list)))

    if sum(flip_list) > 0:
        with arcpy.da.UpdateCursor(line,["Flip"]) as cursor:
            i = 0
            for row in cursor:
                row[0] = flip_list[i]
                cursor.updateRow(row)
                i += 1 
        del row, cursor

        arcpy.MakeFeatureLayer_management(line, "lyrLines")
        arcpy.SelectLayerByAttribute_management("lyrLines", "NEW_SELECTION", '"Flip" > 0')
        arcpy.AddMessage("The number of fliped lines is: " + str(sum(flip_list)))
        arcpy.FlipLine_edit("lyrLines")  ##Need to change to lyrLines
        arcpy.SelectLayerByAttribute_management("lyrLines", "CLEAR_SELECTION")

    arcpy.DeleteField_management (line, "Flip")
    arcpy.Delete_management(line3d)
'''

#------------------------------------------------------------------------------------------------------------
# This is the main function for the ice thickness calculation based on Volta program. This function can be call by the overall model
#------------------------------------------------------------------------------------------------------------
def Ice_Thickness_Volta (flowline, dem, in_outline, ice_density, slope_limit, min_slope, point_res_check, point_res, shear_stress_test, shear_stress_value, final_points,
                     interpolate_check, cellsize_interpolate_check, cellsize_interpolate_user_spec, raster_out):

    cellsize = arcpy.GetRasterProperties_management(dem,"CELLSIZEX")
    cellsize_float = float(cellsize.getOutput(0)) # use float cell size

    arcpy.env.outputCoordinateSystem = dem

    #output_list = []
    gravity = 9.81
    effective_slope_limit = slope_limit

    ##Check and Flip flowline direction
    arcpy.AddMessage("Checking flowline direction....")
    print("Checking flowline direction....")
    Check_If_Flip_Line_Direction (flowline, dem)

    exist_fields_list_centrelines = [f.name for f in arcpy.ListFields(flowline)] #List of current field names in outline layer
    line_fields = ["line_id"]
    for field in line_fields:
        if field not in exist_fields_list_centrelines:
            arcpy.AddField_management(flowline, field, "LONG")
    arcpy.CalculateField_management(flowline,"line_id",str("!"+str(arcpy.Describe(flowline).OIDFieldName)+"!"),"PYTHON_9.3")

    #Make sure that the outline is selected with the overlap of flowlines; added on 2/27/2023
    outline = "in_memory\\outline"
    arcpy.SpatialJoin_analysis(in_outline, flowline, outline,"JOIN_ONE_TO_ONE", "KEEP_COMMON","", "INTERSECT")
  
    exist_fields_list_outlines = [f.name for f in arcpy.ListFields(outline)] #List of current field names in outline layer
    outline_fields = ["stress_pa", "outline_id"]
    for field in outline_fields:
        if field not in exist_fields_list_outlines:
            arcpy.AddField_management(outline, field, "LONG")
    arcpy.CalculateField_management(outline,"outline_id",str("!"+str(arcpy.Describe(outline).OIDFieldName)+"!"),"PYTHON_9.3")

    if "va_thick" not in exist_fields_list_outlines:
         arcpy.AddField_management(outline, "va_thick", "DOUBLE")

    flow_id_list = []
    with arcpy.da.SearchCursor(flowline, ["line_id"]) as cursor:
        for row in cursor:
            flow_id_list.append(row[0])
    del row, cursor

    with arcpy.da.UpdateCursor(outline, ["Shape@Area","va_thick"]) as cursor:
        for row in cursor:
            area_km2 = row[0]/1000.0/1000.0
            ####To calculate the average thickness based on the glacier area; Should be able to find the original equation
            if area_km2 <= 25.0:
                row[1]= (0.0435*area_km2**1.23)/area_km2*1000 ###Where this equation come from????
            else:
                row[1]= (0.0540*area_km2**1.20)/area_km2*1000  ###Where this equation come from????
            cursor.updateRow(row)
    del row, cursor

    order_dict = {}   ##this is the global variable
    max_length_dict = {}

    start_time = time.time()    
    polyArray = arcpy.da.FeatureClassToNumPyArray(outline,('outline_id'))##,'stress_pa'))
    outlineID = np.array([item[0] for item in polyArray])
    #uniqueiceID = np.unique(iceID)
    ss_list = []
    for gid in outlineID:
        query = "outline_id = " + str(gid)
        arcpy.AddMessage("Processing #" + str(gid) +"/" + str(len(outlineID)) + " of reconstructed glaciers...")                                                                                       
        arcpy.Select_analysis (outline, "in_memory\\single_outline", query)
        arcpy.SpatialJoin_analysis(flowline, "in_memory\\single_outline", "in_memory\\subset_flowline", "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "HAVE_THEIR_CENTER_IN")

        flowlineArray = arcpy.da.FeatureClassToNumPyArray("in_memory\\subset_flowline",("SHAPE@LENGTH"))
        lengths = np.array([item[0] for item in flowlineArray])
        max_length = max(lengths)

        max_length_dict[gid] = max_length
        if shear_stress_test == "true":
            shear_stress_value = shear_stress_calculation("in_memory\\subset_flowline", "in_memory\\single_outline", dem, 50000, 200000)
            arcpy.AddMessage("derived shear stress: " + str(shear_stress_value))
        else:
            arcpy.AddMessage("Use the user-specified shear stress for glacier outline "+str(gid)+": "+str(int(shear_stress_value))+" Pa")

        ss_list.append(int(shear_stress_value))        

     ##Add shearstress value to the outline  
    with arcpy.da.UpdateCursor(outline, ["stress_pa"]) as cursor:
        i = 0
        for row in cursor:
            row[0] = ss_list[i]
            cursor.updateRow(row)
            i += 1
    process_time = time.time() - start_time
    arcpy.AddMessage("the processing time is: " + str(process_time))
    '''
    start_time = time.time()    
    with arcpy.da.UpdateCursor(outline, ["outline_id","stress_pa"]) as cursor:
        for row in cursor:
            outline_query = "outline_id = "+str(row[0])
            single_outline = arcpy.Select_analysis(outline, "in_memory\\single_outline", outline_query)  ###use the select analysis tool to replace all select by attribute and copy feature
            #subset_flowline = arcpy.Clip_analysis(flowline, single_outline, "in_memory\\subset_flowlines") ### Use clip analysis to replace the select by location
            subset_flowline = arcpy.SpatialJoin_analysis(flowline, single_outline,"in_memory\\subset_flowline","JOIN_ONE_TO_ONE", "KEEP_COMMON","", "HAVE_THEIR_CENTER_IN")

            lengthArr = arcpy.da.FeatureClassToNumPyArray(subset_flowline, ('SHAPE@LENGTH'))
            lengths = np.array([item[0] for item in lengthArr])
            max_length = max(lengths)
            ##Get the flowline length
            #max_length_dict[row[0]] = max_length_outline
            max_length_dict[row[0]] = max_length
            if shear_stress_test == "true":
                shear_stress_value = shear_stress_calculation(subset_flowline, single_outline, dem, 50000, 200000) ##revised based onupdated function
                arcpy.AddMessage("The shear stress derived for glacier outline "+str(row[0])+": "+str(int(shear_stress_value))+" Pa")
            else:
                arcpy.AddMessage("Use the user-specified shear stress for glacier outline "+str(row[0])+": "+str(int(shear_stress_value))+" Pa")
            row[1] = int(shear_stress_value)
            cursor.updateRow(row)
    del row, cursor
    process_time = time.time() - start_time
    arcpy.AddMessage("the processing time is: " + str(process_time))
    '''
    flowline_joined = arcpy.SpatialJoin_analysis(flowline, outline,"in_memory\\flowline_joined","JOIN_ONE_TO_ONE", "KEEP_COMMON","", "HAVE_THEIR_CENTER_IN")

    flowline_layer_split_vertex = arcpy.SplitLine_management(flowline, "in_memory\\flowline_layer_split_vertex")

    #the loop for each flowline
    with arcpy.da.SearchCursor(flowline_joined, ["line_id","outline_id", "stress_pa","SHAPE@LENGTH"]) as cursor:
        counter = 0
        for row in cursor:
            fl_length = row[3]
            arcpy.AddMessage("Calculating ice thickness on centreline "+str(row[0])) ##+ ' of ' +str(len(flow_id_list)))
            #print("Calculating ice thickness on centreline "+str(row[0])) ##+ ' of ' +str(len(flow_id_list)))
            shape_factor_list = []
            yield_stress = row[2]
            points_glac = "in_memory\\points_glac"+str(row[0])
            flowline_id = row[0]
            outline_id = row[1]
            #arcpy.AddMessage("outline_id: " + str(outline_id))
            flowline_query = "line_id = "+str(flowline_id)
            outline_query = "outline_id = "+str(outline_id)

            single_flowline = arcpy.Select_analysis(flowline, "in_memory\\single_flowline", flowline_query) ## use the select analysis tool to replace all select by attribute and copy feature
            single_outline = arcpy.Select_analysis(outline, "in_memory\\single_outline", outline_query)  ###use the select analysis tool to replace all select by attribute and copy feature

            with arcpy.da.SearchCursor(single_outline, "va_thick") as cursor2:
                for row2 in cursor2:
                    mean_ice_depth = row2[0]
            del row2, cursor2
            
            search_dist_va = mean_ice_depth*15
            if search_dist_va > 2500:
                search_dist_va = 2500
            search_dist_va_round = int(round(search_dist_va/point_res)*point_res)
            if search_dist_va_round > fl_length*0.075:
                search_dist_va_round = int(round((fl_length*0.075)/point_res)*point_res)
            points_on_flowline = points_on_line_withID(single_flowline, point_res)
            points_on_flowline_altitude = ExtractValuesToPoints(points_on_flowline, dem, "in_memory\\points_altitude")                
            split_flowline = arcpy.SplitLineAtPoint_management(single_flowline, points_on_flowline_altitude, "in_memory\\split_flowline", "1 meters")

            perpendiculars = create_perpendiculars_line_sections(split_flowline, 1000000)

            clipped_perpendiculars = arcpy.Clip_analysis(perpendiculars, single_outline, "in_memory\\clipped_perpendiculars")
            simplified_perpendiculars = arcpy.SimplifyLine_cartography(clipped_perpendiculars,"in_memory\\simplified_perpendiculars","POINT_REMOVE", 1,"","NO_KEEP","NO_CHECK")
            split_simplified_perpendiculars = arcpy.SplitLine_management(simplified_perpendiculars,"in_memory\\split_simplified_perpendiculars")
            final_perpendiculars = arcpy.Erase_analysis(split_simplified_perpendiculars, split_flowline, "in_memory\\final_perpendiculars")
            

            arcpy.AddField_management(final_perpendiculars,"orig_perp_ID", "LONG")
            arcpy.AddField_management(final_perpendiculars, "EW", "DOUBLE")
            arcpy.AddField_management(final_perpendiculars, "FW", "DOUBLE")
            arcpy.AddField_management(final_perpendiculars,"inter", "SHORT")
            arcpy.CalculateField_management(final_perpendiculars,"orig_perp_ID",str("!"+str(arcpy.Describe(final_perpendiculars).OIDFieldName)+"!"),"PYTHON_9.3")
            arcpy.CalculateField_management(final_perpendiculars,"FW","!SHAPE.LENGTH!","PYTHON_9.3")
            slope_reclass_string = "0"+" "+str(effective_slope_limit)+" "+"1"+";"+str(effective_slope_limit)+" "+"90"+" "+"NODATA"

            dem_glac_clip = arcpy.Clip_management(dem, "#", "in_memory\\dem_glac_clip",single_outline,"","ClippingGeometry")
            #slope_raster = arcpy.Slope_3d(dem_glac_clip,"in_memory\\slope_raster", "DEGREE")
            slope_raster = Slope(dem_glac_clip)
            reclass_slope_raster = arcpy.Reclassify_3d(slope_raster,"Value",slope_reclass_string,"in_memory\\recls_slp","NODATA")
            reclass_slope_polygons = arcpy.RasterToPolygon_conversion(reclass_slope_raster,"in_memory\\reclass_slope_polygons")
            final_perpendiculars_EW = arcpy.Clip_analysis(final_perpendiculars,reclass_slope_polygons, "in_memory\\final_perpendiculars_EW")
            arcpy.CalculateField_management(final_perpendiculars_EW,"EW","!SHAPE.LENGTH!","PYTHON_9.3")

            final_perpendiculars_EW_layer = arcpy.MakeFeatureLayer_management(final_perpendiculars_EW, "in_memory\\final_perpendiculars_EW_layer")
            final_perpendiculars_layer = arcpy.MakeFeatureLayer_management(final_perpendiculars, "in_memory\\final_perpendiculars_layer")

            if ArcGISPro == 0: ##For ArcMAP
                arcpy.AddJoin_management(final_perpendiculars_layer, "orig_perp_ID", final_perpendiculars_EW_layer, "orig_perp_ID")
                arcpy.CalculateField_management(final_perpendiculars_layer,str("final_perpendiculars")+".EW","["+str("final_perpendiculars_EW")+".EW]")
                arcpy.RemoveJoin_management(final_perpendiculars_layer)
            else:
                arr=arcpy.da.FeatureClassToNumPyArray(final_perpendiculars_EW_layer, ["orig_perp_ID", "EW"])
                perpIdLst = [item[0] for item in arr]
                EWsLst = [item[1] for item in arr]
                with arcpy.da.UpdateCursor(final_perpendiculars_layer, ["orig_perp_ID", "EW"]) as perpcursor:
                    for perprow in perpcursor:
                        try:
                            idx = perpIdLst.index(perprow[0])
                            perprow[1] = EWsLst[idx]
                            cursor.updateRow(perprow)
                        except:
                            pass
                del perpcursor, perprow

            '''
            arr=arcpy.da.FeatureClassToNumPyArray(final_perpendiculars_EW, ["orig_perp_ID", "EW"])
            perpIdLst = [item[0] for item in arr]
            EWsLst = [item[1] for item in arr]
            #with arcpy.da.UpdateCursor(final_perpendiculars_layer, ["orig_perp_ID", "EW"]) as perpcursor:
            with arcpy.da.UpdateCursor(final_perpendiculars, ["orig_perp_ID", "EW"]) as perpcursor:
                for perprow in perpcursor:
                    try:
                        idx = perpIdLst.index(perprow[0])
                        perprow[1] = EWsLst[idx]
                        cursor.updateRow(perprow)
                    except:
                        pass
            del perpcursor, perprow
            '''
            
            cross_test = arcpy.Erase_analysis(flowline_layer_split_vertex, single_flowline, "in_memory\\cross_test") ## use erase function to replace the above processes

            arcpy.SelectLayerByLocation_management(final_perpendiculars_layer, 'intersect', cross_test)
            arcpy.CalculateField_management(final_perpendiculars_layer,"inter",1)
            arcpy.SelectLayerByAttribute_management(final_perpendiculars_layer,"CLEAR_SELECTION")

            arcpy.SpatialJoin_analysis(points_on_flowline_altitude, final_perpendiculars, points_glac, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")

            arcpy.AddField_management(points_glac, "f_type", "TEXT","","",1)
            fields = ["FW","EW","distance","RASTERVALU","max_dist", "ew_thi", "degrees", "sm_thi","Thickness","f_type","inter","p_g_alpha","poly_id","shape_f","bed_ele","shape_ori","order_fl","flow_id","seq"]
            exist_fields_list = [f.name for f in arcpy.ListFields(points_glac)] #List of current field names
            for field in fields:
                if field not in exist_fields_list:  #Add required fields if they do not already exist
                    arcpy.AddField_management(points_glac, field, "DOUBLE")

            #Loop #1 to find the max_distance
            max_dist_list  = []
            ele_list = []

            with arcpy.da.SearchCursor(points_glac, ["flow_id","distance","RASTERVALU"]) as cursor3:
                for row3 in cursor3:  
                    if row3[0] == flowline_id:
                        max_dist_list.append(row3[1])
                        ele_list.append(row3[2])
            del row3, cursor3
            
            max_dist = max(max_dist_list)
            min_dist = min(max_dist_list)
            min_ele = min(ele_list)

            #Loop #2 to assign max distance
            with arcpy.da.UpdateCursor(points_glac, ["flow_id","max_dist"]) as cursor4:
                for row4 in cursor4:
                    if row4[0] == flowline_id:
                        row4[1] = max_dist
                    cursor4.updateRow(row4)  
            del row4, cursor4

            #Loop #3 to create a dictionary for distance and DEM elavation (RasterVALU extracted from DEM, should be ice surface elevation)
            #error_dict = {}    
            search_dist = search_dist_va_round 
            #arcpy.AddMessage(str(search_dist))            
            ele_dict = {}
            with arcpy.da.SearchCursor(points_glac,fields) as cursor5:
                for row5 in cursor5:
                    index = int(row5[2])
                    value = int(row5[3])
                    ele_dict[index] = value       
            del row5, cursor5
            
            if min_dist > 0: ##if the minimum distance is not zero, add the zero and min-elevation into the dict
                index = 0
                value = int(min_ele)
                ele_dict[index] = value 
                
            #Loop #4 to derive all parameters
            
            with arcpy.da.UpdateCursor(points_glac,fields) as cursor6:
                for row6 in cursor6:
                    #row[16] = order_dict[row[17]]  #?????? error
                    row6[16] = flowline_id  #Revised by Yingkui Li 11052019 ##May need to check if it is correct order_dict is a local variable in shearstress, how to pass the value out to here????
                    row6[12] = outline_id
                    start_distance = int(row6[2]) - (search_dist)
                    end_distance = int(row6[2]) + (search_dist)
                    if start_distance >= 0 and end_distance <= row6[4]:
                        diff_distance = end_distance - start_distance
                        start = int(start_distance)
                        start_ele = ele_dict[start]
                        end = int(end_distance)
                        #arcpy.AddMessage(str(end))
                        end_ele = ele_dict[end]
                        diff_ele = end_ele - start_ele
                        if diff_distance != 0:
                            m = diff_ele/diff_distance
                        else:
                            m = 0
                        radians = math.atan(m)                
                        degrees = abs(radians*(180.0/math.pi))
                        row6[6] = degrees
                        if degrees <= min_slope:            #can be replace by the max of degrees and min_slope
                            degrees = min_slope
                        rad_limit = degrees*(math.pi/180.0)
                        pgtana = ice_density*gravity*math.tan(rad_limit)
                        ice_depth_sm = yield_stress/pgtana ###Benn 2010 excel paper 11/22/2020
                        row6[7] = ice_depth_sm
                        row6[11] = pgtana
                        if row6[1] != None: ##row6[1] > 0:                              
                            row6[5] = (0.9*(row6[1]/2.0)*ice_depth_sm)/(0.9*(row6[1]/2.0)-ice_depth_sm) # Equation 4 in the paper
                        else: ## Added for ArcGIS Pro version
                            row6[5] = 0
                        row6[8] = row6[5]
                        row6[9] = "W"
                        cursor6.updateRow(row6)  ####Why update row here, can update at the end of this loop.
                        if row6[5] <= 0 or row6[1] == None or row6[10] ==1 or row6[0] > 0.66*(max_length_dict[row6[12]]):
                            row6[9] = "A"
                        if row6[10] == 1:
                            row6[15] = 1
                        if row6[1] == None:
                            row6[15] = 2
                        if row6[5] <= 0:
                            row6[15] = 3
                        if row6[0] > 0.66*(max_length_dict[row6[12]]): 
                            row6[15] = 4
                        if row6[9] == "W":
                            shape_factor = yield_stress/(row6[5]*pgtana )
                            row6[13] = shape_factor
                            row6[15] = shape_factor
                            if shape_factor <= 0.445:
                                row6[9] = "A"
                        if row6[9] == "W":
                            shape_factor_list.append(shape_factor)
                        if row6[9] == "A":
                            shape_factor_list.append(0.8)
                        #update row just one time
                        cursor6.updateRow(row6)   
                    else: ##this is the palce to delete points: Need to modify it later.
                        cursor6.deleteRow()    ####delete the row does not meet the conditions
            del row6, cursor6

            ##determine the mean shape factor
            if len(shape_factor_list) == 0:
                mean_shape_factor = 0.8
            else:
                mean_shape_factor = np.mean(shape_factor_list)
            if mean_shape_factor < 0 or mean_shape_factor > 0.95:
                mean_shape_factor = 0.8      ####why set shapefactor to 0.8????


            #Loop #5 adjust the thickness based on the shape factor
            with arcpy.da.UpdateCursor(points_glac,fields) as cursor7:
                for row7 in cursor7:
                    if row7[9] == "A":
                        row7[8] = (yield_stress/(mean_shape_factor*row7[11]))  ###can simply adjust the row [8] by divide the shape factor, row[11] is not necessary
                        row7[13] = mean_shape_factor  
                    row7[14] = row7[3] - row7[8]
                    cursor7.updateRow(row7)   #just one update is enough
            del row7, cursor7

            #delete the extra fields created                    
            del_fields = ["Join_Count","TARGET_FID","distance","flow_id","orig_perp_ID","inter","max_dist","ew_thi","sm_thi","p_g_alpha","shape_ori"]
            exist_fields_list = [f.name for f in arcpy.ListFields(points_glac)] #List of current field names
            for field in del_fields:
                if field in exist_fields_list:  #delete required fields 
                    arcpy.DeleteField_management(points_glac, field)

            #the first time, copy feature, others, append feature
            if counter < 1:
                initial_fl_points = arcpy.CopyFeatures_management(points_glac, "in_memory\\initial_fl_points")
            else:
                merged_points = arcpy.Append_management(points_glac, initial_fl_points)

            counter = counter + 1
    #arcpy.AddMessage("Counter is: " + str(counter))
    del row, cursor
    
    if (counter == 1): #Handle just have one glacier outline and center line issue; need to define merge_points
       merged_points = initial_fl_points

    identical_table = "in_memory\\identical_table"
    points_layer = "in_memory\\points_layer"
    arcpy.FindIdentical_management(merged_points, identical_table, ["Shape"], output_record_option="ONLY_DUPLICATES")
    arcpy.MakeFeatureLayer_management(merged_points, points_layer)
    if ArcGISPro == 0:
        arcpy.AddJoin_management(points_layer,arcpy.Describe(points_layer).OIDFieldName,identical_table, "IN_FID")
        arcpy.CalculateField_management(points_layer,"initial_fl_points.seq","!identical_table.FEAT_SEQ!","PYTHON_9.3")
        arcpy.RemoveJoin_management (points_layer)
    else:
        JoinID = arcpy.Describe(merged_points).OIDFieldName
        arr=arcpy.da.FeatureClassToNumPyArray(identical_table, ["IN_FID", "FEAT_SEQ"])
        IdLst = [item[0] for item in arr]
        SEQsLst = [item[1] for item in arr]
        with arcpy.da.UpdateCursor(merged_points, [JoinID, "seq"]) as pntcursor:
            for pntrow in pntcursor:
                try:
                    idx = IdLst.index(pntrow[0])
                    pntrow[1] = SEQsLst[idx]
                    cursor.updateRow(pntrow)
                except:
                    pass
        del pntcursor, pntrow


    seq_list = []
    with arcpy.da.SearchCursor(merged_points, ["seq","order_fl"]) as cursor:
        for row in cursor:
            if row[0] != None:
                if row[0] >= 1:
                    seq_list.append(row[0])
    del row, cursor

    for value in seq_list:
        order_list = []

        with arcpy.da.SearchCursor(merged_points, ["seq","order_fl"]) as cursor:
            for row in cursor:
                if row[0] == value:
                    order_list.append(row[1])
        del row, cursor
        '''
        #Check if this can reduce the number of points along the flowline: Yingkui 2/9/23
        with arcpy.da.UpdateCursor(merged_points, ["seq","order_fl"]) as cursor:   ###delete some points??????
            for row in cursor:
                if row[0] == value:
                    if row[1] != min(order_list):
                        cursor.deleteRow() ##This is another place to delete points
        del row, cursor
        '''
    arcpy.DeleteField_management(merged_points, "seq")
    arcpy.CopyFeatures_management(merged_points, final_points)

    #arcpy.AddMessage("Centreline point thickness calculation complete")

    if (interpolate_check == "true" and len(raster_out) > 0):  #make sure the output raster has been assigned
        arcpy.AddMessage("Interpolating ice thickness raster...")
        print("Interpolating ice thickness raster...")

        dissolve_outlines = arcpy.Dissolve_management(outline, "in_memory\\dissolve_outlines")
        singlepart_outlines = arcpy.MultipartToSinglepart_management(dissolve_outlines, "in_memory\\singlepart_outlines")
        outline_lines_in = arcpy.PolygonToLine_management(singlepart_outlines, "in_memory\\outlines_line_in")
        arcpy.AddField_management(outline_lines_in, "contour", "SHORT")
        arcpy.CalculateField_management(outline_lines_in,"contour",0)
        if (cellsize_interpolate_check == "true" and len(cellsize_interpolate_user_spec) > 0) :
            #cellsize_interp = int(cellsize_interpolate_user_spec)
            cellsize_interp = float(cellsize_interpolate_user_spec)
        else:
            cellsize_interp = cellsize_float

        interpolated_ice_depth_original = TopoToRaster([TopoPointElevation([[merged_points, 'Thickness']]), TopoContour([[outline_lines_in, 'contour']]), TopoBoundary ([singlepart_outlines])], cellsize_interp, "", '20', '0', '#', 'NO_ENFORCE', "SPOT", '1', '#', '1', '0', '0', '200') ##,"","","","",0.1)
        if ArcGISPro == 0:
            arcpy.CopyRaster_management (interpolated_ice_depth_original, raster_out)
        else:
            interpolated_ice_depth_original.save(raster_out)

###Start the main program
if __name__ == '__main__':

    arcpy.Delete_management("in_memory")

    flowline = arcpy.GetParameterAsText(0)
    dem = arcpy.GetParameterAsText(1)

    cellsize = arcpy.GetRasterProperties_management(dem,"CELLSIZEX")
    cellsize_float = float(cellsize.getOutput(0)) # use float cell size

    outline = arcpy.GetParameterAsText(2)
    ice_density = int(arcpy.GetParameterAsText(3))
    slope_limit = arcpy.GetParameterAsText(4)
    min_slope = float(arcpy.GetParameterAsText(5))
    point_res_check = arcpy.GetParameterAsText(6)

    if point_res_check == "true":
        point_res = int(cellsize_float+0.5) ##cellsize_int  (old value)   
    else:
        point_res = float(arcpy.GetParameterAsText(7))

    shear_stress_test = arcpy.GetParameterAsText(8)
    shear_stress_value = arcpy.GetParameterAsText(9)
    final_points = arcpy.GetParameterAsText(10)
    #interpolate_check = arcpy.GetParameterAsText(11)
    #cellsize_interpolate_check = arcpy.GetParameterAsText(12)
    #cellsize_interpolate_user_spec = arcpy.GetParameterAsText(13)

    raster_out = arcpy.GetParameterAsText(11)

    ##Setup the default values fice_densityor other parameters of the ice thickness function
    interpolate_check = "true"
    cellsize_interpolate_check = "true"
    cellsize_interpolate_user_spec = ""



    ##make sure the projection of the glacier outlines is the same with the UTM and the same with the DEM projection 
    spatial_flowline = arcpy.Describe(flowline).spatialReference
    spatial_outline = arcpy.Describe(outline).spatialReference
    spatial_ref_dem = arcpy.Describe(dem).spatialReference

    '''
    if "UTM" in spatial_ref_dem.name:
        arcpy.AddMessage("The DEM projection is: " + spatial_ref_dem.name)
    else:
        arcpy.AddMessage("The DEM projection is not UTM. Please re-project the DEM to a UTM projection for the analysis!")
        exit()   
    if "UTM" in spatial_flowline.name:
        arcpy.AddMessage("The flowline projection is: " + spatial_flowline.name)
    else:
        arcpy.AddMessage("The flowline projection is not UTM. Please re-project the flowline to a UTM projection for the analysis!")
        exit()   

    if "UTM" in spatial_outline.name:
        arcpy.AddMessage("The outline projection is: " + spatial_outline.name)
    else:
        arcpy.AddMessage("The outline projection is not UTM. Please re-project the outline to a UTM projection for the analysis!")
        exit()   
    '''

    if spatial_ref_dem.linearUnitName == "Meter":
        arcpy.AddMessage("The DEM projection is: " + spatial_ref_dem.name)
        print("The DEM projection is: " + spatial_ref_dem.name)
    else:
        arcpy.AddMessage("The unit of the DEM projection is not in meter. Please re-project the DEM to a projected coordinate system for the analysis!")
        print("The DEM projection is: " + spatial_ref_dem.name)
        exit()   

    if spatial_flowline.linearUnitName == "Meter":
        arcpy.AddMessage("The flowline projection is: " + spatial_flowline.name)
        print("The DEM projection is: " + spatial_ref_dem.name)
    else:
        arcpy.AddMessage("The unit of the flowline projection is not in meter. Please re-project it to a projected coordinate system for the analysis!")
        print("The DEM projection is: " + spatial_ref_dem.name)
        exit()   


    #if "UTM" in spatial_outline.name:
    if spatial_outline.linearUnitName == "Meter":
        arcpy.AddMessage("The outline projection is: " + spatial_outline.name)
        print("The DEM projection is: " + spatial_ref_dem.name)
    else:
        arcpy.AddMessage("The unit of the outline projection is not in meter. Please re-project it to a projected coordinate system for the analysis!")
        print("The DEM projection is: " + spatial_ref_dem.name)
        exit()   

    #arcpy.AddMessage(str(spatial_ref_dem.PCSCode))
    #arcpy.AddMessage(str(spatial_flowline.PCSCode))
    #arcpy.AddMessage(str(spatial_outline.PCSCode))
    if (spatial_ref_dem.PCSCode == spatial_outline.PCSCode) and (spatial_ref_dem.PCSCode == spatial_flowline.PCSCode):
    
        arcpy.AddMessage("DEM, flowline, and outlines all have the same projected coordinate system: " + spatial_ref_dem.name)
        print("The DEM projection is: " + spatial_ref_dem.name)
    else:
        arcpy.AddMessage("DEM ,flowline,or outlines have different map projections. Please re-project the datasets to the same projection!")
        print("The DEM projection is: " + spatial_ref_dem.name)
        exit()   

    ##Need to resample the DEM to a DEM with integer resolution
    ##Maybe it is not necessary 5/7/2024
    '''
    indem = "in_memory\\indem"
    if cellsize_float > int (cellsize_float):
        arcpy.AddMessage("Resample DEM to an integer resolution DEM")
        new_cell_size = int(cellsize_float+0.5)
        #arcpy.AddMessage(str(new_cell_size))
        arcpy.Resample_management(dem, indem, str(new_cell_size), "NEAREST")
        #arcpy.CopyRaster_management(indem, "d:\\temp\\indem.tif")
    else:
        arcpy.CopyRaster_management(dem, indem)
    
    int_dem = Int(indem)
    '''

    #Ice_Thickness_Volta (flowline, int_dem, outline, ice_density, slope_limit, min_slope, point_res_check, point_res, shear_stress_test, shear_stress_value, final_points,
    #                     interpolate_check, cellsize_interpolate_check, cellsize_interpolate_user_spec, raster_out)

    arcpy.AddMessage("Volta with float cell resolution")

    Ice_Thickness_Volta (flowline, Raster(dem), outline, ice_density, slope_limit, min_slope, point_res_check, point_res, shear_stress_test, shear_stress_value, final_points,
                         interpolate_check, cellsize_interpolate_check, cellsize_interpolate_user_spec, raster_out)

    arcpy.Delete_management("in_memory")
                            

