#-------------------------------------------------------------------------------
# Name: Volta Centerline
# Purpose: This tool derives centerlines for modern glaciers based on the method
# proposed by James and Carrivick (2016) in VOLTA. 
# The codes here are modified based on the original codes in VOLTA. The modification includes:
# 1) Check if the long axis apporach is suitable for the centerline delineation, if not, use the
#    highest and lowest points instead
# 2) extend the centerlines to the lowest point of the glacier outlines
# 3) The centerline network is revised to remove the duplicated parts in main valleys
# 4) Fixed runtime errors
# 
# Author:      Yingkui Li
# Created:     08/30/2020 to 12/21/2020
# finilized:   06/24/2022
# Copyright:   (c) Yingkui Li 2022
# Department of Geography, University of Tennessee
# Knoxville, TN 37996, USA
#-------------------------------------------------------------------------------

#import SharedFunctions  
from SharedFunctions import *  

def new_branch_bak(input_flowline, dem, branch_outline, input_outline, cellsize_float, TributaryRatio, TributarySourceThreshold, glacier_outline):  
    #clipped_dem_raster = arcpy.sa.ExtractByMask(dem,input_outline)
    clipped_dem_raster = ExtractByMask(dem,input_outline)
    #up_list = []

    arcpy.env.extent = clipped_dem_raster
    #arcpy.env.SnapRaster = clipped_dem_raster
    #global new_branch_count
    new_branch_count = 0
    dem_numpy = arcpy.RasterToNumPyArray(clipped_dem_raster,"","","",0)
    d8_structure = ndimage.generate_binary_structure(2, 2)
    upstream_area_original = 1
    glacier_numpy = dem_numpy.copy()
    glacier_numpy[glacier_numpy != 0] = 1
    with arcpy.da.SearchCursor(input_flowline, ["SHAPE@LENGTH"]) as cursor:
        for row in cursor:
            flowline_length = row[0]
    del row, cursor

    flowline_point_res = flowline_length/100.0  ## Need to make a smallest distance for this, if glacier is smaller,  it is not needed to divide 100
    #flowline_point_res = max(flowline_point_res, cellsize_float*3) ##set DEM resolution (30 meter) as the minimum distance by Yingkui Li

    flowline_points = points_on_line(input_flowline, flowline_point_res)
    #Added Debug line
    #arcpy.CopyFeatures_management(flowline_points, "C:/temp/point.shp") # pased
    branch_points_alt = ExtractValuesToPoints(flowline_points, dem, "in_memory\\branch_points_altitude")
    #arcpy.CopyFeatures_management(branch_points_alt, "C:/temp/pointalt.shp") #passed
    branch_points_alt_layer = arcpy.MakeFeatureLayer_management(branch_points_alt, "in_memory\\branch_points_layer")
    #with arcpy.da.SearchCursor(branch_points_alt_layer, [arcpy.Describe(branch_points_alt_layer).OIDFieldName, "RASTERVALU"]) as cursor:
    FcID = arcpy.Describe(branch_points_alt).OIDFieldName
    with arcpy.da.SearchCursor(branch_points_alt, [FcID, "RASTERVALU"]) as cursor:
        for row in cursor:
            if row[0] > 1:
                tempRs = arcpy.env.scratchFolder + "\\r" + str(row[0])
                #tempRs = flow_line_path + "\\r" + str(row[0])
                #if arcpy.Exists (tempRs):
                #    arcpy.Delete_management(tempRs)
                if new_branch_count == 0:
                    #query = arcpy.Describe(branch_points_alt_layer).OIDFieldName+" = "+str(row[0])

                    query = FcID +" = "+str(row[0])
                    arcpy.SelectLayerByAttribute_management(branch_points_alt_layer, "NEW_SELECTION", query)
                    #arcpy.Select_analysis(branch_points_alt,"in_memory\\sel_branch_point",query)

                    arcpy.PointToRaster_conversion(branch_points_alt_layer, "RASTERVALU", tempRs, "","", cellsize_float)
                    #raster_flowline_general1 = arcpy.PointToRaster_conversion("in_memory\\sel_branch_point", "RASTERVALU", tempRs, "","", cellsize_float)
                    #arcpy.PointToRaster_conversion("in_memory\\sel_branch_point", "RASTERVALU", tempRs, "","", cellsize_float)

                    #raster_flowline_general = arcpy.Raster (raster_flowline_general1)
                    #arcpy.AddMessage("Pass = " + str(row[0]))
                    #flowline_numpy = arcpy.RasterToNumPyArray(raster_flowline_general,"","","",0)
                    flowline_numpy = arcpy.RasterToNumPyArray(tempRs,"","","",0)
                    #arcpy.AddMessage("Pass2 " + str(row[0]))

                    original_cell = np.amax(flowline_numpy)
                    itemindex = np.where(flowline_numpy == original_cell)
                    row = itemindex[0]
                    col = itemindex[1]
                    dem_copy = dem_numpy.copy()
                    dem_copy[dem_copy < original_cell] = 0
                    dem_copy[dem_copy >= original_cell] = 1
                    labeled_array, num_features = scipy.ndimage.measurements.label(dem_copy, structure=d8_structure)
                    zone = labeled_array[row,col]
                    zone_cells = np.argwhere(labeled_array == zone)
                    total_cells = int(len(zone_cells))
                    #upstream_area = total_cells*cellsize_int*cellsize_int
                    upstream_area = total_cells*cellsize_float*cellsize_float
                    upstream_area_change = upstream_area - upstream_area_original
                    percentage_change = (float(upstream_area)/float(upstream_area_original))*100
                    #up_list.append(upstream_area)
                    upstream_area_original = upstream_area
                    #delete the tempRs data
                    arcpy.Delete_management(tempRs)
                    
                    if upstream_area < TributarySourceThreshold or percentage_change < (100.0 + TributaryRatio*100) or percentage_change >500: 
                        newarray = labeled_array.copy()
                        newarray[newarray != zone] = 0
                        newarray[newarray == zone] = 1
                        glacier_numpy[newarray == 1] = 0
                        #arcpy.AddMessage("Return empty!!")
                    else:
                        search_branch = 1
                        left_co_ord = arcpy.GetRasterProperties_management(clipped_dem_raster,"LEFT")
                        bottom_co_ord = arcpy.GetRasterProperties_management(clipped_dem_raster,"BOTTOM")
                        left_co_ord_float = float(left_co_ord.getOutput(0))
                        bottom_co_ord_float = float(bottom_co_ord.getOutput(0))
                        bottom_left_point = arcpy.Point(left_co_ord_float,bottom_co_ord_float)
                        new_branch_raster = Int(arcpy.NumPyArrayToRaster(glacier_numpy, bottom_left_point, cellsize_float, cellsize_float, 0))
                        new_branch_outline_unclipped = arcpy.RasterToPolygon_conversion(new_branch_raster, "in_memory\\new_branch_outline_unclipped", "SIMPLIFY")
                        new_branch_outline = arcpy.Clip_analysis(new_branch_outline_unclipped, glacier_outline, branch_outline)
                        new_branch_count = 1
                        break
                        #arcpy.AddMessage("Return new_branch_outline!!")
                        #return new_branch_outline, new_branch_count
    del row, cursor

    if new_branch_count > 0:
        #arcpy.AddMessage("Return new_branch_outline!!")
        return new_branch_outline, new_branch_count
    else:
        arcpy.AddMessage("Return no new branch_outline!!")
        return None, new_branch_count
        


#------------------------------------------------------------------------------------------------------------
# This is the main function to derive the centerlines from glacier outlines and DEM.
# This function can be call by the overall model
#------------------------------------------------------------------------------------------------------------
def Centerline_Volta (glacier_outlines, dem, TributaryRatio, tributary_threshold, flow_line):
    cellsize = arcpy.GetRasterProperties_management(dem,"CELLSIZEX")
    cellsize_float = float(cellsize.getOutput(0))
    arcpy.env.outputCoordinateSystem = dem
    arcpy.env.snapRaster = dem

    flow_line_name = "flowlinecp"


    ###Start the main program
    arcpy.env.extent = glacier_outlines

    ##First make sure the outline does not contain any spurious polygons
    minarea = cellsize_float * cellsize_float * 5.0 ##five size of the pixels
    with arcpy.da.UpdateCursor(glacier_outlines, ['SHAPE@AREA']) as cursor:
        for row in cursor:
            area = float(row[0])
            if area < minarea: ##This step reomve the small spurious ploygons as well
                #arcpy.AddMessage("delete spurious polygons..." )
                cursor.deleteRow()
    del cursor, row

    FcID = arcpy.Describe(glacier_outlines).OIDFieldName

    arr=arcpy.da.FeatureClassToNumPyArray(glacier_outlines, FcID)
    FIds = np.array([item[0] for item in arr])

    progress_counter = 1

    glacier_outline = "in_memory\\layer_test_hope" ##Single outline for each loop
    flow_line_final = arcpy.CreateFeatureclass_management(arcpy.env.scratchGDB, flow_line_name, "POLYLINE")
    #flow_line_final_individual = arcpy.CreateFeatureclass_management(arcpy.env.scratchGDB, flow_line_name, "POLYLINE")  ##Maybe not necessary??, so the output is flow_line_final now

    for i in range(len(FIds)):
        flow_line_final_individual = arcpy.CreateFeatureclass_management("in_memory","flow_line_final_individual", "POLYLINE")
        arcpy.AddMessage("Generating centreline(s) for glacier "+str(progress_counter)+" of "+str(len(FIds)))
        query = FcID +" = "+str(FIds[i])

        arcpy.Select_analysis(glacier_outlines,glacier_outline,query)
        outline_area = 0
        with arcpy.da.SearchCursor(glacier_outline, ["SHAPE@AREA"]) as cursor:
            for row in cursor:
                outline_area = row[0]
        del row, cursor
        
        sourcearea_threshold = outline_area * tributary_threshold

        ##extract by mask the DEM so that the analysis is only for the extracted DEM
        arcpy.Buffer_analysis(glacier_outline, "in_memory\\glacier_outline_buf", "100 Meter")
        clipped_dem = ExtractByMask(dem,"in_memory\\glacier_outline_buf")
        #clipped_dem = dem
        
        ##Create the axis
        axis, axispoint = create_axis(glacier_outline, clipped_dem)
        extended_axis = extend_line(axis, 0.25)  ###Extend the line bigger enough so that the perp lines can cover all glacier coverage

        #get the length of the axis and use it to determine the the step of the perpendicular distance
        with arcpy.da.SearchCursor(axis, ["SHAPE@LENGTH"]) as cursor:
            for row in cursor:
                axis_length = row[0]
        del row, cursor
        
        perp_space = min(max(int(axis_length/3.0/cellsize_float), 1),  3) * cellsize_float
        perpendiculars = create_perpendiculars(extended_axis, perp_space, 100000) ##Change to 100 km long


        central_points_with_alt = central_points_on_perps(perpendiculars, glacier_outline, axispoint, clipped_dem)

        initial_flowline = flowline(central_points_with_alt, clipped_dem, "in_memory\\initial_flow_line",  glacier_outline, 1000)
        #arcpy.CopyFeatures_management(initial_flowline, "d:\\temp\\initial_flowline.shp")

        new_branch_poly, new_branch_count  = new_branch(initial_flowline, clipped_dem, "in_memory\\new_outline", glacier_outline, TributaryRatio, sourcearea_threshold)
        #new_branch_poly, new_branch_count  = new_branch(initial_flowline, clipped_dem, "in_memory\\new_outline", glacier_outline, cellsize_float, TributaryRatio, sourcearea_threshold, glacier_outline)
        #arcpy.CopyFeatures_management(new_branch_poly, "d:\\temp\\new_branch_poly.shp")
        #arcpy.AddMessage("the number of new branch is:" +str(new_branch_count))
        if new_branch_count == 0:
            arcpy.AddMessage("Add one line!")
            arcpy.Append_management(initial_flowline, flow_line_final, "NO_TEST")

        ##Back to the old codes
        while new_branch_count == 1:
            arcpy.AddMessage("Multiple Branches Detected")
            secondary_axis, axispoint = create_axis(new_branch_poly, clipped_dem)
            secondary_axis_count = arcpy.GetCount_management(secondary_axis)          #check outline is a single feature - count"
            secondary_axis_result = int(secondary_axis_count.getOutput(0))      #check outline is a single feature - get result"            
            if (secondary_axis_result > 0):
                secondary_extended_axis = extend_line(secondary_axis, 0.25)
                secondary_perpendiculars = create_perpendiculars(secondary_extended_axis, perp_space, 100000) 
                secondary_central_points_with_alt = central_points_on_perps(secondary_perpendiculars, new_branch_poly, axispoint, clipped_dem)
                secondary_flowline = flowline(secondary_central_points_with_alt, clipped_dem, "in_memory\\secondary_flow_line", new_branch_poly, 1000)
                secondary_flowline_copied = arcpy.CopyFeatures_management(secondary_flowline, "in_memory\\secondary_flow_line_copied") ##For some reason, this is necessary!!!

                corrected_poly, new_branch_count = new_branch(secondary_flowline, clipped_dem, "in_memory\\corrected_poly1", glacier_outline, TributaryRatio, sourcearea_threshold)
                #corrected_poly, new_branch_count = new_branch(secondary_flowline, clipped_dem, "in_memory\\corrected_poly1", glacier_outline, cellsize_float, TributaryRatio, sourcearea_threshold, glacier_outline)
                if new_branch_count == 1:
                    corrected_axis, correctedaxispoint = create_axis(corrected_poly, clipped_dem)
                    corrected_extended_axis = extend_line(corrected_axis, 0.25)
                    corrected_perpendiculars = create_perpendiculars(corrected_extended_axis, perp_space, 100000) 
                    corrected_central_points_with_alt = central_points_on_perps(corrected_perpendiculars, corrected_poly, axispoint, dem)
                    corrected_flowline = flowline(corrected_central_points_with_alt, clipped_dem, "in_memory\\corrected_additional_flowline", corrected_poly, 1000)
                    corrected_flowline_copied = arcpy.CopyFeatures_management(corrected_flowline, "in_memory\\corrected_flowline_copied")

                    flowline_intersect_points = arcpy.Intersect_analysis([corrected_flowline_copied, secondary_flowline_copied],"in_memory\\flow_intersect_points","","","POINT")
                    secondary_fl_split = arcpy.SplitLineAtPoint_management(secondary_flowline_copied, flowline_intersect_points, "in_memory\\secondary_fl_split", 5)
                    #secondary_fl_split_layer = arcpy.MakeFeatureLayer_management(secondary_fl_split, "in_memory\\secondary_fl_split_layer")
                    '''
                    max_midpoint = 0.0
                    selected_secondary_fl = "in_memory\\selected_secondary_fl"
                    #with arcpy.da.UpdateCursor(secondary_fl_split_layer, arcpy.Describe(secondary_fl_split_layer).OIDFieldName) as cursor:
                    with arcpy.da.UpdateCursor(secondary_fl_split, arcpy.Describe(secondary_fl_split).OIDFieldName) as cursor:
                        for row in cursor:
                            query = arcpy.Describe(secondary_fl_split).OIDFieldName+"="+str(row[0])
                            #arcpy.SelectLayerByAttribute_management(secondary_fl_split_layer, "NEW_SELECTION", query)
                            arcpy.Select_analysis(secondary_fl_split, selected_secondary_fl,query)
                            midpoint_seecondary_fl_split = arcpy.FeatureVerticesToPoints_management(selected_secondary_fl, "in_memory\\midpoint_seecondary_fl_split", "MID")
                            midpoint_seecondary_fl_split_ele = ExtractValuesToPoints(midpoint_seecondary_fl_split, dem, "in_memory\\midpoint_seecondary_fl_split_elevation")
                            with arcpy.da.UpdateCursor(midpoint_seecondary_fl_split_ele, "RASTERVALU") as point_cursor:
                                for point_row in point_cursor:
                                    elevation = point_row[0]
                            del point_row, point_cursor
                            if elevation > max_midpoint:
                                max_midpoint = elevation
                            else:
                                cursor.deleteRow()
                            #arcpy.SelectLayerByAttribute_management(secondary_fl_split_layer, "CLEAR_SELECTION")
                    del row, cursor
                    '''  
                    midpoint_secondary_fl_split = arcpy.FeatureVerticesToPoints_management(secondary_fl_split, "in_memory\\midpoint_secondary_fl_split", "MID")
                    midpoint_secondary_fl_split_ele = ExtractValuesToPoints(midpoint_secondary_fl_split, clipped_dem, "in_memory\\midpoint_secondary_fl_split_ele")
                    max_midpoint = 0.0
                    FcID1 = arcpy.Describe(midpoint_secondary_fl_split_ele).OIDFieldName
                    DeleteID = []
                    with arcpy.da.SearchCursor(midpoint_secondary_fl_split_ele, ["RASTERVALU", FcID1]) as point_cursor:
                        for point_row in point_cursor:
                            elevation = point_row[0]
                            if elevation > max_midpoint:  ##Record the highest elevation section and its Fc ID
                                max_midpoint = elevation
                            else:
                                DeleteID.append(point_row[1])
                    del point_row, point_cursor
                    
                    FcID2 = arcpy.Describe(secondary_fl_split).OIDFieldName
                    selected_secondary_fl = "in_memory\\selected_secondary_fl"
                    with arcpy.da.UpdateCursor(secondary_fl_split, FcID1) as point_cursor:
                        for point_row in point_cursor:
                            if point_row[0] in DeleteID:
                                point_cursor.deleteRow()
                    del point_row, point_cursor

                    new_branch_poly, new_branch_count = new_branch(secondary_flowline, clipped_dem, new_branch_poly, new_branch_poly, TributaryRatio, sourcearea_threshold)
                    #new_branch_poly, new_branch_count = new_branch(secondary_flowline, clipped_dem, new_branch_poly, new_branch_poly, cellsize_float, TributaryRatio, sourcearea_threshold, glacier_outline)
                    #arcpy.Append_management([secondary_fl_split_layer,corrected_flowline_copied], flow_line_final_individual, "NO_TEST")
                    arcpy.Append_management([secondary_fl_split,corrected_flowline_copied], flow_line_final_individual, "NO_TEST")
                    #new_branch_count = 0 ##Force to stop if????
                    ##Good to here!!! 12/31/2022
                    '''
                    ##The following steps just want to make the flowline from the start to the lowest point; it is not necessary if do not want it
                    centreline_count_glacier = arcpy.GetCount_management(flow_line_final_individual) #COUNT HOW MANY centrelines IN glacier
                    centreline_count_result = int(centreline_count_glacier.getOutput(0)) #GET INT RESULT
                    split_points = arcpy.FeatureVerticesToPoints_management(flow_line_final_individual, "in_memory\\flowline_split_points", "BOTH_ENDS")
                    split_flolwine_glacier_ori = arcpy.SplitLineAtPoint_management (flow_line_final_individual, split_points, "in_memory\\split_flolwine_glacier_ori", 1)
                    #split_flolwine_glacier = arcpy.SplitLineAtPoint_management (flow_line_final_individual, split_points, "in_memory\\split_flolwine_glacier", 1)

                    split_flolwine_glacier = arcpy.MakeFeatureLayer_management(split_flolwine_glacier_ori, "in_memory\\split_flolwine_glacier")

                    start_ele_list = []
                    start_ele_dict = {}
                    with arcpy.da.SearchCursor(split_flolwine_glacier, [arcpy.Describe(split_flolwine_glacier).OIDFieldName]) as cursor:
                        for row in cursor:
                            flowline_id = row[0]
                            flowline_query = arcpy.Describe(split_flolwine_glacier).OIDFieldName+" = "+str(flowline_id)
                            arcpy.SelectLayerByAttribute_management(split_flolwine_glacier,"NEW_SELECTION",flowline_query)
                            #arcpy.Select_analysis(split_flolwine_glacier, selected_secondary_fl,flowline_query)
                            line_startpoint = arcpy.FeatureVerticesToPoints_management(split_flolwine_glacier, "in_memory\\line_in_startpoint", "START")
                            line_endpoint = arcpy.FeatureVerticesToPoints_management(split_flolwine_glacier, "in_memory\\line_in_endpoint", "END")
                            startpoint_ele_point = ExtractValuesToPoints(line_startpoint, dem,"in_memory\\startpoint_point_ele") #Get elevation of each point
                            endpoint_ele_point = ExtractValuesToPoints(line_endpoint, dem,"in_memory\\endpoint_point_ele") #Get elevation of each point
                            with arcpy.da.SearchCursor(startpoint_ele_point, "RASTERVALU") as cursor2:
                                for row2 in cursor2:
                                    startpoint_ele = row2[0]
                            del row2, cursor2
                            with arcpy.da.SearchCursor(endpoint_ele_point, "RASTERVALU") as cursor3:
                                for row3 in cursor3:
                                    endpoint_ele = row3[0]
                            del row3, cursor3
                            if startpoint_ele > endpoint_ele:
                                arcpy.FlipLine_edit(split_flolwine_glacier)
                                startpoint_checked = line_endpoint
                                endpoint_checked = line_startpoint
                                startpoint_ele_checked = endpoint_ele
                                endpoint_ele_checked = startpoint_ele
                            else:
                                startpoint_checked = line_startpoint
                                endpoint_checked = line_endpoint
                                startpoint_ele_checked = startpoint_ele
                                endpoint_ele_checked = endpoint_ele
                            start_ele_list.append(startpoint_ele_checked)
                            start_ele_dict[row[0]] = startpoint_ele_checked
                    del row, cursor
                    min_startpoint_ele = min(start_ele_list)
                    arcpy.SelectLayerByAttribute_management(split_flolwine_glacier,"CLEAR_SELECTION")

                    dangle_points = arcpy.FeatureVerticesToPoints_management(split_flolwine_glacier, "in_memory\\dangle_points", "DANGLE")
                    start_points = arcpy.FeatureVerticesToPoints_management(split_flolwine_glacier, "in_memory\\start_points", "START")
                    dangle_points_layer = arcpy.MakeFeatureLayer_management(dangle_points, "in_memory\\dangle_points_layer")
                    start_points_layer = arcpy.MakeFeatureLayer_management(start_points, "in_memory\\start_points_layer")
                    arcpy.SelectLayerByLocation_management(dangle_points_layer, "ARE_IDENTICAL_TO", start_points_layer)
                    arcpy.SelectLayerByLocation_management(dangle_points_layer, "","","","SWITCH_SELECTION")
                    dangle_count_result = arcpy.GetCount_management(dangle_points_layer)
                    dangle_count = int(dangle_count_result.getOutput(0))
                    arcpy.SelectLayerByLocation_management(split_flolwine_glacier, "BOUNDARY_TOUCHES", dangle_points_layer)
                    with arcpy.da.SearchCursor(split_flolwine_glacier, [arcpy.Describe(split_flolwine_glacier).OIDFieldName]) as cursor:
                        for row in cursor:
                            start_ele = start_ele_dict[row[0]]
                            dangle_query = arcpy.Describe(split_flolwine_glacier).OIDFieldName+" = "+str(row[0])
                            arcpy.SelectLayerByAttribute_management(split_flolwine_glacier,"NEW_SELECTION",dangle_query)
                            start_segment = arcpy.CopyFeatures_management(split_flolwine_glacier, "in_memory\\single_dangle")
                            single_flowline = arcpy.CopyFeatures_management(split_flolwine_glacier, "in_memory\\flowline")
                            while start_ele > min_startpoint_ele:
                                arcpy.SelectLayerByLocation_management(split_flolwine_glacier, "INTERSECT", start_segment)
                                with arcpy.da.SearchCursor(split_flolwine_glacier, [arcpy.Describe(split_flolwine_glacier).OIDFieldName]) as cursor:
                                    for row in cursor:
                                        seg_ele = start_ele_dict[row[0]]
                                        if seg_ele < start_ele:
                                            seg_select_query = arcpy.Describe(split_flolwine_glacier).OIDFieldName+" ="+str(row[0])
                                            arcpy.SelectLayerByAttribute_management(split_flolwine_glacier,"NEW_SELECTION",seg_select_query)
                                            arcpy.Append_management(split_flolwine_glacier, single_flowline, "NO_TEST")
                                            start_segment = arcpy.CopyFeatures_management(split_flolwine_glacier, "in_memory\\new_start")
                                            start_ele = seg_ele
                                #del row2, cursor2
                            single_flowline_dissolve = arcpy.Dissolve_management(single_flowline, "in_memory\\single_flowline_dissolve")
                            arcpy.Append_management(single_flowline_dissolve, flow_line_final, "NO_TEST")
                    del row, cursor
                    '''
                    arcpy.Append_management(flow_line_final_individual, flow_line_final, "NO_TEST")
                else:
                    arcpy.Append_management(secondary_flowline, flow_line_final, "NO_TEST")
            else:
                new_branch_count = 0
                arcpy.Append_management(initial_flowline, flow_line_final, "NO_TEST")
            

        '''
        ##New codes
        while new_branch_count > 0:
            arcpy.AddMessage("Multiple Branches Detected...")
            secondary_axis, axispoint = create_axis(new_branch_poly, clipped_dem)
            secondary_axis_count = arcpy.GetCount_management(secondary_axis)          
            secondary_axis_result = int(secondary_axis_count.getOutput(0))
            arcpy.AddMessage(str(secondary_axis_result))
            if (secondary_axis_result > 0):
                secondary_extended_axis = extend_line(secondary_axis, 0.25)
                secondary_perpendiculars = create_perpendiculars(secondary_extended_axis, perp_space, 100000) 
                secondary_central_points_with_alt = central_points_on_perps(secondary_perpendiculars, new_branch_poly, axispoint, clipped_dem)
                secondary_flowline = flowline(secondary_central_points_with_alt, clipped_dem, "in_memory\\secondary_flow_line", new_branch_poly, 1000)
                secondary_flowline_copied = arcpy.CopyFeatures_management(secondary_flowline, "in_memory\\secondary_flow_line_copied") ##For some reason, this is necessary!!!
                arcpy.CopyFeatures_management(secondary_flowline, "d:\\temp\\secondary_flowline.shp")
                ####Now, the initial flowline needs to be fixed by rerun the flowline function
                corrected_poly, new_branch_count = new_branch(secondary_flowline_copied, clipped_dem, "in_memory\\corrected_poly", glacier_outline, TributaryRatio, sourcearea_threshold) ## use the glacier outline herer to fix the inital flowline
                arcpy.AddMessage(str(new_branch_count))
                if (new_branch_count > 0 ):
                    arcpy.AddMessage(str(new_branch_count))
                    arcpy.CopyFeatures_management(corrected_poly, "d:\\temp\\corrected_poly.shp") 
                    corrected_axis, axispoint = create_axis(corrected_poly, clipped_dem)
                    corrected_extended_axis = extend_line(corrected_axis, 0.25)
                    corrected_perpendiculars = create_perpendiculars(corrected_extended_axis, perp_space, 100000)
                    corrected_central_points_with_alt = central_points_on_perps(corrected_perpendiculars, corrected_poly, axispoint, clipped_dem)
                    corrected_flowline = flowline(corrected_central_points_with_alt, clipped_dem, "in_memory\\corrected_additional_flowline", corrected_poly, 1000)
                    corrected_flowline_copied = arcpy.CopyFeatures_management(corrected_flowline, "in_memory\\corrected_flowline_copied")
                    flowline_intersect_points = arcpy.Intersect_analysis([corrected_flowline, secondary_flowline_copied],"in_memory\\flow_intersect_points","","","POINT")
                    secondary_fl_split = arcpy.SplitLineAtPoint_management(secondary_flowline_copied, flowline_intersect_points, "in_memory\\secondary_fl_split", 5)

                    midpoint_secondary_fl_split = arcpy.FeatureVerticesToPoints_management(secondary_fl_split, "in_memory\\midpoint_secondary_fl_split", "MID")
                    midpoint_secondary_fl_split_ele = ExtractValuesToPoints(midpoint_secondary_fl_split, clipped_dem, "in_memory\\midpoint_secondary_fl_split_ele")
                    max_midpoint = 0.0
                    FcID1 = arcpy.Describe(midpoint_secondary_fl_split_ele).OIDFieldName
                    SelectID = []
                    with arcpy.da.SearchCursor(midpoint_secondary_fl_split_ele, ["RASTERVALU", FcID1]) as point_cursor:
                        for point_row in point_cursor:
                            elevation = point_row[0]
                            if elevation > max_midpoint:  ##Record the highest elevation section and its Fc ID
                                max_midpoint = elevation
                                SelectID.append(point_row[1])
                    del point_row, point_cursor
                    arcpy.AddMessage(SelectID)
                    
                    arcpy.Append_management(corrected_flowline, flow_line_cp, "NO_TEST")
                    FcID2 = arcpy.Describe(secondary_fl_split).OIDFieldName
                    selected_secondary_fl = "in_memory\\selected_secondary_fl"
                    for k in len(SelectID):
                        query = FcID2 +"="+str(SelectID[k])
                        arcpy.Select_analysis(secondary_fl_split, selected_secondary_fl,query)
                        arcpy.Append_management(selected_secondary_fl, flow_line_cp, "NO_TEST")

                    #arcpy.Append_management(selected_secondary_fl, corrected_flowline, "NO_TEST") ##Add the selected secondary to corrected flowline
                    ##update the new_branch_poly and check if more loops are needed
                    new_branch_poly, new_branch_count = new_branch(secondary_flowline_copied, clipped_dem, "in_memory\\corrected_poly_cp", new_branch_poly, TributaryRatio, sourcearea_threshold) ##Continue to test if there is other new branches
                #new_branch_poly, new_branch_count = new_branch(corrected_flowline, clipped_dem, "in_memory\\corrected_poly_cp", new_branch_poly, TributaryRatio, sourcearea_threshold) ##Continue to test if there is other new branches

                else:
                    arcpy.AddMessage("Add thirdry line!")
                    #arcpy.Append_management(initial_flowline, flow_line_cp, "NO_TEST")
                    #arcpy.Append_management(corrected_flowline, flow_line_cp, "NO_TEST")
                    arcpy.Append_management(secondary_flowline_copied, flow_line_cp, "NO_TEST")
                
            else:
                new_branch_count = 0
                arcpy.AddMessage("Add secondary line!")
                arcpy.Append_management(initial_flowline, flow_line_cp, "NO_TEST")
            '''
        progress_counter = progress_counter + 1

    ##Integrate the flowline to make sure lines are connected if there are multiple flowlines
    #if progress_counter > 2:
    #    try:
    #        arcpy.Integrate_management(flow_line_final, 10)
    #    except:
    #        pass

    ##Need to clip the flowline using the ouline polygon, so that the flowline is 100% within the glaciers
    arcpy.Clip_analysis(flow_line_final, glacier_outlines, flow_line)

    arcpy.Delete_management (flow_line_final)
        
    ##The following part is to check and flip the line; This part also line the flowline when they merged, so that all lines run from the beginning to the end
    Check_If_Flip_Line_Direction(flow_line, dem)

#------------------------------------------------------------------------------------------------------------
# Start the main program
#------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # Script arguments
    glacier_outlines = arcpy.GetParameterAsText(0)
    dem = arcpy.GetParameterAsText(1)
    TributaryRatio = float(arcpy.GetParameter(2))
    tributary_threshold = float(arcpy.GetParameter(3))
    flow_line = arcpy.GetParameterAsText(4)

    arcpy.Delete_management("in_memory")
    ##make sure the projection of the glacier outlines is the same with the UTM and the same with the DEM projection 
    spatial_ref_outlines = arcpy.Describe(glacier_outlines).spatialReference
    spatial_ref_dem = arcpy.Describe(dem).spatialReference

    if "UTM" in spatial_ref_dem.name:
        arcpy.AddMessage("The DEM projection is: " + spatial_ref_dem.name)
    else:
        arcpy.AddMessage("The DEM projection is not UTM. Please re-project the DEM to a UTM projection for the analysis!")
        exit()   

    if "UTM" in spatial_ref_outlines.name:
        arcpy.AddMessage("The outline projection is: " + spatial_ref_outlines.name)
    else:
        arcpy.AddMessage("The outline projection is not UTM. Please re-project the outlines to a UTM projection for the analysis!")
        exit()   

    Centerline_Volta (glacier_outlines, dem, TributaryRatio, tributary_threshold, flow_line)

    arcpy.Delete_management("in_memory")
   
    

 





