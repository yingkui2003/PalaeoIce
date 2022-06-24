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
                arcpy.AddMessage("delete spurious polygons..." )
                cursor.deleteRow()
    del cursor, row

    FcID = arcpy.Describe(glacier_outlines).OIDFieldName

    arr=arcpy.da.FeatureClassToNumPyArray(glacier_outlines, FcID)
    FIds = np.array([item[0] for item in arr])

    progress_counter = 1

    glacier_outline = "in_memory\\layer_test_hope" ##Single outline for each loop
    flow_line_cp = arcpy.CreateFeatureclass_management(arcpy.env.scratchGDB, flow_line_name, "POLYLINE")  ##Maybe not necessary??, so the output is flow_line_final now

    for i in range(len(FIds)):
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

        new_branch_poly, new_branch_count  = new_branch(initial_flowline, clipped_dem, "in_memory\\new_outline", glacier_outline, TributaryRatio, sourcearea_threshold)

        if new_branch_count == 0:
            arcpy.Append_management(initial_flowline, flow_line_cp, "NO_TEST")
        while new_branch_count > 0:
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

                ####Now, the initial flowline needs to be fixed by rerun the flowline function
                corrected_poly, new_branch_count = new_branch(secondary_flowline_copied, clipped_dem, "in_memory\\corrected_poly", glacier_outline, TributaryRatio, sourcearea_threshold) ## use the glacier outline herer to fix the inital flowline
                if new_branch_count > 0:
                    corrected_axis, axispoint = create_axis(corrected_poly, clipped_dem)
                    corrected_extended_axis = extend_line(corrected_axis, 0.25)
                    corrected_perpendiculars = create_perpendiculars(corrected_extended_axis, perp_space, 100000)
                    corrected_central_points_with_alt = central_points_on_perps(corrected_perpendiculars, corrected_poly, axispoint, clipped_dem)
                    corrected_flowline = flowline(corrected_central_points_with_alt, clipped_dem, "in_memory\\corrected_additional_flowline", corrected_poly, 1000)

                    flowline_intersect_points = arcpy.Intersect_analysis([corrected_flowline, secondary_flowline_copied],"in_memory\\flow_intersect_points","","","POINT")
                    secondary_fl_split = arcpy.SplitLineAtPoint_management(secondary_flowline_copied, flowline_intersect_points, "in_memory\\secondary_fl_split", 5)

                    midpoint_secondary_fl_split = arcpy.FeatureVerticesToPoints_management(secondary_fl_split, "in_memory\\midpoint_secondary_fl_split", "MID")
                    midpoint_secondary_fl_split_ele = ExtractValuesToPoints(midpoint_secondary_fl_split, clipped_dem, "in_memory\\midpoint_secondary_fl_split_ele")
                    max_midpoint = 0.0
                    FcID1 = arcpy.Describe(midpoint_secondary_fl_split_ele).OIDFieldName
                    with arcpy.da.SearchCursor(midpoint_secondary_fl_split_ele, ["RASTERVALU", FcID1]) as point_cursor:
                        for point_row in point_cursor:
                            elevation = point_row[0]
                            if elevation > max_midpoint:  ##Record the highest elevation section and its Fc ID
                                max_midpoint = elevation
                                SelectID = point_row[1]
                    del point_row, point_cursor
                    
                    FcID2 = arcpy.Describe(secondary_fl_split).OIDFieldName
                    query = FcID2 +"="+str(SelectID)
                    selected_secondary_fl = "in_memory\\selected_secondary_fl"
                    arcpy.Select_analysis(secondary_fl_split, selected_secondary_fl,query)

                    arcpy.Append_management([selected_secondary_fl,corrected_flowline], flow_line_cp, "NO_TEST")
                    ##update the new_branch_poly and check if more loops are needed
                    new_branch_poly, new_branch_count = new_branch(secondary_flowline_copied, clipped_dem, new_branch_poly, new_branch_poly, TributaryRatio, sourcearea_threshold) ##Continue to test if there is other new branches

                else:
                    arcpy.Append_management(initial_flowline, flow_line_cp, "NO_TEST")
            
            else:
                new_branch_count = 0
                arcpy.Append_management(initial_flowline, flow_line_cp, "NO_TEST")
        
        progress_counter = progress_counter + 1

    ##Integrate the flowline to make sure lines are connected if there are multiple flowlines
    if progress_counter > 2:
        try:
            arcpy.Integrate_management(flow_line_cp, 10)
        except:
            pass

    ##Need to clip the flowline using the ouline polygon, so that the flowline is 100% within the glaciers
    arcpy.Clip_analysis(flow_line_cp, glacier_outlines, flow_line)

    arcpy.Delete_management (flow_line_cp)
        
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
   
    

 





