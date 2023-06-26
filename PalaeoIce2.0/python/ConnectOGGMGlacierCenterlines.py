# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# ConnectOGGMGlacierCenterlines.py
# Created on: 2022-11-23
# Purpose: This tool connect the centerlines derived from OGGM to generate connected flowlines for modern glaciers 
# 
# Author:      Yingkui Li
# Created:     2022-11-23
# Copyright:   (c) Yingkui Li 2022
# Department of Geography, University of Tennessee
# Knoxville, TN 37996, USA
#-------------------------------------------------------------------------------
#import SharedFunctions  
from SharedFunctions import *  

def Connect_OGGM_Centerline (glacier_outlines, OGGMcenterlines, dem, search_dis, outflowlines):

    #Step 1: Add a field to Glacier outline
    outlines = "in_memory\\outlines"
    arcpy.CopyFeatures_management(glacier_outlines, outlines)
    fieldlist = []
    ListFields=arcpy.ListFields(outlines)
    for x in ListFields:
        fieldlist.append(x.baseName)
    #arcpy.AddMessage(fieldlist)
    GlaicerID = "GlaicerID"
    if GlaicerID in fieldlist:  
        pass
    else:
        #add fieds to attribute tables to be populated with calculated values
        arcpy.AddField_management(outlines, GlaicerID, "LONG", 10)

    arcpy.CalculateField_management(outlines,GlaicerID, str("!"+str(arcpy.Describe(outlines).OIDFieldName)+"!"),"PYTHON_9.3")

    ##Step 1.1: Add the elevation of the lowest point for each flowline into flowlines
    Check_If_Flip_Line_Direction(OGGMcenterlines, dem)
    flowlineStartNodes = "in_memory\\flowlineStartNodes"
    #arcpy.FeatureVerticesToPoints_management(OGGMcenterlines, flowlineStartNodes, "START")
    ##Extract the elevation for each point
    flowlineStartNodes_with_ele = "in_memory\\flowlineStartNodes_with_ele"
    #ExtractValuesToPoints(flowlineStartNodes, dem, flowlineStartNodes_with_ele)
    flowlines_with_start_ele = "in_memory\\flowlines_with_start_ele"
    #arcpy.SpatialJoin_analysis(OGGMcenterlines, flowlineStartNodes_with_ele, flowlines_with_start_ele, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "0.1 Meter", "#")

    ##Step 2: Spatial join the glacier ID to flowlines
    flowlineswithGlacerID = "in_memory\\flowlineswithGlacerID"
    #arcpy.SpatialJoin_analysis(flowlines_with_start_ele, outlines, flowlineswithGlacerID, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "HAVE_THEIR_CENTER_IN", "#", "#")
    ##Step 3 loop for the flowline with the same GlacierID
    #lineArray = arcpy.da.FeatureClassToNumPyArray(flowlineswithGlacerID, GlaicerID)
    #GIDlist = np.array([item[0] for item in lineArray])
    #uniqueGIDArr = np.unique(GIDlist)    
    #arcpy.AddMessage(uniqueGIDArr)
    sel_flowlines = "in_memory\\sel_flowlines"
    lowest_flowline = "in_memory\\lowest_flowline"
    lowest_flowline_points = "in_memory\\lowest_flowline_points"
    otherline_start_nodes = "in_memory\\otherline_start_nodes"
    sel_outline = "in_memory\\sel_outline"

    new_flowlines = arcpy.CreateFeatureclass_management("in_memory", "new_flowlines","POLYLINE","","","",OGGMcenterlines)

    FcID = arcpy.Describe(outlines).OIDFieldName
    polygonArray = arcpy.da.FeatureClassToNumPyArray(outlines, FcID)
    GIDlist = np.array([item[0] for item in polygonArray])
    
    #uniqueGIDArr = np.unique(GIDlist)    

    for i in range(len(GIDlist)):
        #query = GlaicerID +" = "+str(GIDlist[i])
        query = FcID +" = "+str(GIDlist[i])

        arcpy.AddMessage("Generating centreline(s) for glacier "+str(GIDlist[i])+" of "+str(len(GIDlist)))
        arcpy.Select_analysis(outlines,sel_outline,query)
        arcpy.Clip_analysis(OGGMcenterlines, sel_outline, sel_flowlines, "0.001 Meter")

        ##Need to delete the small lines that cross the outline boundary
        with arcpy.da.UpdateCursor(sel_flowlines, ['SHAPE@LENGTH']) as cursor:
            for row in cursor:
                if row[0] < 50: ##if the line is less than 50 m, delete
                    cursor.deleteRow()
        del row, cursor
         

        #arcpy.Select_analysis(flowlineswithGlacerID,sel_flowlines_with_GlacerID,query)
        #arcpy.Select_analysis(outlines,sel_outline_GlacerID,query)
        
        ## find the lowest flowline
        arcpy.FeatureVerticesToPoints_management(sel_flowlines, flowlineStartNodes, "START")
        ExtractValuesToPoints(flowlineStartNodes, dem, flowlineStartNodes_with_ele)
        arcpy.SpatialJoin_analysis(sel_flowlines, flowlineStartNodes_with_ele, flowlines_with_start_ele, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "0.1 Meter", "#")
        arcpy.CopyFeatures_management(flowlines_with_start_ele, "d:\\temp\\flowlines_with_start_ele.shp")
        
        flowlineArr = arcpy.da.FeatureClassToNumPyArray(flowlines_with_start_ele, 'RASTERVALU')
        heights = np.array([item[0] for item in flowlineArr])
        ##If only one height value (one flowline), do nothing and return
        if len(heights) < 2:
            arcpy.AddMessage("Only one centerline in the glacier outlines")
            ##Clip the flowlines using the outline
            #arcpy.Clip_analysis(sel_flowlines_with_GlacerID, sel_outline_GlacerID, "in_memory\\clipped_sel_flowline", "0.001 Meter")
            arcpy.Append_management(sel_flowlines, new_flowlines, "NO_TEST")
        else:
            minheight = min(heights)
            #arcpy.AddMessage(str(minheight))
            ##Copy the line with minheight as a seperated flowline
            query = 'RASTERVALU' +" = "+str(minheight)
            arcpy.Select_analysis(flowlines_with_start_ele,lowest_flowline,query)

            ##remove the lowest flowline
            #with arcpy.da.UpdateCursor(flowlines_with_start_ele, ['RASTERVALU']) as cursor:
            #    for row in cursor:
            #        if row[0] == minheight:
            #            cursor.deleteRow()
            #del row, cursor
            
            ##Do the connection with the lowest flowline
            arcpy.FeatureVerticesToPoints_management(lowest_flowline, lowest_flowline_points, "ALL") ## check if the centerline can be extended to the nearest vertices

            #arcpy.FeatureVerticesToPoints_management(flowlines_with_start_ele, otherline_start_nodes, "START") ## check if the centerline can be extended to the nearest vertices
            arcpy.FeatureVerticesToPoints_management(flowlines_with_start_ele, otherline_start_nodes, "BOTH_ENDS") ## check if the centerline can be extended to the nearest vertices

            bLoop = True
            numLoop = 1
            while (bLoop):
                search_radius = str(numLoop* search_dis) + " Meters"
                #arcpy.AddMessage(search_radius)
                if numLoop* search_dis < 500: ##if larger than 100 m, try using the end points
                    arcpy.Near_analysis(otherline_start_nodes, lowest_flowline_points, search_radius, "LOCATION", "")
                else: ##if larger than 100 m, try using the end points
                    arcpy.FeatureVerticesToPoints_management(flowlines_with_start_ele, otherline_start_nodes, "End")
                    arcpy.Near_analysis(otherline_start_nodes, lowest_flowline_points, search_radius, "LOCATION", "")
                    
                leftover = 0
                x1 = []
                y1 = []
                x2 = []
                y2 = []
                NearID = []
                fields = ['NEAR_FID', 'SHAPE@XY', 'NEAR_X', 'NEAR_Y', 'NEAR_DIST']
                with arcpy.da.SearchCursor(otherline_start_nodes, fields)as cursor:
                    for row in cursor:
                        if row[0] >= 0:
                            if row[4] > 0: ##only connect when near distance is > 0
                                NearID.append(row[0])
                                x1.append(row[1][0])
                                y1.append(row[1][1])
                                x2.append(row[2])
                                y2.append(row[3])
                        else:
                            #arcpy.AddMessage("Cannot find the near ID!!")
                            leftover += 1
                del cursor, row
                number_connected = 0
                if len(NearID) > 0:
                    arcpy.AddMessage("Create connection lines")
                    connectionline = arcpy.CreateFeatureclass_management("in_memory", "connectionline","POLYLINE","","","",dem)
                    new_line_cursor = arcpy.da.InsertCursor(connectionline, ('SHAPE@'))
                    for i in range(len(NearID)):
                        #arcpy.AddMessage("Add line: " + str(i))
                        array = arcpy.Array([arcpy.Point(x2[i],y2[i]),arcpy.Point(x1[i], y1[i])])
                        polyline = arcpy.Polyline(array)
                        new_line_cursor.insertRow([polyline])
                    del new_line_cursor
                    number_connected = len(NearID)

                    arcpy.Erase_analysis(flowlines_with_start_ele, lowest_flowline, "in_memory\\flowlines_without_lowestflowlines")

                    arcpy.SpatialJoin_analysis("in_memory\\flowlines_without_lowestflowlines", connectionline, "in_memory\\connectedflowlines", "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
                    #arcpy.CopyFeatures_management(flowlines_with_start_ele, "d:\\temp\\flowlines_with_start_ele.shp")

                    arcpy.Append_management(connectionline, "in_memory\\connectedflowlines", "NO_TEST")

                    arcpy.Dissolve_management("in_memory\\connectedflowlines", "in_memory\\dissolve_connectedflowlines", "","", "SINGLE_PART")
                    #arcpy.CopyFeatures_management("in_memory\\dissolve_connectedflowlines", "d:\\temp\\dissolve_connectedflowlines.shp")
                    #arcpy.CopyFeatures_management(lowest_flowline, "d:\\temp\\lowest_flowline.shp")

                    #arcpy.AddMessage("Append")    
                    arcpy.Append_management("in_memory\\dissolve_connectedflowlines", lowest_flowline, "NO_TEST")

                if leftover > 0:
                    if number_connected < 1: ##Need to increase the search distance
                        numLoop += 1
                        #arcpy.AddMessage("numLoop is " + str(numLoop))
                    else: ##Need to update the points for lowese_flowlines
                        arcpy.FeatureVerticesToPoints_management(lowest_flowline, lowest_flowline_points, "ALL") ## check if the centerline can be extended to the nearest vertices
                        #arcpy.FeatureVerticesToPoints_management(sel_flowlines_with_GlacerID, otherline_start_nodes, "START") ## check if the centerline can be extended to the nearest vertices
                else:
                    bLoop = False ## stop the loop
            
            ##Clip the lowest flowline using the outline
            #arcpy.Clip_analysis(lowest_flowline, sel_outline_GlacerID, "in_memory\\clipped_sel_flowline", "0.001 Meter")
            arcpy.Append_management(lowest_flowline, new_flowlines, "NO_TEST")

    #arcpy.cartography.SmoothLine(new_flowlines, outflowlines, "PAEK", 200)

    arcpy.CopyFeatures_management(new_flowlines, outflowlines)
    return


#------------------------------------------------------------------------------------------------------------
# Start the main program
#------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # Script arguments
    glacier_outlines = arcpy.GetParameterAsText(0)
    OGGMcenterlines = arcpy.GetParameterAsText(1)
    dem = arcpy.GetParameterAsText(2)
    searchDis = float(arcpy.GetParameter(3))
    flow_line = arcpy.GetParameterAsText(4)

    arcpy.Delete_management("in_memory")

    OGGMcenterlines_clip = "in_memory\\OGGMcenterlines_clip"
    #arcpy.CopyFeatures_management(OGGMcenterlines, OGGMcenterlines_cp)
    outflowlines = "in_memory\\outflowlines"

    ##Use the convex hull of glacier outlines to clip OGGMcenterlines
    #arcpy.MinimumBoundingGeometry_management(glacier_outlines, "in_memory\\mbg", "CONVEX_HULL", "NONE", "","MBG_FIELDS")
    arcpy.Clip_analysis(OGGMcenterlines, glacier_outlines, OGGMcenterlines_clip)


    spatial_ref_outlines = arcpy.Describe(glacier_outlines).spatialReference
    spatial_ref_centerlines = arcpy.Describe(OGGMcenterlines_clip).spatialReference

    #OGGMcenterlines_prj = "in_memory\\OGGMcenterlines_prj"
    OGGMcenterlines_prj = arcpy.env.scratchGDB + "\\OGGMcenterlines_prj"
    if spatial_ref_centerlines.name == spatial_ref_outlines.name:
        arcpy.AddMessage("Both outline and centerline projection are: " + spatial_ref_centerlines.name)
        arcpy.CopyFeatures_management(OGGMcenterlines_clip, OGGMcenterlines_prj)
    else:
        arcpy.AddMessage("The centerline projection is not the same with the outline projection!")
        arcpy.AddMessage("The glacier outline projection is: " + spatial_ref_outlines.name)
        arcpy.Project_management(OGGMcenterlines_clip, OGGMcenterlines_prj, glacier_outlines)

    ##Remove the small lines generated by clipping
    ##Need to delete the small lines that cross the outline boundary
    #arcpy.CopyFeatures_management(OGGMcenterlines_prj, "c:\\test3\\OGGMcenterlines_prj.shp")
    #with arcpy.da.UpdateCursor(OGGMcenterlines_prj, ['SHAPE@LENGTH']) as cursor:
    #    for row in cursor:
    #        if row[0] < 50: ##if the line is less than 50 m, delete
    #            arcpy.AddMessage("Delete small lines")
    #            cursor.deleteRow()
    #del row, cursor

    
    
    Connect_OGGM_Centerline (glacier_outlines, OGGMcenterlines_prj, dem, searchDis, flow_line)
    '''    
    ##make sure the projection of the glacier outlines is the same with the UTM and the same with the DEM projection 
    spatial_ref_outlines = arcpy.Describe(glacier_outlines).spatialReference
    spatial_ref_dem = arcpy.Describe(dem).spatialReference

    ##Stream network analysis
    cellsize = arcpy.GetRasterProperties_management(dem,"CELLSIZEX")
    cellsize_int = int(cellsize.getOutput(0))


    FcID = arcpy.Describe(glacier_outlines).OIDFieldName

    arr=arcpy.da.FeatureClassToNumPyArray(glacier_outlines, FcID)
    FIds = np.array([item[0] for item in arr])

    progress_counter = 1

    glacier_outline = "in_memory\\glacier_outline" ##Single outline for each loop
    #flow_line_name = "flowlinecp"
    flow_line_cp = "in_memory\\flow_line_cp"
    #flow_line_cp = arcpy.CreateFeatureclass_management(arcpy.env.scratchGDB, flow_line_name, "POLYLINE")  ##Maybe not necessary??, so the output is flow_line_final now

    for i in range(len(FIds)):
        arcpy.AddMessage("Generating centreline(s) for glacier "+str(progress_counter)+" of "+str(len(FIds)))
        query = FcID +" = "+str(FIds[i])

        arcpy.Select_analysis(glacier_outlines,glacier_outline,query)

        ##extract by mask the DEM so that the analysis is only for the extracted DEM
        #arcpy.Buffer_analysis(glacier_outline, "in_memory\\glacier_outline_buf", "100 Meter")
        #clipped_dem = ExtractByMask(dem,"in_memory\\glacier_outline_buf")
        clipped_dem = ExtractByMask(dem,glacier_outline)

        centerline = Centerline_for_single_outline (glacier_outline, clipped_dem, StreamThreshold, TributaryThreshold)

        if progress_counter == 1:
            arcpy.CopyFeatures_management(centerline, flow_line_cp)
        else:
            arcpy.Append_management(centerline, flow_line_cp, "NO_TEST")

        progress_counter += 1

    ##Integrate the flowline to make sure lines are connected if there are multiple flowlines
    if progress_counter > 2:
        try:
            arcpy.Integrate_management(flow_line_cp, 10)
        except:
            pass
        

    Check_If_Flip_Line_Direction (flow_line_cp, dem)

    #arcpy.CopyFeatures_management(flow_line_cp, flow_line)
    
    ##Smooth the line
    lineSmooth(flow_line_cp, "in_memory\\final_flow_line", "Max_Max", cellsize_int*10)

    arcpy.AddMessage ("Merge and Add Glacier ID...")

    #arcpy.CopyFeatures_management("in_memory\\final_flow_line", flow_line)    
    Merge_and_Add_GlacierID_by_Topology ("in_memory\\final_flow_line", "Max_Max", "GlacierID", "MergeID", flow_line)
    '''    
        
    ##Delete intermidiate data
    arcpy.Delete_management("in_memory") ### Empty the in_memory
