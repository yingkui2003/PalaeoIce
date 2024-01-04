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

def Connect_OGGM_Centerline (glacier_outlines, OGGMcenterlines, indem, search_dis, outflowlines):

    ##Step 0: convert the DEM to intergal DEM
    dem = Int(indem)
    #Step 1: Add a field to Glacier outline
    outlines = arcpy.env.scratchGDB + "\\outlines"
    sel_centerlines = arcpy.env.scratchGDB + "\\sel_centerlines"
    sel_flowlines = arcpy.env.scratchGDB + "\\sel_flowlines"
    longest_flowline = arcpy.env.scratchGDB + "\\longest_flowline"
    longest_flowline_points = arcpy.env.scratchGDB + "\\longest_flowline_points"
    otherline_start_nodes = arcpy.env.scratchGDB + "\\otherline_start_nodes"
    sel_outline = arcpy.env.scratchGDB + "\\sel_outline"
    flowlines_without_longestflowlines = arcpy.env.scratchGDB + "\\flowlines_without_longestflowlines"
    connectedflowlines = arcpy.env.scratchGDB + "\\connectedflowlines"
    dissolve_connectedflowlines = arcpy.env.scratchGDB + "\\dissolve_connectedflowlines"
    centerlines_in_outline = arcpy.env.scratchGDB + "\\centerlines_in_outline"

    arcpy.AddMessage("Make a copy of glaicier outlines...")
    print("Make a copy of glaicier outlines...")

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

    new_flowlines = arcpy.CreateFeatureclass_management(arcpy.env.scratchGDB, "new_flowlines","POLYLINE","","","",OGGMcenterlines)

    FcID = arcpy.Describe(outlines).OIDFieldName
    polygonArray = arcpy.da.FeatureClassToNumPyArray(outlines, FcID)
    GIDlist = np.array([item[0] for item in polygonArray])
    
    #uniqueGIDArr = np.unique(GIDlist)    
    arcpy.AddMessage("Check flowline direction...")
    print("Check flowline direction...")

    Check_If_Flip_Line_Direction(OGGMcenterlines, dem)
    
    for i in range(len(GIDlist)):
        #query = GlaicerID +" = "+str(GIDlist[i])
        query = FcID +" = "+str(GIDlist[i])

        arcpy.AddMessage("Generating centreline(s) for glacier "+str(GIDlist[i])+" of "+str(len(GIDlist)))
        print("Generating centreline(s) for glacier "+str(GIDlist[i])+" of "+str(len(GIDlist)))

        arcpy.Select_analysis(outlines,sel_outline,query)

        ##First use spatial join to select the centerlines within the sel_outline
        arcpy.SpatialJoin_analysis(OGGMcenterlines, sel_outline, centerlines_in_outline, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "HAVE_THEIR_CENTER_IN", None, "#")
        ##Then, clip the centerlines just within the outline       
        arcpy.Clip_analysis(centerlines_in_outline, sel_outline, sel_centerlines, "0.001 Meter")

        ##Need to make sure that there is centerlines

        LineID = arcpy.Describe(sel_centerlines).OIDFieldName
        lineArray = arcpy.da.FeatureClassToNumPyArray(sel_centerlines, 'SHAPE@LENGTH')
        lengthArr = np.array([item[0] for item in lineArray])
        if (len(lengthArr) == 0 or np.max(lengthArr) < 50):
            arcpy.AddMessage("There is no centerline for this glacier, move to the next one")
            continue ##move to the next loop
        
        ##Need to delete the small lines that cross the outline boundary
        with arcpy.da.UpdateCursor(sel_centerlines, ['SHAPE@LENGTH']) as cursor:
            for row in cursor:
                if row[0] < 50: ##if the line is less than 50 m, delete
                    arcpy.AddMessage("Delete small lines")
                    cursor.deleteRow()
        del row, cursor
         
        arcpy.cartography.SmoothLine(sel_centerlines, sel_flowlines, "PAEK", 200)

        #arcpy.Select_analysis(flowlineswithGlacerID,sel_flowlines_with_GlacerID,query)
        #arcpy.Select_analysis(outlines,sel_outline_GlacerID,query)
        
        ## find the lowest flowline
        #arcpy.FeatureVerticesToPoints_management(sel_flowlines, flowlineStartNodes, "START")
        #ExtractValuesToPoints(flowlineStartNodes, dem, flowlineStartNodes_with_ele)
        #arcpy.SpatialJoin_analysis(sel_flowlines, flowlineStartNodes_with_ele, flowlines_with_start_ele, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "0.1 Meter", "#")
        #arcpy.CopyFeatures_management(flowlines_with_start_ele, "d:\\temp\\flowlines_with_start_ele.shp")
        arcpy.AddField_management(sel_flowlines, "LineLength", "LONG", 20)
        with arcpy.da.UpdateCursor(sel_flowlines, ['SHAPE@LENGTH', 'LineLength']) as cursor:
            for row in cursor:
                row[1] = int(float(row[0]))
                cursor.updateRow(row)
        del row, cursor

        #arcpy.CopyFeatures_management(sel_flowlines, "d:\\temp\\flowlines_with_start_ele.shp")
        
        
        flowlineArr = arcpy.da.FeatureClassToNumPyArray(sel_flowlines, 'LineLength')  ###??? should just find the longest one as the lowest one in the fitst step!!!!
        #heights = np.array([item[0] for item in flowlineArr])
        lengths = np.array([item[0] for item in flowlineArr])
        #arcpy.AddMessage(lengths)
        ##If only one height value (one flowline), do nothing and return
        if len(lengths) < 2:
            arcpy.AddMessage("Only one centerline in the glacier outlines")
            print("Only one centerline in the glacier outlines")
            ##Clip the flowlines using the outline
            #arcpy.Clip_analysis(sel_flowlines_with_GlacerID, sel_outline_GlacerID, arcpy.env.scratchGDB + "\\clipped_sel_flowline", "0.001 Meter")
            arcpy.Append_management(sel_flowlines, new_flowlines, "NO_TEST")
        else:
            #minheight = min(heights)
            maxlength = max(lengths)
            #arcpy.AddMessage(str(minheight))
            ##Copy the line with minheight as a seperated flowline
            #query = 'RASTERVALU' +" = "+str(minheight)
            query = 'LineLength' +" = "+str(maxlength)
            arcpy.Select_analysis(sel_flowlines,longest_flowline,query)

            ##remove the lowest flowline
            #with arcpy.da.UpdateCursor(sel_flowlines, ['LineLength']) as cursor:
            #    for row in cursor:
            #        if row[0] == maxlength:
            #            cursor.deleteRow()
            #del row, cursor
            
            ##Do the connection with the lowest flowline
            #arcpy.FeatureVerticesToPoints_management(longest_flowline, longest_flowline_points, "ALL") ## check if the centerline can be extended to the nearest vertices
            #arcpy.CopyFeatures_management(longest_flowline_points, "d:\\temp\\longest_flowline_points.shp")

            #arcpy.FeatureVerticesToPoints_management(flowlines_with_start_ele, otherline_start_nodes, "START") ## check if the centerline can be extended to the nearest vertices

            #arcpy.FeatureVerticesToPoints_management(sel_flowlines, otherline_start_nodes, "BOTH_ENDS") ## check if the centerline can be extended to the nearest vertices
            arcpy.FeatureVerticesToPoints_management(sel_flowlines, otherline_start_nodes, "START") ## check if the centerline can be extended to the nearest vertices

            bLoop = True
            numLoop = 1
            while (bLoop):
                search_radius = str(numLoop* search_dis) + " Meters"
                #arcpy.AddMessage(search_radius)
                if numLoop* search_dis < 300: ##if larger than 500 m, try using the end points
                    arcpy.Near_analysis(otherline_start_nodes, longest_flowline, search_radius, "LOCATION", "")
                else: ##if larger than 500 m, try using the end points
                    arcpy.AddMessage("use the other end...")
                    print("use the other end...")
                    arcpy.FeatureVerticesToPoints_management(sel_flowlines, otherline_start_nodes, "END")
                    arcpy.Near_analysis(otherline_start_nodes, longest_flowline, search_radius, "LOCATION", "")

                #arcpy.CopyFeatures_management(otherline_start_nodes, "d:\\temp\\otherline_start_nodes.shp")
                    
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
                    arcpy.AddMessage("Create connection lines ... ") 
                    print("Create connection lines ... ") 
                    connectionline = arcpy.CreateFeatureclass_management(arcpy.env.scratchGDB, "connectionline","POLYLINE","","","",dem)
                    new_line_cursor = arcpy.da.InsertCursor(connectionline, ('SHAPE@'))
                    number_connected = 0
                    for i in range(len(NearID)):
                        #arcpy.AddMessage("Add line: " + str(i))
                        dist = math.sqrt((math.pow((x2[i] -x1[i]),2) + math.pow((y2[i] -y1[i]),2)))
                        #arcpy.AddMessage("Dist: " + str(dist))
                        if dist > 1: ##if distance is larger than 1 m to prevent zero connection line case
                            array = arcpy.Array([arcpy.Point(x2[i],y2[i]),arcpy.Point(x1[i], y1[i])])
                            polyline = arcpy.Polyline(array)

                            #arcpy.AddMessage("Add one connection line...")
                            number_connected += 1
                            new_line_cursor.insertRow([polyline])
                    del new_line_cursor
                    #number_connected = len(NearID)

                    arcpy.Erase_analysis(sel_flowlines, longest_flowline, flowlines_without_longestflowlines)

                    arcpy.SpatialJoin_analysis(flowlines_without_longestflowlines, connectionline, connectedflowlines, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
                    #arcpy.CopyFeatures_management(flowlines_with_start_ele, "d:\\temp\\flowlines_with_start_ele.shp")

                    arcpy.Append_management(connectionline, connectedflowlines, "NO_TEST")
                    #arcpy.CopyFeatures_management(connectedflowlines, "d:\\temp\\connectedflowlines.shp")

                    arcpy.Dissolve_management(connectedflowlines, dissolve_connectedflowlines, "","", "SINGLE_PART")
                    #arcpy.CopyFeatures_management(arcpy.env.scratchGDB + "\\dissolve_connectedflowlines", "d:\\temp\\dissolve_connectedflowlines.shp")
                    #arcpy.CopyFeatures_management(longest_flowline, "d:\\temp\\longest_flowline.shp")

                    #arcpy.AddMessage("Append")    
                    arcpy.Append_management(dissolve_connectedflowlines, longest_flowline, "NO_TEST")

                    arcpy.Delete_management(connectionline)

                if leftover > 0:
                    if number_connected < 1: ##Need to increase the search distance
                        numLoop += 1
                        #arcpy.AddMessage("numLoop is " + str(numLoop))
                        #if numLoop > 3:
                        #    bLoop = False
                    else: ##Need to update the points for lowese_flowlines
                        #arcpy.FeatureVerticesToPoints_management(longest_flowline, longest_flowline_points, "ALL") ## check if the centerline can be extended to the nearest vertices
                        pass
                        #arcpy.FeatureVerticesToPoints_management(sel_flowlines_with_GlacerID, otherline_start_nodes, "START") ## check if the centerline can be extended to the nearest vertices
                else:
                    bLoop = False ## stop the loop
            
            ##Clip the lowest flowline using the outline
            #arcpy.Clip_analysis(longest_flowline, sel_outline_GlacerID, arcpy.env.scratchGDB + "\\clipped_sel_flowline", "0.001 Meter")
            arcpy.Append_management(longest_flowline, new_flowlines, "NO_TEST")


    arcpy.CopyFeatures_management(new_flowlines, outflowlines)

    ##Delete intermidiate files
    arcpy.Delete_management(centerlines_in_outline)
    arcpy.Delete_management(outlines)
    arcpy.Delete_management(sel_flowlines)
    arcpy.Delete_management(sel_centerlines)
    arcpy.Delete_management(longest_flowline)
    arcpy.Delete_management(longest_flowline_points)
    arcpy.Delete_management(otherline_start_nodes)
    arcpy.Delete_management(sel_outline)
    arcpy.Delete_management(flowlines_without_longestflowlines)
    arcpy.Delete_management(new_flowlines)
    arcpy.Delete_management(connectedflowlines)
    arcpy.Delete_management(dissolve_connectedflowlines)

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

    OGGMcenterlines_clip = arcpy.env.scratchGDB + "\\OGGMcenterlines_clip"

    #arcpy.Clip_analysis(OGGMcenterlines, glacier_outlines, OGGMcenterlines_clip)
    arcpy.SpatialJoin_analysis(OGGMcenterlines, glacier_outlines, OGGMcenterlines_clip, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "HAVE_THEIR_CENTER_IN", None, "#")


    spatial_ref_outlines = arcpy.Describe(glacier_outlines).spatialReference
    spatial_ref_centerlines = arcpy.Describe(OGGMcenterlines_clip).spatialReference

    #OGGMcenterlines_prj = arcpy.env.scratchGDB + "\\OGGMcenterlines_prj"
    OGGMcenterlines_prj = arcpy.env.scratchGDB + "\\OGGMcenterlines_prj"
    if spatial_ref_centerlines.name == spatial_ref_outlines.name:
        arcpy.AddMessage("Both outline and centerline projection are: " + spatial_ref_centerlines.name)
        arcpy.CopyFeatures_management(OGGMcenterlines_clip, OGGMcenterlines_prj)
    else:
        arcpy.AddMessage("The centerline projection is not the same with the outline projection!")
        arcpy.AddMessage("The glacier outline projection is: " + spatial_ref_outlines.name)
        arcpy.Project_management(OGGMcenterlines_clip, OGGMcenterlines_prj, glacier_outlines)

    Connect_OGGM_Centerline (glacier_outlines, OGGMcenterlines_prj, dem, searchDis, flow_line)
        
    ##Delete intermidiate data
    arcpy.Delete_management(OGGMcenterlines_clip) ### Empty the in_memory
    arcpy.Delete_management(OGGMcenterlines_prj) ### Empty the in_memory
