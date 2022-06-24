#-------------------------------------------------------------------------------
# Name: Combine Flowlines with Centerlines
# Purpose: This tool combines the flowlines generated from stream network with centerlines derived from modern glacier outlines.
# The flowlines derived from modern glacier surface may be not correct becasue the topography in the middle of glaciers may be higher
# than the margin of glaciers espacilly in the ablation zone. This tool replaces the flowline within the glacier outlines with the
# centerlines derived from glacier outlines (based on Volta, other programs, or manually delineated centerlines). 
# 
# Author:      Yingkui Li
# Created:     08/30/2020 to 12/16/2020
# Copyright:   (c) Yingkui Li 2020
# Department of Geography, University of Tennessee
# Knoxville, TN 37996, USA
#-------------------------------------------------------------------------------

#import SharedFunctions  
from SharedFunctions import *  ## 

#------------------------------------------------------------------------------------------------------------
# This fuction is the main program to combine flowlines and centerlines.
#------------------------------------------------------------------------------------------------------------
def Combine_Flowlines_with_Centerlines (flowlineinput, centerlineinput, outlinepolysinput, inputWS, dem, SearchDis, combineflowline):

    centerlineinWatershed = "in_memory\\centerlineinWatershed"
    centerlineToNodes = "in_memory\\centerlineToNodes"
    flowlineToNodes = "in_memory\\flowlineToNodes"
    dissolvedflowline = "in_memory\\dissolvedflowline"
    centerline = "in_memory\\centerline"
    flowlinecopy = "in_memory\\flowlinecopy"
    flowline = "in_memory\\flowline"
    outlinepolys = "in_memory\\outlinepolys"

    ####Flow direction and accumulation analysis
    cellsize = arcpy.GetRasterProperties_management(dem,"CELLSIZEX")
    cellsize_float = float(cellsize.getOutput(0)) # use float cell size

    arcpy.CopyFeatures_management(centerlineinput, centerline)
    arcpy.CopyFeatures_management(flowlineinput, flowline)
    arcpy.CopyFeatures_management(outlinepolysinput, outlinepolys)

    ##Only consider the centerline within the watershed
    arcpy.Clip_analysis(centerline, inputWS, centerlineinWatershed)
    centerlinecountResult= arcpy.GetCount_management(centerlineinWatershed)
    centerlinecount = int(centerlinecountResult.getOutput(0))

    if centerlinecount == 0:
        arcpy.AddMessage ("no glacier centerline in the watershed of the flowlines")
        arcpy.Delete_management("in_memory") ### Empty the in_memory
        exit()

    ##Make sure centerline is from low to high elevations
    Check_If_Flip_Line_Direction(centerline, dem)
    Check_If_Flip_Line_Direction(flowline, dem)

    points=[]
    lines=[]
    #Add Max_Max field to centerline and assign the value with shape@length
    arcpy.AddField_management(centerline, "Max_Max", "LONG")
    with arcpy.da.UpdateCursor(centerline, ["SHAPE@", "SHAPE@LENGTH", "Max_Max"]) as flows:
        for flow in flows:
            flow[2] = int(flow[1])
            points.append(((flow[0].firstPoint).X, (flow[0].firstPoint).Y))
            lines.append(flow[0])
            flows.updateRow(flow)
    del flow, flows

    ##Just get the main flow start points (only one start points) for the flowline
    centerlineToNodes = arcpy.CreateFeatureclass_management("in_memory", "centerlineToNodes","POINT","","","",centerline)
    new_point_cursor = arcpy.da.InsertCursor(centerlineToNodes, ('SHAPE@'))


    array=np.array(points)

    ##Get the lowest flowline
    for i in range(len(array)):
        nTouches=0
        punto=arcpy.Point(array[i][0],array[i][1])
        for a in range (len(array)):
            if a != i: ##don't touch itself
                touches=(punto.touches(lines[a]) or punto.within(lines[a]))
                if touches==True: ##if the start point touches others
                    nTouches = 1
                    break  ##break the loop
        if nTouches==0: ##this is for lowest points
            new_point_cursor.insertRow([punto])
    del new_point_cursor


    ##Add Merge ID to the flowline, so that the splited flowline can be merged back later
    arcpy.AddField_management(flowline, "MergeID", "Long", "", "", "", "", "", "", "")
    arcpy.CalculateField_management(flowline,"MergeID",str("!"+str(arcpy.Describe(flowline).OIDFieldName)+"!"),"PYTHON_9.3")

    
    split_flowline = line_split_by_intersect(flowline)

    ##added by Yingkui Li on 6/2/2021
    ##only consider the flowline sections that is intersected with the glacier outline
    ##Save the other flowlines that are not intersected with the outlinepolys to a seperated file
    flowline_intersected = "in_memory\\flowline_intersected"
    flowline_leftover = "in_memory\\flowline_leftover"
    
    outlinepoly_layer = arcpy.MakeFeatureLayer_management(outlinepolys, "in_memory\\outlinepoly_layer")
    flowline_layer = arcpy.MakeFeatureLayer_management(split_flowline, "in_memory\\stream_layer")


    arcpy.SelectLayerByLocation_management(flowline_layer,"INTERSECT", outlinepoly_layer,"1 METERS","NEW_SELECTION","")
    arcpy.CopyFeatures_management(flowline_layer, flowline_intersected)

    arcpy.SelectLayerByAttribute_management (flowline_layer, "SWITCH_SELECTION")
    arcpy.CopyFeatures_management(flowline_layer, flowline_leftover)
    
    
    ###remove the flowline within the glacier outlines by select by erase
    feature_count_result = arcpy.GetCount_management(flowline_intersected)          
    feature_count = int(feature_count_result.getOutput(0))
    arcpy.AddMessage("Intersected flowlines: " + str(feature_count))
    if feature_count < 1: ##No feature selected, do nothing, simply return the flowline input                                       
        arcpy.CopyFeatures_management(flowlineinput, combineflowline)
        return

    #flowlinecopy, bEmpty = Remove_Flowline_in_Polygon_by_Erase(flowline_intersected, outlinepolys, 100)
    flowlinecopy, bEmpty = Remove_Flowline_in_Polygon_by_Erase(flowline_intersected, outlinepolys, cellsize_float)
    
    if bEmpty == True:
        ##copy the centerline as the combined flowline
        arcpy.AddMessage("Copy the centerline as flowlines")
        arcpy.CopyFeatures_management(centerline, combineflowline)
        return
    
    ##Find the end points for selected flowline
    #arcpy.FeatureVerticesToPoints_management(flowlinecopy, flowlineToNodes, "END")
    arcpy.FeatureVerticesToPoints_management(flowlinecopy, flowlineToNodes, "ALL") ## check if the centerline can be extended to the nearest vertices

    ###Start to process the connectionline
    search_radius = str(SearchDis) + " Meters"
    arcpy.Near_analysis(centerlineToNodes, flowlineToNodes, search_radius, "LOCATION", "")


    nodecountresult = arcpy.GetCount_management(centerlineToNodes)
    nodecount = int(nodecountresult.getOutput(0))
    arcpy.AddMessage("The number of near IDs is: " + str(nodecount))
    if nodecount == 0: ##If no nearID
        arcpy.AddMessage("No nearID is selected, exit the program!!")
        exit()
    else:
        #if nodecount > 0:
        ##Get the points from the centerlinetonodes and flowlinetonodes
        x1 = []
        y1 = []
        x2 = []
        y2 = []
        NearID = []
        fields = ['NEAR_FID', 'SHAPE@XY', 'NEAR_X', 'NEAR_Y']
        with arcpy.da.SearchCursor(centerlineToNodes, fields)as cursor:
            for row in cursor:
                if row[0] >= 0:
                    NearID.append(row[0])
                    x1.append(row[1][0])
                    y1.append(row[1][1])
                    x2.append(row[2])
                    y2.append(row[3])
        del cursor, row
        ##Create the connectionline
        arcpy.AddMessage("Near ID number is: " + str(len(NearID)))
        if len(NearID) > 0:
            arcpy.AddMessage("Create connection lines")
            connectionline = arcpy.CreateFeatureclass_management("in_memory", "connectionline","POLYLINE","","","",dem)
            new_line_cursor = arcpy.da.InsertCursor(connectionline, ('SHAPE@'))
            for i in range(len(NearID)):
                arcpy.AddMessage("Create connection lines: " + str(i+1))
                array = arcpy.Array([arcpy.Point(x2[i],y2[i]),arcpy.Point(x1[i], y1[i])])
                polyline = arcpy.Polyline(array)
                new_line_cursor.insertRow([polyline])
            del new_line_cursor

            ##Need to add the mergeID to connectionline by spatial join with the conncectionline
            connectionlinecopy = "in_memory\\connectionlinecopy"      
            arcpy.SpatialJoin_analysis(connectionline, flowlinecopy, connectionlinecopy, "JOIN_ONE_TO_ONE", "KEEP_ALL", '#', "INTERSECT", "1 Meters", "#")

            arcpy.Append_management(connectionlinecopy, flowlinecopy, "NO_TEST")

    ##Need to add the mergeID to centerline by spatial join with the conncectionline
    centerlinecopy = "in_memory\\centerlinecopy"      
    arcpy.SpatialJoin_analysis(centerline, connectionlinecopy, centerlinecopy, "JOIN_ONE_TO_ONE", "KEEP_ALL", '#', "INTERSECT", "1 Meters", "#")
    
    arcpy.Append_management(centerlinecopy, flowlinecopy, "NO_TEST")

    ##only select the flowlinecopy that is intersected/or touched with the connection line
    connected_flowline = "in_memory\\connected_flowline"
    arcpy.SpatialJoin_analysis(flowlinecopy, connectionlinecopy, connected_flowline, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
    
    ##Append the connected_flowline back to the flowline_leftover
    arcpy.Append_management(connected_flowline, flowline_leftover, "NO_TEST")
    
    arcpy.Dissolve_management(flowline_leftover, dissolvedflowline, "MergeID", "#", 'SINGLE_PART', 'UNSPLIT_LINES')

    #Handle the issue of the original dissolve with the same connection point issue
    dissolvedlines = MergeUndissolvedFlowlines(dissolvedflowline, "MergeID")

    try:
        arcpy.Integrate_management(dissolvedlines, cellsize_float)
    except:
        pass
    dissolvedlines2 = "in_memory\\dissolvedlines2"
    ##Dissolve again to make sure the line connected
    arcpy.Dissolve_management(dissolvedlines, dissolvedlines2, "MergeID", "#", 'SINGLE_PART', 'UNSPLIT_LINES')

    ###remove the flowline within the glacier outlines by select by attribute
    exist_fields = [f.name for f in arcpy.ListFields(dissolvedlines2)] #List of current field names in outline layer
    #arcpy.AddMessage(exist_fields)
    field = ["GlacierID"]
    if field in exist_fields:
        arcpy.DeleteField_management(dissolvedlines2,field) #Make sure to delete GlacierID. It will be added later in the paleoice reconstruction

    ##Remove big turns along the flowlines before smoothing it
    newline = flowline_remove_bigturn(dissolvedlines2, 120, cellsize_float)
    arcpy.cartography.SmoothLine(newline, "in_memory\\smmothline", "PAEK", 200)

    
    arcpy.CopyFeatures_management("in_memory\\smmothline", combineflowline)
    ##Do the intergrate again in case the smmothline does not break the connection of tributaries and the main flowline
    try:
        arcpy.Integrate_management(combineflowline, cellsize_float)
    except:
        pass
    
##Main program
if __name__ == '__main__':
    # Script arguments
    flowlineinput = arcpy.GetParameterAsText(0)          ##Input flowline
    centerlineinput = arcpy.GetParameterAsText(1)        ##Input Centerline
    outlinepolys = arcpy.GetParameterAsText(2)      ##Input Glacier bnd
    inputWS = arcpy.GetParameterAsText(3)           ## Input watershed
    dem = arcpy.GetParameterAsText(4)               ##Input DEM
    SearchDis = arcpy.GetParameter(5)                                 ## Search distance for line match
    combineflowline = arcpy.GetParameterAsText(6)   ##Output combined flowline


    arcpy.Delete_management("in_memory") ### Empty the in_memory

    Combine_Flowlines_with_Centerlines (flowlineinput, centerlineinput, outlinepolys, inputWS, dem, SearchDis, combineflowline)

    arcpy.Delete_management("in_memory") ### Empty the in_memory

