# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# GlacierCenterlineNew.py
# Created on: 2022-11-13 15:40:52.00000
# Purpose: This tool derives centerlines for modern glaciers based on ????? 
# 
# Author:      Yingkui Li
# Created:     2022-11-13 15:40:52.00000
# Copyright:   (c) Yingkui Li 2022
# Department of Geography, University of Tennessee
# Knoxville, TN 37996, USA
#-------------------------------------------------------------------------------
#import SharedFunctions  
from SharedFunctions import *  

#---------------------------------------------------------------------------------------------------------------
# This function calculates the distance between two points
#--------------------------------------------------------------------------------------------------------------- 
def Dist(x1,y1,x2,y2):
    return math.sqrt(math.pow(math.fabs(x1-x2),2)+math.pow(math.fabs(y1-y2),2))

def delete_multiple_element(list_object, indices):
    indices = sorted(indices, reverse=True)
    for idx in indices:
        if idx < len(list_object):
            list_object.pop(idx)

            
##This is for bigturn in 2D
def get_lowest_and_localmax_points(line, dem, max_angle, min_dis):
    arcpy.FeatureVerticesToPoints_management(line, "in_memory\\line_points", 'All')
    ExtractValuesToPoints("in_memory\\line_points", dem, "in_memory\\line_points_with_alt", "INTERPOLATE", "VALUE_ONLY")
    ###Create the new line after removing the outlier points
    spatialref=arcpy.Describe(line).spatialReference
    field = 'ORIG_FID' ##The first two fields are FID and Geometry

    pointarray = arcpy.da.FeatureClassToNumPyArray("in_memory\\line_points_with_alt", ('SHAPE@X', 'SHAPE@Y','RASTERVALU',field))
    line_ids = np.array([item[3] for item in pointarray])
    unique_line_ids = np.unique(line_ids)

    t_points = []
    t_angles   = []
    t_FID    = []
    t_dis     = []
    t_ele     = []
    
    for fid in unique_line_ids:
        arr = pointarray[line_ids == fid]
        pntx = []
        pnty = []
        pntz = []
        cumLen = []
        dis = 0
        cumLen.append(dis) ##set the first distance to zero
        for i in range(len(arr)):
            pntx.append(arr[i][0])
            pnty.append(arr[i][1])
            pntz.append(arr[i][2])
            if i>0:
                point_dis = Dist(arr[i][0], arr[i][1], arr[i-1][0], arr[i-1][1])
                dis += point_dis
                cumLen.append(dis) ##set the first distance to zero
                
            
        ##Add the second point at the end in order to calculate the value for the end or start point 
        pntx.append(arr[1][0])
        pnty.append(arr[1][1])
        pntz.append(arr[1][2])

        ##NEED to first filter the points and only keep the local min and local maximum points
        #arcpy.AddMessage(str(len(pntx)))
        removelist = []
        for i in range(len(pntx)):
            if i <(len(pntx)-1) and i > 0:#if it is not first or last point of all
                #arcpy.AddMessage(str(i))
                #arcpy.AddMessage(str(pntz[i-1]))
                #arcpy.AddMessage(str(pntz[i]))
                #arcpy.AddMessage(str(pntz[i+1]))
                
                if (((pntz[i] < max(pntz[i-1], pntz[i+1])) and ((pntz[i] > min(pntz[i-1], pntz[i+1]))))):
                    #arcpy.AddMessage("Remove one!")
                    removelist.append(i)
                    #if i == (len(pntx)-2): ##if deleting the last points of the poylgon, also delete the first points because they are the same
                    #    removelist.append(0)
                        
        #arcpy.AddMessage(str(len(pntx))) 
        if len(removelist) > 0:
            delete_multiple_element(pntx, removelist)
            delete_multiple_element(pnty, removelist)
            delete_multiple_element(pntz, removelist)
            delete_multiple_element(cumLen, removelist)
        
        points = np.array(list(zip(pntx,pnty)))

        #arcpy.AddMessage(str(len(removelist))) 
        #arcpy.AddMessage(str(len(points))) 
        
        ##derive the updated angles
        end = points[2:]
        start = points[0:-2]
        center = points[1:-1]
        ba = start-center
        bc = end-center
        ca = end - start
        dot = []
        for i in range(len(bc)):
            dot.append(np.dot(ba[i], bc[i]))

        cosine_angle = np.array(dot) / (np.linalg.norm(ba, axis=-1) * np.linalg.norm(bc, axis=-1))
        #angles = np.degrees(np.arccos(cosine_angle))
        angles = np.degrees(np.arccos(np.maximum(np.minimum(1,cosine_angle), -1)))
        #Derive the distance to the line, so that to determine the convex slope
        dist = np.cross(ca, -ba, axis=-1) / np.linalg.norm(ca, axis=-1)

        tmp_points = []
        tmp_angles   = []
        tmp_FID    = []
        tmp_dis     = []
        tmp_ele     = []
        tmp_cumdis  = []
        ##Need to keep the minimum and maximum elevation points
        minZ = min(pntz)
        maxZ = max(pntz)
        #arcpy.AddMessage("MinZ is " + str(minZ))
        #arcpy.AddMessage("MaxZ is " + str(maxZ))
        #arcpy.AddMessage(str(len(angles)))
        #arcpy.AddMessage(str(len(cumLen)))
        #arcpy.AddMessage(cumLen)
        
        for i in range(len(points)):
            if i <(len(points)-1) and i > 0:#if it is not first or last point of all
                pntangle = angles[i-1] ##get the angle
                pntdist = dist[i-1]    ##get the direction of the angle
                #arcpy.AddMessage(str(i))
                #arcpy.AddMessage(str(pntz[i-1]))
                #arcpy.AddMessage(str(pntz[i]))
                #arcpy.AddMessage(str(pntz[i+1]))

                if (abs(pntangle) < max_angle and pntz[i] > max(pntz[i-1],pntz[i+1])) or (int(pntz[i]) == int(minZ)):# or (int(pntz[i]) == int(maxZ)) : ##The max point can remove
                    #arcpy.AddMessage("Add one!")
                    tmp_points.append(points[i])
                    tmp_angles.append(pntangle)
                    tmp_dis.append(pntdist)
                    tmp_ele.append(pntz[i])
                    tmp_FID.append(fid)
                    tmp_cumdis.append(cumLen[i])
        ##Check the turning points using the cumLen
        #arcpy.AddMessage(tmp_cumdis)         
        if len(tmp_points) > 1:
            ids = []
            for i in range(1, len(tmp_points)):
                #point_dis = Dist(t_points[i][0], t_points[i][1], t_points[i-1][0], t_points[i-1][1])
                point_dis = tmp_cumdis[i] - tmp_cumdis[i-1]
                #arcpy.AddMessage(str(point_dis) )         
                if (point_dis < min_dis): ##if the two point distance are two close, only pick the largest elevation one
                    if (tmp_ele[i] < tmp_ele[i-1]):
                        if tmp_ele[i] > minZ: ## then remove i 
                            index = i
                        elif tmp_ele[i-1] < maxZ:
                            index = i-1
                    else:
                        if tmp_ele[i-1] > minZ: ## then remove i 
                            index = i-1
                        elif tmp_ele[i] < maxZ:
                            index = i
                    ids.append(index)
            #remove the value at index
            #arcpy.AddMessage(ids) 
            
            if len(ids) > 0:
                delete_multiple_element(tmp_points, ids)
                delete_multiple_element(tmp_angles, ids)
                delete_multiple_element(tmp_dis, ids)
                delete_multiple_element(tmp_ele, ids)
                delete_multiple_element(tmp_FID, ids)
                 
        ##Add the cleaned turning points to         
        t_points.extend(tmp_points)
        t_angles.extend(tmp_angles)
        t_dis.extend(tmp_dis)
        t_ele.extend(tmp_ele)
        t_FID.extend(tmp_FID)
                    
    ##Create the point file
    t_points_FC = arcpy.CreateFeatureclass_management("in_memory", "t_points_FC","POINT","", "DISABLED", "DISABLED", spatialref)
    arcpy.AddField_management(t_points_FC, "LineID", "Long")
    arcpy.AddField_management(t_points_FC, "Angle", "DOUBLE")
    arcpy.AddField_management(t_points_FC, "Dist", "DOUBLE")
    arcpy.AddField_management(t_points_FC, "Elevation", "DOUBLE")
    #arcpy.AddField_management(t_points_FC, "CumLen", "DOUBLE")
        
    CursorPnts = arcpy.da.InsertCursor(t_points_FC, ['SHAPE@', 'LineID', 'Angle', 'Dist', 'Elevation'])
    #LineIDlist.append(lineId)
    for i in range(len(t_points)):
        Pnt = arcpy.Point(t_points[i][0], t_points[i][1])
        CursorPnts.insertRow((Pnt, t_FID[i], t_angles[i], t_dis[i], t_ele[i])) 
    del CursorPnts ##Delete cursor object

    ##Get the lowest point and other points
    query = "Elevation < " + str(min(t_ele) + 1)

    lowestpnt = "in_memory\\lowestpnt"
    arcpy.Select_analysis (t_points_FC, lowestpnt, query)

    ##erase the lowest point from t_points_FC
    local_max_pnts = "in_memory\\local_max_pnts"
    arcpy.Erase_analysis(t_points_FC, lowestpnt, local_max_pnts)
    
    return lowestpnt, local_max_pnts


def clean_untouched_lines(inlines, inpoints):

    linearray = arcpy.da.FeatureClassToNumPyArray(inlines, ('OID@'))
    lineIDs = np.array([item[0] for item in linearray])
    if (len(lineIDs)) == 0:
        #arcpy.AddMessage("No lines!!")
        return inlines
    
    fromnodes=[]
    
    with arcpy.da.SearchCursor(inlines, ["SHAPE@"]) as flows: 
        for flow in flows:
            fromnodes.append(((flow[0].firstPoint).X, (flow[0].firstPoint).Y))
    del flow, flows

    inpointX = []
    inpointY = []
    with arcpy.da.SearchCursor(inpoints, ["SHAPE@XY"]) as cursor: 
        for row in cursor:
            inpointX.append(row[0][0])
            inpointY.append(row[0][1])
    del row, cursor


    array=np.array(fromnodes)
    keep_flag = []
    from_touch = []
    for i in range(len(array)):
        Touch=0
        for a in range(len(inpointX)):
            point_dis = Dist(array[i][0], array[i][1], inpointX[a], inpointY[a])
            #arcpy.AddMessage(str(point_dis))
            if point_dis < 100: ##the distance of the start line is < 100 from the inpoints
                Touch=1
                break
        keep_flag.append(Touch)

    #arcpy.AddMessage(keep_flag)
    with arcpy.da.UpdateCursor(inlines, ["OID@"]) as cursor:
        i = 0
        for row in cursor:
            if keep_flag[i] == 0:
                #arcpy.AddMessage("Delete one line away from the lowest point!")
                cursor.deleteRow()
            i += 1
    del row, cursor

    return inlines

#------------------------------------------------------------------------------------------------------------
# This function check each line in the line feature and make sure the line is from low elevation to high
# elevation (Glacier flowline needs from low to high elevation in order to reconstruct the paleo ice thickness).
# It is revised from the codes by Pellitero et al.(2016) in GlaRe.
#------------------------------------------------------------------------------------------------------------
def Check_Sink_And_Flip(inline, dem):

    fromnodes=[]
    tonodes=[]
    lines=[]
    start_ele = []
    end_ele = []
    ##First loop to get the startpoints, lines and the glacierID list
    with arcpy.da.SearchCursor(inline, ["SHAPE@"]) as flows: 
        for flow in flows:
            fromnodes.append(((flow[0].firstPoint).X, (flow[0].firstPoint).Y))
            coord= str((flow[0].firstPoint).X)+" "+str((flow[0].firstPoint).Y)
            Cellvalue = arcpy.GetCellValue_management(dem, coord)
            start_ele.append(int(Cellvalue.getOutput(0)))
            tonodes.append(((flow[0].lastPoint).X, (flow[0].lastPoint).Y))
            coord= str((flow[0].lastPoint).X)+" "+str((flow[0].lastPoint).Y)
            Cellvalue = arcpy.GetCellValue_management(dem, coord)
            end_ele.append(int(Cellvalue.getOutput(0)))
            lines.append(flow[0])
    del flow, flows


    ##Check if the elevation of the centerline network without sinks
    array=np.array(fromnodes)
    flip = [0]*len(array)
    loop = 0
    while loop < 5:
        for i in range(len(array)):
            near_ele = []
            punto=arcpy.Point(array[i][0],array[i][1])
            for a in range (len(array)):
                touches=(punto.touches(lines[a]) or punto.within(lines[a]))
                if touches==True: ##if the start point touches others
                    #arcpy.AddMessage("touched")
                    near_ele.append(end_ele[a])
                    #near_ele.append(start_ele[a])
            #arcpy.AddMessage(near_ele)
            if (start_ele[i] < min(near_ele)) and (len(near_ele)> 1):
                new_ele = []
                for k in range(len(near_ele)):
                    if near_ele[k] > start_ele[i]:
                        new_ele.append(near_ele[k])
                #arcpy.AddMessage(new_ele)
                if end_ele[i] == min(new_ele):
                    arcpy.AddMessage("The line ID needs to flip: " + str(i))
                    flip[i] = 1
                    start_ele[i] = min(new_ele) + 1 ##increase the start_ele value to the minimum surrounding elevation + 1
            loop += 1
            
    #arcpy.AddMessage(flip)            
    if max(flip) > 0:
        arcpy.AddField_management(inline, "Flip", "Long", "", "", "", "", "", "", "")
        with arcpy.da.UpdateCursor(inline, "Flip") as cursor:
            i = 0
            for row in cursor:
                row[0] = flip[i]
                cursor.updateRow(row)
                i += 1
        del cursor, row
        
        arcpy.MakeFeatureLayer_management(inline, "lyrLines")
        arcpy.SelectLayerByAttribute_management("lyrLines", "NEW_SELECTION", '"Flip" > 0')
        #arcpy.AddMessage("Count of selected features is " + str(nFlip))
        arcpy.FlipLine_edit("lyrLines")  ##Need to change to lyrLines
        arcpy.SelectLayerByAttribute_management("lyrLines", "CLEAR_SELECTION")
    arcpy.DeleteField_management (inline, "Flip")


def upstreamLength(sel_tonode, fromnodes, tonodes, lines, lengths, flags, section_Length):
    cumLength = section_Length
    new_tonodes = []
    new_tonodes.append(sel_tonode)

    while (len(new_tonodes) > 0):
        nTouches=0
        tonode = new_tonodes[0]
        punto=arcpy.Point(tonode[0],tonode[1])
        new_tonodes.pop(0) ##remove this tonode
        for a in range (len(flags)):
            if flags[a] == 0: ##the line has not been used
                touches=(punto.touches(lines[a]) or punto.within(lines[a]))
                if touches==True:
                    #arcpy.AddMessage("Touch!")
                    nTouches += 1
                    flags[a] = 1
                    cumLength += lengths[a]
                    new_tonodes.append(tonodes[a])
                        
        #arcpy.AddMessage(str(cumLength))     
    #arcpy.AddMessage(flags)     
    return cumLength 
    
#---------------------------------------------------------------------------------------------------------------
## The function to clean extrlines based on from and to nodes
## if only one to node and no corresponding from node, except for the highest facc section, marking for deletion
## The same processes are iterated to remove all extra lines
##This is more efficient than another clean line method based on the intersect of the to node points
#---------------------------------------------------------------------------------------------------------------
def UpdateUpstreamLength(inline,dem, outline, field): ##The line should be checked with directions
    #Check_If_Flip_Line_Direction(inline, dem)
    Check_Sink_And_Flip(inline, dem)

    fromnodes=[]
    tonodes=[]
    lines=[]
    lengths = []
    start_ele = []
    ##First loop to get the startpoints, lines and the glacierID list
    with arcpy.da.SearchCursor(inline, ["SHAPE@", "SHAPE@LENGTH"]) as flows: 
        for flow in flows:
            fromnodes.append(((flow[0].firstPoint).X, (flow[0].firstPoint).Y))
            coord= str((flow[0].firstPoint).X)+" "+str((flow[0].firstPoint).Y)
            Cellvalue = arcpy.GetCellValue_management(dem, coord)
            start_ele.append(Cellvalue.getOutput(0))
            tonodes.append(((flow[0].lastPoint).X, (flow[0].lastPoint).Y))
            lines.append(flow[0])
            lengths.append(flow[1])
    del flow, flows
    
    
    start_elearr = np.array(start_ele) ##04 ##Still have problems
    #arcpy.AddMessage(start_elearr)
    #end_elearr = np.array(end_ele) ##04

    ids = np.argsort(start_elearr) ##05 sort the ids based on the elevation from the lowest to highest
    #arcpy.AddMessage(ids)

    cumLengths = [0]*len(ids)
    array=np.array(tonodes)
    for i in range(len(ids)):
        index = ids[i]
        sel_tonode = array[index]
        section_Length = lengths[index]
        flags = [0] * len(ids)
        flags[index] = 1 ##make sure not touch itself
        cumLength = upstreamLength(sel_tonode, fromnodes, tonodes, lines, lengths, flags, section_Length)
        #arcpy.AddMessage(str(cumLength)) 
        cumLengths[index] = cumLength

    ##Delete the line marked for deletion
    arcpy.AddField_management(inline, field, "DOUBLE")
    
    with arcpy.da.UpdateCursor(inline, field) as cursor:
        i = 0
        for row in cursor:
            #if (from_touch[i] == 0) and (to_touch[i] == 0):
            #    cursor.deleteRow()
            #else:
            row[0] = cumLengths[i]
            cursor.updateRow(row)
            i += 1
    del cursor, row

    arcpy.CopyFeatures_management(inline, outline)

    return outline    

##This is the method based on Kienholz
def get_localmax(outline, dem, radius, min_cumdis, cellsize_int):
    outline_line = "in_memory\\outline_line"
    outline_buf = "in_memory\\outline_buf"
    lowestpnt = "in_memory\\lowestpnt"
    highestpnt = "in_memory\\highestpnt"
    lowestpnt_buf = "in_memory\\lowestpnt_buf"
    outline_erased = "in_memory\\outline_erased"
    simplified_outline = "in_memory\\simplified_outline"
    outmostline = "in_memory\\outmostline"

    
    arcpy.PolygonToLine_management(outline, outline_line)

    arcpy.CopyFeatures_management(outline_line, outmostline)

    linearray = arcpy.da.FeatureClassToNumPyArray(outmostline, ('SHAPE@LENGTH'))
    Length_list = np.array([item[0] for item in linearray])
    max_length = max(Length_list)

    with arcpy.da.UpdateCursor(outmostline, 'SHAPE@LENGTH') as cursor:  ##outline_cp is the outmost simplified outlines
        for row in cursor:
            if row[0] < max_length:
                cursor.deleteRow()
    del row, cursor

    ##Find the lowest and highest points of only the major glacier outlines, not the holes
    LineID = arcpy.Describe(outmostline).OIDFieldName

    arcpy.cartography.SimplifyLine(outline_line, simplified_outline, "POINT_REMOVE", int(cellsize_int/2))
    '''
    outZonalmin = ZonalStatistics(outmostline, LineID, dem, "MINIMUM")
    OutConMin = Con(dem == outZonalmin,1)
    arcpy.RasterToPoint_conversion(OutConMin, lowestpnt, "VALUE")

    outZonalmax = ZonalStatistics(outmostline, LineID, dem, "Maximum")
    OutConMax = Con(dem == outZonalmax,1)
    arcpy.RasterToPoint_conversion(OutConMax, highestpnt, "VALUE")
    '''
    #Get the lower 1/3 elevation of the DEM as the elevation threshold
    std_value = dem.standardDeviation
    mean_value = dem.mean #float(arcpy.GetRasterProperties_management(dem,"MINIMUM").getOutput(0))
 

    ele_threshold = int(mean_value - std_value/2.0) ##This elevation is about the lower 30% of the elevation if the DEM is a normal distrbution 

    arcpy.FeatureVerticesToPoints_management(simplified_outline, "in_memory\\line_points", 'All')
    ExtractValuesToPoints("in_memory\\line_points", dem, "in_memory\\line_points_with_alt", "INTERPOLATE", "VALUE_ONLY")
    ###Create the new line after removing the outlier points
    spatialref=arcpy.Describe(outline).spatialReference
    field = 'ORIG_FID' ##The first two fields are FID and Geometry

    pointarray = arcpy.da.FeatureClassToNumPyArray("in_memory\\line_points_with_alt", ('SHAPE@X', 'SHAPE@Y','RASTERVALU',field))
    line_ids = np.array([item[3] for item in pointarray])
    unique_line_ids = np.unique(line_ids)

    t_points = []
    t_angles   = []
    t_FID    = []
    t_cumLen     = []
    t_ele     = []
    
    for fid in unique_line_ids:
        arr = pointarray[line_ids == fid]
        pntx = []
        pnty = []
        pntz = []
        cumLen = []
        dis = 0
        cumLen.append(dis) ##set the first distance to zero
        for i in range(len(arr)):
            pntx.append(arr[i][0])
            pnty.append(arr[i][1])
            pntz.append(int(arr[i][2]))
            if i>0:
                point_dis = Dist(arr[i][0], arr[i][1], arr[i-1][0], arr[i-1][1])
                dis += point_dis
                cumLen.append(dis) ##set the first distance to zero
                
            
        ##Add the second point at the end in order to calculate the value for the end or start point 
        pntx.append(arr[1][0])
        pnty.append(arr[1][1])
        pntz.append(int(arr[1][2]))


        ##NEED to first filter the points and only keep the local min and local maximum points
        #removelist = []
        #removelist.append(0) ##remove the first point because it will be judged at the end
        localmax_X = []
        localmax_Y = []
        localmax_Z = []
        localmax_CumLen = []
        localmax_FID = []
        
        for i in range(len(pntx)):
            if i <(len(pntx)-1) and i > 0:#if it is not first or last point of all
                #if (((pntz[i] < max(pntz[i-1], pntz[i+1])) and ((pntz[i] > min(pntz[i-1], pntz[i+1]))))):
                if (pntz[i] > max(pntz[i-1], pntz[i+1], ele_threshold)):
                    localmax_X.append(pntx[i])
                    localmax_Y.append(pnty[i])
                    localmax_Z.append(pntz[i])
                    localmax_CumLen.append(cumLen[i])
                    localmax_FID.append(fid)

        ##Check the turning points using the cumLen
        if len(localmax_Z) > 1:
            ids = []
            for i in range(1, len(localmax_Z)):
                #point_dis = Dist(t_points[i][0], t_points[i][1], t_points[i-1][0], t_points[i-1][1])
                point_dis = localmax_CumLen[i] - localmax_CumLen[i-1]
                #arcpy.AddMessage(str(point_dis) )         
                if (point_dis < min_cumdis): ##if the two point distance are two close, only pick the largest elevation one
                    if (localmax_Z[i] < localmax_Z[i-1]):
                        index = i
                    else:
                        index = i-1
                    ids.append(index)
           
            if len(ids) > 0:
                delete_multiple_element(localmax_X, ids)
                delete_multiple_element(localmax_Y, ids)
                delete_multiple_element(localmax_Z, ids)
                delete_multiple_element(localmax_CumLen, ids)
                delete_multiple_element(localmax_FID, ids)

        points = np.array(list(zip(localmax_X,localmax_Y)))      
        ##Add the cleaned turning points to         
        t_points.extend(points)
        t_ele.extend(localmax_Z)
        t_cumLen.extend(localmax_CumLen)
        t_FID.extend(localmax_FID)
        
    ##Create the point file
    t_points_FC = arcpy.CreateFeatureclass_management("in_memory", "t_points_FC","POINT","", "DISABLED", "DISABLED", spatialref)
    arcpy.AddField_management(t_points_FC, "LineID", "Long")
    arcpy.AddField_management(t_points_FC, "Elevation", "Long")
    arcpy.AddField_management(t_points_FC, "CumLen", "DOUBLE")
        
    CursorPnts = arcpy.da.InsertCursor(t_points_FC, ['SHAPE@', 'LineID', 'Elevation', 'CumLen'])
    #LineIDlist.append(lineId)
    for i in range(len(t_points)):
        Pnt = arcpy.Point(t_points[i][0], t_points[i][1])
        CursorPnts.insertRow((Pnt, t_FID[i], int(t_ele[i]), t_cumLen[i])) 
    del CursorPnts ##Delete cursor object

    #arcpy.CopyFeatures_management(t_points_FC, "d:\\temp\\left_t_points_FC.shp")
    ##use the radius to check the local max
    final_points = "in_memory\\final_points"
    #arcpy.AddMessage("Start the buffer erase...")
    if len(t_ele) > 1:
        ZArr = np.array(t_ele) 
        ids = ZArr.argsort()[::-1] ##05 sort the ids based on the elevation from largest to lowest
        point_buf = "in_memory\\point_buf"
        #arcpy.AddMessage(str(len(ids)))
        #arcpy.AddMessage(str(radius))
        for i in range(len(ids)): 
            elev = int(ZArr[ids[i]])
            #arcpy.AddMessage("The elevation is: " +str(elev))
            query = "Elevation = " + str(elev)
            arcpy.Select_analysis (t_points_FC, "in_memory\\single_point", query)
            pntarray = arcpy.da.FeatureClassToNumPyArray("in_memory\\single_point", 'OID@')
            #arcpy.AddMessage("The number of points is: " +str(len(pntarray)))
            if len(pntarray) > 0:
                if i < 1:
                    #arcpy.AddMessage("copy the highest point")
                    arcpy.CopyFeatures_management("in_memory\\single_point", final_points)
                else:
                    #arcpy.AddMessage("Append the local max point")
                    arcpy.Append_management("in_memory\\single_point", final_points,"NO_TEST")
                if i < len(ids)-1: ##for the last point, do nothing  
                    ##erase the final points from all_points
                    arcpy.Erase_analysis(t_points_FC, final_points, "in_memory\\left_points")
                    arcpy.Buffer_analysis("in_memory\\single_point", "in_memory\\single_point_buf", (str(radius)+ " Meter"))
                    #arcpy.CopyFeatures_management("in_memory\\single_point_buf", "d:\\temp\\single_point_buf.shp")
                    arcpy.Clip_analysis(outline, "in_memory\\single_point_buf", "in_memory\\clipped_outline")
                    arcpy.MultipartToSinglepart_management("in_memory\\clipped_outline", "in_memory\\clipped_outline_singleParts")
                    ##spatial the single part that intersect with the selected point
                    arcpy.SpatialJoin_analysis("in_memory\\clipped_outline_singleParts", "in_memory\\single_point", "in_memory\\singleParts_spatialjoin", "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
                    ##erase the left_over points           
                    arcpy.Erase_analysis("in_memory\\left_points", "in_memory\\singleParts_spatialjoin", t_points_FC)
                    #arcpy.CopyFeatures_management(t_points_FC, "d:\\temp\\left_t_points_FC.shp")
    else:
        arcpy.CopyFeatures_management(t_points_FC, final_points)

    return final_points

#------------------------------------------------------------------------------------------------------------
# This function creates equally spaced points along a line from start to end.
# It is revised from the codes in Volta.
#------------------------------------------------------------------------------------------------------------
def points_on_line_withID(line, resolution):
    new_points = arcpy.CreateFeatureclass_management("in_memory", "points","POINT", line)
    arcpy.AddField_management(new_points, "distance", "DOUBLE")
    arcpy.AddField_management(new_points, "line_id", "LONG")
    new_points_cursor = arcpy.da.InsertCursor(new_points, ('SHAPE@',"distance","line_id"))
    FcID = arcpy.Describe(line).OIDFieldName
    with arcpy.da.SearchCursor(line, ["SHAPE@LENGTH", "SHAPE@", FcID]) as cursor:
        for row in cursor:
            start = 0
            length = row[0]
            while start < length:
                new_point = row[1].positionAlongLine(start)
                new_points_cursor.insertRow((new_point,start,row[2]))
                start = start + resolution          
    del row, cursor
    del new_points_cursor

    return new_points

#------------------------------------------------------------------------------------------------------------
# This function creates equally spaced points along a line from start to end.
# It is revised from the codes in Volta.
#------------------------------------------------------------------------------------------------------------
def compute_flowline_upslope_index (flowline, dem, resolution):
    flowline_points = points_on_line_withID(flowline, resolution)
    flowline_points_with_alt = "in_memory\\flowline_points_with_alt"
    ExtractValuesToPoints(flowline_points, dem, flowline_points_with_alt, "INTERPOLATE")

    pointarray = arcpy.da.FeatureClassToNumPyArray(flowline_points_with_alt, ('line_id','RASTERVALU'))
    line_ids = np.array([item[0] for item in pointarray])
    ele_arr = np.array([item[1] for item in pointarray])
    unique_line_ids = np.unique(line_ids)
    m_list = []
    lineid_list = []
    #arcpy.AddMessage("The number of centerlines is: " + str(line_ids))
    
    for line_id in unique_line_ids:
        lineid_list.append(line_id)
        line_ele_arr = ele_arr[line_ids == line_id]
        sum_Zup = 0
        m = 0
        ##For small glaicers with the length < resolution, there will be only one start points
        if len(line_ele_arr) > 1:
            start_eles = line_ele_arr[0:-1]
            end_eles = line_ele_arr[1:]
            ele_diffs = end_eles - start_eles

            #arcpy.AddMessage(ele_diffs)
            sum_Zup = sum(ele_diffs[ele_diffs > 10]) ## if the upslope is < 10 m out of 100 m, it should be OK
            #m = 0
            #arcpy.AddMessage("Sum_zUp is: " + str(sum_Zup))
            if sum_Zup > 0:
                ##Determine the max number of up points
                index = np.where(ele_diffs > 10)[0]
                #arcpy.AddMessage(index)
                max_num_Zup = 1 ##at lease one point
                if len(index) > 1:
                    index_diff = ((index[1:] - index[0:-1]) == 1)
                    #arcpy.AddMessage(index_diff)
                    
                    consecutive_Zup_number = np.diff(np.where(np.concatenate(([index_diff[0]], index_diff[:-1] != index_diff[1:], [True])))[0])[::2]
                    #arcpy.AddMessage(consecutive_Zup_number)
                    if sum(consecutive_Zup_number) > 0:
                        max_num_Zup = max(consecutive_Zup_number) + 1 ##Add 1 becasue the number of difference is from number + 1 points
                    else:
                        max_num_Zup = 1
                    #arcpy.AddMessage("max_num_Zup: " + str(max_num_Zup))
            
                m = 0.1 * max_num_Zup + 0.01 * pow(sum_Zup, 0.7)
                #arcpy.AddMessage("m = " + str(m))
        m_list.append(m)
        
    #arcpy.AddMessage(m_list)
    return m_list, lineid_list


def cost_path_centerline (normalreverseDis, normalDEM, dem, lowestpnt, local_max, cellsize_int, resolution):

    #outBkLinkRaster = "in_memory\\outbklink"
    outBkLinkRaster = arcpy.env.scratchGDB + "\\outbklink"
    all_centerlines = "in_memory\\all_centerlines"
    local_max_snap = "in_memory\\local_max_snap"
    #arcpy.CopyFeatures_management(local_max, local_max_cp)
    
    a = 4.25
    b0 = 3.5
    del_bmax = 0.5
    PntID = arcpy.Describe(lowestpnt).OIDFieldName
    PntID2 = arcpy.Describe(local_max).OIDFieldName

    #outSnapLowestPnt = SnapPourPoint(lowestpnt, costRaster, cellsize_int, PntID)
    
    outSnapLowestPnt = SnapPourPoint(lowestpnt, (1 - normalDEM), cellsize_int*3, PntID) ##Increase the snap distance to 3 times of cellsize
    #outSnapLowestPnt.save("c:\\test\\outSnapLowestPnt.tif")

    
    #arcpy.CopyFeatures_management(local_max, "c:\\test\\local_max_cp.shp")
    #outSnaplocalMax = SnapPourPoint(local_max_cp, costRaster, cellsize_int, PntID2)
    outSnaplocalMax = SnapPourPoint(local_max, normalDEM, cellsize_int*3, PntID2) ##Increase the snap distance to 3 times of cellsize to the highest elevation
    #outSnaplocalMax.save("c:\\test\\outSnaplocalMax.tif")
    arcpy.RasterToPoint_conversion(outSnaplocalMax, local_max_snap, "VALUE")
    #arcpy.CopyFeatures_management(local_max_snap, "c:\\test\\local_max_snap.shp")
    
    costpath_flowlines = arcpy.CreateFeatureclass_management("in_memory", "costpath_flowlines","POLYLINE","","","",local_max)
    costDis = Power(normalreverseDis * 1000, a)
    #costDis.save("c:\\test\\costDis.tif")
    
    for i in range(6):
        delta_b = i * 0.1
        b = b0 + i * 0.1
        #arcpy.AddMessage("b increases to: " + str(b))
        # Derive the cost raster
        costRaster = costDis + Power(normalDEM * 3000, b)
        #costRaster.save("c:\\test\\costRaster.tif")

        outCostDistance = CostDistance(outSnapLowestPnt, costRaster, "#", outBkLinkRaster)
        #arcpy.AddMessage("Moximum cost distance is: " + str(outCostDistance.maximum))
        if str(outCostDistance.maximum) == "None":
            #arcpy.AddMessage("Change the snapLowestPnt")
            outSnapLowestPnt2 = SnapPourPoint(lowestpnt, (1-normalreverseDis), cellsize_int*3, PntID)
            #outSnapLowestPnt2.save("d:\\temp\\outSnapLowestPnt2.tif")
            outCostDistance = CostDistance(outSnapLowestPnt2, costRaster, "#", outBkLinkRaster)
            
        #outCostDistance.save("c:\\test\\outCostDistance.tif")
        #arcpy.CopyRaster_management(outBkLinkRaster, "c:\\test\\outBkLinkRaster.tif")
        #outBkLinkRaster.save("d:\\temp\\outBkLinkRaster.tif")

        #CostPathAsPolyline(outSnaplocalMax, outCostDistance, outBkLinkRaster, all_centerlines, "EACH_CELL", "#", "INPUT_RANGE")
        CostPathAsPolyline(local_max_snap, outCostDistance, outBkLinkRaster, all_centerlines, "EACH_CELL", "#", "#")

        ##check if there is line created
        linearr  = arcpy.da.FeatureClassToNumPyArray(all_centerlines, ('OID@'))
        #arcpy.AddMessage("Number of costpath_centerlines: " + str(len(linearr)))
        if len(linearr)==0:
            ##use the maximum outCostDistance Raster
            arcpy.AddMessage("Snap the highest outcostdistance point")
            PntID3 = arcpy.Describe(local_max_snap).OIDFieldName
            #outSnaplocalMax = SnapPourPoint(local_max_snap, normalreverseDis, cellsize_int*3, PntID3)
            outSnaplocalMax = SnapPourPoint(local_max_snap, outCostDistance, cellsize_int*3, PntID3)
            #outSnaplocalMax.save("d:\\temp\\outSnaplocalMax.tif")
            CostPathAsPolyline(outSnaplocalMax, outCostDistance, outBkLinkRaster, all_centerlines, "EACH_CELL", "#", "#")

        linearr  = arcpy.da.FeatureClassToNumPyArray(all_centerlines, ('OID@'))
        #arcpy.AddMessage("Number of costpath_centerlines: " + str(len(linearr)))

        if len(linearr)==0:
            arcpy.AddMessage("The glacier outline is too small to generate the centerline(s)!")
            break

        #arcpy.CopyFeatures_management(all_centerlines, "c:\\test\\all_centerlines2.shp")

        ##Need to clean up the costpath created centerlines with the length of < the resolution, otherwise there is potential errors
        #linearr  = arcpy.da.FeatureClassToNumPyArray(all_centerlines, ('OID@'))
        if len(linearr)==1: ##If only having one centerline, return the line directly
            #arcpy.AddMessage("Just one line, return directly with the line")
            arcpy.Append_management(all_centerlines, costpath_flowlines, "NO_TEST")
            arcpy.Delete_management(outBkLinkRaster)
            return costpath_flowlines
            
        with arcpy.da.UpdateCursor(all_centerlines, ['SHAPE@LENGTH']) as cursor:
            for row in cursor:
                if row[0] < resolution:
                    cursor.deleteRow()
        del row, cursor    
        
         
        lineID = arcpy.Describe(all_centerlines).OIDFieldName

        m_list, lineid_list = compute_flowline_upslope_index (all_centerlines, dem, resolution)
        #arcpy.AddMessage(m_list)
        #arcpy.AddMessage(str(len(m_list)))

        leftover = len(m_list)
        for j in range(len(m_list)):
            m = m_list[j]
            if m <= delta_b:
                ##Copy the correcponding line to output line and remove it from local_lax
                query = lineID + " = " + str(lineid_list[j])
                arcpy.Select_analysis (all_centerlines, "in_memory\\selected_centerline", query)
                arcpy.Append_management("in_memory\\selected_centerline", costpath_flowlines, "NO_TEST")
                #arcpy.AddMessage("Added one line!")

        #arcpy.CopyFeatures_management(costpath_flowlines, "c:\\test\\costpath_flowlines.shp") 
        ##use the selected line to remove the local max point
        linearr  = arcpy.da.FeatureClassToNumPyArray(costpath_flowlines, ('OID@'))
        #arcpy.AddMessage("Number of costpath_centerlines: " + str(len(linearr)))
        if len(linearr) == len(m_list):
            #arcpy.AddMessage("The same number of centerlines with the M-list")
            break
        #arcpy.AddMessage(str(len(linearr)))    
        if len(linearr) > 0:
            arcpy.SpatialJoin_analysis(local_max_snap, costpath_flowlines, "in_memory\\local_max_spatialjoin", "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "100 Meters", "#")
            pntarr  = arcpy.da.FeatureClassToNumPyArray("in_memory\\local_max_spatialjoin", ('OID@'))
            #arcpy.AddMessage("Selected local max point: " + str(len(pntarr)))
            arcpy.Erase_analysis(local_max_snap, "in_memory\\local_max_spatialjoin", "in_memory\\leftover_localmax")
            ##Need to update localMax
            arcpy.CopyFeatures_management("in_memory\\leftover_localmax", local_max_snap)
            pntarr2  = arcpy.da.FeatureClassToNumPyArray(local_max_snap, ('OID@'))
            #arcpy.AddMessage("Number of leftover local max points: " + str(len(pntarr2)))
            leftover = len(pntarr2)

        if i == 5: ##the last loop, add all centerlines to costpath_flowlines
            arcpy.Append_management(all_centerlines, costpath_flowlines, "NO_TEST")
            
        if leftover == 0: ##if all local max points are deleted; all lines are corrected, then break the loop
            #arcpy.AddMessage("All lines are corrected!")
            break

    #arcpy.Delete_management(outCostDistance)
    arcpy.Delete_management(outBkLinkRaster)
    #arcpy.AddMessage("Return costpath_flowline")
    return costpath_flowlines

def Centerline_for_single_outline (glacier_outline, dem, cellsize_int, radius):
    
    outline_line = "in_memory\\outline_line"
    outline_buf = "in_memory\\outline_buf"
    lowestpnt = "in_memory\\lowestpnt"
    highestpnt = "in_memory\\highestpnt"
    simplified_outline = "in_memory\\simplified_outline"
    outline_cp = "in_memory\\outline_cp"

    ##buf the sel_outline a little bit to clip the DEM
    arcpy.Buffer_analysis(glacier_outline, outline_buf, (str(int(cellsize_int/2))+ " Meter"))
  
    
    # Step 1: get the lowest point and the local maximum points
    arcpy.PolygonToLine_management(glacier_outline, outline_line)

    arcpy.CopyFeatures_management(outline_line, outline_cp)
    linearray = arcpy.da.FeatureClassToNumPyArray(outline_cp, ('SHAPE@LENGTH'))
    Length_list = np.array([item[0] for item in linearray])
    max_length = max(Length_list)
    #arcpy.AddMessage(str(max_length))
    #arcpy.AddMessage(str(len(Length_list)))

    with arcpy.da.UpdateCursor(outline_cp, 'SHAPE@LENGTH') as cursor:  ##outline_cp is the outmost simplified outlines
        for row in cursor:
            if row[0] < max_length:
                cursor.deleteRow()
    del row, cursor

    ##Find the lowest and highest points of only the major glacier outlines, not the holes
    LineID = arcpy.Describe(outline_cp).OIDFieldName

    arcpy.cartography.SimplifyLine(outline_cp, simplified_outline, "POINT_REMOVE", int(cellsize_int/2))

    outZonalmin = ZonalStatistics(outline_cp, LineID, dem, "MINIMUM")
    OutConMin = Con(dem == outZonalmin,1)
    arcpy.RasterToPoint_conversion(OutConMin, lowestpnt, "VALUE")
    #arcpy.CopyFeatures_management(lowestpnt, "d:\\temp\\lowestpnt.shp") 

    outZonalmax = ZonalStatistics(outline_cp, LineID, dem, "Maximum")
    OutConMax = Con(dem == outZonalmax,1)
    arcpy.RasterToPoint_conversion(OutConMax, highestpnt, "VALUE")

    local_max = get_localmax(glacier_outline, dem, radius, 500, cellsize_int)
    #arcpy.CopyFeatures_management(local_max, "d:\\temp\\local_max.shp") 
    
    pointarray = arcpy.da.FeatureClassToNumPyArray(local_max, ('SHAPE@LENGTH'))
    point_list = np.array([item[0] for item in pointarray])
    #arcpy.AddMessage(str(len(point_list)))
    if len(point_list) == 0: ## if no point, add the highest one
        arcpy.Append_management(highestpnt, local_max, "NO_TEST") ##only add the lowest point
    
    ##Step 2: create the cost raster
    # Process: Euclidean Distance
    tempEnvironment0 = arcpy.env.mask
    #arcpy.env.mask = outline_buf
    arcpy.env.mask = glacier_outline
    #distancerst = "in_memory\\distancerst"
    outEucDistance = EucDistance(outline_line, "#", dem)
    #outEucDistance.save("c:\\test\\outEucDistance.tif")

    max_value = dem.maximum 
    min_value = dem.minimum

    normalDEM = (dem - min_value) / (max_value - min_value)

    #normalDEM.save("d:\\temp\\normalDEM.tif")

    ##Normalize the reversed distance raster
    max_dis = outEucDistance.maximum 
    min_dis = outEucDistance.minimum 

    normalreverseDis = (max_dis - outEucDistance) / (max_dis - min_dis)

    #arcpy.Delete_management(outEucDistance)
    ##reset the mask environment
    arcpy.env.mask = tempEnvironment0

    all_centerlines = cost_path_centerline (normalreverseDis, normalDEM, dem, lowestpnt, local_max, cellsize_int, 100)

    #arcpy.CopyFeatures_management(all_centerlines, "d:\\temp\\all_centerlines.shp")
    
    ##Method 01:
    #Simply dissolve the all centerlines and then use upstream method to merge the centerlines
    #arcpy.Dissolve_management(all_centerlines, "in_memory\\dissolve_centerlines", "#", "#", 'SINGLE_PART', '#')
    #Check_If_Flip_Line_Direction(all_centerlines, dem)
    #final_centerlines = "in_memory\\final_centerlines"
    #UpdateUpstreamLength("in_memory\\dissolve_centerlines",dem, final_centerlines, "CumLen")
    
    ##Method 02 based on the published paper of Kienholz 2014
    ##Add Length into the attribute table

    final_centerlines = "in_memory\\final_centerlines"

    lineArray = arcpy.da.FeatureClassToNumPyArray(all_centerlines,['OID@'])
    if len(lineArray) ==0:
        return final_centerlines
    
    arcpy.AddField_management(all_centerlines, "Len", "Long")
    with arcpy.da.UpdateCursor(all_centerlines, ["Len", 'SHAPE@LENGTH']) as cursor:
        for row in cursor:
            row[0] = int(row[1])
            cursor.updateRow(row)
    del row, cursor    
    #arcpy.CopyFeatures_management(all_centerlines, "d:\\temp\\all_centerlines.shp")
    
    ##Step 3: using buffer to erase the centerlines within a certain distance
    lineArray = arcpy.da.FeatureClassToNumPyArray(all_centerlines,['Len'])
    Length_list = np.array([item[0] for item in lineArray])
    #arcpy.AddMessage(Length_list)
    lenarr = np.array(Length_list) 
    #arcpy.AddMessage(lenarr)
    ids = lenarr.argsort()[::-1] ##05 sort the ids based on the facc from largest to lowest
    centerline_buf = "in_memory\\centerline_buf"
    if len(ids) > 1: ##at least two lines
        for i in range(len(ids)): ##Don't need to buffer the last line
            length = lenarr[ids[i]]
            #arcpy.AddMessage("The length is: " +str(length))
            query = "Len = " + str(length)
            #arcpy.AddMessage(query)
            arcpy.Select_analysis (all_centerlines, "in_memory\\single_centerline", query)
            if i < 1:
                arcpy.CopyFeatures_management("in_memory\\single_centerline", final_centerlines)
            else:
                arcpy.Append_management("in_memory\\single_centerline", final_centerlines,"NO_TEST")
            ##The last line will not do the following
            if i < len(ids)-1:
                ##erase the final centerlines from all_centerlines
                #The radius calcualted based on the paper is too big; set three times of the cell size here!!!
                arcpy.Buffer_analysis(final_centerlines, centerline_buf, (str(int(cellsize_int*3))+ " Meter"))
                #arcpy.Buffer_analysis(final_centerlines, centerline_buf, (str(radius)+ " Meter"))
                arcpy.Erase_analysis(all_centerlines, centerline_buf, "in_memory\\left_centerline")
                #arcpy.AddMessage("Erase_analysis")
                arcpy.MultipartToSinglepart_management("in_memory\\left_centerline", "in_memory\\left_centerline_singleParts")
                ##select the line close to the local max
                arcpy.SpatialJoin_analysis("in_memory\\left_centerline_singleParts", local_max, "in_memory\\left_centerline_spatialjoin", "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "100 Meters", "#")
                
                #arcpy.CopyFeatures_management("in_memory\\left_centerline_spatialjoin", "d:\\temp\\left_centerline_spatialjoin.shp")
                ##Need to update all centerlines
                arcpy.CopyFeatures_management("in_memory\\left_centerline_spatialjoin", all_centerlines)


        ##Delete the small length lines < radius
        #arcpy.AddMessage("Delete the small length lines < radius")
        lineArray2 = arcpy.da.FeatureClassToNumPyArray(final_centerlines,['SHAPE@LENGTH'])
        len_list = np.array([item[0] for item in lineArray2]) 
        max_len = max(len_list)        
        #arcpy.AddMessage("The max_len is: " + str(max_len))
        if len(len_list) > 1:  
            #number_deleted = 0        
            with arcpy.da.UpdateCursor(final_centerlines, 'SHAPE@LENGTH') as cursor:  ##outline_cp is the outmost simplified outlines
                for row in cursor:
                    if row[0] < min(max_len, radius):
                        #arcpy.AddMessage("Delete one")
                        #number_deleted += 1
                        cursor.deleteRow()
            del row, cursor
    else:
        arcpy.CopyFeatures_management(all_centerlines, final_centerlines)
    

    return final_centerlines

def Connect_centerline (outline, centerlines, dem, search_dis): ##for a single outline

    Check_If_Flip_Line_Direction(centerlines, dem)

    #arcpy.CopyFeatures_management(centerlines, "d:\\temp\\centerlines.shp")

    flowlineStartNodes = "in_memory\\flowlineStartNodes"
    flowlineStartNodes_with_ele = "in_memory\\flowlineStartNodes_with_ele"
    flowlines_with_start_ele = "in_memory\\flowlines_with_start_ele"
    sel_flowlines = "in_memory\\sel_flowlines"
    lowest_flowline = "in_memory\\lowest_flowline"
    lowest_flowline_points = "in_memory\\lowest_flowline_points"
    otherline_start_nodes = "in_memory\\otherline_start_nodes"

    connected_flowlines = arcpy.CreateFeatureclass_management("in_memory", "connected_flowlines","POLYLINE","","","",centerlines)

    arcpy.Clip_analysis(centerlines, outline, sel_flowlines)
    arcpy.MultipartToSinglepart_management(sel_flowlines, "in_memory\\sel_flowlines_singleParts")
    with arcpy.da.UpdateCursor("in_memory\\sel_flowlines_singleParts", ['SHAPE@LENGTH']) as cursor:
        for row in cursor:
            if row[0] < 50: ##if the line is less than 50 m, delete
                cursor.deleteRow()
    del row, cursor
     
    ## find the lowest flowline
    arcpy.FeatureVerticesToPoints_management("in_memory\\sel_flowlines_singleParts", flowlineStartNodes, "START")
    ExtractValuesToPoints(flowlineStartNodes, dem, flowlineStartNodes_with_ele)
    arcpy.SpatialJoin_analysis("in_memory\\sel_flowlines_singleParts", flowlineStartNodes_with_ele, flowlines_with_start_ele, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "0.1 Meter", "#")

    #arcpy.CopyFeatures_management(flowlines_with_start_ele, "d:\\temp\\flowlines_with_start_ele.shp")

    flowlineArr = arcpy.da.FeatureClassToNumPyArray(flowlines_with_start_ele, 'RASTERVALU')
    heights = np.array([item[0] for item in flowlineArr])
    ##If only one height value (one flowline), do nothing and return
    if len(heights) < 2:
        #arcpy.AddMessage("Only one centerline in the glacier outlines")
        arcpy.Append_management("in_memory\\sel_flowlines_singleParts", connected_flowlines, "NO_TEST")
        return connected_flowlines
    else:
        minheight = min(heights)
        #arcpy.AddMessage(str(minheight))
        ##Copy the line with minheight as a seperated flowline
        query = 'RASTERVALU' +" = " + str(minheight)
        #arcpy.AddMessage(query)
        arcpy.Select_analysis(flowlines_with_start_ele,lowest_flowline,query)

        ##remove the lowest flowline
        arcpy.Erase_analysis(flowlines_with_start_ele, lowest_flowline, "in_memory\\flowlines_without_lowestflowlines")

        ##Check if there is flowlines_without_lowestflowlines
        flowlineArr2 = arcpy.da.FeatureClassToNumPyArray("in_memory\\flowlines_without_lowestflowlines", 'OID@')
        if (len(flowlineArr2)== 0): 
            #arcpy.AddMessage("Not other flowlines left except the lowest one")
            ##Only keep the first one of the lowest flowlines
            with arcpy.da.UpdateCursor(lowest_flowline, ['OID@']) as cursor:
                i = 0
                for row in cursor:
                    if i > 0: ##
                        cursor.deleteRow()
                    i += 1
            del row, cursor

            arcpy.Append_management(lowest_flowline, connected_flowlines, "NO_TEST")
            return connected_flowlines
        
        ##Do the connection with the lowest flowline
        arcpy.FeatureVerticesToPoints_management(lowest_flowline, lowest_flowline_points, "ALL") ## check if the centerline can be extended to the nearest vertices
        #arcpy.CopyFeatures_management(lowest_flowline_points, "d:\\temp\\lowest_flowline_points.shp")

        #arcpy.FeatureVerticesToPoints_management("in_memory\\flowlines_without_lowestflowlines", otherline_start_nodes, "START") ## check if the centerline can be extended to the nearest vertices
        arcpy.FeatureVerticesToPoints_management("in_memory\\flowlines_without_lowestflowlines", otherline_start_nodes, "BOTH_ENDS") ## check if the centerline can be extended to the nearest vertices
        #arcpy.CopyFeatures_management(otherline_start_nodes, "d:\\temp\\otherline_start_nodes.shp")

        bLoop = True
        numLoop = 1
        k = 0
        while (bLoop):
            search_radius = str(numLoop * search_dis + 1) + " Meters"
            #arcpy.AddMessage(search_radius)
            #arcpy.Near_analysis(otherline_start_nodes, lowest_flowline_points, search_radius, "LOCATION", "")
            arcpy.Near_analysis(otherline_start_nodes, lowest_flowline, search_radius, "LOCATION", "")
            #arcpy.CopyFeatures_management(otherline_start_nodes, "d:\\temp\\otherline_start_nodes.shp")
            leftover = 0
            x1 = []
            y1 = []
            x2 = []
            y2 = []
            NearID = []
            #leftID = 0
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
            #if (leftID > len(NearID)):
            #    leftover += 1
            number_connected = 0
            if len(NearID) > 0:
                arcpy.AddMessage("Create connection lines")
                connectionline = arcpy.CreateFeatureclass_management("in_memory", "connectionline","POLYLINE","","","",centerlines)
                new_line_cursor = arcpy.da.InsertCursor(connectionline, ('SHAPE@'))
                for j in range(len(NearID)):
                    #arcpy.AddMessage("Add line: " + str(j))
                    array = arcpy.Array([arcpy.Point(x2[j],y2[j]),arcpy.Point(x1[j], y1[j])])
                    polyline = arcpy.Polyline(array)
                    new_line_cursor.insertRow([polyline])
                del new_line_cursor
                number_connected = len(NearID)

                arcpy.Erase_analysis(flowlines_with_start_ele, lowest_flowline, "in_memory\\flowlines_without_lowestflowlines")

                arcpy.SpatialJoin_analysis("in_memory\\flowlines_without_lowestflowlines", connectionline, "in_memory\\connectedflowlines", "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")

                arcpy.AddField_management("in_memory\\connectedflowlines", "LineID", "Long")
                arcpy.arcpy.CalculateField_management("in_memory\\connectedflowlines","LineID",str("!"+str(arcpy.Describe("in_memory\\connectedflowlines").OIDFieldName)+"!"),"PYTHON_9.3")
                #arcpy.CopyFeatures_management("in_memory\\connectedflowlines", "d:\\temp\\connectedflowlines.shp")

                arcpy.SpatialJoin_analysis(connectionline, connectionline, "in_memory\\connectedline_spatialjoin", "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
                arcpy.AddField_management("in_memory\\connectedline_spatialjoin", "LineID", "Long")
                #connectlineArr = arcpy.da.FeatureClassToNumPyArray("in_memory\\connectedline_spatialjoin", 'OID@')
                #arcpy.AddMessage(str(len(connectlineArr)))
                with arcpy.da.UpdateCursor("in_memory\\connectedline_spatialjoin", ["TARGET_FID", "LineID"]) as cursor:
                    for row in cursor:
                        row[1] = row[0]
                        cursor.updateRow(row)
                del cursor, row
                #arcpy.CopyFeatures_management("in_memory\\connectedline_spatialjoin", "d:\\temp\\connectedline_spatialjoin.shp")
                arcpy.Append_management("in_memory\\connectedline_spatialjoin", "in_memory\\connectedflowlines", "NO_TEST")
                #arcpy.CopyFeatures_management("in_memory\\connectedflowlines", "d:\\temp\\connectedflowlines2.shp")

                arcpy.Dissolve_management("in_memory\\connectedflowlines", "in_memory\\dissolve_connectedflowlines", "LineID","", "SINGLE_PART")
                ##Need to smooth the line here

                arcpy.Append_management("in_memory\\dissolve_connectedflowlines", lowest_flowline, "NO_TEST")

            
            #arcpy.AddMessage("leftover is " + str(leftover))
            if leftover > 0:
                if number_connected < 1: ##Need to increase the search distance
                    numLoop += 1
                    arcpy.AddMessage("numLoop is " + str(numLoop))
                else: ##Need to update the points for lowese_flowlines
                    arcpy.FeatureVerticesToPoints_management(lowest_flowline, lowest_flowline_points, "ALL") ## check if the centerline can be extended to the nearest vertices
            else:
                bLoop = False ## stop the loop
            k += 1
            if k > 5: ##if more than five loops quit
                arcpy.AddMessage("There is more than five loops, quite the loop")
                bLoop = False
        arcpy.Append_management(lowest_flowline, connected_flowlines, "NO_TEST")

    return connected_flowlines

#------------------------------------------------------------------------------------------------------------
# This function splits the flowline to seperated lines by the intersection points. The flowlines for paleo glacier
# reconstruction are lines from the source to end. The tools splits the flowlines back to the inital lines derived
# from the stream network.
#------------------------------------------------------------------------------------------------------------
def line_split_by_intersect(lines):
    ##get the start and end points of the lines
    arcpy.FeatureVerticesToPoints_management(lines, "in_memory\\line_points", 'BOTH_ENDS')
    arcpy.DeleteIdentical_management("in_memory\\line_points", "Shape")
    arcpy.SplitLineAtPoint_management(lines, "in_memory\\line_points", "in_memory\\split_line", "1 Meters")

    return "in_memory\\split_line"

#------------------------------------------------------------------------------------------------------------
# This function creates a set of points along the line based on a specified distance
# The code is revised from Volta centerline.
#------------------------------------------------------------------------------------------------------------
def points_on_line(line, resolution):   ###Need to add the end point into the output point feature. In this way, extend line is not necessary??!!!
    new_points = arcpy.CreateFeatureclass_management("in_memory", "points","POINT", line)
    new_points_cursor = arcpy.da.InsertCursor(new_points, ('SHAPE@'))
    start = 0
    with arcpy.da.SearchCursor(line, ["SHAPE@LENGTH", "SHAPE@"]) as cursor:
        for row in cursor:
            length = row[0]
            while start < length:
                new_point = row[1].positionAlongLine(start)
                new_points_cursor.insertRow((new_point,))
                start = start + resolution
    del row, cursor

    del new_points_cursor

    return new_points


#---------------------------------------------------------------------------------------------------------------
def cleandangleline(inputline, outline): ##The line should be checked with directions

    inline = "in_memory\\inline"
    arcpy.CopyFeatures_management(inputline, inline)

    fromnodes=[]
    tonodes=[]
    lines=[]
    ##First loop to get the startpoints, lines and the glacierID list
    with arcpy.da.SearchCursor(inline, ["SHAPE@"]) as flows: 
        for flow in flows:
            fromnodes.append(((flow[0].firstPoint).X, (flow[0].firstPoint).Y))
            tonodes.append(((flow[0].lastPoint).X, (flow[0].lastPoint).Y))
            lines.append(flow[0])
    del flow, flows

    array=np.array(tonodes)
    to_touch = []
    for i in range(len(array)):
        nTouches=0
        punto=arcpy.Point(array[i][0],array[i][1])
        for a in range (len(array)):
            if a != i: ##don't touch itself
                touches=(punto.touches(lines[a]) or punto.within(lines[a]))
                if touches==True: ##if the start point touches others
                    nTouches += 1
        to_touch.append(nTouches)
    #arcpy.AddMessage(to_touch)

    array=np.array(fromnodes)
    from_touch = []
    for i in range(len(array)):
        nTouches=0
        punto=arcpy.Point(array[i][0],array[i][1])
        for a in range (len(array)):
            if a != i: ##don't touch itself
                touches=(punto.touches(lines[a]) or punto.within(lines[a]))
                if touches==True: ##if the start point touches others
                    nTouches += 1
        from_touch.append(nTouches)

    #arcpy.AddMessage(from_touch)

    if len(from_touch) > 1: ##only process when there are more than one lines
        with arcpy.da.UpdateCursor(inline, ["SHAPE@LENGTH", "SHAPE@"]) as cursor:
            i = 0
            for row in cursor:
                if (from_touch[i] == 0) and (to_touch[i] == 0): 
                    arcpy.AddMessage("Delete one line away from the lowest point!")
                    cursor.deleteRow()
                i += 1
        del cursor, row

    arcpy.CopyFeatures_management(inline, outline)

    return outline    

   
def Centerline_costPath(glacier_outlines, dem, flow_lines):
    arcpy.Delete_management("in_memory")

    ##make sure the projection of the glacier outlines is the same with the UTM and the same with the DEM projection 
    #spatial_flowline = arcpy.Describe(flowline).spatialReference
    spatial_outline = arcpy.Describe(glacier_outlines).spatialReference
    spatial_ref_dem = arcpy.Describe(dem).spatialReference

    if "UTM" in spatial_ref_dem.name:
        arcpy.AddMessage("The DEM projection is: " + spatial_ref_dem.name)
    else:
        arcpy.AddMessage("The DEM projection is not UTM. Please re-project the DEM to a UTM projection for the analysis!")
        exit()   

    if "UTM" in spatial_outline.name:
        arcpy.AddMessage("The outline projection is: " + spatial_outline.name)
    else:
        arcpy.AddMessage("The outline projection is not UTM. Please re-project the outline to a UTM projection for the analysis!")
        exit()   

    if spatial_outline.name != spatial_ref_dem.name:
        arcpy.AddMessage("The DEM and galcier outlines have different UTM projections. Please re-project them to the same UTM projection!")
        exit()   

    
    arcpy.env.snapRaster = dem

    ##Get the area of the outline
    FcID = arcpy.Describe(glacier_outlines).OIDFieldName
    polyarray = arcpy.da.FeatureClassToNumPyArray(glacier_outlines, ('SHAPE@AREA', FcID))
    area_list = np.array([item[0] for item in polyarray])
    FIds = np.array([item[1] for item in polyarray])

    radius_list = area_list * 2e-6 + 500
    smooth_list = area_list * 2e-6 + 200
    for i in range(len(radius_list)):
        if radius_list[i] > 1000:
            radius_list[i] = 1000
        if smooth_list[i] > 400:
            smooth_list[i] = 400

    ##Stream network analysis
    cellsize = arcpy.GetRasterProperties_management(dem,"CELLSIZEX")
    cellsize_int = int(cellsize.getOutput(0))

    #arcpy.AddMessage(str(cellsize_int))


    glacier_outline = "in_memory\\glacier_outline" ##Single outline for each loop
    #flow_line_cp = "in_memory\\flow_line_cp"
    Fc_Added = 0
    for i in range(len(FIds)):
        #try:
        arcpy.AddMessage("Generating centreline(s) for glacier "+str(i+1)+" of "+str(len(FIds)))
        query = FcID +" = "+str(FIds[i])

        arcpy.Select_analysis(glacier_outlines,glacier_outline,query)


        ##buf the sel_outline a little bit to clip the DEM
        arcpy.Buffer_analysis(glacier_outline, "in_memory\\outline_buf", (str(cellsize_int*3)+ " Meter"))

        
        clipped_dem = ExtractByMask(dem,"in_memory\\outline_buf")

        radius = radius_list[i]
        smooth_dis = smooth_list[i]

        try:
            centerline = Centerline_for_single_outline (glacier_outline, clipped_dem, cellsize_int, radius)
            #arcpy.CopyFeatures_management(centerline, "d:\\temp\\centerline.shp")
            
            lineArr = arcpy.da.FeatureClassToNumPyArray(centerline, 'OID@')
            
            if len(lineArr) > 1:
                #arcpy.AddMessage("start connecting lines")
                connected_flowline = Connect_centerline (glacier_outline, centerline, clipped_dem, cellsize_int*3)
                ##smooth centerlines
                ##Add Merge ID to the flowline, so that the splited flowline can be merged back later
                #arcpy.AddMessage("Smooth!!!")
                ##First clip the connected flowlines using the glacier outlines
                ## and then delete the small segments created by clip
                ##Need to find the smallest line length first
                linearray = arcpy.da.FeatureClassToNumPyArray(connected_flowline, ('SHAPE@LENGTH'))
                length_list1 = np.array([item[0] for item in linearray])
                min_length = min(length_list1)
                line_count1 = len(length_list1)
                #arcpy.AddMessage(str(line_count1))
                #arcpy.AddMessage(str(min_length))

                arcpy.Clip_analysis(connected_flowline, glacier_outline, "in_memory\\clipped_flowlines")
                arcpy.MultipartToSinglepart_management("in_memory\\clipped_flowlines", "in_memory\\clipped_flowlines_singleParts")
                linearray = arcpy.da.FeatureClassToNumPyArray("in_memory\\clipped_flowlines_singleParts", ('SHAPE@LENGTH'))
                length_list2 = np.array([item[0] for item in linearray])
                line_count2 = len(length_list2)
                #arcpy.AddMessage(str(line_count2))
                
                if line_count2 > line_count1:
                    number_delete = line_count2 - line_count1
                    length_list2.sort() ##sort length from smallest to largest
                    del_length = length_list2[0:number_delete]
                    del_length = [ int(x) for x in del_length ]
                    #arcpy.AddMessage(del_length)

                    with arcpy.da.UpdateCursor("in_memory\\clipped_flowlines_singleParts", ['SHAPE@LENGTH']) as cursor:
                        for row in cursor:
                            if int(row[0]) in del_length: 
                                cursor.deleteRow()
                    del row, cursor
                else: ##Still need to delete small flowlines
                    #arcpy.AddMessage("Remove small segments")
                    with arcpy.da.UpdateCursor("in_memory\\clipped_flowlines_singleParts", ['SHAPE@LENGTH']) as cursor:
                        for row in cursor:
                            #arcpy.AddMessage(str(row[0]))
                            if row[0] < min(150,min_length): ##less than 100 m
                                cursor.deleteRow()

                arcpy.AddField_management("in_memory\\clipped_flowlines_singleParts", "MergeID", "Long", "", "", "", "", "", "", "")
                arcpy.arcpy.CalculateField_management("in_memory\\clipped_flowlines_singleParts","MergeID",str("!"+str(arcpy.Describe("in_memory\\clipped_flowlines_singleParts").OIDFieldName)+"!"),"PYTHON_9.3")
                split_flowline = line_split_by_intersect("in_memory\\clipped_flowlines_singleParts")

                #clean the standalone lines without connection with other lines
                cleandangleline(split_flowline, "in_memory\\clipped_flowlines_cleaned")
                #smooth lines
                #arcpy.PolygonToLine_management(glacier_outline, "in_memory\\glacier_outline_line")
                arcpy.cartography.SmoothLine("in_memory\\clipped_flowlines_cleaned", "in_memory\\smooth_flowline", "PAEK", smooth_dis) ##, "#", "#", "in_memory\\glacier_outline_line")
                arcpy.Dissolve_management("in_memory\\smooth_flowline", "in_memory\\smooth_flowline_dissove", "MergeID", "#", 'SINGLE_PART', '#')
                arcpy.DeleteField_management("in_memory\\smooth_flowline_dissove","MergeID")
            
            else:
                arcpy.AddMessage("just one line")
                arcpy.cartography.SmoothLine(centerline, "in_memory\\smooth_flowline_dissove", "PAEK", smooth_dis) ##, "#", "#", "in_memory\\glacier_outline_line")
                
            
            if Fc_Added == 0:
                arcpy.CopyFeatures_management("in_memory\\smooth_flowline_dissove", flow_lines)
            else:
                arcpy.Append_management("in_memory\\smooth_flowline_dissove", flow_lines, "NO_TEST")

            Fc_Added += 1
            
        except:
            arcpy.AddMessage("No centerline derived from this glacier outline")
            #arcpy.AddMessage("There is an error in the centerline delineation! move to the next one")
            pass

           
    arcpy.Delete_management("in_memory") ### Empty the in_memory    


#------------------------------------------------------------------------------------------------------------
# Start the main program
#------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # Script arguments
    glacier_outlines = arcpy.GetParameterAsText(0)
    dem = arcpy.GetParameterAsText(1)
    #StreamThresholdKM2 = float(arcpy.GetParameter(2))
    #TributaryThresholdKM2 = float(arcpy.GetParameter(3))
    flow_lines = arcpy.GetParameterAsText(2)

    Centerline_costPath(glacier_outlines, dem, flow_lines)
    '''
    arcpy.Delete_management("in_memory")

    ##make sure the projection of the glacier outlines is the same with the UTM and the same with the DEM projection 
    #spatial_flowline = arcpy.Describe(flowline).spatialReference
    spatial_outline = arcpy.Describe(glacier_outlines).spatialReference
    spatial_ref_dem = arcpy.Describe(dem).spatialReference

    if "UTM" in spatial_ref_dem.name:
        arcpy.AddMessage("The DEM projection is: " + spatial_ref_dem.name)
    else:
        arcpy.AddMessage("The DEM projection is not UTM. Please re-project the DEM to a UTM projection for the analysis!")
        exit()   

    if "UTM" in spatial_outline.name:
        arcpy.AddMessage("The outline projection is: " + spatial_outline.name)
    else:
        arcpy.AddMessage("The outline projection is not UTM. Please re-project the outline to a UTM projection for the analysis!")
        exit()   

    if spatial_outline.name != spatial_ref_dem.name:
        arcpy.AddMessage("The DEM and galcier outlines have different UTM projections. Please re-project them to the same UTM projection!")
        exit()   

    
    arcpy.env.snapRaster = dem

    ##Get the area of the outline
    FcID = arcpy.Describe(glacier_outlines).OIDFieldName
    polyarray = arcpy.da.FeatureClassToNumPyArray(glacier_outlines, ('SHAPE@AREA', FcID))
    area_list = np.array([item[0] for item in polyarray])
    FIds = np.array([item[1] for item in polyarray])

    radius_list = area_list * 2e-6 + 500
    smooth_list = area_list * 2e-6 + 200
    for i in range(len(radius_list)):
        if radius_list[i] > 1000:
            radius_list[i] = 1000
        if smooth_list[i] > 400:
            smooth_list[i] = 400

    ##Stream network analysis
    cellsize = arcpy.GetRasterProperties_management(dem,"CELLSIZEX")
    cellsize_int = int(cellsize.getOutput(0))


    glacier_outline = "in_memory\\glacier_outline" ##Single outline for each loop
    #flow_line_cp = "in_memory\\flow_line_cp"

    for i in range(len(FIds)):
        #try:
        arcpy.AddMessage("Generating centreline(s) for glacier "+str(i+1)+" of "+str(len(FIds)))
        query = FcID +" = "+str(FIds[i])

        arcpy.Select_analysis(glacier_outlines,glacier_outline,query)


        ##buf the sel_outline a little bit to clip the DEM
        arcpy.Buffer_analysis(glacier_outline, "in_memory\\outline_buf", (str(cellsize_int*3)+ " Meter"))

        
        clipped_dem = ExtractByMask(dem,"in_memory\\outline_buf")

        radius = radius_list[i]
        smooth_dis = smooth_list[i]

        centerline = Centerline_for_single_outline (glacier_outline, clipped_dem, cellsize_int, radius)
        #arcpy.CopyFeatures_management(centerline, "d:\\temp\\centerline.shp")
        
        lineArr = arcpy.da.FeatureClassToNumPyArray(centerline, 'OID@')
        if len(lineArr) > 1:
            #arcpy.AddMessage("start connecting lines")
            connected_flowline = Connect_centerline (glacier_outline, centerline, clipped_dem, cellsize_int*3)
            ##smooth centerlines
            ##Add Merge ID to the flowline, so that the splited flowline can be merged back later
            #arcpy.AddMessage("Smooth!!!")
            ##First clip the connected flowlines using the glacier outlines
            ## and then delete the small segments created by clip
            ##Need to find the smallest line length first
            linearray = arcpy.da.FeatureClassToNumPyArray(connected_flowline, ('SHAPE@LENGTH'))
            length_list1 = np.array([item[0] for item in linearray])
            min_length = min(length_list1)
            line_count1 = len(length_list1)
            #arcpy.AddMessage(str(line_count1))
            #arcpy.AddMessage(str(min_length))

            arcpy.Clip_analysis(connected_flowline, glacier_outline, "in_memory\\clipped_flowlines")
            arcpy.MultipartToSinglepart_management("in_memory\\clipped_flowlines", "in_memory\\clipped_flowlines_singleParts")
            linearray = arcpy.da.FeatureClassToNumPyArray("in_memory\\clipped_flowlines_singleParts", ('SHAPE@LENGTH'))
            length_list2 = np.array([item[0] for item in linearray])
            line_count2 = len(length_list2)
            #arcpy.AddMessage(str(line_count2))
            
            if line_count2 > line_count1:
                number_delete = line_count2 - line_count1
                length_list2.sort() ##sort length from smallest to largest
                del_length = length_list2[0:number_delete]
                del_length = [ int(x) for x in del_length ]
                #arcpy.AddMessage(del_length)

                with arcpy.da.UpdateCursor("in_memory\\clipped_flowlines_singleParts", ['SHAPE@LENGTH']) as cursor:
                    for row in cursor:
                        if int(row[0]) in del_length: 
                            cursor.deleteRow()
                del row, cursor
            else: ##Still need to delete small flowlines
                #arcpy.AddMessage("Remove small segments")
                with arcpy.da.UpdateCursor("in_memory\\clipped_flowlines_singleParts", ['SHAPE@LENGTH']) as cursor:
                    for row in cursor:
                        #arcpy.AddMessage(str(row[0]))
                        if row[0] < min(150,min_length): ##less than 100 m
                            cursor.deleteRow()

            arcpy.AddField_management("in_memory\\clipped_flowlines_singleParts", "MergeID", "Long", "", "", "", "", "", "", "")
            arcpy.arcpy.CalculateField_management("in_memory\\clipped_flowlines_singleParts","MergeID",str("!"+str(arcpy.Describe("in_memory\\clipped_flowlines_singleParts").OIDFieldName)+"!"),"PYTHON_9.3")
            split_flowline = line_split_by_intersect("in_memory\\clipped_flowlines_singleParts")

            #clean the standalone lines without connection with other lines
            cleandangleline(split_flowline, "in_memory\\clipped_flowlines_cleaned")
            #smooth lines
            #arcpy.PolygonToLine_management(glacier_outline, "in_memory\\glacier_outline_line")
            arcpy.cartography.SmoothLine("in_memory\\clipped_flowlines_cleaned", "in_memory\\smooth_flowline", "PAEK", smooth_dis) ##, "#", "#", "in_memory\\glacier_outline_line")
            arcpy.Dissolve_management("in_memory\\smooth_flowline", "in_memory\\smooth_flowline_dissove", "MergeID", "#", 'SINGLE_PART', '#')
            arcpy.DeleteField_management("in_memory\\smooth_flowline_dissove","MergeID")
        else:
            arcpy.AddMessage("just one line")
            arcpy.cartography.SmoothLine(centerline, "in_memory\\smooth_flowline_dissove", "PAEK", smooth_dis) ##, "#", "#", "in_memory\\glacier_outline_line")

        if i == 0:
            arcpy.CopyFeatures_management("in_memory\\smooth_flowline_dissove", flow_lines)
        else:
            arcpy.Append_management("in_memory\\smooth_flowline_dissove", flow_lines, "NO_TEST")
        #except:
        #    arcpy.AddMessage("There is an error in the centerline delineation! move to the next one")
        #    pass

           
    #arcpy.CopyFeatures_management(connected_flowline, flow_lines)

    ##Delete intermidiate data
    arcpy.Delete_management("in_memory") ### Empty the in_memory
    '''
