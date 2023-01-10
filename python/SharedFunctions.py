#-------------------------------------------------------------------------------
# Name: SharedFunctions.py
# Purpose: This python file include all standard alone functions used for palaeoglacier reconstruction
# The code is modified by GlaRe (Pellitero et al. 2016), VOLTA (James and Carrivick, 2016), and AutoCirque (Li and Zhao 2022)
# This file will be import for other python files associated with each tool
#
# Author: Dr. Yingkui Li
# Department of Geography, University of Tennessee
# Created:     09/01/2020
# Finalized:   06/24/2022 
# Copyright:   (c) Yingkui Li 2022
#-------------------------------------------------------------------------------
from __future__ import division
import arcpy
from arcpy import env
from arcpy.sa import *
import numpy as np
import math
import os, sys
import scipy
from scipy.spatial import cKDTree as KDTree
from scipy import ndimage
import arcpy.cartography as CA

arcpy.env.overwriteOutput = True
arcpy.env.XYTolerance= "0.01 Meters"

##Check the python version to determine ArcGIS or ArcGIS Pro
ArcGISPro = 0
arcpy.AddMessage("The current python version is: " + str(sys.version_info[0]))
if sys.version_info[0] == 2:  ##For ArcGIS 10, need to check the 3D and Spatial Extensions
    try:
        if arcpy.CheckExtension("Spatial")=="Available":
            arcpy.CheckOutExtension("Spatial")
        else:
            raise Exception ("not extension available")
    except:
        raise Exception ("unable to check out extension")

    try:
        if arcpy.CheckExtension("3D")=="Available":
            arcpy.CheckOutExtension("3D")
        else:
            raise Exception ("not extension available")
    except:
        raise Exception ("unable to check out extension")
elif sys.version_info[0] == 3:  ##For ArcGIS Pro
    ArcGISPro = 1
else:
    raise Exception("Must be using Python 2.x or 3.x")
    exit()   

#------------------------------------------------------------------------------------------------------------
# The following functions are for the transition of dictinary data type operations between different python versions
#------------------------------------------------------------------------------------------------------------
def dict_listvalues(d):
    if ArcGISPro == 1: ##for python 3, ArcGIS Pro
        return list(d.values())
    else: #for python 2; ArcMap
        return d.values()

def dict_listitems(d):
    if ArcGISPro == 1: ##for python 3, ArcGIS Pro
        return list(d.items())
    else: #for python 2; ArcMap
        return d.items()

def dict_itervalues(d):
    if ArcGISPro == 1: ##for python 3, ArcGIS Pro
        return iter(d.values())
    else: #for python 2; ArcMap
        return d.itervalues()

def dict_iteritems(d):
    if ArcGISPro == 1: ##for python 3, ArcGIS Pro
        return iter(d.items())
    else: #for python 2; ArcMap
        return d.iteritems()


#------------------------------------------------------------------------------------------------------------
# This function calculates the ice thickness of each point along the flowline based on the Excel flowline model
# introduced by the Benn and Houlton (2010). It is revised from the codes by Pellitero et al.(2016) in GlaRe.
#------------------------------------------------------------------------------------------------------------
def Ice_Thickness_Calculation(flowPnts):
    Array = arcpy.da.TableToNumPyArray(flowPnts,('RASTERVALU','OFID','NEAR_FID','POINT_X','POINT_Y',"SSTRESS", "ffactor","ice","thick"))
    TerminusDist= 0
    for i in range(len(Array)):
        if TerminusDist == 0:
            Array[i][7]=Array[i][0]
            Array[i][8]=Array [i][7]-Array[i][0]
            TerminusDist=1
        elif Array[i][1] == Array[i-1][1]:
            d=math.sqrt(math.pow((Array[i][3]-Array[i-1][3]),2)+math.pow((Array[i][4]-Array[i-1][4]),2))#euclidean distance
            b = -((Array[i-1][0]+Array[i][0]))
            c = Array[i-1][7]*(Array[i][0]-((Array[i-1][7]-Array[i-1][0])))-((2*d*(((Array[i][5]+Array[i-1][5])/2)/Array[i][6]))/(900*9.81))
            Array[i][7]=(-b+math.pow((math.pow(b,2))-(4*c),0.5))/2
            Array[i][8]=Array [i][7]-Array[i][0]
        else:
            Array[i][7]=Array[(Array[i][2]-1)][7]
            if 'nan' in str(Array[i][7]).lower():
                Array[i][7] = Array[i][0]
            if Array[i][8] != 0:
                Array[i][8]=Array [i][7]-Array[i][0]
            else:
                Array[i][7] = Array[i][0]  ##reset the ice surface to bedDEM

    #Save and update the new calculations to flowPnts
    count=0
    with arcpy.da.UpdateCursor(flowPnts,("ice","thick")) as cursor:
        for row in cursor:
            row[0]=Array[count][7]
            row[1]=Array[count][8]
            cursor.updateRow(row)
            count+=1
    #delete cursors
    del row, cursor

#------------------------------------------------------------------------------------------------------------
# This function check each line in the line feature and make sure the line is from low elevation to high
# elevation (Glacier flowline needs from low to high elevation in order to reconstruct the paleo ice thickness).
# It is revised from the codes by Pellitero et al.(2016) in GlaRe.
#------------------------------------------------------------------------------------------------------------
def Check_If_Flip_Line_Direction(line, dem):
    arcpy.AddField_management(line, "Flip", "Long", "", "", "", "", "", "", "")
    updCursor = arcpy.da.UpdateCursor(line, ["Shape@", "Flip"])
    nFlip = 0
    i = 0
    for curLine in updCursor:
        Startpoint2 = curLine[0].firstPoint
        Startpoint = arcpy.Geometry('Point',Startpoint2)
        coord= str(Startpoint2.X)+" "+str(Startpoint2.Y)
        Cellvalue=arcpy.GetCellValue_management(dem, coord)
        Startpoint.Z=Cellvalue.getOutput(0)
        Lastpoint2 = curLine[0].lastPoint
        Lastpoint = arcpy.Geometry('Point',Lastpoint2)
        coord= str(Lastpoint2.X)+" "+str(Lastpoint2.Y)
        Cellvalue=arcpy.GetCellValue_management(dem, coord)
        Lastpoint.Z=Cellvalue.getOutput(0)
        if Startpoint.Z >= Lastpoint.Z:  ##Flip = True use equal in case the start and end point are the same
            nFlip = nFlip + 1
            curLine[1] = 1
        else:  ##Flip = False
            curLine[1] = 0
        updCursor.updateRow(curLine)
        i += 1 
        
    if nFlip > 0:
        arcpy.MakeFeatureLayer_management(line, "lyrLines")
        arcpy.SelectLayerByAttribute_management("lyrLines", "NEW_SELECTION", '"Flip" > 0')
        #arcpy.AddMessage("Count of selected features is " + str(nFlip))
        arcpy.FlipLine_edit("lyrLines")  ##Need to change to lyrLines
        arcpy.SelectLayerByAttribute_management("lyrLines", "CLEAR_SELECTION")

    arcpy.DeleteField_management (line, "Flip")
    
    del updCursor
    if i>0:
        del curLine


#------------------------------------------------------------------------------------------------------------
# This function check the connectivity of the flowlines and the order of each flowline for plaeo ice reconstruction
# based on the elevation of it start point from the lowest to highest order.
# This function is revised from the codes by Pellitero et al.(2016) in GlaRe.
#------------------------------------------------------------------------------------------------------------
def Check_Flowline_Connectivity(flowline, Distance):
    arcpy.AddMessage ("Checking connectivity and order between flowlines...")
    points=[]
    lines=[]
    distances=[]
    bad=False
    with arcpy.da.SearchCursor(flowline, ["SHAPE@","OID@"]) as flows: ##This function needs to revise further to only check the flowline with the same GlacierID
        for flow in flows:
            points.append(((flow[0].firstPoint).X, (flow[0].firstPoint).Y))
            lines.append(flow[0])

        array=np.array(points)

        for i in range(len(array)):
            if i==0:
                pass
            else:
                warning=0
                punto=arcpy.Point(array[i][0],array[i][1])
                for a in range (len(array)):
                    #if a < i:
                    if a != i:
                        distances.append(math.sqrt(math.pow((array[i][0]-array[a][0]),2)+math.pow((array[i][1]-array[a][1]),2)))
                        touches = (punto.touches(lines[a]) or punto.within(lines[a]))
                        if touches==True:
                            warning=1
                            break
                        #else:
                        #    pass
                if warning==0:
                    bad=True
                    #arcpy.AddWarning( "The flowline number %r is likely to provoke the \"Flowline Ice Thickness\" tool to crash, please check increasing numbering and snap connectivity with its parental flowline" %(i))
                #else:
                #    pass

    del flow, flows
    
    #if bad==False:
    #    arcpy.AddMessage("The start points of flowlines are correctly numbered and snapped to their parental flowline")
    #else:
    #    arcpy.AddWarning( "The flowline number %r is likely to provoke the \"Flowline Ice Thickness\" tool to crash, please check increasing numbering and snap connectivity with its parental flowline" %(i))
        

    #arcpy.AddMessage ("Checking distance between flowlines first points")
    if Distance >= min(distances):
        bad = False
        #arcpy.AddWarning("This flowlines file is only guaranteed to work in the \"Flowline Ice Thickness\" tool with an interval distance lower than %r (in your local units)" %int(min(distances)))
    else:
        bad = True
        #arcpy.AddMessage("Flowline firstpoints are correctly spaced")

    return bad

#------------------------------------------------------------------------------------------------------------
# This function derives the F factor based on a set of XYZ points of a cross section. This function also include
# a boundary thickness parameter to account for the case that if the ice thickness is too small, it will not
# affect the shear stress of the flowline. 
# This function is revised from the codes by Pellitero et al.(2016) in GlaRe.
#------------------------------------------------------------------------------------------------------------
def Calculate_ffactor (pointx, pointy, pointz, bndthickness):
    distan3D=[]
    distan2D=[]
    avalue=[]
    width = 0
    if len(pointz) > 0:
        altmax= max(pointz) + bndthickness ### Add bndthickness to consider the ice boundary 
        altmin=min(pointz) 
        H = altmax-altmin 
    else:
        H = 0
    if len(pointx) > 1:
        for i in range(len(pointx)):
            if i!=0:
                distance2D=math.sqrt((math.pow(((pointx[i])-(pointx[i-1])),2)) +math.pow(((pointy[i])-(pointy[i-1])),2))
                distan2D.append(distance2D)
                distance3D=math.sqrt(math.pow((distance2D),2) + math.pow(((pointz[i])-(pointz[i-1])),2))
                distan3D.append(distance3D)
                aval=0.5 * (distance2D*((altmax-(pointz[i-1])+(altmax-(pointz[i])))))
                avalue.append(aval)
        width = sum(distan2D)
        P = sum(distan3D)
        A = sum(avalue)
    else:
        P = 0
        A = 0
    ##Calculate the F factor
    if P*H > 0:
        ffactor = A/(P*H)
    else:
        #ffactor = 1
        ffactor = 0.8
    if ffactor < 0.445:
        #arcpy.AddMessage("The ffactor is lower than 0.445, therefore it will be setted as 0.0445")
        #ffactor = 0.445
        ffactor = 0.8  ##Based on the Volta paper, the maximum F factor is 0.8 and when the value is < 0.445, use 0.8. This is just to check if the output is close to the vota thickness
        ##then, after the whole factors calcualted, using the average F factor to replace all 0.8.    
    return ffactor, width

#------------------------------------------------------------------------------------------------------------
# This function derives the F factor based on a set of XYZ points of a cross section. This function also include
# a boundary thickness parameter to account for the case that if the ice thickness is too small, it will not
# affect the shear stress of the flowline. The method to derive the F factor is based on the method proposed
# by Li et al (2012) for a polyfit of the assumed parabolic cross section.
#------------------------------------------------------------------------------------------------------------
def Calculate_ffactor_by_polyfit (pointx, pointy, pointz, bndthickness):
    distan2D=[]
    width = 0
    if len(pointz) > 2:
        altmax= max(pointz) + bndthickness ### Add bndthickness to consider the ice boundary
        altmin=min(pointz) 
        h = altmax-altmin
        cumdis = 0
        for i in range(len(pointz)):
            if i ==0:
                distan2D.append(0)
            else:
                dist=math.sqrt((math.pow(((pointx[i])-(pointx[i-1])),2)) +math.pow(((pointy[i])-(pointy[i-1])),2))
                cumdis += dist 
                distan2D.append(cumdis)

        #arcpy.AddMessage(distan2D)
        width = max(distan2D) ##get the width of the cross section
        polyfit = np.polyfit(distan2D,pointz, 2)
        a = polyfit[0]
        if (a > 0):
            r = math.sqrt(1.0/(h * a))
            ffactor = 0.9 * r / (1+ 0.9 * r)
        else:
            #ffactor = 1.0
            ffactor = 0.8 ##Based Volta the max value is 0.8; This is just for the comparision with Volta; need to check with the real dataset
    else:
        #ffactor = 1.0
        ffactor = 0.8

    if ffactor < 0.445:
        #arcpy.AddMessage("The ffactor is lower than 0.445, therefore it will be setted as 0.0445")
        #ffactor = 0.445
        ffactor = 0.8 ##Based Volta the max value is 0.8; This is just for the comparision with Volta; need to check with the real dataset

    return ffactor, width

#------------------------------------------------------------------------------------------------------------
# This function derives the shear stress of a glacier based on the primary flowline and ice surface elevations.
# To make sure to extract the valid elevations for the start and end points of the flowline, this tool applied a
# small buffer around these two points (2 cellsize) and extract the minimum and maximum elevations of the buffer.
#
# This function is revised based on the codes by James and Carrivick (2016) in VOLTA.
# Specifically, the contour step is changed to 300 m (1000ft) as proposed by the original paper (1986)
# and counted the missing part from the highest contour to the highest elevation of the flowline.
#------------------------------------------------------------------------------------------------------------
def shear_stress_calculation(mainflowline, outline, icedem, min_ss, max_ss): 
    #arcpy.AddMessage("Start shear stress calculation...")
    cellsize = arcpy.GetRasterProperties_management(icedem,"CELLSIZEX")
    cellsize_value = float(cellsize.getOutput(0)) ### should change to float ?????

    primary_flowline = arcpy.CopyFeatures_management(mainflowline, "in_memory\\primary_flowline")

    ##Get the start point and end point and make a small buffer around each point to make sure to get the elevation for each point
    arcpy.FeatureVerticesToPoints_management(primary_flowline, "in_memory\\startpnt", "START")
    arcpy.Buffer_analysis("in_memory\\startpnt", "in_memory\\startpntbuf", (str(cellsize_value*2)+ " Meter"))
    extractDEM = ExtractByMask(icedem, "in_memory\\startpntbuf")
    try:
        MinZresult = arcpy.GetRasterProperties_management(extractDEM,"MINIMUM")
    except:
        MinZresult = arcpy.GetRasterProperties_management(icedem,"MINIMUM")
    startz = float(MinZresult.getOutput(0))

    arcpy.FeatureVerticesToPoints_management(primary_flowline, "in_memory\\endpnt", "END")
    arcpy.Buffer_analysis("in_memory\\endpnt", "in_memory\\endpntbuf", (str(cellsize_value)+ " Meter"))
    extractDEM = ExtractByMask(icedem, "in_memory\\endpntbuf")
    try:
        MaxZresult = arcpy.GetRasterProperties_management(extractDEM,"MAXIMUM")
    except:
        MaxZresult = arcpy.GetRasterProperties_management(icedem,"MAXIMUM")
    endz = float(MaxZresult.getOutput(0))
    ##Get the flowline length
    with arcpy.da.SearchCursor(primary_flowline,["SHAPE@LENGTH"]) as cursor:
        for row in cursor:
            flowlinelength = row[0]
    del row, cursor

    if startz > endz:  ##This is not necessary if the direction of the primary flowline is checked
        arcpy.FlipLine_edit(primary_flowline)
        startz_original = startz
        endz_original = endz
        startz = endz_original ##reverse the value
        endz = startz_original

    Zdiff = float(endz - startz)
    
    arcpy.AddField_management(primary_flowline, "zero", "FLOAT")
    arcpy.AddField_management(primary_flowline, "max_dist", "FLOAT")
    #arcpy.CalculateField_management(primary_flowline,"zero",0)
    #arcpy.CalculateField_management(primary_flowline,"max_dist",flowlinelength)  ### All these can be simplified by cursor
    with arcpy.da.UpdateCursor(primary_flowline, ["zero","max_dist"]) as cursor:
        for row in cursor:
            row[0] = 0
            row[1] = flowlinelength
            cursor.updateRow(row)
    del row, cursor

    dem_glac_clip = "in_memory\\dem_glac_clip"
    arcpy.Clip_management(icedem, "#", dem_glac_clip, outline, "ClippingGeometry") ##Can also use Extract by mask tool
    dem_glac_clip_square = "in_memory\\dem_glac_clip_square"
    arcpy.Clip_management(icedem, "#", dem_glac_clip_square,outline) ##Extract the dem by the box of outline

    route_table = "in_memory\\route_table"
    intersect_point_ss = "in_memory\\intersect_point_ss"
    flowline_route = "in_memory\\flowline_route"
    contour_shear_stress = "in_memory\\contour_shear_stress"  
    contour_shear_stress_dissolve =  "in_memory\\contour_shear_stress_dissolve"   

    raster = arcpy.Raster(dem_glac_clip)
    array = arcpy.RasterToNumPyArray(raster,"","","",0)
    count = np.count_nonzero(array)

    if Zdiff >= 300.0: ### the orginal paper was 300; so change to 300 m
        contour_interval = 300.0
        arcpy.CreateRoutes_lr(primary_flowline, "zero", flowline_route, "TWO_FIELDS", "zero","max_dist")
        shear_list = []
        contour = startz
        oldsegmentlength = 0
        cum_segment_length = 0
        while contour < endz:
            start_contour = contour
            end_contour = contour + contour_interval
            area = (((( (start_contour < array) & (array < end_contour) ))).sum())*cellsize_value*cellsize_value
            segment_length = 0
            
            if end_contour < endz:
                #arcpy.AddMessage("Attention: end contour is less than endZ")
                ContourList(dem_glac_clip_square, contour_shear_stress, end_contour)
                arcpy.Dissolve_management(contour_shear_stress, contour_shear_stress_dissolve)
                arcpy.Intersect_analysis ([contour_shear_stress_dissolve, primary_flowline], intersect_point_ss, "","","POINT")
                if ArcGISPro == 0: ##for ArcMap
                    arcpy.LocateFeaturesAlongRoutes_lr(intersect_point_ss, flowline_route , "zero", 2, route_table)
                else: ##For ArcGIS Pro
                    try:
                        arcpy.LocateFeaturesAlongRoutes_lr(intersect_point_ss, flowline_route , "zero", 2, route_table, "RID; Point; MEAS")
                    except:
                        arcpy.LocateFeaturesAlongRoutes_lr(intersect_point_ss, flowline_route , "zero", 2, route_table)
                ##check if the route_table has value
                segment_length_list = []
                if (arcpy.management.GetCount(route_table)[0]) > "0":
                    with arcpy.da.SearchCursor(route_table, ["MEAS"]) as cursor:
                        for row in cursor: ##How about there are more than one record, sum or average???
                            segment_length_list.append(row[0])
                    cum_segment_length = max(segment_length_list)
                    segment_length = cum_segment_length - oldsegmentlength
                    oldsegmentlength = cum_segment_length
                else:
                    #arcpy.AddMessage("Attention: error!!!")
                    cum_segment_length = 0
            else: ##The last peice to the highest Z
                #arcpy.AddMessage("WARNING: end contour is larger than endZ!!!")
                segment_length = flowlinelength - cum_segment_length
            #arcpy.AddMessage("segment_length is:" +str(segment_length) )
            if segment_length > 0:
                angle_rad = math.atan(contour_interval/segment_length)
                area_cos = area/(math.cos(angle_rad))
                shear_list.append(area_cos)
                
            contour += contour_interval
        
        sum_a_cos = sum(shear_list)
        shear_stress = 27000* sum_a_cos**0.106  ### this equation was from Glacier volume estimation on Cascade volcanoes 1986 paper in Annual of Glaciology
        if shear_stress < min_ss:
            shear_stress = min_ss
            #arcpy.AddMessage("WARNING: The calculated shear stress is smaller than the specified minimum value, using the specified minimum value instead")
        if shear_stress > max_ss: 
            shear_stress = max_ss
            #arcpy.AddMessage("WARNING: The calculated shear stress is larger than the specified maximum value, using the specified maximum value instead")
            

    else:
        #arcpy.AddMessage("WARNING: glacier altitudinal extent < 300 m. using the default value instead")
        shear_stress = min_ss

    return shear_stress

#------------------------------------------------------------------------------------------------------------
# This fuction regroups flowlines to individual GlacierID, but not disslove the flowlines. It is assumed that
# each flowline is already connected from the top source the lowest points or the confluence points of another
# flowline. The flowline direction has to be from low to high for each flowline.
#------------------------------------------------------------------------------------------------------------
def Add_GlacierID_by_Touches (flowlines, field, outflowline): 

    arcpy.CopyFeatures_management(flowlines, outflowline)
    arcpy.AddField_management(outflowline, field, "Long", "", "", "", "", "", "", "") ##Add GID to the flowline layer
    arcpy.CalculateField_management(outflowline, field, "-1", '#', '#') ##Assign the default value as -1

    points=[]
    lines=[]
    glaciers = []
    GlacierID = []
    ##First loop to get the startpoints, lines and the glacierID list
    with arcpy.da.SearchCursor(outflowline, ["SHAPE@",field]) as flows:
        for flow in flows:
            points.append(((flow[0].firstPoint).X, (flow[0].firstPoint).Y))
            lines.append(flow[0])
            GlacierID.append(flow[1]) ##default are all -1 for each line
    del flow, flows
    
    array=np.array(points)

    ##Second loop: to get the lowest flowline and assign seperated glacier ID
    iceId = 0
    idlist = []
    for i in range(len(array)):
        nTouches=0
        punto=arcpy.Point(array[i][0],array[i][1])
        for a in range (len(array)):
            if a != i: ##don't touch itself
                touches=(punto.touches(lines[a]) or punto.within(lines[a]))
                if touches==True: ##if the start point touches others
                    #arcpy.AddMessage("touched")
                    nTouches = 1
                    break  ##break the loop
        if nTouches==0: ##this is for glaicer iD
            #arcpy.AddMessage("not touched")
            glaciers.append(lines[i])
            idlist.append(iceId)
            GlacierID[i] = iceId
            iceId += 1

    ##Loop 3: to test if eachline touches the glaciers
    leftover = len(array) - len(glaciers)
    while leftover > 0:
        for i in range(len(array)):
            nTouches=0
            punto=arcpy.Point(array[i][0],array[i][1])
            for a in range (len(glaciers)): 
                touches=(punto.touches(glaciers[a]) or punto.within(glaciers[a]))
                if touches==True: ##if the start point touches others
                    if GlacierID[i] == -1: ##if this line has not been assigned
                        GlacierID[i] = idlist[a]
                        glaciers.append(lines[i])
                        idlist.append(idlist[a])
            ##update leftover
            leftover = len(array) - len(glaciers)
            #arcpy.AddMessage("Leftover is: " + str(leftover))

    ##Finally add GlacierID to the outflowline
    i = 0
    with arcpy.da.UpdateCursor(outflowline, field) as cursor:
        for row in cursor:
            row[0] = GlacierID[i]
            cursor.updateRow(row)
            i += 1
    del row, cursor

    return outflowline


#------------------------------------------------------------------------------------------------------------
# This fuction creates perpendicular lines along the flowlines and then create a set of points along these
# perpendicular lines. This tool aslo extracts the elevation of each points. The max_width is used to limit
# the extent of these points. The division between the perp lines in different valleys are determined by the
# EU Allocation tool in ArcGIS.
#------------------------------------------------------------------------------------------------------------
def cross_section_points(flowlinepoints, flowline, watershed, beddem, max_width, resolution): 

    cellsize = arcpy.GetRasterProperties_management(beddem,"CELLSIZEX")
    cellsize_float = float(cellsize.getOutput(0)) # use float cell size
    spatialref=arcpy.Describe(flowlinepoints).spatialReference
    
    width = max_width/2.0  ##only use the half the width for each side of the flowline

    oldextent = arcpy.env.extent
    arcpy.Buffer_analysis(watershed, "in_memory\\watershedbuf", (str(cellsize_float)+ " Meter"))
    arcpy.env.extent = "in_memory\\watershedbuf"
    
    ##use the buffer so that the extract elevations of the points along the watershed boundary are not zero
    try:
        error = False
        extractDEM = ExtractByMask(beddem, "in_memory\\watershedbuf") 
    except:
        error = True

    ##Step 1: data preparation
    flowlinepointscp = arcpy.CopyFeatures_management(flowlinepoints, "in_memory\\flowlinepointscp")
    arcpy.AddField_management(flowlinepointscp, 'PntID', 'Long', 6) ##OFID is not FID, but it is related to FID
    arcpy.CalculateField_management(flowlinepointscp,"PntID",str("!"+str(arcpy.Describe(flowlinepointscp).OIDFieldName)+"!"),"PYTHON_9.3")

    arcpy.AddField_management(flowline, 'SegmentID', 'Long', 6) ##OFID is not FID, but it is related to FID
    arcpy.CalculateField_management(flowline,"SegmentID",str("!"+str(arcpy.Describe(flowline).OIDFieldName)+"!"),"PYTHON_9.3")
    flowlinepointscopy = "in_memory\\flowlinepointscopy"
    arcpy.SpatialJoin_analysis(flowlinepointscp, flowline, flowlinepointscopy, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")

    arr=arcpy.da.FeatureClassToNumPyArray(flowlinepointscopy, ('SHAPE@X', 'SHAPE@Y', 'PntID', 'OFID', 'SegmentID'))
    segment_ids = np.array([item[4] for item in arr])
    unique_segment_ids = np.unique(segment_ids)

    ##Step 2: create perpendicular lines along the flowpoints
    distance = 10000
        
    perpendiculars = arcpy.CreateFeatureclass_management("in_memory", "perpendiculars","POLYLINE", "","","", spatialref)
    arcpy.AddField_management(perpendiculars, 'FlowlineID', 'Long', 6) ##OFID is not FID, but it is related to FID
    arcpy.AddField_management(perpendiculars, 'SegmentID', 'Long', 6) ##OFID is not FID, but it is related to FID
    arcpy.AddField_management(perpendiculars, 'FlowPntID', 'Long', 6) ##OFID is not FID, but it is related to FID
    exclude_PntID = []
    new_line_cursor = arcpy.da.InsertCursor(perpendiculars, ('SHAPE@', 'FlowPntID' , 'FlowlineID', 'SegmentID'))
    for row in range(len(arr)):
        if row > 0: ##not the first point
            startx = arr[row][0]
            starty = arr[row][1]
            creat_perp = 1
            if row <(len(arr)-1):
                if arr[row][3]==arr[row+1][3] and arr[row][3]==arr[row-1][3]: ##inside points
                    endx = arr[row+1][0]
                    endy = arr[row+1][1]
                elif arr[row][3]==arr[row+1][3]: ##start tributary points
                    exclude_PntID.append(arr[row][2])
                    creat_perp = 0
                else: ##end of the tributary
                    endx = arr[row-1][0]
                    endy = arr[row-1][1]
            else:
                endx = arr[row-1][0]
                endy = arr[row-1][1]

            if creat_perp > 0:
                if starty==endy or startx==endx:
                    if starty == endy:
                        y1 = starty + distance
                        y2 = starty - distance
                        x1 = startx
                        x2 = startx
                    if startx == endx:
                        y1 = starty
                        y2 = starty 
                        x1 = startx + distance
                        x2 = startx - distance     
                else:
                    m = ((starty - endy)/(startx - endx)) #get the slope of the line
                    negativereciprocal = -1*((startx - endx)/(starty - endy))    #get the negative reciprocal
                    if m > 0:
                        if m >= 1:
                            y1 = negativereciprocal*(distance)+ starty
                            y2 = negativereciprocal*(-distance) + starty
                            x1 = startx + distance
                            x2 = startx - distance
                        if m < 1:
                            y1 = starty + distance
                            y2 = starty - distance
                            x1 = (distance/negativereciprocal) + startx
                            x2 = (-distance/negativereciprocal)+ startx           
                    if m < 0:
                        if m >= -1:
                            y1 = starty + distance
                            y2 = starty - distance
                            x1 = (distance/negativereciprocal) + startx
                            x2 = (-distance/negativereciprocal)+ startx     
                        if m < -1:
                            y1 = negativereciprocal*(distance)+ starty
                            y2 = negativereciprocal*(-distance) + starty
                            x1 = startx + distance
                            x2 = startx - distance
                array = arcpy.Array([arcpy.Point(x1,y1),arcpy.Point(x2, y2)])
                polyline = arcpy.Polyline(array)
                lineID = arr[row][3]
                pntID = arr[row][2]
                segID = arr[row][4]
                new_line_cursor.insertRow([polyline, pntID, lineID, segID])
                
    del new_line_cursor

    ##Step 3: Use EU Alloation to cut perp lines
    ##divide watershed by different flowlines
    ##Delete the exclude_points
    if len(exclude_PntID) > 0:
        with arcpy.da.UpdateCursor(flowlinepointscopy, 'PntID') as cursor:
            for row in cursor:
                if row in exclude_PntID:
                    cursor.deleteRow()
        del row, cursor

    # Process: Euclidean Allocation
    eucAllocate = EucAllocation(flowlinepointscopy, "", "", cellsize_float,'SegmentID', "", "")
    
    arcpy.RasterToPolygon_conversion(eucAllocate, "in_memory\\eucAllocate", "SIMPLIFY", "VALUE")
    arcpy.Clip_analysis ("in_memory\\eucAllocate", watershed, "in_memory\\clipAllocation")

    selAllocation = "in_memory\\selAllocation"
    sel_perp = "in_memory\\sel_perp"
    clipedlines = "in_memory\\clipedlines"

    for segment_id in unique_segment_ids:
        query = "gridcode = " + str(segment_id) 
        arcpy.Select_analysis ("in_memory\\clipAllocation", selAllocation, query)
        query = "SegmentID = "+ str(segment_id) 
        arcpy.Select_analysis (perpendiculars, sel_perp, query)

        arcpy.Clip_analysis (sel_perp, selAllocation, "in_memory\\clip_perp")

        if segment_id == unique_segment_ids[0]:
            arcpy.CopyFeatures_management("in_memory\\clip_perp", clipedlines)
        else:
            arcpy.Append_management("in_memory\\clip_perp", clipedlines, "NO_TEST" )
        
    
    ##Multipart to singleparts
    singlepartlines = arcpy.MultipartToSinglepart_management(clipedlines, "in_memory\\singlepartlines")

    ##Step 4: Clean perp lines, so that only the part intersect the flowpoints are preserved
    lines = []
    flowlineids = []
    flowpntids = []
    keep = []
    with arcpy.da.SearchCursor(singlepartlines, ['SHAPE@', 'FlowlineID', 'FlowPntID']) as perps:
        i = 0
        for perp in perps:
            lines.append(perp[0])
            flowlineids.append(perp[1]) 
            flowpntids.append(perp[2])
            keep.append(0)
            i += 1
    if i > 0:
        del perp
    del perps

    flowpntid_arr = np.array(flowpntids)

    with arcpy.da.SearchCursor(flowlinepointscopy, ['SHAPE@', 'OFID', 'OID@']) as pnts:
        for pnt in pnts:
            index = np.where(flowpntid_arr == pnt[2])
            indexlen = len(index[0])
            if indexlen > 0:
                if indexlen > 1: ##more than two lines for a point, only select the lines intersected with the point
                    for a in range (indexlen):
                        value = index[0][a]
                        within = pnt[0].within(lines[value])
                        if within==True: 
                            keep[value] = 1
                            break
                else:
                    value = index[0][0]
                    keep[value] = 1
    del pnts, pnt

    ##Delete the lines with keep = 0
    with arcpy.da.UpdateCursor(singlepartlines, 'SHAPE@') as perps:
        i = 0
        for perp in perps:
            if keep[i] == 0:
                perps.deleteRow()
            i += 1
    if i>0:
        del perp
    del perps

    ##Step 5: intersection of the perp lines
    ##Maybe this step is not necessary, instead, using a maximum width to limit the extremely long cross sections
    '''
    ##Solve the intersection issue of the perp lines
    splitted_perps = "in_memory\\splitted_perps"
    arcpy.PolygonToLine_management("in_memory\\clipAllocation", "in_memory\\clipAlloc_to_lines")
    arcpy.Append_management("in_memory\\clipAlloc_to_lines", singlepartlines, 'NO_TEST')
    arcpy.FeatureToPolygon_management(singlepartlines, "in_memory\\perppolys", "0.05 Meters")
    arcpy.PolygonToLine_management("in_memory\\perppolys", "in_memory\\perppolys_to_lines")
    arcpy.SplitLine_management("in_memory\\perppolys_to_lines", "in_memory\\splitline_perppolys")
    arcpy.SpatialJoin_analysis("in_memory\\splitline_perppolys", flowlinepointscopy, "in_memory\\final_perps", "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
    arcpy.SplitLineAtPoint_management("in_memory\\final_perps", flowlinepointscopy, splitted_perps, '1 Meters')

    #arcpy.CopyFeatures_management(splitted_perps, "c:\\test\\singlepartlines3.shp")
    '''    
    ##Step 6: Create a set of points along the perp lines on both side of the flowlinepointscopy use the steps of resolution; also add the maxmium width constrain
    ##the following three lines to replace step 5
    splitted_perps = "in_memory\\splitted_perps"
    arcpy.SpatialJoin_analysis(singlepartlines, flowlinepointscopy, "in_memory\\final_perps", "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
    arcpy.SplitLineAtPoint_management("in_memory\\final_perps", flowlinepointscopy, splitted_perps, '1 Meters')

    linearray = arcpy.da.FeatureClassToNumPyArray(splitted_perps, ('PntID', 'OFID'))  
    pnt_ids = np.array([item[0] for item in linearray])
    unique_pnt_ids = np.unique(pnt_ids)

    perp_points = arcpy.CreateFeatureclass_management("in_memory", "perp_points","POINT", "","","", spatialref)
    arcpy.AddField_management(perp_points, 'PntID', 'Long', 6) ##PntID related to the flowline point ID
    arcpy.AddField_management(perp_points, 'SortID', 'Long', 6) ## SortID: right from 0 to max, left from -min to -1, so that it can be used to derive the F factor
    arcpy.AddField_management(perp_points, 'PointZ', 'Double', 6, 2)
    
    sel_splitted_perps = "in_memory\\sel_splitted_perps"
    sel_flowline_pnt = "in_memory\\sel_flowline_pnt"
    perp_points_cursor = arcpy.da.InsertCursor(perp_points, ('SHAPE@', 'pntID', "SortID"))
    
    for pnt_id in unique_pnt_ids:
        query = "PntID = " + str(pnt_id)
        arcpy.Select_analysis (splitted_perps, sel_splitted_perps, query)
        arcpy.Select_analysis (flowlinepointscopy, sel_flowline_pnt, query)
        pntGeometry = arcpy.CopyFeatures_management(sel_flowline_pnt, arcpy.Geometry())
        pointx = pntGeometry[0].firstPoint.X
        pointy = pntGeometry[0].firstPoint.Y

        ##determine Left and right lines
        with arcpy.da.SearchCursor(sel_splitted_perps, ["SHAPE@","SHAPE@LENGTH"]) as cursor:
            for line in cursor:
                pntx = line[0].firstPoint.X
                pnty = line[0].firstPoint.Y
                length = line[1]
                effect_length = min(width,length)
                if (abs(pointx - pntx) < 0.01 and abs(pointy - pnty) < 0.01): ##start from the flowline point
                    start = 0
                    i = 0
                    while start < effect_length:
                        new_point = line[0].positionAlongLine(start)
                        perp_points_cursor.insertRow([new_point, pnt_id, i])
                        start = start + resolution
                        i += 1
                else:  #end at the flowline point
                    num_pnts = int(effect_length/resolution)
                    if length < width:
                        leftover = length - resolution * num_pnts
                    else:
                        leftover = length - width 
                    for i in range(num_pnts):
                        start = leftover + i * resolution
                        new_point = line[0].positionAlongLine(start)
                        perp_points_cursor.insertRow([new_point, pnt_id, (i-num_pnts)])
        del line, cursor            
    del perp_points_cursor

    ##Step 7: Extract Z values to points
    perp_points_with_Z = "in_memory\\perp_points_with_Z"
    if error == False:
        ExtractValuesToPoints(perp_points, extractDEM, perp_points_with_Z, "INTERPOLATE", "VALUE_ONLY")
    else:
        ExtractValuesToPoints(perp_points, beddem, perp_points_with_Z, "INTERPOLATE", "VALUE_ONLY")

    arcpy.CalculateField_management(perp_points_with_Z, 'PointZ', str("!"+"RASTERVALU"+"!"), 'PYTHON_9.3')
    arcpy.DeleteField_management(perp_points_with_Z, 'RASTERVALU')

    arcpy.env.extent = oldextent

    return perp_points_with_Z

#------------------------------------------------------------------------------------------------------------
# This function create a small rectangle polygon perpendicular to the flowline. The small polygon will be used
# to determine the outlet of the watershed, so that the watershed boundary will be ended at the start/end points
# of the flowline
#------------------------------------------------------------------------------------------------------------
def Startpoint_perpendiculars(line, step, distance):
    perpendicular_poly = arcpy.CreateFeatureclass_management("in_memory", "perps","POLYGON","","","",line)
    new_poly_cursor = arcpy.da.InsertCursor(perpendicular_poly, ('SHAPE@'))
    with arcpy.da.SearchCursor(line, "SHAPE@") as cursor:
        for row in cursor:
            firstPnt = row[0].firstPoint
            startx = firstPnt.X
            starty = firstPnt.Y

            endPnt = row[0].positionAlongLine(step).getPart()
            endx = endPnt.X
            endy = endPnt.Y

            if starty==endy or startx==endx:
                if starty == endy:
                    sy1 = starty + distance
                    sy2 = starty - distance
                    sx1 = startx
                    sx2 = startx
                    ey1 = endy + distance
                    ey2 = endy - distance
                    ex1 = endx
                    ex2 = endx
                if startx == endx:
                    sy1 = starty
                    sy2 = starty 
                    sx1 = startx + distance
                    sx2 = startx - distance     
                    ey1 = endy
                    ey2 = endy 
                    ex1 = endx + distance
                    ex2 = endx - distance     
            else:
                m = ((starty - endy)/(startx - endx)) #get the slope of the line
                negativereciprocal = -1*((startx - endx)/(starty - endy))    #get the negative reciprocal
                if m > 0:
                    if m >= 1:
                        sy1 = negativereciprocal*(distance)+ starty
                        sy2 = negativereciprocal*(-distance) + starty
                        sx1 = startx + distance
                        sx2 = startx - distance
                        ey1 = negativereciprocal*(distance)+ endy
                        ey2 = negativereciprocal*(-distance) + endy
                        ex1 = endx + distance
                        ex2 = endx - distance
                    if m < 1:
                        sy1 = starty + distance
                        sy2 = starty - distance
                        sx1 = (distance/negativereciprocal) + startx
                        sx2 = (-distance/negativereciprocal)+ startx           
                        ey1 = endy + distance
                        ey2 = endy - distance
                        ex1 = (distance/negativereciprocal) + endx
                        ex2 = (-distance/negativereciprocal)+ endx           
                if m < 0:
                    if m >= -1:
                        sy1 = starty + distance
                        sy2 = starty - distance
                        sx1 = (distance/negativereciprocal) + startx
                        sx2 = (-distance/negativereciprocal)+ startx     
                        ey1 = endy + distance
                        ey2 = endy - distance
                        ex1 = (distance/negativereciprocal) + endx
                        ex2 = (-distance/negativereciprocal)+ endx     
                    if m < -1:
                        sy1 = negativereciprocal*(distance)+ starty
                        sy2 = negativereciprocal*(-distance) + starty
                        sx1 = startx + distance
                        sx2 = startx - distance
                        ey1 = negativereciprocal*(distance)+ endy
                        ey2 = negativereciprocal*(-distance) + endy
                        ex1 = endx + distance
                        ex2 = endx - distance

            array = arcpy.Array([arcpy.Point(sx1,sy1),arcpy.Point(sx2, sy2),arcpy.Point(ex2, ey2), arcpy.Point(ex1, ey1)])
            polygon = arcpy.Polygon(array)
            new_poly_cursor.insertRow([polygon])

    del row, new_poly_cursor, cursor  

    return perpendicular_poly


#------------------------------------------------------------------------------------------------------------
# This function calculates ice surface elevation for the points along the flowlines and assign the ice surface elevation
# of the these points to the closest cross section points. Then, the ice surface elevations of the the cross section
# points will be used to interpret the ice surface raster based on the Topo to Raster tool and determine the boundary
# of the reconstructed glaicer.
#------------------------------------------------------------------------------------------------------------
def IceSurfaceCalculation_with_crosssection_pnts (flowline_points, beddem, watershed, cross_section_pnts, bFeatureComparison, tagetFc, field, outicesurface, outicepoly):
    arcpy.env.extent = beddem
    arcpy.env.cellSize = beddem
    
    ##Use the cross section points that already created to interpret the ice surfaces
    clip_CS_points = "in_memory\\clip_CS_points"
    arcpy.Clip_analysis (cross_section_pnts, watershed, clip_CS_points)
    arcpy.Near_analysis (clip_CS_points, flowline_points)#identify nearest flowline point
    
    pointarray = arcpy.da.FeatureClassToNumPyArray(flowline_points,["OID@",field])#create array of ice values and populate a list with them
    fidval = np.array([item[0] for item in pointarray])
    iceval = np.array([item[1] for item in pointarray])

    ###Add field into the cross section points
    exist_fields = [f.name for f in arcpy.ListFields(clip_CS_points)] #List of current field names in outline layer
    if field not in exist_fields:
        arcpy.AddField_management(clip_CS_points, field, "FLOAT",10,6) #field for ice value

    with arcpy.da.UpdateCursor(clip_CS_points,("NEAR_FID", field)) as cursor:   #populate ice field with value from the nearest flowline point
        i = 0
        for row in cursor:
            fid = row[0]
            idx_result = np.where(fidval == fid)
            idx = idx_result[0][0]
            row[1]=iceval[idx]
            cursor.updateRow(row)
            i += 1
    if i>0:   
        del row
    del cursor

    if bFeatureComparison:
        #convert target feature to points
        targetpoints = "in_memory\\targetpoints"
        arcpy.FeatureVerticesToPoints_management(tagetFc, targetpoints, 'ALL')
        targetpoints3D = ExtractValuesToPoints(targetpoints, beddem, "in_memory\\targetpoints3D", "INTERPOLATE")
        inPointElevations = TopoPointElevation([[clip_CS_points,field],[targetpoints3D, 'RASTERVALU'] ])
    else:
        inPointElevations = TopoPointElevation([[clip_CS_points,field]])

    palaeoShape= TopoToRaster([inPointElevations])

    dod = palaeoShape - beddem
    Palaeo1 = Con(dod > 0, 1)

    Palaeoice=ExtractByMask(Palaeo1,watershed)
    PaleoGlacierPoly = "in_memory\\palepice" ##For some reasons, the inmemory does not work here!!
    arcpy.RasterToPolygon_conversion(Palaeoice, PaleoGlacierPoly, "NO_SIMPLIFY", "VALUE")

    ##Use spatialjoin to extract the polygon corresponding to the flowline points, this step remove the small polygons that not intersected with the flowline points
    arcpy.SpatialJoin_analysis(PaleoGlacierPoly, flowline_points, outicepoly, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "0.1 Meters", "#")

    ##Check if outicepoly area is zero
    if (arcpy.management.GetCount(outicepoly)[0]) > "0":
        paleoicesurface = ExtractByMask(palaeoShape,outicepoly)
    else:
        paleoicesurface = ExtractByMask(palaeoShape,watershed)
        arcpy.CopyFeatures_management(watershed, outicepoly)

    paleoicesurface.save(outicesurface)

    return outicesurface, outicepoly

#------------------------------------------------------------------------------------------------------------
# This fuction adjusts the shape (F) factor based on the cross section points. The F factor is determined by the
# cross section area, perimeter and the maximum height of the cross section. 
#------------------------------------------------------------------------------------------------------------
def AdjustFfactor_with_cross_section_pnts (flowlinepoints, cross_section_pnts, icefield, icepoly):
    #clip cross section points by icepoly
    clip_CS_points = "in_memory\\clip_CS_points"
    arcpy.Clip_analysis (cross_section_pnts, icepoly, clip_CS_points)

    ##add the ice surface elevation into the cross section points
    arcpy.Near_analysis (clip_CS_points, flowlinepoints)#identify nearest flowline point
    
    pointarray = arcpy.da.FeatureClassToNumPyArray(flowlinepoints,["OID@",icefield])#create array of ice values and populate a list with them
    fidval = np.array([item[0] for item in pointarray])
    iceval = np.array([item[1] for item in pointarray])

    ###Add field into the cross section points
    exist_fields = [f.name for f in arcpy.ListFields(clip_CS_points)] #List of current field names in outline layer
    if icefield not in exist_fields:
        arcpy.AddField_management(clip_CS_points, icefield, "FLOAT",10,6) #field for ice value

    with arcpy.da.UpdateCursor(clip_CS_points,("NEAR_FID", icefield)) as cursor:   #populate ice field with value from the nearest flowline point
        i = 0
        for row in cursor:
            fid = row[0]
            idx_result = np.where(fidval == fid)
            idx = idx_result[0][0]
            row[1]=iceval[idx]
            cursor.updateRow(row)
            i += 1
    if i>0:   
        del row
    del cursor


    ffactorlist = []
    CSpntarray = arcpy.da.FeatureClassToNumPyArray(clip_CS_points, ('PntID', 'SortID', 'SHAPE@X', 'SHAPE@Y', 'PointZ', icefield))
    pnt_ids = np.array([item[0] for item in CSpntarray])
    unique_pnt_ids = np.unique(pnt_ids)

    for pntid in unique_pnt_ids:
        selected_array = CSpntarray[pnt_ids == pntid]
        sort_ids = np.array([item[1] for item in selected_array])

        i = np.argsort(sort_ids)
        pointx = np.array([item[2] for item in selected_array])
        pointx = pointx[i[::-1]].tolist()

        pointy = np.array([item[3] for item in selected_array])
        pointy = pointy[i[::-1]].tolist()

        pointz = np.array([item[4] for item in selected_array])
        pointz = pointz[i[::-1]].tolist()

        ffactor, width = Calculate_ffactor (pointx, pointy, pointz, 0) ##Use a thickness , if icethickness less than this value will not contribute to the the F factor calculation
        
        #ffactor, width = Calculate_ffactor (pointx, pointy, pointz, 0)
        ffactorlist.append(ffactor)

    ##update ffactors in the flowlinepoints based on the treatment of Volta
    mean_ffactor = np.mean(ffactorlist)
    if mean_ffactor < 0 or mean_ffactor > 0.95:
        mean_ffactor = 0.8
    
    with arcpy.da.UpdateCursor(flowlinepoints,("OID@","ffactor")) as curs:
        for row in curs:
            pid = row[0]
            if pid in unique_pnt_ids:
                idx_result = np.where(unique_pnt_ids == pid)
                idx = idx_result[0][0]
                factor = ffactorlist[idx]
                if factor == 0.8:
                    row[1] = mean_ffactor
                else:
                    row[1] = factor
            curs.updateRow(row)

    del row, curs

    return

#------------------------------------------------------------------------------------------------------------
# This fuction adjusts the F factor based on the cross section points. The F factor is determined by the
# polyfit og the cross section based on Li et al (2012). 
#------------------------------------------------------------------------------------------------------------
def AdjustFfactor_Ployfit_with_cross_section_pnts (flowlinepoints, cross_section_pnts,icefield, icepoly):
    #clip cross section points by icepoly
    clip_CS_points = "in_memory\\clip_CS_points"
    arcpy.Clip_analysis (cross_section_pnts, icepoly, clip_CS_points)

    ##add the ice surface elevation into the cross section points
    arcpy.Near_analysis (clip_CS_points, flowlinepoints)#identify nearest flowline point
    
    pointarray = arcpy.da.FeatureClassToNumPyArray(flowlinepoints,["OID@",icefield])#create array of ice values and populate a list with them
    fidval = np.array([item[0] for item in pointarray])
    iceval = np.array([item[1] for item in pointarray])

    ###Add field into the cross section points
    exist_fields = [f.name for f in arcpy.ListFields(clip_CS_points)] #List of current field names in outline layer
    if icefield not in exist_fields:
        arcpy.AddField_management(clip_CS_points, icefield, "FLOAT",10,6) #field for ice value

    with arcpy.da.UpdateCursor(clip_CS_points,("NEAR_FID", icefield)) as cursor:   #populate ice field with value from the nearest flowline point
        i = 0
        for row in cursor:
            fid = row[0]
            idx_result = np.where(fidval == fid)
            idx = idx_result[0][0]
            row[1]=iceval[idx]
            cursor.updateRow(row)
            i += 1
    if i>0:   
        del row
    del cursor
    

    ffactorlist = []

    CSpntarray = arcpy.da.FeatureClassToNumPyArray(clip_CS_points, ('PntID', 'SortID', 'SHAPE@X', 'SHAPE@Y', 'PointZ', icefield))
    pnt_ids = np.array([item[0] for item in CSpntarray])
    unique_pnt_ids = np.unique(pnt_ids)

    for pntid in unique_pnt_ids:
        selected_array = CSpntarray[pnt_ids == pntid]
        sort_ids = np.array([item[1] for item in selected_array])

        i = np.argsort(sort_ids)
        pointx = np.array([item[2] for item in selected_array])
        pointx = pointx[i[::-1]].tolist()

        pointy = np.array([item[3] for item in selected_array])
        pointy = pointy[i[::-1]].tolist()

        pointz = np.array([item[4] for item in selected_array])
        pointz = pointz[i[::-1]].tolist()


        ffactor, width = Calculate_ffactor_by_polyfit (pointx, pointz, pointz, 0)
        ffactorlist.append(ffactor)

    ##update ffactors in the flowlinepoints based on the treatment of Volta
    mean_ffactor = np.mean(ffactorlist)
    if mean_ffactor < 0 or mean_ffactor > 0.95:
        mean_ffactor = 0.8
    
    with arcpy.da.UpdateCursor(flowlinepoints,("OID@","ffactor")) as curs:
        for row in curs:
            pid = row[0]
            if pid in unique_pnt_ids:
                idx_result = np.where(unique_pnt_ids == pid)
                idx = idx_result[0][0]
                factor = ffactorlist[idx]
                if factor == 0.8:
                    row[1] = mean_ffactor
                else:
                    row[1] = factor
            curs.updateRow(row)

    del row, curs

    return

#------------------------------------------------------------------------------------------------------------
# This fuction calculates the distance from the reconstructed ice boundary to target linear features. The function
# is based on Li et al (2008) and only determine the mean distance. This function aslo derives a percentage of the
# target features outside the ice boundary and use it to determine the direction of the distance (more inside or
# outside of the ice boundary)
#------------------------------------------------------------------------------------------------------------
def Distance_to_Target_Features (icepoly, TargetFeatures, field, watershed, cellsize):
    ##Convert icepoly to line features
    arcpy.PolygonToLine_management(icepoly, "in_memory\\icepolyline")
    oldextent = arcpy.env.extent
    arcpy.env.extent = watershed
    outEucDistance = EucDistance("in_memory\\icepolyline", "#", cellsize, "#")
    ##Zonal statistics to table
    outZSaT = ZonalStatisticsAsTable(TargetFeatures, field, outEucDistance, "in_memory\\zonalmean", "#", "ALL")
    ##Get the mean distance value from the table
    with arcpy.da.SearchCursor(outZSaT, ["MEAN", "SUM"]) as cursor:
        for row in cursor:
            mean_dist = row[0]
            sum_dist = row[1]
    del row, cursor

    ##just the distance outside of the icepoly
    outEucDistance = EucDistance(icepoly, "#", cellsize, "#")
    outZSaT = ZonalStatisticsAsTable(TargetFeatures, field, outEucDistance, "in_memory\\zonalsum", "#", "SUM")
    ##Get the sum distance outside of the ice poly from the table
    with arcpy.da.SearchCursor(outZSaT, ["SUM"]) as cursor:
        for row in cursor:
            outside_sum_dist = row[0]
    del row, cursor

    outside_percent = min(1.0, outside_sum_dist / sum_dist)

    arcpy.env.extent = oldextent

    return mean_dist, outside_percent ##Return the average distance and outside percentage of the target features
    
#------------------------------------------------------------------------------------------------------------
# This fuction smooths the flowlines based on the size of the watershed (flow accumulation)
#------------------------------------------------------------------------------------------------------------
def lineSmooth(inline, outline, smoothfield, cellsize):
    LineID = arcpy.Describe(inline).OIDFieldName
    countResult = arcpy.GetCount_management(inline)
    count = int(countResult.getOutput(0))
    cellarea = float(cellsize) * float(cellsize)
    for i in range (count):
        query = LineID +" = "+ str(i+1)  ##THe FID starts from 1 for the in_memory data; Which is different from shape (FID from 0) and feature class in the geodatabase (objectID from 1)
        #arcpy.AddMessage(query)
        arcpy.Select_analysis(inline, "in_memory\\linesection", query)
        #Determine the smooth tolerance
        areacount = 0
        with arcpy.da.SearchCursor("in_memory\\linesection", smoothfield) as cursor:
            for row in cursor:
                areacount = row[0]
        del cursor, row

        tolerance = int(areacount * cellarea * 2 / 1.0e6) + 200 ##The function provided by Kienholz et al. (2014) and  James & Carrivick (2016)
        if tolerance > 1000:
            tolerance = 1000
        #arcpy.AddMessage(str(tolerance))
        arcpy.cartography.SmoothLine("in_memory\\linesection", "in_memory\\linesectsmooth", "PAEK", tolerance)

        if i == 0: ##The first loop
            arcpy.CopyFeatures_management("in_memory\\linesectsmooth", outline)
        else:
            arcpy.Append_management("in_memory\\linesectsmooth", outline, "NO_TEST")    
    return outline

#------------------------------------------------------------------------------------------------------------
# The function cleans extrlines based on from and to nodes. If only one to node and no corresponding from node, 
# except for the highest facc section, marking for deletion. The same processes are iterated to remove all extra 
# lines. This is more efficient than other clean line methods based on the intersect of the to node points.
#------------------------------------------------------------------------------------------------------------
def cleanextralineswithtopology(inline,outline, field):
    bflag = 1
    while bflag:
        bflag = 0
        lineArray = arcpy.da.FeatureClassToNumPyArray(inline,['OID@','from_node','to_node', field])
        fromnode = np.array([item[1] for item in lineArray])
        tonode = np.array([item[2] for item in lineArray])
        facc = np.array([item[3] for item in lineArray])
        uniquetonode = np.unique(tonode)
        maxfacc = max(facc)

        lineid = [] ##Record the id for deletion
        for i in range(len(lineArray)):
            linetonode = lineArray[i][2]
            nodecount = np.count_nonzero(uniquetonode == linetonode)
            if nodecount == 1 and not (linetonode in fromnode) and lineArray[i][3] < maxfacc: ###only one tonode except for the highest facc section
                #print "mark for deletion"
                lineid.append(lineArray[i][0])
                bflag = 1

        ##Delete the line marked for deletion
        with arcpy.da.UpdateCursor(inline, "OID@") as cursor:
            for row in cursor:
                if int(row[0]) in lineid:
                    cursor.deleteRow()     
        del cursor, row

    arcpy.CopyFeatures_management(inline, outline)
    return outline

#------------------------------------------------------------------------------------------------------------
# This fuction regroups flowlines to individual GlacierID and dissolve the flowline sections from the top 
# to the lowest points or the confluence points of another flowline. The flowline direction has to be from 
# low to high for each flowline.
#------------------------------------------------------------------------------------------------------------
def Merge_and_Add_GlacierID_by_Topology (flowlines, FaccField, GlacierID, MergeID, outflowline): 

    flowlinecopy = "in_memory\\flowlinecopy"
    arcpy.CopyFeatures_management(flowlines, flowlinecopy)
    arcpy.AddField_management(flowlinecopy, GlacierID, "Long") #Add GID to the flowline layer
    arcpy.AddField_management(flowlinecopy, MergeID, "Long") #Add MergeID to the flowline layer
    #Calculate Field does not work in ArcGIS Pro, using the following code instead
    with arcpy.da.UpdateCursor(flowlinecopy, [GlacierID, MergeID]) as cursor:
        for row in cursor:
            row[0] = -1
            row[1] = -1
            cursor.updateRow(row)
    del row, cursor

    points=[]
    lines=[]
    facc = []
    glaciers = []
    glacierID = []
    mergeid = []
    mergeused = []
    ##First loop to get the startpoints, lines and the glacierID list
    with arcpy.da.SearchCursor(flowlinecopy, ["SHAPE@", FaccField, GlacierID, MergeID]) as flows:
        for flow in flows:
            points.append(((flow[0].firstPoint).X, (flow[0].firstPoint).Y))
            lines.append(flow[0])
            facc.append(flow[1])
            glacierID.append(flow[2]) ##default are all -1 for each line
            mergeid.append(flow[3])
            mergeused.append(0)
    del flow, flows
    
    lenarr = np.array(facc)
    ids = lenarr.argsort()[::-1]

    array=np.array(points)
    ##Second loop: to get the lowest flowline and assign seperated glacier ID
    iceId = 0
    idlist = []
    mergidlist = []
    for i in range(len(array)):
        nTouches=0
        punto=arcpy.Point(array[i][0],array[i][1])
        for a in range (len(array)):
            if a != i: ##don't touch itself
                touches=(punto.touches(lines[a]) or punto.within(lines[a]))
                if touches==True: ##if the start point touches others
                    #arcpy.AddMessage("touched")
                    nTouches = 1
                    break  ##break the loop
        if nTouches==0: ##this is for glaicer iD
            #arcpy.AddMessage("not touched")
            glaciers.append(lines[i])
            idlist.append(iceId)
            glacierID[i] = iceId
            mergeid[i] = iceId
            mergidlist.append(iceId)
            iceId += 1

    ##Loop 3: to test if eachline touches the glaciers
    leftover = len(array) - len(glaciers)
    while leftover > 0:
        #arcpy.AddMessage("leftover: " + str(leftover))
        for i in range(len(array)):
            nTouches=0
            lineid = ids[i]
            if mergeid[lineid] == -1:
                punto=arcpy.Point(array[lineid][0],array[lineid][1])
                for a in range (len(glaciers)): 
                    touches = punto.touches(glaciers[a])
                    if touches==True: ##if the start point touches others
                        if glacierID[lineid] == -1: ##if this line has not been assigned
                            glacierID[lineid] = idlist[a]
                            glaciers.append(lines[lineid])
                            idlist.append(idlist[a])
                            if (mergeused[a] == 0):
                                mergeid[lineid] = mergidlist[a]
                                mergidlist.append(mergidlist[a])
                                mergeused[a] = 1
                            else:
                                ##Start a new mergeID
                                maxmergid = max(mergidlist)
                                mergeid[lineid] = maxmergid + 1
                                mergidlist.append(maxmergid + 1)
                    else:
                        within = punto.within(glaciers[a])
                        if within==True: ##if the start point touches others
                            if glacierID[lineid] == -1: ##if this line has not been assigned
                                glacierID[lineid] = idlist[a]
                                glaciers.append(lines[lineid])
                                idlist.append(idlist[a])
                                ##start a new mergeid with the max mergid + 1
                                maxmergid = max(mergidlist)
                                mergeid[lineid] = maxmergid + 1
                                mergidlist.append(maxmergid + 1)
                        
                ##update leftover
            leftover = len(array) - len(glaciers)

    ##Finally add GlacierID to the outflowline
    i = 0
    with arcpy.da.UpdateCursor(flowlinecopy, [GlacierID, MergeID]) as cursor:
        for row in cursor:
            row[0] = glacierID[i]
            row[1] = mergeid[i]
            cursor.updateRow(row)
            i += 1
    del row, cursor

    ##Disslove the flowline
    #field_treatment = "'" + FaccField +" SUM; " + GlacierID +" FIRST'"
    field_treatment = FaccField +" SUM; " + GlacierID +" FIRST"
    arcpy.Dissolve_management(flowlinecopy, outflowline, MergeID, field_treatment, 'SINGLE_PART', 'UNSPLIT_LINES')

    exist_fields = [f.name for f in arcpy.ListFields(outflowline)] #List of current field names in outline layer
    #arcpy.AddMessage(exist_fields)
    for field in exist_fields:
        if field[0:6] == "FIRST_":
            FirstGID = field
        elif field[0:4] == "SUM_":
            MaxFcc = field
    
    new_fields = [FaccField, GlacierID]
    for field in new_fields:
        if field not in exist_fields:
            arcpy.AddField_management(outflowline, field, "LONG")
    
    fields = [FaccField, GlacierID, MaxFcc, FirstGID]
    #arcpy.AddMessage(fields)
    with arcpy.da.UpdateCursor(outflowline, fields) as cursor:
        for row in cursor:
            row[0] = row[2]
            row[1] = row[3]
            cursor.updateRow(row)
    del row, cursor

    arcpy.DeleteField_management(outflowline,[FirstGID, MaxFcc, MergeID])

    return outflowline

#------------------------------------------------------------------------------------------------------------
# This function remove the flowline within the boundary of the modern glaciers based on a buffer distance.
#------------------------------------------------------------------------------------------------------------
def Remove_Flowline_in_Polygon_by_Erase(inline, polygons, buffer_dis):
    #buffer polygon a little bit, so that the centerline can be extended smoothly and erase the flowline in the polygon completely
    Str_buffer_dis = str(buffer_dis) + " Meters"
    polybuf = "in_memory\\polybuf"
    eraseoutline = "in_memory\\eraseoutline"
    outline = "in_memory\\outline"
    ##Determine the minimum length of the inline
    bEmpty = False
    minlength = 100000
    with arcpy.da.SearchCursor(inline, 'SHAPE@LENGTH') as cursor:
        for row in cursor:
            if row[0] < minlength:
                minlength = row[0]
    del cursor, row
    minlength = min(minlength, buffer_dis) ##Keep the minlength less than the buffer_distance
    #arcpy.AddMessage (str(minlength))

    arcpy.Buffer_analysis(polygons, polybuf, Str_buffer_dis)
    arcpy.Erase_analysis(inline, polybuf, eraseoutline)
    ##Need to get rid of small line sections and clean the folowline## MAX_MAX is the maximum flow accumulation value created by the flowline code
    arcpy.MultipartToSinglepart_management(eraseoutline, outline)
    ##Any line with the length smaller than the minlength is created by the erase and should be removed for the cleanning propose.
    with arcpy.da.UpdateCursor(outline, 'SHAPE@LENGTH') as cursor:
        i = 0
        for row in cursor:
            if row[0] < minlength:
                cursor.deleteRow()
            i += 1 
    del cursor

    if i>0:
        del row
    else:
        bEmpty = True
        return outline, bEmpty 

    ##Remove the lines with start point touch/within the polybuf
    remove = []
    lines = []
    points = []
    ##First loop to get the startpoints, lines and the glacierID list
    with arcpy.da.SearchCursor(outline, ["SHAPE@", "ORIG_FID"]) as flows:
        i = 0
        for flow in flows:
            points.append(((flow[0].firstPoint).X, (flow[0].firstPoint).Y))##Start points
            remove.append(0)
            i += 1
    del flows
    if i>0:
        del flow
    else:
        bEmpty = True
        return outline, bEmpty 
        
    ##second: loop to get polygonbuf outlines
    single_polybuf= arcpy.MultipartToSinglepart_management(polybuf, "in_memory\\single_polybuf")
    polybuf_outline = arcpy.PolygonToLine_management(single_polybuf, "in_memory\\polybuf_outline")
    
    with arcpy.da.SearchCursor(polybuf_outline, ["SHAPE@"]) as cursor:
        for row in cursor:
            lines.append(row[0])
    del row, cursor

    array=np.array(points)

    for i in range(len(array)):
        nTouches=0
        punto=arcpy.Point(array[i][0],array[i][1])
        for a in range (len(lines)):
            touches=(punto.touches(lines[a]) or punto.within(lines[a]))
            if touches==True: ##if the start point touches others
                remove[i] = 1
                break  ##break the loop because only need one touch

    ##Remove the lines corresponding to each point
    with arcpy.da.UpdateCursor(outline, ["SHAPE@"]) as cursor:
        i = 0
        j = 0
        for row in cursor:
            if remove[i] > 0:
                #arcpy.AddMessage("remove")
                cursor.deleteRow()
                j += 1
            i += 1
    del cursor

    if i>0:
        del row
    if (i == 0) or (i == j):
        bEmpty = True
        return outline, bEmpty 

    return outline, bEmpty

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
# This function calcuates the 2D distance of two points
#------------------------------------------------------------------------------------------------------------
def distance2points(x1, y1, x2, y2):
    import math
    return math.sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1))

#------------------------------------------------------------------------------------------------------------
# This function derives the angle between two lines based on the length of the three corresponding points
#------------------------------------------------------------------------------------------------------------
def angle(a, b, c):
    import math
    try:
        cosc = (a*a + b*b -c*c)/(2*a*b)
        return math.acos(cosc) / 3.14159 * 180
    except:  
        #arcpy.AddMessage("math error ignored!!!")
        return 180

#------------------------------------------------------------------------------------------------------------
# This function insert new features into an existing feature class (polygon, polyline or multipoint) based on
# a NumPy array. This function is from by online free code.
#------------------------------------------------------------------------------------------------------------
def numpy_array_to_features(in_fc, in_array, geom_fields, id_field):

    # Establish minimum number of x,y pairs to create proper geometry
    min_xy_dict = {'Polygon': 3, 'Polyline': 2, 'Multipoint': 1}
    min_xy_pairs = min_xy_dict[arcpy.Describe(in_fc).shapeType]
 
    if isinstance(geom_fields, list) and len(geom_fields) == 1:
        # Can't access a single field via a list later, extract the
        # only value
        geom_fields = geom_fields[0]
 
    with arcpy.da.InsertCursor(in_fc, ['SHAPE@', id_field]) as cursor:
        unique_array = np.unique(in_array[id_field])  # unique ids
 
        # Iterate through unique sets, get array that matches unique
        # value, convert coordinates to a list and insert via cursor.
        for unique_value in unique_array:
            a = in_array[in_array[id_field] == unique_value]
            if len(a) >= min_xy_pairs:  # skip if not enough x,y pairs
                cursor.insertRow([a[geom_fields].tolist(), unique_value])
            else:
                pass  # skip if not enough x,y pairs
    del cursor
    
    return

#------------------------------------------------------------------------------------------------------------
# This function smooths the flowline by adjusting the big turns.
#------------------------------------------------------------------------------------------------------------
def flowline_remove_bigturn(flowline, max_angle, cellsize):
    
    arcpy.SimplifyLine_cartography(flowline, "in_memory\\simply_line", 'POINT_REMOVE', (str(cellsize) + ' Meters'))
    arcpy.FeatureVerticesToPoints_management("in_memory\\simply_line", "in_memory\\flowline_points", 'All')

    ###Create the new line after removing the outlier points
    spatialref=arcpy.Describe(flowline).spatialReference
    new_line = arcpy.CreateFeatureclass_management("in_memory", "new_line","POLYLINE", flowline,"","", spatialref)
    arcpy.AddField_management(new_line, "ORIG_FID", "LONG")
    exist_fields = [f.name for f in arcpy.ListFields(flowline)] #List of current field names in outline layer
    fields = exist_fields[2:] ##The first two fields are FID and Geometry
    
    linearray = arcpy.da.FeatureClassToNumPyArray("in_memory\\simply_line", fields)

    pointarray = arcpy.da.FeatureClassToNumPyArray("in_memory\\flowline_points", ('SHAPE@X', 'SHAPE@Y','ORIG_FID'))
    line_ids = np.array([item[2] for item in pointarray])
    unique_line_ids = np.unique(line_ids)

    for fid in unique_line_ids:
        arr = pointarray[line_ids == fid]
        ##Methd 2: move this point until the angle is larger than the max_angle
        for row in range(len(arr)):
            if row <(len(arr)-1) and row > 0:#if it is not first or last point of all
                x1 = float(arr[row-1][0])
                y1 = float(arr[row-1][1])
                x = float(arr[row][0])
                y = float(arr[row][1])
                x2 = float(arr[row+1][0])
                y2 = float(arr[row+1][1])
                length1 = distance2points(x1, y1, x, y)
                length2 = distance2points(x2, y2, x, y)
                length  = distance2points(x1, y1, x2, y2)
                pntangle = angle(length1, length2, length)
                if pntangle < max_angle:
                    midx = (x1 + x2)/2
                    midy = (y1 + y2)/2
                    for i in range(5):
                        newx = x + (midx - x) * (i+1) / 5
                        newy = y + (midy - y) * (i+1) / 5
                        length1 = distance2points(x1, y1, newx, newy)
                        length2 = distance2points(x2, y2, newx, newy)
                        pntangle = angle(length1, length2, length)
                        if pntangle > max_angle:
                            arr[row][0] = newx
                            arr[row][1] = newy
                            break
     
        numpy_array_to_features(new_line, arr, ['SHAPE@X', 'SHAPE@Y'], 'ORIG_FID')

    ##Assign field to the new_line
    arcpy.DeleteField_management(new_line, 'ORIG_FID')

    with arcpy.da.UpdateCursor(new_line, fields) as cursor:
        i = 0
        for row in cursor:
            for j in range(len(fields)):
                row[j] = linearray[i][j]
            cursor.updateRow(row)
            i += 1
    del cursor, row

    return new_line


#------------------------------------------------------------------------------------------------------------
# This function creates the perpendiculars lines based on resolution (step) and the maximum distance.
# The code is revised from Volta centerline.
#------------------------------------------------------------------------------------------------------------
def create_perpendiculars(line, resolution, distance):
    perpendiculars = arcpy.CreateFeatureclass_management("in_memory", "perps","POLYLINE","","","",line)
    new_line_cursor = arcpy.da.InsertCursor(perpendiculars, ('SHAPE@'))
    startxlist = []
    startylist = []
    endxlist  = []
    endylist = []
    with arcpy.da.SearchCursor(line, ["SHAPE@LENGTH", "SHAPE@"]) as cursor:
        for row in cursor:
            length = row[0]
            ##Add the first point
            firstPnt = row[1].firstPoint
            start = 1   ###Start from 1 because the first point already recorded
            while start < length:
                startxlist.append(firstPnt.X)
                startylist.append(firstPnt.Y)
                lastPnt = row[1].positionAlongLine(start).getPart()
                endxlist.append(lastPnt.X)
                endylist.append(lastPnt.Y)
                firstPnt = lastPnt
                start = start + resolution
            ##It is not necessary becasue the line will be extended
    
    ##Create the perpendiculars lines
    for i in range(len(startxlist)):
        startx = startxlist[i]
        starty = startylist[i]
        endx = endxlist[i]
        endy = endylist[i]
        midx = (startx+endx)/2.0
        midy = (starty+endy)/2.0
        if starty==endy or startx==endx:
            if starty == endy:
                y1 = starty + distance
                y2 = starty - distance
                x1 = startx
                x2 = startx
            if startx == endx:
                y1 = starty
                y2 = starty 
                x1 = startx + distance
                x2 = startx - distance     
        else:
            m = ((starty - endy)/(startx - endx)) #get the slope of the line
            negativereciprocal = -1*((startx - endx)/(starty - endy))    #get the negative reciprocal
            if m > 0:
                if m >= 1:
                    y1 = negativereciprocal*(distance)+ starty
                    y2 = negativereciprocal*(-distance) + starty
                    x1 = startx + distance
                    x2 = startx - distance
                if m < 1:
                    y1 = starty + distance
                    y2 = starty - distance
                    x1 = (distance/negativereciprocal) + startx
                    x2 = (-distance/negativereciprocal)+ startx           
            if m < 0:
                if m >= -1:
                    y1 = starty + distance
                    y2 = starty - distance
                    x1 = (distance/negativereciprocal) + startx
                    x2 = (-distance/negativereciprocal)+ startx     
                if m < -1:
                    y1 = negativereciprocal*(distance)+ starty
                    y2 = negativereciprocal*(-distance) + starty
                    x1 = startx + distance
                    x2 = startx - distance
        array = arcpy.Array([arcpy.Point(x1,y1),arcpy.Point(x2, y2)])
        polyline = arcpy.Polyline(array)
        new_line_cursor.insertRow([polyline])
    del row, new_line_cursor, cursor  
    return perpendiculars

#------------------------------------------------------------------------------------------------------------
# This function creates the perpendicular line along each line using a distance.
# It is revised from the codes by James et al.(2016) in Volta.
#------------------------------------------------------------------------------------------------------------
def create_perpendiculars_line_sections(in_lines, distance):
    perpendiculars = arcpy.CreateFeatureclass_management("in_memory", "perpendiculars","POLYLINE")
    new_line_cursor = arcpy.da.InsertCursor(perpendiculars, ('SHAPE@'))
    with arcpy.da.SearchCursor(in_lines, ["SHAPE@"]) as cursor:                                   
        for row in cursor:
            startx = row[0].firstPoint.X
            starty = row[0].firstPoint.Y
            endx = row[0].lastPoint.X
            endy = row[0].lastPoint.Y
            midx = row[0].centroid.X
            midy = row[0].centroid.Y
            if starty==endy or startx==endx:
                if starty == endy:
                    y1 = starty + distance
                    y2 = starty - distance
                    x1 = startx
                    x2 = startx
                if startx == endx:
                    y1 = starty
                    y2 = starty 
                    x1 = startx + distance
                    x2 = startx - distance     
            else:
                m = ((starty - endy)/(startx - endx)) #get the slope of the line
                negativereciprocal = -1*((startx - endx)/(starty - endy))    #get the negative reciprocal
                if m > 0:
                    if m >= 1:
                        y1 = negativereciprocal*(distance)+ starty
                        y2 = negativereciprocal*(-distance) + starty
                        x1 = startx + distance
                        x2 = startx - distance
                    if m < 1:
                        y1 = starty + distance
                        y2 = starty - distance
                        x1 = (distance/negativereciprocal) + startx
                        x2 = (-distance/negativereciprocal)+ startx           
                if m < 0:
                    if m >= -1:
                        y1 = starty + distance
                        y2 = starty - distance
                        x1 = (distance/negativereciprocal) + startx
                        x2 = (-distance/negativereciprocal)+ startx     
                    if m < -1:
                        y1 = negativereciprocal*(distance)+ starty
                        y2 = negativereciprocal*(-distance) + starty
                        x1 = startx + distance
                        x2 = startx - distance
            #arcpy.AddMessage(str(x1) + "; " + str(y1) + " | " + str(x2)+"; " + str(y2))
            array = arcpy.Array([arcpy.Point(x1,y1),arcpy.Point(x2, y2)])
            polyline = arcpy.Polyline(array)
            new_line_cursor.insertRow([polyline])
    del row, cursor
    del new_line_cursor
    
    return perpendiculars

#------------------------------------------------------------------------------------------------------------
# This function creates the axis based on the minmum box of the polygon
# The code is revised from Volta centerline.
#------------------------------------------------------------------------------------------------------------
def create_axis(outline, dem):
    outline_feature_count = arcpy.GetCount_management(outline)          #check outline is a single feature - count"
    outline_feature_count_result = int(outline_feature_count.getOutput(0))      #check outline is a single feature - get result"
    #arcpy.AddMessage(str(outline_feature_count_result))
    if outline_feature_count_result > 1:                                       #check outline is a single feature"
        outline = arcpy.Dissolve_management(outline, "in_memory\\dissolved_outline")
    mbg = arcpy.MinimumBoundingGeometry_management(outline, "in_memory\\mbg", "CONVEX_HULL", "NONE", "","MBG_FIELDS") #create minimum bounding geometry, convex hull method"
    axis = arcpy.XYToLine_management(mbg, "in_memory\\axis_out","MBG_APodX1", "MBG_APodY1", "MBG_APodX2", "MBG_APodY2") # Create long axis from fields in mbg

    ##check if the elevations of the two axis points are biger enough > 80% of the elevation difference
    outline_line = arcpy.PolygonToLine_management(outline, "in_memory\\outline_line")

    ##create a one cellsize buffer around the outline
    cellsize = arcpy.GetRasterProperties_management(dem,"CELLSIZEX")
    cellsize_float = float(cellsize.getOutput(0))
    arcpy.Buffer_analysis(outline_line, "in_memory\\outlinebuf2", (str(cellsize_float)+ " Meter"))

    clipped_dem = ExtractByMask(dem,"in_memory\\outlinebuf2")

    #outDEMrange = ZonalStatistics(outline, arcpy.Describe(outline).OIDFieldName, clipped_dem, "RANGE")
    EleRange =  clipped_dem.maximum - clipped_dem.minimum
    #arcpy.AddMessage("EleRange is: " + str(EleRange))
    
    outaxisRange = ZonalStatistics(axis, arcpy.Describe(axis).OIDFieldName, clipped_dem, "RANGE")
    axisRange = outaxisRange.maximum
    #arcpy.AddMessage("axisRange is: " + str(axisRange))
    
    if axisRange > (0.7 * EleRange):
        arcpy.AddMessage ("Determine the centerline using the long axis of the convex_hull")
        pass
    else: ##calculate the axis based on the highest and lowest points
        arcpy.AddMessage ("Determine the centerline from the highest and lowest points")
        arcpy.Delete_management(axis)

        OutHighestEle = Con(clipped_dem == clipped_dem.maximum,1)  ##Determine the highest flowaccumuation part
        arcpy.RasterToPoint_conversion(OutHighestEle, "in_memory\\highestPoint")
        with arcpy.da.SearchCursor("in_memory\\highestPoint", "SHAPE@XY") as cursor:
            for row in cursor:
                highpointX = row[0][0]
                highpointY = row[0][1]
                break  ##only consider the first point
        del row, cursor

        OutLowestEle = Con(clipped_dem == clipped_dem.minimum,1)  ##Determine the highest flowaccumuation part
        arcpy.RasterToPoint_conversion(OutLowestEle, "in_memory\\lowestPoint")
        with arcpy.da.SearchCursor("in_memory\\lowestPoint", "SHAPE@XY") as cursor:
            for row in cursor:
                lowpointX = row[0][0]
                lowpointY = row[0][1]
                break  ##only consider the first point
        del row, cursor

        axis = arcpy.CreateFeatureclass_management("in_memory", "axis","POLYLINE", "", "", "", outline)
        new_line_cursor = arcpy.da.InsertCursor(axis, ('SHAPE@'))
        array = arcpy.Array([arcpy.Point(lowpointX,lowpointY),arcpy.Point(highpointX, highpointY)])
        polyline = arcpy.Polyline(array)
        new_line_cursor.insertRow([polyline])

    ###Save the two end points of this axis line
    axispoint = arcpy.FeatureVerticesToPoints_management(axis, "in_memory\\axispoint", "BOTH_ENDS") # Export the both end of the axis

    return axis, axispoint

#------------------------------------------------------------------------------------------------------------
# This function extends the orginal line be an extend percentage
# The code is revised from Volta centerline.
#------------------------------------------------------------------------------------------------------------
def extend_line(original_line, extend_percentage):  ###Only valid to extend the line based on the start and end points, and extend the both ends
    arcpy.env.extent = "MAXOF"
    with arcpy.da.SearchCursor(original_line, ["SHAPE@", "SHAPE@LENGTH"]) as cursor:                                         #set up search cursor to access geometry of line
        for row in cursor:
            extend_length_axis = row[1]*extend_percentage
            firstpoint = row[0].firstPoint                                                          #access firstpoint co-ords
            lastpoint = row[0].lastPoint                                                            #access lastpoint co-ords
            if firstpoint.X - lastpoint.X == 0:                                                     #test to see if line is vertical (infinite slope so calculating slope would result in error)
                new_firstpoint_X = firstpoint.X                                                     #Line is vertical so new x co-ordinates will be the same as original
                new_lastpoint_X = lastpoint.X                                                       #Line is vertical so new x co-ordinates will be the same as original
                if firstpoint.Y > lastpoint.Y:                                                      #Test to see if firstpoint above lastpoint
                    new_firstpoint_Y = firstpoint.Y + extend_length_axis                            #Firstpoint is above, so add the extend length
                    new_lastpoint_Y = lastpoint.Y - extend_length_axis                              #Lastpoint is below, so minus the extend length
                else:                                                                               #test if firstpoint is below
                    new_firstpoint_Y = firstpoint.Y - extend_length_axis                            #firstpoint is below, so minus the extend length
                    new_lastpoint_Y = lastpoint.Y + extend_length_axis                              #lastpoint is above so add extend length    
            else:                                                                                   #Line is not vertical, so slope can be calculated
                axis_slope = ((firstpoint.Y - lastpoint.Y)/(firstpoint.X - lastpoint.X))            #Calculate slope (change in y/change in x)
                if axis_slope <0:                                                                   #Test if slope is negative
                    axis_slope = axis_slope*-1                                                      #Invert negative slopes to make calculations simpler
                axis_slope_radians = math.atan(axis_slope)                                          #Convert slope to radians
                changex = float(math.cos(axis_slope_radians))*float(extend_length_axis)             #Calculate amount to change X co-ord by (cos of slope * extend length)
                changey = float(math.sin(axis_slope_radians))*float(extend_length_axis)             #Calculate amount to change Y co-ord by (sin of slope * extend length)
                if firstpoint.X > lastpoint.X:                                                      #Test if firstpoint X co-ord is greater than lastpoint (ie to the right)
                    new_firstpoint_X = firstpoint.X + changex                                       #Firstpoint X  co-ord is greater (to the right) of lastpoint, so add changeX
                    new_lastpoint_X = lastpoint.X - changex                                         #Firstpoint X co-ord is smaller (to the left) of lastpoint, so minus changeX
                else:                                                                               #Test if firstpoint X co-ord is smaller (left) of lastpoint
                    new_firstpoint_X = firstpoint.X - changex                                       #Firstpoint X co-ord is smaller (left) of lastpoint, so minus changex
                    new_lastpoint_X = lastpoint.X + changex                                         #Lastpoint X co-ord is greater (right) of firstpoint, so add changex
                if firstpoint.Y >= lastpoint.Y:                                                     #Test if firstpoint X co-ord is greater (above) lastpoint. Equal implies a horizontal line
                    new_firstpoint_Y = firstpoint.Y + changey                                       #Firstpoint is above, so add change Y
                    new_lastpoint_Y = lastpoint.Y - changey                                         #Lastpoint is below so minus change Y
                else:                                                                               #Test if lastpoint is above firstpoint
                    new_firstpoint_Y = firstpoint.Y - changey                                       #firstpoint is below, so minus change Y
                    new_lastpoint_Y = lastpoint.Y + changey                                         #lastpoint is above so add change Y
            pointfc = arcpy.CreateFeatureclass_management("in_memory","extend_line_pt","POINT","","","",original_line)             #create new feature class (in memory) for new xy points
            point_cursor = arcpy.da.InsertCursor(pointfc, ["SHAPE@XY"])                                         #set up insert cursor to populate feature class
            new_firstpoint = (new_firstpoint_X, new_firstpoint_Y)                                               #gather firstpoint co-ords
            new_lastpoint =(new_lastpoint_X, new_lastpoint_Y)                                                   #gather lastpoint co-ords
            point_cursor.insertRow([new_firstpoint])                                                            #populate firstpoint
            point_cursor.insertRow([new_lastpoint])                                                             #populate lastpoint
            extended_line = arcpy.PointsToLine_management(pointfc, "in_memory\\extended_line_line")                  #draw line (new extended axis) between first and last points

    del row, cursor
    del point_cursor
    
    return extended_line

#------------------------------------------------------------------------------------------------------------
# This function creates a set of points along the line based on a specified distance
# The code is revised from Volta centerline.
#------------------------------------------------------------------------------------------------------------
def points_on_line(line, resolution):   ###Need to add the end point into the output point feature. In this way, extend line is not necessary??!!!
    new_points = arcpy.CreateFeatureclass_management("in_memory", "points","POINT")
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

#------------------------------------------------------------------------------------------------------------
# This function creates equally spaced points along a line from start to end.
# It is revised from the codes in Volta.
#------------------------------------------------------------------------------------------------------------
def points_on_line_withID(line, resolution):
    new_points = arcpy.CreateFeatureclass_management("in_memory", "points","POINT")
    arcpy.AddField_management(new_points, "distance", "DOUBLE")
    arcpy.AddField_management(new_points, "flow_id", "LONG")
    new_points_cursor = arcpy.da.InsertCursor(new_points, ('SHAPE@',"distance","flow_id")) 
    with arcpy.da.SearchCursor(line, ["SHAPE@LENGTH", "SHAPE@", "line_id"]) as cursor:
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
# This function creates the center points of a set of perpentcular lines and extract the elevation for these points.
# The code is revised from Volta centerline.
#------------------------------------------------------------------------------------------------------------
def central_points_on_perps(perps, glacier_outline, axispoint, dem):  ## Need to add the lowest point, so that the flowline can go to the boundary!!!
    clipped_axis_perps = arcpy.Clip_analysis(perps, glacier_outline, "in_memory\\clipped_axis_perps")
    split_axis_perp_lines = arcpy.SplitLine_management(clipped_axis_perps, "in_memory\\split_axis_perps")
    central_points_no_alt = arcpy.FeatureToPoint_management (split_axis_perp_lines, "in_memory\\central_no_alt")
    '''
    ##Need to add the lowest points into the central points_no_alt or just extent the
    axispoint_alt = ExtractValuesToPoints(axispoint, dem, "in_memory\\axispoint_alt")
    max_min_list = []
    with arcpy.da.SearchCursor(axispoint_alt, "RASTERVALU") as cursor:  ##Need to write a code to combine central points with alt with this part using cursor
        for row in cursor:
            max_min_list.append(row[0])  ##This just get the elevation values for these central points, should be able to combined with the previous function
    max_central_point_remove = max(max_min_list)

    with arcpy.da.UpdateCursor(axispoint_alt, "RASTERVALU") as cursor:        ##DELETE MAX AND MIN POUNTS TO STOP CODE FAILING CLOSE TO MARGIN
        for row in cursor:                                                              ###!!This is the problem that the flowline does not extend to the boundary
            if row[0] == max_central_point_remove:
                cursor.deleteRow()  ##Just remove the max elevation and keep the lowest elevation

    arcpy.DeleteField_management(axispoint_alt, "RASTERVALU")
    '''
    arcpy.Append_management(axispoint, central_points_no_alt, "NO_TEST")

    #arcpy.Append_management(axispoint, central_points_no_alt, "NO_TEST")
    ##Extract the elevation for each point
    central_points_with_alt = ExtractValuesToPoints(central_points_no_alt, dem, "in_memory\\central_points_with_alt")

    return central_points_with_alt

#------------------------------------------------------------------------------------------------------------
# This function creates the flowlines based on the all center points derived from the perpendcular lines.
# The code is revised from Volta centerline.
#------------------------------------------------------------------------------------------------------------
def flowline (central_points_with_alt, dem, flow_line_output, flowline_glacier_outline, smooth_tolerance):
    #inaccessible_lowpoint = 0
    glacier_outline_line = arcpy.FeatureToLine_management(flowline_glacier_outline, "in_memory\\glacier_outline_line")
    dissolve_outline_line = arcpy.Dissolve_management(glacier_outline_line, "in_memory\\dissolve_outline_line")
    outline_line_geom = arcpy.CopyFeatures_management(dissolve_outline_line, arcpy.Geometry())

    ##create a 0.5 cellsize buffer around the outline
    cellsize = arcpy.GetRasterProperties_management(dem,"CELLSIZEX")
    cellsize_float = float(cellsize.getOutput(0))
    outline_buf = arcpy.Buffer_analysis(dissolve_outline_line, "in_memory\\outline_buf", (str(cellsize_float)+ " Meter"))
    outline_buf_geom = arcpy.CopyFeatures_management(outline_buf, arcpy.Geometry())

    
    fields = ["RASTERVALU", "SHAPE@XY"] #list of fields used in central points: row[0] = elevation, row[1] = co-ords
    central_point_dict = {} #create empty dictionary for all central points
    flowline_points_dict = {}
    
    max_min_list = []
    with arcpy.da.SearchCursor(central_points_with_alt, "RASTERVALU") as cursor:  ##Need to write a code to combine central points with alt with this part using cursor
        for row in cursor:
            max_min_list.append(row[0])  ##This just get the elevation values for these central points, should be able to combined with the previous function
    max_central_point_remove = max(max_min_list)

    with arcpy.da.UpdateCursor(central_points_with_alt, "RASTERVALU") as cursor:        ##DELETE MAX AND MIN POUNTS TO STOP CODE FAILING CLOSE TO MARGIN
        for row in cursor:                                                              ###!!This is the problem that the flowline does not extend to the boundary
            if row[0] == max_central_point_remove:
                cursor.deleteRow()  ##Just remove the max elevation and keep the lowest elevation
    
    with arcpy.da.SearchCursor(central_points_with_alt, fields) as cursor:  ##I think this is wrong!! The first is Elevation value, the second is point location (x,y)
        for row in cursor:
            if row[0] != -9999:  #Excludes any 'no data' points
                index = row[1] # use cood as index  ##I think the value should assigned reversely
                value = row[0] # use Elevation as value
                central_point_dict[row[1]] = row[0]  #populate dictionary

    #maximum_central_point_alt = max(central_point_dict.itervalues()) #get max elevation value
    maximum_central_point_alt = max(dict_itervalues(central_point_dict)) #get max elevation value
    #min_central_point_alt = min(central_point_dict.itervalues()) #get min elevation value
    min_central_point_alt = min(dict_itervalues(central_point_dict)) #get min elevation value
    max_co_ords_list = [] #set up empty list for high point(s)
    #for index, value in central_point_dict.iteritems():
    for index, value in dict_iteritems(central_point_dict):
        if value == maximum_central_point_alt:
            max_co_ords_list.append(index)

    min_co_ords_list = [] #set up empty list for low point(s)
    #for index, value in central_point_dict.iteritems():
    for index, value in dict_iteritems(central_point_dict):
        if value == min_central_point_alt:
            min_co_ords_list.append(index)

    end_co_ords = min_co_ords_list[0] ##use the first lowest points as the end coordinate

    if len(max_co_ords_list) > 1: #test if more than 1 high point
        #arcpy.AddMessage("WARNING multiple high points, using first point")
        pass

    start_ele = maximum_central_point_alt #set initial starting elevation to maximum 
    start_co_ords = max_co_ords_list[0] #set initial starting co-ords to that of maximum
    flowline_points_dict = {} # Dictionary to hold initial flowline points in
    line_segment_list = [] # list to hold all line segments in 
    lowpoint_stop = 0 #boolean to continue or stop
    while (lowpoint_stop == 0):
        with arcpy.da.SearchCursor(central_points_with_alt, fields) as cursor:
            nearest_neighbour_list = [] #create new empty list for each loop
            for row in cursor:
                if row[0] < start_ele and row[0] != -9999: #find all lower points (excluding nodata)
                    nearest_neighbour_list.append(row[1]) #make list of x co ord and y co ord ###change to start_co_ords
                    if start_co_ords not in nearest_neighbour_list:
                        nearest_neighbour_list.append(start_co_ords) #add starting point to list  ##This is not necessary?? Already added in previous statement; likely pair the point with the highest point toghther???
            len_nearest_neighbour_list = len(nearest_neighbour_list) #get length of list
        if len_nearest_neighbour_list == 0:
            lowpoint_stop = 1 #stop as no more lower nearest neighbours
        else:
            index_val = nearest_neighbour_list.index(start_co_ords)
            ArrayOfPoints = np.array(nearest_neighbour_list) #convert list to numpy array
            KDTreeOfPoints = KDTree(ArrayOfPoints)
            closest_neighbour = 2 #closest neighbour = 2 as 1 = original)
            return_neighbour = 1 #closets neighbour
            distances, indicies = KDTreeOfPoints.query(ArrayOfPoints[index_val],closest_neighbour) #closest neighbour = 2 as 1 = original)
            NearestNeighbor = (nearest_neighbour_list[indicies[return_neighbour]])
            pointlist = []                  #create new list to hold start and end point for segment
            pointlist.append([start_co_ords,NearestNeighbor]) #append start and end point
            features = []
            for feature in pointlist:
                features.append(arcpy.Polyline(arcpy.Array([arcpy.Point(*coords) for coords in feature]),flowline_glacier_outline))
            for segment in features:
                intersect = segment.crosses(outline_line_geom[0])
                if intersect == 0: ##if the section does not intersection the glacier outline
                    line_segment_list.append(segment)
                else:
                    while intersect == 1: ##if the section intersect with the outline
                        if NearestNeighbor != end_co_ords: ##This condition is checked to make sure the line does not cross the outline
                            closest_neighbour = closest_neighbour + 1
                            return_neighbour = return_neighbour + 1 ##Take the next nearest value
                            distances, indicies = KDTreeOfPoints.query(ArrayOfPoints[index_val],closest_neighbour) #Check the next point
                            NearestNeighbor = (nearest_neighbour_list[indicies[return_neighbour]])
                            pointlist_intersect = []                  #create new list to hold start and end point for segment
                            pointlist_intersect.append([start_co_ords,NearestNeighbor]) #append start and end point
                            features_intersect = []
                            for feature in pointlist_intersect:
                                features_intersect.append(arcpy.Polyline(arcpy.Array([arcpy.Point(*coords) for coords in feature]),flowline_glacier_outline))
                            for segment_intersect in features_intersect:
                                intersect = segment_intersect.crosses(outline_line_geom[0])
                                if intersect == 0:
                                    line_segment_list.append(segment_intersect)
                        else: ##This lowest point is OK
                            intersect = 0   
            start_co_ords = NearestNeighbor  #reset start co-ords to that of new nearest neighbour
            start_ele = central_point_dict[NearestNeighbor] #reset start elevation to that of new nearest neighbour
            closest_neighbour = 2 #closest neighbour = 2 as 1 = original)
            return_neighbour = 1 #closest neighbour
    
    #arcpy.AddMessage(len(line_segment_list))
    segment_flowline = arcpy.CopyFeatures_management(line_segment_list, "in_memory\\segment_flowline")

    ##test if the last line segment cross the glaacier outline
    last_line_segment = line_segment_list[-1]
    if last_line_segment.crosses(outline_buf_geom[0]) == 0: ##did not intersect the glacier outline buffer
        #arcpy.AddMessage("The last line segment does not cross the glacier outline")
        ##find the lowest points
        ##create a 0.5 cellsize buffer around the outline
        #cellsize = arcpy.GetRasterProperties_management(dem,"CELLSIZEX")
        #cellsize_float = float(cellsize.getOutput(0))
        #arcpy.Buffer_analysis(dissolve_outline_line, "in_memory\\tempbuf", (str(int(cellsize_float/2))+ " Meter"))
        clipped_dem = ExtractByMask(dem,outline_buf)
        ##get the lowest point
        lowestEle = Con(clipped_dem == clipped_dem.minimum,1)  
        arcpy.RasterToPoint_conversion(lowestEle, "in_memory\\lowestElePoint")
        arcpy.CopyFeatures_management(last_line_segment, "in_memory\\last_line_segment")
        arcpy.FeatureVerticesToPoints_management("in_memory\\last_line_segment", "in_memory\\last_line_points", "BOTH_ENDS")
        arcpy.Append_management("in_memory\\lowestElePoint", "in_memory\\last_line_points", "NO_TEST" )
        arcpy.PointsToLine_management("in_memory\\last_line_points", "in_memory\\last_line")
        arcpy.Append_management("in_memory\\last_line", segment_flowline, "NO_TEST" )
        
    dissolved_flowline = arcpy.Dissolve_management(segment_flowline, flow_line_output)

    ###Need to make sure that the flowline_output extended to the lowest points of the glacier outlines
    #arcpy.FeatureVerticesToPoints_management(inFeatures, outFeatureClass, "MID")

    
    smooth_flowline = CA.SmoothLine(dissolved_flowline, "in_memory\\flow_line", "PAEK", smooth_tolerance)
    
    return smooth_flowline


#------------------------------------------------------------------------------------------------------------
# This function detect new branch based on the tributary ratio and tributary source threshold.
# a new branch will be returned if there is a new branch.
# The code is revised from Volta centerline.
#------------------------------------------------------------------------------------------------------------
def new_branch(input_flowline, dem, branch_outline, input_outline, TributaryRatio, TributarySourceThreshold):  
    #arcpy.AddMessage("New branch Check...")
    cellsize = arcpy.GetRasterProperties_management(dem,"CELLSIZEX")
    cellsize_float = float(cellsize.getOutput(0))

    clipped_dem_raster = ExtractByMask(dem,input_outline)  
    arcpy.env.extent = clipped_dem_raster
    new_branch_count = 0
    dem_numpy = arcpy.RasterToNumPyArray(clipped_dem_raster,"","","",0)
    d8_structure = ndimage.generate_binary_structure(2, 2)
    upstream_area_original = 1
    glacier_numpy = dem_numpy.copy()
    glacier_numpy[glacier_numpy != 0] = 1  ##So galcier_numpy is a 0/1 array

    with arcpy.da.SearchCursor(input_flowline, ["SHAPE@LENGTH"]) as cursor: ##Just get the flowline length
        for row in cursor:
            flowline_length = row[0]
    del row, cursor
    
    flowline_point_res = flowline_length/100.0  ## Need to make a smallest distance for this, if glacier is smaller,  it is not needed to divide 100
    flowline_point_res = max(flowline_point_res, cellsize_float*3) ##set DEM resolution (30 meter) as the minimum distance by Yingkui Li
    flowline_points = points_on_line(input_flowline, flowline_point_res) ##Divide the flowline into 100 points along the flowline
    branch_points_alt = ExtractValuesToPoints(flowline_points, dem, "in_memory\\branch_points_altitude")

    FcID = arcpy.Describe(branch_points_alt).OIDFieldName
    pntcount_result = arcpy.GetCount_management(branch_points_alt)
    pntcount = int(pntcount_result.getOutput(0))
    #arcpy.AddMessage("point count is:" + str(pntcount))
    for i in range(1, pntcount):
        tempRs = arcpy.env.scratchFolder + "\\r" + str(i) ###Create a temp raster file. In-memory doesnot work for raster???
        if new_branch_count == 0:
            query = FcID+" = "+str(i+1)
            arcpy.Select_analysis(branch_points_alt, "in_memory\\select_branch_points", query) ###Just do a simply select analysis
            arcpy.PointToRaster_conversion("in_memory\\select_branch_points", "RASTERVALU", tempRs, "","", cellsize_float)
            flowline_numpy = arcpy.RasterToNumPyArray(tempRs,"","","",0)
            #delete the tempRs data
            arcpy.Delete_management(tempRs)
            
            original_cell = np.amax(flowline_numpy)
            itemindex = np.where(flowline_numpy == original_cell)
            row = itemindex[0]
            col = itemindex[1]
            dem_copy = dem_numpy.copy()
            dem_copy[dem_copy < original_cell] = 0
            dem_copy[dem_copy >= original_cell] = 1
            labeled_array, num_features = scipy.ndimage.measurements.label(dem_copy, structure=d8_structure) ## All these can be simplified using the array comparsion function and then count how many values are higher than the original elevation
            zone = labeled_array[row,col]   ##The scipy.ndimage.measurements function consider the connectivity issue of the cells, so cannot simply replaced by greater than, unless using region group
            zone_cells = np.argwhere(labeled_array == zone) ##Find the same zone with the corresponding point, ##How about using a flow accumulation along a cross section (horizontal surface cut the DEM/facc) of the point
            total_cells = int(len(zone_cells))                 ## It seems that the scipy.ndimage method is still the best
            upstream_area = total_cells*cellsize_float*cellsize_float
            #arcpy.AddMessage(str(upstream_area_original))   
            #if (upstream_area_original > 0):   # in memory always starts from 1
            increase_ratio = float(i)/float(i+1) ##This is assume a natural increase in area becasue the point moves down
            percentage_change = (float(upstream_area* increase_ratio)/float(upstream_area_original+1))
            upstream_area_original = upstream_area
                
            if upstream_area < TributarySourceThreshold or percentage_change < (1.0 + TributaryRatio) or percentage_change >5.0: ##>5.0 is for the start point
                newarray = labeled_array.copy()
                newarray[newarray != zone] = 0
                newarray[newarray == zone] = 1
                glacier_numpy[newarray == 1] = 0     ##So in this case, a new branch will not be created???, but change glacier_numpy array
            else:
                left_co_ord = arcpy.GetRasterProperties_management(clipped_dem_raster,"LEFT")
                bottom_co_ord = arcpy.GetRasterProperties_management(clipped_dem_raster,"BOTTOM")
                left_co_ord_float = float(left_co_ord.getOutput(0))
                bottom_co_ord_float = float(bottom_co_ord.getOutput(0))
                bottom_left_point = arcpy.Point(left_co_ord_float,bottom_co_ord_float)
                new_branch_raster = Int(arcpy.NumPyArrayToRaster(glacier_numpy, bottom_left_point, cellsize_float, cellsize_float, 0))
                new_branch_outline_unclipped = arcpy.RasterToPolygon_conversion(new_branch_raster, "in_memory\\new_branch_outline_unclipped", "SIMPLIFY")

                ##Make sure to remove spurious polygons
                branchcount = arcpy.GetCount_management(new_branch_outline_unclipped) 
                branchcount_result = int(branchcount.getOutput(0))
                #arcpy.AddMessage(str(branchcount_result))
                if branchcount_result > 0:
                    if branchcount_result > 1: ##if more than one polygon, just use the largest one
                        maxArea = 0
                        with arcpy.da.SearchCursor(new_branch_outline_unclipped, "SHAPE@AREA") as cursor:        
                            for row in cursor:                                                             
                                if row[0] > maxArea:
                                    maxArea = row[0]   ##Just remove the max elevation and keep the lowest elevation
                        del cursor, row
                        with arcpy.da.UpdateCursor(new_branch_outline_unclipped, "SHAPE@AREA") as cursor:        
                            for row in cursor:                                                              
                                if row[0] < maxArea:
                                    arcpy.AddMessage("Delete spurous polygon!")
                                    cursor.deleteRow()  
                        del cursor, row
                    #try:
                    new_branch_outline = arcpy.Clip_analysis(new_branch_outline_unclipped, input_outline, branch_outline)
                    #    #arcpy.AddMessage("an process correctly")
                    #    new_branch_count = 1
                    #except:
                    #    new_branch_outline = branch_outline
                    #    arcpy.AddMessage("an error happens")
                    #    new_branch_count = 0
                    new_branch_count = 1
                    break

    if new_branch_count > 0:
        #arcpy.AddMessage("Return new_branch_outline!!")
        return new_branch_outline, new_branch_count
    else:
        #arcpy.AddMessage("Return no new branch_outline!!")
        return None, new_branch_count

    #return branch_outline, new_branch_count

#------------------------------------------------------------------------------------------------------------
# This function creates parallels points based on a series of points.
# It is revised from the codes by in Volta.
#------------------------------------------------------------------------------------------------------------
def create_parallels_points(pntP, field, distance):
    spatial_reference = arcpy.Describe(pntP).spatialReference 
    parallels = arcpy.CreateFeatureclass_management("in_memory", "parallels","POINT", pntP, "DISABLED", "DISABLED", spatial_reference)
    #Get the field values for all points
    fieldValueList = []
    with arcpy.da.SearchCursor(pntP,field) as cursor:
        for row in cursor:
            fieldValueList.append(row[0])
    del row, cursor
    
    #Create a polyline feature based on the points
    inline = arcpy.PointsToLine_management(pntP, "in_memory\\inline")
    splitline = arcpy.SplitLine_management(inline, "in_memory\\splitline")

    new_pnt_cursor = arcpy.da.InsertCursor(parallels, ('SHAPE@', field))

    with arcpy.da.SearchCursor(splitline, ["SHAPE@"]) as cursor:
        i = 0
        for row in cursor:
            startx = row[0].firstPoint.X
            starty = row[0].firstPoint.Y
            endx = row[0].lastPoint.X
            endy = row[0].lastPoint.Y
            midx = row[0].centroid.X
            midy = row[0].centroid.Y
            if starty==endy or startx==endx:
                if starty == endy:
                    y1 = starty + distance
                    y2 = starty - distance
                    x1 = startx
                    x2 = startx
                if startx == endx:
                    y1 = starty
                    y2 = starty 
                    x1 = startx + distance
                    x2 = startx - distance     
            else:
                m = ((starty - endy)/(startx - endx)) #get the slope of the line
                negativereciprocal = -1*((startx - endx)/(starty - endy))    #get the negative reciprocal
                if m > 0:
                    if m >= 1:
                        y1 = negativereciprocal*(distance)+ starty
                        y2 = negativereciprocal*(-distance) + starty
                        x1 = startx + distance
                        x2 = startx - distance
                    if m < 1:
                        y1 = starty + distance
                        y2 = starty - distance
                        x1 = (distance/negativereciprocal) + startx
                        x2 = (-distance/negativereciprocal)+ startx           
                if m < 0:
                    if m >= -1:
                        y1 = starty + distance
                        y2 = starty - distance
                        x1 = (distance/negativereciprocal) + startx
                        x2 = (-distance/negativereciprocal)+ startx     
                    if m < -1:
                        y1 = negativereciprocal*(distance)+ starty
                        y2 = negativereciprocal*(-distance) + starty
                        x1 = startx + distance
                        x2 = startx - distance
            pnt1 = arcpy.Point(x1,y1)
            new_pnt_cursor.insertRow([pnt1, fieldValueList[i]])
            pnt2 = arcpy.Point(x2, y2)
            new_pnt_cursor.insertRow([pnt2, fieldValueList[i]])
            i = i + 1
    del row, cursor

    del new_pnt_cursor

    return parallels
 

#------------------------------------------------------------------------------------------------------------
# This function is to merge undissolved flowlines in combining the flowlines with centerlines
# Added on 5/26/2022 by Yingkui Li
#------------------------------------------------------------------------------------------------------------
def MergeUndissolvedFlowlines(lines, field):

    dissolvedlines = "in_memory\\dissolvedlines"   
    lines_layer = arcpy.MakeFeatureLayer_management(lines, "in_memory\\lines_layer")

    arr=arcpy.da.FeatureClassToNumPyArray(lines, field)
    merge_ids = np.array([item[0] for item in arr])
    #unique_merge_ids = np.unique(merge_ids)
    unique_merge_ids, counts = np.unique(merge_ids, return_counts=True)

    idx_result = np.where(counts > 2) ##only find the lines with at least three lines with the same MergeID
    idx_arr = idx_result[0]
    if len(idx_arr) > 0:
        ChangeID = []
        for i in range(len(idx_arr)):
            idx = idx_arr[i]
            value_id = unique_merge_ids[idx]

            query = field + " = " + str(value_id)

            arcpy.SelectLayerByAttribute_management (lines_layer, "NEW_SELECTION", query)

            startpoints = []
            #endpoints = []
            linelength=[]
            FcID = []
            ##First loop to get the startpoints, lines and the glacierID list
            OBJID = arcpy.Describe(lines).OIDFieldName
            with arcpy.da.SearchCursor(lines_layer, ["SHAPE@", "SHAPE@LENGTH", OBJID]) as flows:
                for flow in flows:
                    startpoints.append(((flow[0].firstPoint).X, (flow[0].firstPoint).Y))
                    #endpoints.append(((flow[0].lastPoint).X, (flow[0].lastPoint).Y))
                    linelength.append(flow[1])
                    FcID.append(flow[2]) ##default are all -1 for each line
            del flow, flows

            startarray = np.array(startpoints)

            for i in range(len(startarray)):
                startpnt = arcpy.Point(startarray[i][0],startarray[i][1])
                
                for a in range (len(startarray)):
                    if a != i: ##don't touch itself
                        comparepoint = arcpy.Point(startarray[a][0],startarray[a][1])
                        touches=(startpnt.equals (comparepoint))
                        if touches==True: ##if the start point touches others
                            if linelength[i] < linelength[a]:
                                ChangeID.append(FcID[i])
                            else:
                                ChangeID.append(FcID[a]) ##set the short line id as -1
                            
        arcpy.SelectLayerByAttribute_management (lines_layer, "CLEAR_SELECTION")
        ChangeIDArr = np.array(ChangeID)
        UniqueChangeIDArr = np.unique(ChangeIDArr)

        if len(UniqueChangeIDArr) > 0:
            with arcpy.da.UpdateCursor(lines, [arcpy.Describe(lines).OIDFieldName, field]) as cursor:
                for row in cursor:
                    flag = np.where(UniqueChangeIDArr == row[0],1,0)
                    if np.sum(flag) > 0:
                        row[1] = -1  ##set the mergeId as -1
                        cursor.updateRow(row)
            del cursor, row

        ##dissolve the lines again based on the revised mergeID
        arcpy.Dissolve_management(lines, dissolvedlines, field, "#", 'SINGLE_PART', 'UNSPLIT_LINES')
    else:
        arcpy.CopyFeatures_management(lines, dissolvedlines)

    arcpy.Delete_management(lines_layer)
    
    return dissolvedlines
