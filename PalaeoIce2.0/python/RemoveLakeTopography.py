#-------------------------------------------------------------------------------
# Name: RemoveLakeTopography.py
# Purpose: This tool is designed to estiamte lake bathmetric topography based on lake outlines and
# specifed depth contours, points or simply a maximum depth in the attrbute table
# Then, the estimated lake thickness is removed from the DEM to generate a bare-earth DEM
#
# Author: Dr. Yingkui Li
# Created:     02/09/2023
# Department of Geography, University of Tennessee
# Knoxville, TN 37996
#-------------------------------------------------------------------------------

# Import arcpy module
import arcpy, sys
from arcpy import env
from arcpy.sa import *
#import numpy
import numpy as np
#import time

arcpy.env.overwriteOutput = True
arcpy.env.XYTolerance= "0.01 Meters"

arcpy.Delete_management("in_memory") ### Empty the in_memory
ArcGISPro = 0
arcpy.AddMessage("The current python version is: " + str(sys.version_info[0]))
if sys.version_info[0] == 2:  ##For ArcGIS 10, need to check the 3D and Spatial Extensions
    try:
        if arcpy.CheckExtension("Spatial")=="Available":
            arcpy.CheckOutExtension("Spatial")
        else:
            raise Exception ("not extension available")
            #print "not extension available"
    except:
        raise Exception ("unable to check out extension")
        #print "unable to check out extension"

    try:
        if arcpy.CheckExtension("3D")=="Available":
            arcpy.CheckOutExtension("3D")
        else:
            raise Exception ("not extension available")
            #print "not extension available"
    except:
        raise Exception ("unable to check out extension")
        #print "unable to check out extension"
elif sys.version_info[0] == 3:  ##For ArcGIS Pro
    ArcGISPro = 1
    #pass ##No need to Check
else:
    raise Exception("Must be using Python 2.x or 3.x")
    exit()   



##main program
inDEM = arcpy.GetParameterAsText(0)
in_lake_outlines = arcpy.GetParameterAsText(1)
in_lake_contours = arcpy.GetParameterAsText(2)
contour_field = arcpy.GetParameterAsText(3)
in_lake_ele_points = arcpy.GetParameterAsText(4)
pnt_field = arcpy.GetParameterAsText(5)

outLakeThickness = arcpy.GetParameterAsText(6)
outAdjDEM = arcpy.GetParameterAsText(7)

arcpy.Delete_management("in_memory")

spatialref=arcpy.Describe(inDEM).spatialReference

cellsize = arcpy.GetRasterProperties_management(inDEM,"CELLSIZEX")
cellsize_int = int(float(cellsize.getOutput(0)))

oldextent = arcpy.env.extent
arcpy.env.snapRaster = inDEM
arcpy.env.cellSize = inDEM

lake_outline = "in_memory\\lake_outline"
outline_lines = "in_memory\\outline_lines"
sel_outline_line = "in_memory\\sel_outline_line"
lakeTckness = arcpy.env.scratchFolder + "\\r" + "lakeTckness"

b_lake_contours = True
if len(in_lake_contours) == 0:
    b_lake_contours = False



FcID = arcpy.Describe(in_lake_outlines).OIDFieldName
lakeArray = arcpy.da.FeatureClassToNumPyArray(in_lake_outlines,FcID)
lakeIDs = np.array([item[0] for item in lakeArray])
uniquelakeID = np.unique(lakeIDs)

arcpy.PolygonToLine_management(in_lake_outlines, outline_lines)

arcpy.AddField_management(outline_lines, "contour", "LONG")
arcpy.CalculateField_management(outline_lines,"contour",0)

b_flag = 0
for lake_id in uniquelakeID:
    try:
        #arcpy.env.extent = oldextent
        query = FcID + " = " + str(lake_id)
        arcpy.AddMessage("Processing Depth Interpretation for Lake #" + str(lake_id))                                                                                       
        arcpy.Select_analysis (in_lake_outlines, lake_outline, query)
        
        arcpy.Buffer_analysis(lake_outline, "in_memory\\lake_outline_buf", "100 Meter", "#", "#", "ALL")
                  
        #arcpy.env.extent = "in_memory\\lake_outline_buf"
        
        arcpy.SpatialJoin_analysis(outline_lines, lake_outline, sel_outline_line,
                                       "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT")
        #arcpy.PolygonToLine_management(lake_outline, outline_lines)
        
        #arcpy.AddField_management(outline_lines, "contour", "LONG")
        #arcpy.CalculateField_management(outline_lines,"contour",0)

        in_boundary = TopoBoundary ([lake_outline])

        ##consider lake contours
        if b_lake_contours:
            arcpy.Clip_analysis (in_lake_contours, lake_outline, "in_memory\\clip_lake_contours")
            lineArray = arcpy.da.FeatureClassToNumPyArray("in_memory\\clip_lake_contours",contour_field)
            if len(lineArray)> 0:
                #fieldlist=[field.name for field in arcpy.ListFields("in_memory\\clip_lake_contours")]
                if contour_field != "contour":
                    arcpy.AddField_management("in_memory\\clip_lake_contours", "contour", "LONG")
                    arcpy.CalculateField_management("in_memory\\clip_lake_contours","contour", str("!"+contour_field+"!"),"PYTHON_9.3")
                   
                arcpy.Append_management("in_memory\\clip_lake_contours", sel_outline_line, "NO_TEST" )

        in_contours = TopoContour([[sel_outline_line, 'contour']])

        ##consider lake elevation points
        arcpy.Clip_analysis (in_lake_ele_points, lake_outline, "in_memory\\clip_lake_points")
        pntArray = arcpy.da.FeatureClassToNumPyArray("in_memory\\clip_lake_points",pnt_field)
        b_Process = True
        #arcpy.AddMessage(str(len(pntArray)))
        
        if len(pntArray)== 0:
            exist_fields = [f.name for f in arcpy.ListFields(lake_outline)] #List of current field names in outline layer
            if pnt_field in exist_fields:
                arcpy.AddMessage("No points within the lake, but has the depth attribute in the lake polygon")
                ##Create the centerpoints for the lake
                arcpy.FeatureToPoint_management(lake_outline, "in_memory\\clip_lake_points", "INSIDE")
            else:
                arcpy.AddMessage("No points within the lake, and also no depth field in the lake polygon!")
                b_Process = False
                
             
        if b_Process == True:
            #in_points = TopoPointElevation([["in_memory\\clip_lake_points", pnt_field]])
            ##Need to interpret more points along the centerlines
            if ArcGISPro == 1:
                arcpy.topographic.PolygonToCenterline(lake_outline, "in_memory\\lakecenterline", None)
            else:
                arcpy.cartography.CollapseDualLinesToCenterline(lake_outline, "in_memory\\lakecenterline", 500) ##Need to check or write a code to extract the centerlines

            ##Interpret the points along the lake centerlines
            arcpy.FeatureVerticesToPoints_management("in_memory\\lakecenterline", "in_memory\\centerline_points", "All")
            arcpy.Near_analysis ("in_memory\\clip_lake_points", "in_memory\\centerline_points")#identify nearest flowline point
            pntArray2 = arcpy.da.FeatureClassToNumPyArray("in_memory\\clip_lake_points",["NEAR_FID",pnt_field])
            nearIDs = np.array([item[0] for item in pntArray2])
            pntFields = np.array([item[1] for item in pntArray2])
            arcpy.AddField_management("in_memory\\centerline_points", pnt_field, "Double")
            FcID = arcpy.Describe("in_memory\\centerline_points").OIDFieldName
            with arcpy.da.UpdateCursor("in_memory\\centerline_points",(FcID, pnt_field)) as cursor:   #populate ice field with value from the nearest flowline point
                for row in cursor:
                    fid = row[0]
                    if fid in nearIDs:
                        idx_result = np.where(nearIDs == fid)
                        idx = idx_result[0][0]
                        row[1]=pntFields[idx]
                        cursor.updateRow(row)
                    else:##delete the other points
                        cursor.deleteRow()  
            del cursor, row
            ##Need to add the start and end points of the lake centerlines to the centerline_points with depth of zero
            arcpy.FeatureVerticesToPoints_management("in_memory\\lakecenterline", "in_memory\\start_end_points", "BOTH_ENDS")
            arcpy.AddField_management("in_memory\\start_end_points", pnt_field, "Double")
            arcpy.CalculateField_management("in_memory\\start_end_points",pnt_field,0)
            arcpy.Append_management("in_memory\\start_end_points", "in_memory\\centerline_points", "NO_TEST" )
            
            ##Split centerline by the cleaned centerline points with the depth value
            arcpy.SplitLineAtPoint_management("in_memory\\lakecenterline", "in_memory\\centerline_points", "in_memory\\split_centerlines", "1 Meters")
            #arcpy.AddField_management("in_memory\\split_centerlines", "Start", "Double")
            #arcpy.AddField_management("in_memory\\split_centerlines", "End", "Double")

            #arcpy.CopyFeatures_management("in_memory\\split_centerlines", "c:\\test\\split_centerlines.shp")
            
            arcpy.FeatureVerticesToPoints_management("in_memory\\split_centerlines", "in_memory\\start_points", "START")
            arcpy.SpatialJoin_analysis("in_memory\\start_points", "in_memory\\centerline_points", "in_memory\\start_pnt_with_depth",
                                       "JOIN_ONE_TO_ONE", "KEEP_ALL", '#', "INTERSECT", "1 Meters", "#")
            StartArray = arcpy.da.FeatureClassToNumPyArray("in_memory\\start_pnt_with_depth", pnt_field)
            Starts = np.array([item[0] for item in StartArray])
            #arcpy.AddMessage(Starts)
            
            arcpy.FeatureVerticesToPoints_management("in_memory\\split_centerlines", "in_memory\\end_points", "END")
            arcpy.SpatialJoin_analysis("in_memory\\end_points", "in_memory\\centerline_points", "in_memory\\end_pnt_with_depth",
                                       "JOIN_ONE_TO_ONE", "KEEP_ALL", '#', "INTERSECT", "1 Meters", "#")
            EndArray = arcpy.da.FeatureClassToNumPyArray("in_memory\\end_pnt_with_depth", pnt_field)
            Ends = np.array([item[0] for item in EndArray])
            #arcpy.AddMessage(Ends)
            
            new_points = arcpy.CreateFeatureclass_management("in_memory", "points","POINT")
            arcpy.AddField_management(new_points, pnt_field, "Double")
            
            new_points_cursor = arcpy.da.InsertCursor(new_points, ('SHAPE@', pnt_field))
            resolution = cellsize_int * 3 ##set the resolution

            with arcpy.da.SearchCursor("in_memory\\split_centerlines", ["SHAPE@LENGTH", "SHAPE@"]) as cursor:
                i = 0
                for row in cursor:
                    start = 0
                    depth = Starts[i]
                    length = row[0]
                    depth_interval = (Ends[i] - Starts[i]) * resolution / length
                    #arcpy.AddMessage(str(length))
                    #arcpy.AddMessage(str(depth_interval))

                    while start < length:
                        new_point = row[1].positionAlongLine(start)
                        new_points_cursor.insertRow((new_point, depth))
                        start += resolution
                        depth += depth_interval
                    i += 1
            del row, cursor

            del new_points_cursor
                    
            #arcpy.CopyFeatures_management(new_points, "c:\\test\\new_points.shp")

            in_points = TopoPointElevation([[new_points, pnt_field]])
            
            interpolated_lake_depth = TopoToRaster([in_points, in_contours, in_boundary])##, cellsize_int, "", '20', '0', '#', 'NO_ENFORCE', "SPOT", '1', '#', '1', '0', '0', '200') ##,"","","","",0.1)
        else: ##Need to consider no lake contour and no points
            interpolated_lake_depth = TopoToRaster([in_contours, in_boundary])##, cellsize_int, "", '20', '0', '#', 'NO_ENFORCE', "SPOT", '1', '#', '1', '0', '0', '200') ##,"","","","",0.1)
            
        outTckness = Con(interpolated_lake_depth > 0, interpolated_lake_depth, 0)

        if b_flag == 0:
            arcpy.CopyRaster_management(outTckness, lakeTckness)
            b_flag = 1
        else:
            arcpy.Mosaic_management(outTckness, lakeTckness, "MEAN","","", "", "", "", "")
    except:
        arcpy.AddMessage("Cannot derive the lake depth for this lake due to the lack of depth information!")
        pass
    
arcpy.env.extent = inDEM

outConTck = Con(Raster(lakeTckness) > 0, Raster(lakeTckness))

outConTck.save(outLakeThickness)


arcpy.CalculateStatistics_management(outConTck)
outFillNull = Con(IsNull(outConTck), 0, outConTck)



outMinus = inDEM - outFillNull

outMinus.save(outAdjDEM)

arcpy.env.extent = oldextent

arcpy.Delete_management("in_memory")
