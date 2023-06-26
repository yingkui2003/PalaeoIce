#-------------------------------------------------------------------------------
# Name: SmoothFlowlineElevations.py
# Purpose: This tool smooths the elevations along the specified flowline section that may be
# affected by the down-cut of rivers under the ice or after glacier retreat. This tool first
# uses a buffer of the flowline section to extract the elevations within the buffer boundary
# lines. Then, interpret the elevation surface within the buffer zone based on the points
# from the buffer boundary lines. Finally, the DEM is adjusted by replacing the interpreted
# elevations within the buffer zone.
#
# Author: Dr. Yingkui Li
# Created:     02/09/2023
# Updated:     
# revised:     
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
inflowlines = arcpy.GetParameterAsText(1)
buffer_dis = arcpy.GetParameter(2)

outDEM = arcpy.GetParameterAsText(3)

arcpy.Delete_management("in_memory")

spatialref=arcpy.Describe(inDEM).spatialReference

cellsize = arcpy.GetRasterProperties_management(inDEM,"CELLSIZEX")
cellsize_int = int(float(cellsize.getOutput(0)))

arcpy.env.extent = inDEM
arcpy.env.snapRaster = inDEM

tmpbuf = "in_memory\\tmpbuf"
str_buffer_dis = str(buffer_dis) + " Meters"
arcpy.Buffer_analysis(inflowlines, tmpbuf, str_buffer_dis)
buf_line = "in_memory\\buf_line"
arcpy.PolygonToLine_management(tmpbuf, buf_line)
arcpy.FeatureVerticesToPoints_management(buf_line, "in_memory\\buf_points", "ALL")
ExtractValuesToPoints("in_memory\\buf_points", inDEM, "in_memory\\points3D")

#arcpy.CopyFeatures_management("in_memory\\points3D", "c:\\test\\points3D.shp")
##Using Kriging to get the elevation surface
kriging = arcpy.env.scratchGDB + "\\kriging"
arcpy.Kriging_3d("in_memory\\points3D", "RASTERVALU", kriging,"Circular", str(cellsize_int), "Variable 5")

##Then use the buffer to extract the elevation within the buffer
extKrig = ExtractByMask(kriging, tmpbuf)
outKrigingCon = Con(IsNull(extKrig), 0, extKrig)

##Convert buffer to a raster
FcID = arcpy.Describe(tmpbuf).OIDFieldName
arcpy.PolygonToRaster_conversion(tmpbuf, FcID, "in_memory\\bufRaster", "#", "#", cellsize_int)
outbufferCon = Con(IsNull("in_memory\\bufRaster"), 1, 0)

outMinus = inDEM * outbufferCon + outKrigingCon
outMinus.save(outDEM)


arcpy.Delete_management("in_memory")
