#-------------------------------------------------------------------------------
# Name: DEMLakeExtraction.py
# Purpose: This tool extracts lakes from a DEM based on a minimum area threshold
# The lake in a DEM should have zero slope within the lake area. Therefore, the slope
# analysis can be used to extract lakes.
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
#import numpy as np
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
min_area_km2 = arcpy.GetParameter(1)
outLakes = arcpy.GetParameterAsText(2)

arcpy.Delete_management("in_memory")

spatialref=arcpy.Describe(inDEM).spatialReference

arcpy.env.extent = inDEM
arcpy.env.snapRaster = inDEM
arcpy.env.cellSize = inDEM

OutSlope = Slope(inDEM)
ConSlope = Con(OutSlope == 0, 1, 0)
outFocalMax = FocalStatistics(ConSlope, "", "MAXIMUM")
OutFlat = Con(outFocalMax > 0, 1)
arcpy.RasterToPolygon_conversion(OutFlat, outLakes)

min_area = min_area_km2 * 1e6
#remove the potential spurious polygons
with arcpy.da.UpdateCursor(outLakes, "SHAPE@AREA") as cursor:        
    for row in cursor:
        if row[0] < min_area:
            cursor.deleteRow()  ##Just remove the max elevation and keep the lowest elevation
del cursor, row

arcpy.Delete_management("in_memory")
