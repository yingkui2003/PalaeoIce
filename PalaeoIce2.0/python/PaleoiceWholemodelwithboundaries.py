#-------------------------------------------------------------------------------
# Name: PaleoiceWholemodelwithboundaries
# Purpose: This tool runs the whole paleoice reconstruction from the basic input of DEM and moraine
# or crosssection features representing the old terminal of glaciers. The user can also specify the
# the paleoglacier polygon boundaries reconstructed by geomorphic evidence. Some other
# parameters can also be specified to run the program, for example, the extent of extant glaicers if
# existing, the tributarr ratio, steps for the ice flowline interpretation, etc.
#
# This tool will first generate the flowline from stream network using the moraine or cross-section
# features and then use the flowline and DEM, as well as target features to reconstruct the exent of
# the paleoglaciers. If modern glaciers exist, the tool will generate centerlines and interpret the ice
# thickness for these glaciers. Then, using the ice thickness to adjust the DEM and combine the
# centerlines with the flowlines. The paleoglacier will then reconstrcuted based on the flowlines and
# adjusted DEM, as well as the target features.
# 
# Author:      Yingkui Li
# Created:     6/6/2022
# Copyright:   (c) Yingkui Li 2022
# Department of Geography, University of Tennessee
# Knoxville, TN 37996, USA
#-------------------------------------------------------------------------------

from SharedFunctions import *  ## 
from volta_centreline import Centerline_Volta
from AdjustDEM import Adjust_DEM_by_Subtracting_Ice_Thickness
from voltathicknessFinal import Ice_Thickness_Volta
from CenterlineFromStream import Flowline_from_Stream_Network
from Paleoicerewithboundaries import PaleoIceReconstructionwithboundary
from CombineFlowlineCenterline import Combine_Flowlines_with_Centerlines

# Script arguments
##Input Datasets
Input_DEM = arcpy.GetParameterAsText(0)
Input_moraines_or_cross_sections = arcpy.GetParameterAsText(1)
Input_Paleoice_boundary = arcpy.GetParameterAsText(2)
Input_Glacier_Outlines = arcpy.GetParameterAsText(3)

##Input arguments
StreamThresholdKM2 = float(arcpy.GetParameter(4))
TributaryThresholdKM2 = float(arcpy.GetParameter(5))
TributaryRatio = float(arcpy.GetParameter(6)) ##Define a tributary ratio
Distance=int(arcpy.GetParameter(7))

CenterlineTributaryRatio = float(arcpy.GetParameter(8))
CenterlineTributaryAreaRatio = float(arcpy.GetParameter(9))
slope_limit = float(arcpy.GetParameter(10))
min_slope = float(arcpy.GetParameter(11))
shearstress = float(arcpy.GetParameter(12))
bFactorPolyfit = arcpy.GetParameter(13)

##Output Datasets
Output_Flowline_Points = arcpy.GetParameterAsText(14)
Output_Ice_Surfaces = arcpy.GetParameterAsText(15)
Output_Ice_Thickness = arcpy.GetParameterAsText(16)


Centerlines = arcpy.env.scratchGDB + "\\Centerlines"
icethick = arcpy.env.scratchGDB + "\\icethick"
adjustDEM = arcpy.env.scratchGDB + "\\adjustDEM"
flowlines = arcpy.env.scratchGDB + "\\flowlines"
watersheds = arcpy.env.scratchGDB + "\\watersheds"
Combined_Flowline = arcpy.env.scratchGDB + "\\Combined_Flowline"
volta_fl_points = arcpy.env.scratchGDB + "\\volta_fl_points"

##Step 1: Generate flowlines from stream network
arcpy.Delete_management("in_memory") ### Empty the in_memory
arcpy.AddMessage("Generate flowlines from stream network...")
Flowline_from_Stream_Network (Input_DEM, Input_moraines_or_cross_sections, StreamThresholdKM2, TributaryThresholdKM2, TributaryRatio, flowlines, watersheds)

##Step 2:
if Input_Glacier_Outlines != "":
    ##Step 2: generate centerlines from glacier outlines
    arcpy.Delete_management("in_memory") ### Empty the in_memory
    arcpy.AddMessage("Generate centerlines from extant glacier outlines...")
    arcpy.Delete_management("in_memory") ### Empty the in_memory
    Centerline_Volta (Input_Glacier_Outlines, Input_DEM, CenterlineTributaryRatio, CenterlineTributaryAreaRatio, Centerlines)
    arcpy.Delete_management("in_memory") ### Empty the in_memory
    arcpy.AddMessage("Interpret ice thickness of extant glaciers...")
    arcpy.Delete_management("in_memory") ### Empty the in_memory
    Ice_Thickness_Volta (Centerlines, Input_DEM, Input_Glacier_Outlines, 900, slope_limit, min_slope, "false", Distance, "true", shearstress, volta_fl_points,
                     "true", "true", "", icethick)
    arcpy.Delete_management("in_memory") ### Empty the in_memory
    arcpy.AddMessage("Remove the ice thickness of extant glacier(s) from topography...")
    Adjust_DEM_by_Subtracting_Ice_Thickness (Input_DEM, icethick, adjustDEM)
    arcpy.Delete_management("in_memory") ### Empty the in_memory
    arcpy.AddMessage("Combine flowlines with centerlines...")
    Combine_Flowlines_with_Centerlines (flowlines, Centerlines, Input_Glacier_Outlines, watersheds, Input_DEM, 700, Combined_Flowline)
    arcpy.Delete_management("in_memory") ### Empty the in_memory
    arcpy.AddMessage("Palaeoglacier reconstruction...")
    ##let the program to recalculate the watershed based on the adjusted DEM
    PaleoIceReconstructionwithboundary(adjustDEM, Combined_Flowline, Distance, Input_Paleoice_boundary, shearstress, 50000, 200000, bFactorPolyfit, Output_Flowline_Points, 
                           Output_Ice_Surfaces, Output_Ice_Thickness)

else:
    arcpy.Delete_management("in_memory") ### Empty the in_memory
    arcpy.AddMessage("Palaeoglacier reconstruction...")
    PaleoIceReconstructionwithboundary(Input_DEM, Combined_Flowline, Distance, Input_Paleoice_boundary, shearstress, 50000, 200000, bFactorPolyfit, Output_Flowline_Points, 
                           Output_Ice_Surfaces, Output_Ice_Thickness)


arcpy.Delete_management("in_memory") ### Empty the in_memory
arcpy.Delete_management(Centerlines)
arcpy.Delete_management(icethick)
arcpy.Delete_management(adjustDEM)
arcpy.Delete_management(flowlines)
arcpy.Delete_management(watersheds)
arcpy.Delete_management(Combined_Flowline)
arcpy.Delete_management(volta_fl_points) 
