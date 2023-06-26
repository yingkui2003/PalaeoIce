# ---------------------------------------------------------------------------
# AdjustDEM.py
# This python code is to subtract the ice thickness from DEM to generate an ice free DEM for paleoglacier reconstruction
# Usage: AdjustDEM <DEM> <Thickness> <AdjustDEM>
# Created by:
# Yingkui Li
# Department of Geography
# University of Tennessee
# Knoxville, TN 37996
# 12/10/2019
# ---------------------------------------------------------------------------

# Import arcpy module
import arcpy
from arcpy.sa import *

#------------------------------------------------------------------------------------------------------------
# This is the main function for adjusting DEM. This function can be call by the overall model
#------------------------------------------------------------------------------------------------------------
def Adjust_DEM_by_Subtracting_Ice_Thickness (DEM, Thickness, AdjDEM):
        
    # set the extent of the raster to DEM
    oldextent = arcpy.env.extent
    arcpy.env.extent = DEM

    # replace nodata to zero from the thickness raster
    outCon = Con(IsNull(Thickness), 0, Thickness)

    # subtract thickness from DEM
    outMinus = DEM - outCon

    # set zero to NoData again for adjusted DEM
    outSetNull = SetNull(outMinus, outMinus, "Value = 0")

    #save the adjusted DEM
    outSetNull.save(AdjDEM)

    # Reset the default extent for raster processing
    arcpy.env.extent = oldextent

#------------------------------------------------------------------------------------------------------------
# Start the main program
#------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # Script arguments
    DEM = arcpy.GetParameterAsText(0)
    Thickness = arcpy.GetParameterAsText(1)
    AdjDEM = arcpy.GetParameterAsText(2)

    Adjust_DEM_by_Subtracting_Ice_Thickness (DEM, Thickness, AdjDEM)
