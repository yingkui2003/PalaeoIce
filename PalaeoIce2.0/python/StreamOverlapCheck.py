#-------------------------------------------------------------------------------
# Name: CenterlineFromStream.py
# Usage: CenterlineFromStream <InputDEM> <InputMoraineorCrossSection> <StreamFc> <SmoothlinePoints> 
# Purpose: This tool create smoothed flowline for glacial valleys based on the ArcGIS hydrological tools
# Author: Dr. Yingkui Li
# Created:     09/21-12/21/2020
# Department of Geography, University of Tennessee
# Knoxville, TN 37996
#-------------------------------------------------------------------------------

#import SharedFunctions  
from SharedFunctions import *  ## 
#------------------------------------------------------------------------------------------------------------
# This function Check if there are potential overlaped flowline sections within the flowline dataset
#------------------------------------------------------------------------------------------------------------
def Check_Flowline_Overlap(inflowline, outflowline):
    arcpy.AddMessage ("Checking and removing flowline overlaps...")
    lines=[]
    ids=[]
    flowline = "in_memory\\flowline"
    arcpy.CopyFeatures_management(inflowline, flowline)
    with arcpy.da.SearchCursor(flowline, ["SHAPE@","OID@"]) as cursor: ##This function needs to revise further to only check the flowline with the same GlacierID
        for row in cursor:
            lines.append(row[0])
            ids.append(row[1])
    del row, cursor

    #arcpy.AddMessage(str(len(lines)))

    removeIds = [] 
    for i in range(len(lines)):
        line = lines[i]
        for j in range (i+1, len(lines)):
            #if j != i: ##don't touch itself
            #within = (line.equals(lines[j]) or line.within(lines[j]))  
            within = line.within(lines[j])  
            equals = line.equals(lines[j])  
            if equals: ##== True: ##if the start point touches others
                #arcpy.AddMessage("Equals")
                #if j not in removeIds:
                removeIds.append(ids[i])
            elif within:## == True:
                #arcpy.AddMessage("Within but not equals")
                removeIds.append(ids[i])
                    
        #arcpy.AddMessage("Finish " + str(i))
    #arcpy.AddMessage(removeIds)
                    
    with arcpy.da.UpdateCursor(flowline, ["OID@"]) as cursor:
        for row in cursor:
            if row[0] in removeIds:
                #arcpy.AddMessage("Delete one")
                cursor.deleteRow()
    del row, cursor

    arcpy.CopyFeatures_management(flowline, outflowline)
    arcpy.Delete_management(flowline)
    return

##Main program
if __name__ == '__main__':
    # Script arguments
    inFlowlines = arcpy.GetParameterAsText(0)
    outFlowlines = arcpy.GetParameterAsText(1)

    Check_Flowline_Overlap(inFlowlines, outFlowlines)

    #Delete intermidiate data
    #arcpy.Delete_management("in_memory") ### Empty the in_memory


   
