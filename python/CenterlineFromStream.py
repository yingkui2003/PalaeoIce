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
# This fuction is the main program to derive flowlines from stream network.
#------------------------------------------------------------------------------------------------------------
def Flowline_from_Stream_Network (InputDEM, InputMoraineorCrossSection, StreamThresholdKM2, TributaryThresholdKM2, TributaryRatio, StreamLine, outWatershed):

    GlacierID = "GlacierID" ##Add a GlacierID for each moriane or cross section

    cellsize = arcpy.GetRasterProperties_management(InputDEM,"CELLSIZEX")
    cellsize_int = int(cellsize.getOutput(0))

    StreamThreshold = int(float(StreamThresholdKM2) * 1e6 / (cellsize_int * cellsize_int))
    TributaryThreshold = int(float(TributaryThresholdKM2) * 1e6 / (cellsize_int * cellsize_int))

    #Calculate Flowdirection
    fillDEM =Fill(InputDEM)  ##Fill the sink first
    fdir = FlowDirection(fillDEM,"NORMAL") ##Flow direction

    #Calculate Flowaccumulation
    facc = FlowAccumulation(fdir) ##Flow accmulation

    FcID = arcpy.Describe(InputMoraineorCrossSection).OIDFieldName

    arr=arcpy.da.FeatureClassToNumPyArray(InputMoraineorCrossSection, FcID)
    FIds = np.array([item[0] for item in arr])
    count = len(FIds)

    if count < 1:
        arcpy.AddMessage("There is no features in moraines or cross sections")
        sys.exit()

    moraineselected = "in_memory\\moraineselected"  ##Set a in_memory file for each moraine feature
    MaxFccTable = "in_memory\\MaxFccTable"
    CleanStream = "in_memory\\CleanStream"
    tmpoutStream = "in_memory\\tmpoutStream"
    TmpStream = "in_memory\\TmpStream"
    smoothline = "in_memory\\smoothline"
    tmpws = "in_memory\\tmpws"
    tmpbuf = "in_memory\\tmpbuf"
    smoothtolerance = cellsize_int * 10



    for imoraine in range (count):
        arcpy.AddMessage("Generating flowline(s) for glacial moraine "+str(imoraine + 1)+" of "+str(count))

        query = FcID +" = "+str(FIds[imoraine])

        arcpy.Select_analysis(InputMoraineorCrossSection, moraineselected, query)
        ##make a small buffer of the cross section to make sure the cross section get the highest fcc
        arcpy.Buffer_analysis(moraineselected, tmpbuf, (str(cellsize_int/2)+ " Meter"))
        bufID = arcpy.Describe(tmpbuf).OIDFieldName
        
        outZonalStatistics = ZonalStatistics(tmpbuf, bufID, facc, "MAXIMUM") #Find the maximum flowaccumulation point on the moriane feature
        OutHighestFcc = Con(facc == outZonalStatistics,facc)  ##Determine the highest flowaccumuation part

        #Snap_Pour Points ## the purpose is to create the same extent raster with only the highest fcc point
        outSnapPour = SnapPourPoint(OutHighestFcc, facc, 0) ## Just create a pourpoint raster with the same extent of the input DEM
        
        #Calculate Watershed
        outWs = Watershed(fdir, outSnapPour)
        #Watershed to outline
        if outWatershed !="":
            arcpy.RasterToPolygon_conversion(outWs, tmpws, "NO_SIMPLIFY", "VALUE")
            if imoraine < 1: ##The first loop
                arcpy.CopyFeatures_management(tmpws, outWatershed)
            else:
                arcpy.Append_management(tmpws, outWatershed, "NO_TEST")        
                    

        # Process: Extract by Mask
        ExtraFcc = ExtractByMask(facc,outWs)
        # Process: Greater Than
        outGreaterThan = Con(ExtraFcc > StreamThreshold, 1,0)  ##Determine the highest flowaccumuation part
        #need to check if outGreaterThan has the 1 values. If not, no stream will be created
        MaxRasterValue = int((arcpy.GetRasterProperties_management(outGreaterThan, "MAXIMUM").getOutput(0)))
        if MaxRasterValue > 0:
            # Process: Stream Link
            outStreamLink = StreamLink(outGreaterThan, fdir)
            
            # Process: Stream to Feature
            StreamToFeature(outStreamLink, fdir, TmpStream, "SIMPLIFY")

            # Process: Zonal Statistics as Table
            ZonalStatisticsAsTable(outStreamLink, "VALUE", ExtraFcc, MaxFccTable, "DATA", "MAXIMUM")

            # Process: Join Field
            arcpy.JoinField_management(TmpStream, "grid_code", MaxFccTable, "Value", "MAX")  ##Join to get the flow accumulation value

            ###This TmpStream already have a to_node in the attibute table, so that it can be used to make the decision
            ##the following is the new to remove and unnecessary lines
            lineArray = arcpy.da.FeatureClassToNumPyArray(TmpStream,['OID@','to_node','MAX'])
            tonode = np.array([item[1] for item in lineArray])
            uniquenode = np.unique(tonode)
            lineid = [] ##Record the id for deletion
            arcpy.AddMessage("checking tributary threshold...")
            for i in range(len(uniquenode)):
                selArr = lineArray[tonode == uniquenode[i]]
                fcclist = []
                if len(selArr) > 1: ##Sometimes having more than two end points
                    for j in range(len(selArr)):
                        fcclist.append(selArr[j][2])

                    numselected = len(fcclist)
                    while numselected > 1:
                        minfcc = min(fcclist)
                        if minfcc < TributaryThreshold:
                            for j in range(len(selArr)):
                                if selArr[j][2] == minfcc: ##Remove the smaller fcc one
                                    lineid.append(selArr[j][0])
                                    fcclist.pop(j)
                                    selArr = np.delete(selArr, j)##Remove this one and loop to delete others
                                    numselected = len(selArr)
                                    break ##Only remove one each time

                        else: ##quit the loop if all minfcc are larger than the theshold?? Try to remove more based on the ratio
                            break

            ##Remove features based on the tributary ratio
            arcpy.AddMessage("checking tributary ratio...")
            for i in range(len(uniquenode)):
                selArr = lineArray[tonode == uniquenode[i]]
                fcclist = []
                if len(selArr) > 1: ##Sometimes having more than two end points
                    for j in range(len(selArr)):
                        fcclist.append(selArr[j][2])

                    sumfcc = sum(fcclist)
                    fccratio = [x / float(sumfcc) for x in fcclist]

                    for j in range(len(fccratio)):
                        if fccratio[j] < float(TributaryRatio): ##Remove the smaller fcc one
                            lineid.append(selArr[j][0])
            
            ##Delete the line marked for deletion
            with arcpy.da.UpdateCursor(TmpStream, "OID@") as cursor:
                for row in cursor:
                    if int(row[0]) in lineid:
                        cursor.deleteRow()     
            del cursor, row				

            ##Dissolve the line file
            arcpy.AddMessage("dissolve and clean lines...")

            ##Clean extralines based on the end points intersection 09/24/2020
            cleanextralineswithtopology(TmpStream,tmpoutStream, 'MAX')  ## clean the extra lines before dissolving

            arcpy.Dissolve_management(tmpoutStream, CleanStream, '#', 'MAX MAX', 'SINGLE_PART', 'UNSPLIT_LINES') 

            ##Remove big turns along the flowlines before smoothing it
            newline = flowline_remove_bigturn(CleanStream, 120, cellsize_int)
            ##Smooth the line
            lineSmooth(newline, smoothline, "Max_Max", cellsize_int)
            
            
            if imoraine < 1: ##The first loop
                arcpy.CopyFeatures_management(smoothline, "in_memory\\outflowline")
            else:
                arcpy.Append_management(smoothline, "in_memory\\outflowline", "NO_TEST")
        else:
            arcpy.AddMessage("No flowline is created for this feature. It seems that the threshold for a stream is too large for this feature.")

    ##make sure there are flowlines created
    countResult = arcpy.GetCount_management("in_memory\\outflowline")
    count = int(countResult.getOutput(0))
    if count < 1:
        arcpy.AddMessage("No flowlines are created for this set of features!!")
        sys.exit()
    ##Merge flowline and add GlacierID
    Check_If_Flip_Line_Direction ("in_memory\\outflowline", fillDEM) ##use fillDEM because the orginal DEM may have problems

    Merge_and_Add_GlacierID_by_Topology ("in_memory\\outflowline", "Max_Max", GlacierID, "MergeID", StreamLine)

##Main program
if __name__ == '__main__':
    # Script arguments
    InputDEM = arcpy.GetParameterAsText(0)
    InputMoraineorCrossSection = arcpy.GetParameterAsText(1)
    StreamThresholdKM2 = arcpy.GetParameter(2)
    TributaryThresholdKM2 = arcpy.GetParameter(3)
    TributaryRatio = arcpy.GetParameter(4) ##Define a tributary ratio
    StreamLine = arcpy.GetParameterAsText(5)
    outWatershed = arcpy.GetParameterAsText(6)

    Flowline_from_Stream_Network (InputDEM, InputMoraineorCrossSection, StreamThresholdKM2, TributaryThresholdKM2, TributaryRatio, StreamLine, outWatershed)

    ##Delete intermidiate data
    arcpy.Delete_management("in_memory") ### Empty the in_memory


   
