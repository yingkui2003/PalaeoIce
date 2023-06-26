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
# This function Check if there are potential overlaped flowline sections within the flowline dataset
#------------------------------------------------------------------------------------------------------------
def Check_Flowline_Overlap(flowline, Distance):
    arcpy.AddMessage ("Checking and removing flowline overlaps...")
    lines=[]
    with arcpy.da.SearchCursor(flowline, ["SHAPE@","OID@"]) as flows: ##This function needs to revise further to only check the flowline with the same GlacierID
        for flow in flows:
            lines.append(flow[0])
    del flow, flows

    removeIds = [] 
    for i in range(len(lines)):
        line = lines[i]
        for j in range (len(lines)):
            if j != i: ##don't touch itself
                within = (line.equals(lines[j]) or line.within(lines[j]))  
                if within == True: ##if the start point touches others
                    arcpy.AddMessage("Within")
                    ##record the i as removed
                    removeIds.append(i)
    arcpy.AddMessage(removeIds)
                    

#------------------------------------------------------------------------------------------------------------
# This fuction is the main program to derive flowlines from stream network.
#------------------------------------------------------------------------------------------------------------
def Flowline_from_Stream_Network (InputDEM, InputMoraineorCrossSection, StreamThresholdKM2, TributaryThresholdKM2, TributaryRatio, StreamLine, outWatershed):

    GlacierID = "GlacierID" ##Add a GlacierID for each moriane or cross section

    cellsize = arcpy.GetRasterProperties_management(InputDEM,"CELLSIZEX")
    cellsize_int = int(cellsize.getOutput(0))
    arcpy.env.snapRaster = InputDEM

    StreamThreshold = int(float(StreamThresholdKM2) * 1e6 / (cellsize_int * cellsize_int))
    TributaryThreshold = int(float(TributaryThresholdKM2) * 1e6 / (cellsize_int * cellsize_int))

    ###Step 1: Stream network
    arcpy.AddMessage("Step 1: Stream extraction...")
 
    #Calculate Flowdirection
    fillDEM =Fill(InputDEM)  ##Fill the sink first
    fdir = FlowDirection(fillDEM,"NORMAL") ##Flow direction

    #Calculate Flowaccumulation
    facc = FlowAccumulation(fdir) ##Flow accmulation

    TmpStream = "in_memory\\TmpStream"
    moraineselected = "in_memory\\moraineselected"  ##Set a in_memory file for each moraine feature
    MaxFccTable = "in_memory\\MaxFccTable"
    CleanStream = "in_memory\\CleanStream"
    StreamClip = "in_memory\\StreamClip"
    Allstreams = "in_memory\\Allstreams"
    tmpoutStream = "in_memory\\tmpoutStream"
    smoothline = "in_memory\\smoothline"
    tmpws = "in_memory\\tmpws"
    tmpbuf = "in_memory\\tmpbuf"
    #smoothtolerance = cellsize_int * 10
    intersect_points = "in_memory\\intersect_points"

    outGreaterThan = Con(facc > StreamThreshold, 1,0)  ##Determine the highest flowaccumuation part
    outStreamLink = StreamLink(outGreaterThan, fdir)
    
    # Process: Stream to Feature
    StreamToFeature(outStreamLink, fdir, TmpStream, "SIMPLIFY")
    #arcpy.CopyFeatures_management(TmpStream, "c:\\test\\TmpStream02172023.shp")


    FcID = arcpy.Describe(InputMoraineorCrossSection).OIDFieldName

    arr=arcpy.da.FeatureClassToNumPyArray(InputMoraineorCrossSection, FcID)
    FIds = np.array([item[0] for item in arr])
    count = len(FIds)

    if count < 1:
        arcpy.AddMessage("There is no features in moraines or cross sections! Quit the program!")
        sys.exit()

    bPolyline = True
    fc_type = arcpy.Describe(InputMoraineorCrossSection).shapeType
    if fc_type == "Polygon": ##quit if not polyline features
        arcpy.AddMessage("Derive the flowlines within the polygons")
        bPolyline = False


    outflowline = arcpy.CreateFeatureclass_management("in_memory", "outflowline","POLYLINE","","","",InputMoraineorCrossSection)
    arcpy.AddField_management(outflowline, "Max_Max", "Long") 

    for imoraine in range (count):
        arcpy.AddMessage("Generating flowline(s) for glacial moraine "+str(imoraine + 1)+" of "+str(count) + " moraine(s)")

        query = FcID +" = "+str(FIds[imoraine])
        arcpy.Select_analysis(InputMoraineorCrossSection, moraineselected, query)

        if bPolyline:
            ##make a small buffer of the cross section to make sure the cross section get the highest fcc
            arcpy.Buffer_analysis(moraineselected, tmpbuf, (str(cellsize_int)+ " Meter"))
            bufID = arcpy.Describe(tmpbuf).OIDFieldName
            
            outZonalStatistics = ZonalStatistics(tmpbuf, bufID, facc, "MAXIMUM") #Find the maximum flowaccumulation point on the moriane feature
        else: ## for polygon input
            outZonalStatistics = ZonalStatistics(moraineselected, FcID, facc, "MAXIMUM")
            
        OutHighestFcc = Con(facc == outZonalStatistics,facc)  ##Determine the highest flowaccumuation part
        
        outSnapPour = SnapPourPoint(OutHighestFcc, facc, 0) ## Just create a pourpoint raster with the same extent of the input DEM
        
        #Calculate Watershed
        outWs = Watershed(fdir, outSnapPour)
        ConOutWs = Con(outWs >= 0, 1)  
        ##Boundary clean
        OutBndCln = BoundaryClean(ConOutWs)

        #arcpy.RasterToPolygon_conversion(ConOutWs, tmpws, "NO_SIMPLIFY", "VALUE")
        arcpy.RasterToPolygon_conversion(OutBndCln, tmpws, "NO_SIMPLIFY", "VALUE")
        #arcpy.RasterToPolygon_conversion(outWs, tmpws, "NO_SIMPLIFY", "VALUE")

        ##Need to process with the moraine bdoundary to extend the watershed if necessary
        
        #bMoraineExtension = True
        #if bMoraineExtension:
        if bPolyline:
            ##erase the moraines by watershed to see if there is some outside of the watershed
            arcpy.Erase_analysis(moraineselected, tmpws, "in_memory\\moraine_outside")
            outsideArr = arcpy.da.FeatureClassToNumPyArray( "in_memory\\moraine_outside", 'OID@')
            arcpy.AddMessage("moraine outside parts: " + str(len(outsideArr)))
            if len(outsideArr) > 0:
                ws_line = "in_memory\\ws_line"
                arcpy.PolygonToLine_management(tmpws, ws_line)
                arcpy.Append_management("in_memory\\moraine_outside", ws_line, "NO_TEST")                    
                arcpy.FeatureToPolygon_management(ws_line, "in_memory\\wsline_poly", "5 Meters", "NO_ATTRIBUTES")
                arcpy.Dissolve_management("in_memory\\wsline_poly", tmpws, '#', '#', 'SINGLE_PART')
        #else:
        #    ##if polygon, append the polygon into the tmpws and and then dissolve
        #    arcpy.Append_management(moraineselected, tmpws, "NO_TEST")
        #    arcpy.Dissolve_management(tmpws, "in_memory\\wsline_poly", '#', '#', 'SINGLE_PART')
        #    arcpy.CopyFeatures_management("in_memory\\wsline_poly", tmpws)
            
        ##delete potential small polygons
        poly_areaArr = arcpy.da.FeatureClassToNumPyArray(tmpws, 'SHAPE@AREA')
        poly_areas = np.array([item[0] for item in poly_areaArr])
        max_area = np.max(poly_areas)
        if len(poly_areas) > 1:
            with arcpy.da.UpdateCursor(tmpws, 'SHAPE@AREA') as cursor:
                for row in cursor:
                    if int(row[0]) < (max_area - 0.5):
                        arcpy.AddMessage("Delete one outline!")
                        cursor.deleteRow()     
            del cursor, row				

        #Get the watershed if required
        if outWatershed !="":
            if imoraine < 1: ##The first loop
                arcpy.CopyFeatures_management(tmpws, outWatershed)
            else:
                arcpy.Append_management(tmpws, outWatershed, "NO_TEST")        


        # Process: Extract by Mask
        try:
            ExtraFcc = ExtractByMask(facc,tmpws)
            
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
                #arcpy.CopyFeatures_management(TmpStream, "c:\\test\\TmpStream02072023.shp")
                ###This TmpStream already have a to_node in the attibute table, so that it can be used to make the decision
                ##the following is the new to remove and unnecessary lines
                lineArray = arcpy.da.FeatureClassToNumPyArray(TmpStream,['OID@','to_node','MAX'])
                tonode = np.array([item[1] for item in lineArray])
                uniquenode = np.unique(tonode)
                lineid = [] ##Record the id for deletion
                #arcpy.AddMessage("Checking tributary threshold...")
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
                #arcpy.AddMessage("Checking tributary ratio...")
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
                #arcpy.AddMessage("dissolve and clean lines...")

                ##Clean extralines based on the end points intersection 09/24/2020
                cleanextralineswithtopology(TmpStream,tmpoutStream, 'MAX')  ## clean the extra lines before dissolving

                arcpy.Dissolve_management(tmpoutStream, CleanStream, '#', 'MAX MAX', 'SINGLE_PART', 'UNSPLIT_LINES') 
            
                streamArr = arcpy.da.FeatureClassToNumPyArray(CleanStream, 'OID@')
                arcpy.AddMessage("The number of streams: " + str(len(streamArr)))
                if len(streamArr) > 0:
                    ##Remove big turns along the flowlines before smoothing it
                    newline = flowline_remove_bigturn(CleanStream, 120, cellsize_int)
                    ##Smooth the line
                    lineSmooth(newline, smoothline, "Max_Max", cellsize_int)
                    
                    arcpy.Append_management(smoothline, outflowline, "NO_TEST")
                else:
                    arcpy.AddMessage("No flowline is created for this feature. It seems that the threshold for a stream is too large!")
            else:
                arcpy.AddMessage("No flowline is created for this feature. It seems that the threshold for a stream is too large!")

        except:
            arcpy.AddMessage("There is an error in extract Facc! move to the next")


        '''
        arcpy.Intersect_analysis([Allstreams, moraineselected], intersect_points, "#", "#", "POINT")
        ##If more than one points, get the point with the highest facc
        
        
        #Snap_Pour Points ## the purpose is to create the same extent raster with only the highest fcc point
        try: 
            outSnapPour = SnapPourPoint(intersect_points, facc, 30) ## Just create a pourpoint raster with the same extent of the input DEM
            
            #Calculate Watershed
            outWs = Watershed(fdir, outSnapPour)
            ConOutWs = Con(outWs >= 0, 1)  
            ##Boundary clean
            OutBndCln = BoundaryClean(ConOutWs)

            #arcpy.RasterToPolygon_conversion(ConOutWs, tmpws, "NO_SIMPLIFY", "VALUE")
            arcpy.RasterToPolygon_conversion(OutBndCln, tmpws, "NO_SIMPLIFY", "VALUE")
            #arcpy.RasterToPolygon_conversion(outWs, tmpws, "NO_SIMPLIFY", "VALUE")

            ##Need to process with the moraine bdoundary to extend the watershed if necessary
            bMoraineExtension = True
            if bMoraineExtension:
                ##erase the moraines by watershed to see if there is some outside of the watershed
                arcpy.Erase_analysis(moraineselected, tmpws, "in_memory\\moraine_outside")
                outsideArr = arcpy.da.FeatureClassToNumPyArray( "in_memory\\moraine_outside", 'OID@')
                arcpy.AddMessage("moraine outside parts: " + str(len(outsideArr)))
                if len(outsideArr) > 0:
                    ws_line = "in_memory\\ws_line"
                    arcpy.PolygonToLine_management(tmpws, ws_line)
                    arcpy.Append_management("in_memory\\moraine_outside", ws_line, "NO_TEST")                    
                    arcpy.FeatureToPolygon_management(ws_line, "in_memory\\wsline_poly", "5 Meters", "NO_ATTRIBUTES")
                    arcpy.Dissolve_management("in_memory\\wsline_poly", tmpws, '#', '#', 'SINGLE_PART')

            ##delete potential small polygons
            poly_areaArr = arcpy.da.FeatureClassToNumPyArray(tmpws, 'SHAPE@AREA')
            poly_areas = np.array([item[0] for item in poly_areaArr])
            max_area = np.max(poly_areas)
            if len(poly_areas) > 1:
                with arcpy.da.UpdateCursor(tmpws, 'SHAPE@AREA') as cursor:
                    for row in cursor:
                        if int(row[0]) < (max_area - 0.5):
                            cursor.deleteRow()     
                del cursor, row				

            #Get the watershed if required
            if outWatershed !="":
                if imoraine < 1: ##The first loop
                    arcpy.CopyFeatures_management(tmpws, outWatershed)
                else:
                    arcpy.Append_management(tmpws, outWatershed, "NO_TEST")        
                        
            #Use the watershed to clip the cleaned streams to derive the flowline
            arcpy.Clip_analysis(TmpStream, tmpws, StreamClip)
            ##Clean extralines based on the end points intersection 09/24/2020
            cleanextralineswithtopology(StreamClip,tmpoutStream, 'MAX')  ## clean the extra lines before dissolving

            arcpy.Dissolve_management(tmpoutStream, CleanStream, '#', 'MAX MAX', 'SINGLE_PART', 'UNSPLIT_LINES') 
            
            streamArr = arcpy.da.FeatureClassToNumPyArray(CleanStream, 'OID@')
            arcpy.AddMessage("The number of streams: " + str(len(streamArr)))
            if len(streamArr) > 0:
                ##Remove big turns along the flowlines before smoothing it
                newline = flowline_remove_bigturn(CleanStream, 120, cellsize_int)
                ##Smooth the line
                lineSmooth(newline, smoothline, "Max_Max", cellsize_int)
                
                arcpy.Append_management(smoothline, outflowline, "NO_TEST")
            else:
                arcpy.AddMessage("No flowline is created for this feature. It seems that the threshold for a stream is too large!")
        except:
            arcpy.AddMessage("No flowline is created for this feature!")
        '''     
    if  bPolyline == False: ##if polygon as the input, clip the flowlines within the polygon
        arcpy.Clip_analysis(outflowline, InputMoraineorCrossSection, "in_memory\\flowline_clip")
        arcpy.CopyFeatures_management("in_memory\\flowline_clip", outflowline)        


    ##make sure there are flowlines created
    countResult = arcpy.GetCount_management(outflowline)
    count = int(countResult.getOutput(0))
    if count < 1:
        arcpy.AddMessage("No flowlines are created for this set of moraine features!!")
        sys.exit()
    ##Merge flowline and add GlacierID
    Check_If_Flip_Line_Direction (outflowline, fillDEM) ##use fillDEM because the orginal DEM may have problems
    #arcpy.CopyFeatures_management(outflowline, "")
    Merge_and_Add_GlacierID_by_Topology (outflowline, "Max_Max", GlacierID, "MergeID", StreamLine)
    #arcpy.CopyFeatures_management(outflowline, StreamLine)

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


   
