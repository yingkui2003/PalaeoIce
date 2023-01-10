#-------------------------------------------------------------------------------
# Name: Paleo Ice Reconstruction
# Purpose: This tool calculates the ice thickness on flowlines based on the Excel flowline model introduced by the Benn and Houlton (2010)
# and the python model, GLaRe, developed by Pellitero et al.(2016). However, both of the above models require manually assign shear stress
# and F factors. It requires a lot of time and efforts to conduct a paleo ice reconstrcution work.
# 
# The new developed tool will automatically adjust shear stress and F factors based on the DEM and target features. For shear stress, the tool
# assumes one value for the whole glacier based on the recommendation of many previous studies. The shear stress is first derived based on a
# revised shearstress code from James et al. (2016) based on the elevation distrbution of ice surface and then adjusted to reach the best fit to
# the assigned target features (if assigned) that represent trimlines and boundaries of the paleo glaciaer.
#
# This tool also automatically adjust the F factor (shape factor) based on the cross sections along the flowlines. There are two options to derive
# the F factor: one is based on the cross section area, ice-contact perimeter, and ice thickness; the other is based on the fit of the polynomial
# function introduced by Li et al (2012). This tool uses a maximum width to prevent the error long cross sections in some sections when
#  the direction is not paralell the valley direction and where the tributary valley joins the main valley. This tool also applied EU
# allocation method to make sure the cross section does not extend to the tributaries (May need to check if it is necessary because the width is
# already used to constrain the exent of the cross section).
# 
# Author:      Yingkui Li
# Created:     08/30/2020 to 12/16/2020
# Copyright:   (c) Yingkui Li 2020
# Department of Geography, University of Tennessee
# Knoxville, TN 37996, USA
#-------------------------------------------------------------------------------
#from __future__ import division
from SharedFunctions import *  ## 

#------------------------------------------------------------------------------------------------------------
# This fuction is the whole process to reconstruct paleoice based on DEM, input flowlines, target features, and default shear stress
#------------------------------------------------------------------------------------------------------------
def PaleoIceReconstruction(BedDEM, inputflowline, Distance, inwatershed, TargetFeatures, shearstress, min_ss, max_ss, bFactorPolyfit, outpoints, outIcePolys, outIceSurfaces, outIceThickness):
    arcpy.env.extent = BedDEM
    arcpy.env.cellSize = BedDEM
    arcpy.env.snapRaster = BedDEM ##setup snap raster
    
    GlacierID = "GlacierID" ##This is an ID field in inputflowline to identify the flowline(s) for each glacier (maybe connected with multiple flowlines)


    t=arcpy.Describe(inputflowline).shapeType
    if not (t == "Polyline"): ##quit if not polyline features
        arcpy.AddError("You did not define a flowline polyline input")
        arcpy.GetMessage(0)
        exit()

    #if TargetFeatures != "":  ##The target features can be point, line and polygon type: Yingkui Li 09/12/2022
    #    targetType = arcpy.Describe(inputflowline).shapeType
    #    if not (targetType == "Polyline"): ##quit if not polyline features
    #        arcpy.AddMessage("The target features are not polyline input! Assuming no target features")
    #        TargetFeatures = ""

    arcpy.AddMessage("Checking the quality of flowline(s)...")
    number=arcpy.GetCount_management(inputflowline)
    if int(number.getOutput(0))==1:
        arcpy.AddMessage("The flowline is good to be used in the toolbox")
    elif int(number.getOutput(0))==0:
        arcpy.AddWarning("There is not any line in the selected file")
        quit()
    else:
        Check_Flowline_Connectivity(inputflowline, Distance)

    #Check and flip lines if necessary    
    arcpy.AddMessage("Checking flowline direction...")
    Check_If_Flip_Line_Direction(inputflowline, BedDEM)


    ####Flow direction and accumulation analysis
    cellsize = arcpy.GetRasterProperties_management(BedDEM,"CELLSIZEX")
    cellsize_float = float(cellsize.getOutput(0)) # use float cell size

    #arcpy.env.extent = BedDEM
    #arcpy.env.cellSize = BedDEM

    burninDEM = BedDEM - Power (cellsize_float / (cellsize_float + EucDistance(inputflowline) ), 2 ) * 10 ##Burn in the DEM to make sure the flow pass through the flowline start points
    ##Start to delineate the watershed
    #Hydro analysis
    fillDEM =Fill(burninDEM)  ##Fill the sink first
    fdir = FlowDirection(fillDEM,"NORMAL") ##Flow direction
        
    facc = FlowAccumulation(fdir) ##Flow accmulation    

    #Make a copy of the flowline
    flowlines = "in_memory\\flowlines"
    flowline = "in_memory\\flowline"
    flowline3dpoints = "in_memory\\flowline3dpoints"
    selflowline3dpoints = "in_memory\\selflowline3dpoints"
    singepoint = "in_memory\\singepoint"
    ws = "in_memory\\ws"
    wsflowpoints = "in_memory\\wsflowpoints"
    icepolys = "in_memory\\icepolys"
    mainflowline = "in_memory\\mainflowline"
    flowlinecopy = "in_memory\\flowlinecopy"
    allicepolys = "in_memory\\allicepolys"
    singeflowpnts = "in_memory\\singeflowpnts"
    icesurs = arcpy.env.scratchFolder + "\\r" + "icesurs" ##the inmemory does not work for raster
    oldsurface = arcpy.env.scratchFolder + "\\r" + "oldsurface" ##the inmemory does not work for raster

    exist_fields = [f.name for f in arcpy.ListFields(inputflowline)] #List of current field names in outline layer
    if GlacierID not in exist_fields:
        arcpy.AddMessage("Assigning Glacier ID...")
        Add_GlacierID_by_Touches (inputflowline, GlacierID, flowlines)
    else:
        arcpy.CopyFeatures_management(inputflowline, flowlines)


    ###The process to ordering the flowlines
    #Obtain the height info for the start of each flowline
    height=[]
    lineid = []
    glaciers = []

    with arcpy.da.SearchCursor(flowlines, ['SHAPE@', GlacierID]) as cursor:
        i = 0
        for row in cursor:
            Startpoint = row[0].firstPoint
            coord= str(Startpoint.X)+" "+str(Startpoint.Y)
            Cellvalue=arcpy.GetCellValue_management(BedDEM, coord)
            Startpoint.Z=Cellvalue.getOutput(0)
            height.append(Startpoint.Z)
            lineid.append(i)
            glaciers.append(row[1])
            i += 1
    del row, cursor
    
    ##Order the line geometries in the list
    arcpy.AddMessage("Ordering flowlines...")
    arcpy.AddField_management(flowlines,"ProcessID","LONG",6)

    order = sorted(range(len(height)), key=lambda k: height[k])  ##order is the ID

    with arcpy.da.UpdateCursor(flowlines, "ProcessID") as cursor: ##Fix the assigning order issue
        i = 0
        for row in cursor:
            row[0] = order.index(i)
            cursor.updateRow(row)
            i += 1
    del row, cursor

    OriginalFID = []
    Points = []
    processPID = []
    glaIds = []
    p = 0

    geometry = arcpy.CopyFeatures_management(flowlines, arcpy.Geometry()) ##Can revise here to include glacier ID as a list?????

    for i in order:
        Length = geometry[i].length
        Length = int(Length)
        try:
            rlist = xrange(0, Length, Distance)
        except: ##python 3 xrange is replaced by range
            rlist = range(0, Length, Distance)
        
        for j in rlist:
            Points.append(geometry[i].positionAlongLine(j))
            OriginalFID.append(i)
            glaIds.append(glaciers[i])
            processPID.append(p)
        p += 1

    PointsFile = "in_memory\\PointsFile"
    arcpy.CopyFeatures_management(Points, PointsFile)

    ##Extract Values and xy coordinates from Elevation Raster
    ExtractValuesToPoints(PointsFile, BedDEM, flowline3dpoints,"INTERPOLATE", "VALUE_ONLY")

    ##Add shear stress value
    arcpy.AddField_management(flowline3dpoints,"SSTRESS","Long",7)
    arcpy.CalculateField_management(flowline3dpoints,"SSTRESS", shearstress)  ####Need to automatic determine the shearstress!!!

    ##Add original FID Field
    arcpy.AddField_management(flowline3dpoints, 'OFID', 'Long', 6) ##OFID is not FID, but it is related to FID
    arcpy.AddField_management(flowline3dpoints, 'ProcessID', 'Long', 6) ##OFID is not FID, but it is related to FID
    arcpy.AddField_management(flowline3dpoints, GlacierID, 'Long', 6) ##OFID is not FID, but it is related to FID

    PointsCursor = arcpy.UpdateCursor(flowline3dpoints, ['OFID','ProcessID', GlacierID])
    i = 0
    for row in PointsCursor:
        row.OFID = OriginalFID[i]
        row.ProcessID = processPID[i]
        row.GlacierID = glaIds[i]
        PointsCursor.updateRow(row)
        i+=1
    del row, PointsCursor  

    arcpy.DeleteField_management(flowline3dpoints,["POINT_Z","POINT_M"])

    fieldlist=[field.name for field in arcpy.ListFields(flowline3dpoints)]
    if not("ffactor" in fieldlist):
        arcpy.AddField_management(flowline3dpoints, "ffactor", "FLOAT", 5, 4)
        arcpy.CalculateField_management(flowline3dpoints,"ffactor", 1.0)                   ###Need to automatically determine the F factors

    if not("ice" in fieldlist):
        arcpy.AddField_management(flowline3dpoints, "ice", "DOUBLE", 6, 2)

    if not("thick" in fieldlist):
        arcpy.AddField_management(flowline3dpoints,"thick", "Double", 6, 2)

    exist_fields = [f.name for f in arcpy.ListFields(flowlines)] #List of current field names in outline layer
    line_fields = ["line_id"]
    for field in line_fields:
        if field not in exist_fields:
            arcpy.AddField_management(flowlines, field, "LONG")
    arcpy.CalculateField_management(flowlines,"line_id",str("!"+str(arcpy.Describe(flowlines).OIDFieldName)+"!"),"PYTHON_9.3")

    ##loop for the flowline of each glacier ID
    lineArray = arcpy.da.FeatureClassToNumPyArray(flowlines,GlacierID)
    iceID = np.array([item[0] for item in lineArray])
    uniqueiceID = np.unique(iceID)

    minArea = cellsize_float * cellsize_float * 5 ##5 pixels 

    for gid in range(len(uniqueiceID)):
        query = GlacierID + " = " + str(uniqueiceID[gid])
        arcpy.AddMessage("Processing #" + str(gid+1) +"/" + str(len(uniqueiceID)) + " of reconstructed glaciers...")                                                                                       
        arcpy.Select_analysis (flowlines, flowline, query)

        ###Select the flowline points corresponding to the flowlines
        arcpy.Select_analysis (flowline3dpoints, selflowline3dpoints, query)

        ###Merge the branches and start to run the ice thickness calculation 
        arcpy.Near_analysis (selflowline3dpoints, selflowline3dpoints, Distance * 2) ##setup a search radius for near analysis (may be problematic to choose the upstream near point; should select the nearest two points and use the average value!! or keep the intersect points as part of the flowpoints)
        arcpy.AddXY_management(selflowline3dpoints)

        ##Get the order of the processing ID
        flowlineArray = arcpy.da.FeatureClassToNumPyArray(flowline,"ProcessID")
        orderID = np.array([item[0] for item in flowlineArray])
        sortID = np.sort(orderID, axis=0) ##sort the uniqueID from small to large

        ##only select the main flowline to delieate the watershed
        queryMainFlowline = "ProcessID = "+ str(sortID[0])  
        arcpy.Select_analysis(flowline, mainflowline, queryMainFlowline)
        arcpy.FeatureVerticesToPoints_management(mainflowline, singepoint, "START") ###the flowline direction is from the downslope to upper slope or need to check the direction

        if inwatershed == "":
            startperppolygon = Startpoint_perpendiculars(mainflowline, int(cellsize_float), 300) ##Just create a short distance perpendicular polygon (need polygon because sometimes line does not cross the highest Facc) from the start point to get the highest Facc point, to avoid the snap pointpoints downstream
            outZonalStatistics = ZonalStatistics(startperppolygon, arcpy.Describe(startperppolygon).OIDFieldName, facc, "MAXIMUM") #Find the maximum flowaccumulation point on the moriane feature
            OutHighestFcc = Con(facc == outZonalStatistics,facc)  ##Determine the highest flowaccumuation part
            #Snap_Pour Points ## the purpose is to create the same extent raster with only the highest fcc point
            outSnapPour = SnapPourPoint(OutHighestFcc, facc, 0) ## Just create a pourpoint raster with the same extent of the input DEM
            outpntWs = Watershed(fdir, outSnapPour)

            ##make sure to have the watershed and iceboundary intersection that can cover all flowline points
            points_count = int(arcpy.GetCount_management(selflowline3dpoints).getOutput(0))
            arcpy.RasterToPolygon_conversion(outpntWs, ws)
            inside_point_count = 0
            i = 1
            while (points_count > inside_point_count): 
                arcpy.Clip_analysis (selflowline3dpoints, ws, "in_memory\\selectedflpoints")
                inside_point_count = int(arcpy.GetCount_management("in_memory\\selectedflpoints").getOutput(0))
                if inside_point_count < points_count:
                    #arcpy.AddMessage("ws does not include all flowline points!!")
                    outSnapPour = SnapPourPoint(singepoint, facc, 100*i) ##each time increase 100 m downstream
                    outpntWs = Watershed(fdir, outSnapPour)
                    arcpy.RasterToPolygon_conversion(outpntWs, ws)

                ##add the maximum loop controls
                if i >= 10:
                    #arcpy.AddMessage("Cannot find the watershed that include all flowline points! use the 1 km buffer of the flowline instead")
                    arcpy.Buffer_analysis(flowline, "in_memory\\flowlinebuf", "1000 Meter", "#", "#", "ALL")
                    ##get the minmum elevation of the startperppolygon
                    extDEM = ExtractByMask(fillDEM, startperppolygon)
                    min_elev = float(arcpy.GetRasterProperties_management(extDEM,"MINIMUM").getOutput(0))
                    ##Extract the DEM within the buffer
                    extDEMbuf = ExtractByMask(fillDEM, "in_memory\\flowlinebuf")
                    ##Con for elevation > min_elev
                    conDEMbuf = Con(extDEMbuf >= min_elev, 1)
                    ##Convert conDEMbuf to polygon
                    arcpy.RasterToPolygon_conversion(conDEMbuf, "in_memory\\conDEMbufpoly")
                    arcpy.MultipartToSinglepart_management("in_memory\\conDEMbufpoly", "in_memory\\singleDEMbufpolys")
                    arcpy.SpatialJoin_analysis("in_memory\\singleDEMbufpolys", flowline, ws, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
                    break

                i += 1

            #remove the potential spurious polygons
            with arcpy.da.UpdateCursor(ws, "SHAPE@AREA") as cursor:        
                for row in cursor:
                    if row[0] < minArea:
                        cursor.deleteRow()  ##Just remove the max elevation and keep the lowest elevation
            del cursor, row

        else:
            ##Use spatial join to select the watershed corresponding to the flowline
            arcpy.SpatialJoin_analysis(inwatershed, flowline, ws, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")

        ##Make sure if there is features used for calibrate the paleoglacier reconstruction
        bFeatureComparison = False
        if TargetFeatures != "":
            ##Clip target features within the watershed ## only consider the features in the watershed
            arcpy.Clip_analysis (TargetFeatures, ws, "in_memory\\clipedFc")
            targetcountResult = arcpy.GetCount_management("in_memory\\clipedFc")
            targetcount = int(targetcountResult.getOutput(0))
            if targetcount > 0:
                bFeatureComparison = True
                arcpy.AddField_management("in_memory\\clipedFc", "Watershed", "Long")
                arcpy.CalculateField_management("in_memory\\clipedFc", "Watershed", 1)
        

        cross_section_pnts = cross_section_points(selflowline3dpoints, flowline, ws, BedDEM, 1600, 60)

        ##Ice thickness calculation along the flowlines
        Ice_Thickness_Calculation (selflowline3dpoints)
        
        ss0 = float(shearstress) ##each glacier only have one shear stress for all flowlines and flow points

        ##Add three lists to record ss, distance to target features, and outside_percentange of the target features
        ss_list = []
        distance_list = []
        
        ss_ratio0 = 0
        bStop = False
        bAdjustSS = False
        Dist2Target = 10000
        icepolyold = "in_memory\\icesurold"
        selflowline3dpointsold = "in_memory\\selflowline3dpointsold"
        icepolyold2 = "in_memory\\icesurold2"
        selflowline3dpointsold2 = "in_memory\\selflowline3dpointsold2"
        nloop = 0
        while (bStop == False) and (nloop < 10):
            icesurface, icepoly = IceSurfaceCalculation_with_crosssection_pnts (selflowline3dpoints, BedDEM, ws, cross_section_pnts, bFeatureComparison, "in_memory\\clipedFc", "ice", "in_memory\\icesurface", "in_memory\\icepoly")

            ##Calculate the dist to Target
            if bFeatureComparison == True:
                distance, outside_percent = Distance_to_Target_Features (icepoly, "in_memory\\clipedFc", "Watershed", ws, cellsize_float)  ##Need more work on this function???
                arcpy.AddMessage("The average distance to target features is: " + str(distance))

                ##Append the ss, distance, and outside_percent for the current loop to the three lists
                ss_list.append (ss0)
                if outside_percent > 0.5:
                    distance_list.append(-distance)
                else:
                    distance_list.append(distance)
                    
                if distance < Dist2Target: 
                    Dist2Target = distance
                    ##copy icesurface and icepoly
                    arcpy.CopyRaster_management(icesurface, oldsurface)
                    arcpy.CopyFeatures_management(icepoly, icepolyold)
                    arcpy.CopyFeatures_management(selflowline3dpoints, selflowline3dpointsold)

                else: ##Stop the ss calculation by glacier size and height distribution!! because the distance already reach to the minimum value; and start to refine the ss
                    arcpy.AddMessage("Reached the minimum distance of " + str(Dist2Target) + " to target features...")

                    arcpy.CopyRaster_management(oldsurface, icesurface)
                    arcpy.CopyFeatures_management(icepolyold, icepoly)
                    arcpy.CopyFeatures_management(selflowline3dpointsold, selflowline3dpoints)

                    if bAdjustSS == True: ##if the adjust SS case, still get even larger distance, stop the loop
                        break

                    ## set bAdjustSS to True
                    bAdjustSS = True 

            if bAdjustSS == False:
                #ss = shear_stress_calculation (mainflowline, icepoly, icesurface)
                ss = shear_stress_calculation(mainflowline, icepoly, icesurface, min_ss, max_ss)
            else:
                if abs(distance_list[-2] - distance_list[-1]) > 0: 
                    m = (ss_list[-2] - ss_list[-1])/(distance_list[-2] - distance_list[-1])
                    ss = ss_list[-1] - m * distance_list[-1]
                else:
                    ss = (ss_list[-2] + ss_list[-1]) / 2
                if ss < min_ss: ##If the shear stress is negative, using the default value
                    ss = min_ss
                if ss > max_ss:
                    ss = max_ss
               
            arcpy.AddMessage("The calculated shear stress is:" + str(ss) + " and the difference with the previous value is " + str(abs(ss-ss0)/ss0))
            ss_ratio = abs(ss - ss0)/ss0
            if abs(ss_ratio - ss_ratio0)< 0.005 or (ss_ratio < 0.01):  ##assuming ss always increase, to stop if decreasing; ##use 1.0% change as the threshold to stop the loop
            #if (abs(ss - ss0)/ss0 < 0.01):  ##assuming ss always increase, to stop if decreasing; ##use 0.01% change as the threshold to stop the loop
                break
                
            ##update ss value
            arcpy.CalculateField_management(selflowline3dpoints,"SSTRESS",ss) ##update ss to the flowline points

            if bFactorPolyfit == True:
                arcpy.AddMessage("Deriving F factor based on the Polyfit of the cross section...")
                cross_section_width = AdjustFfactor_Ployfit_with_cross_section_pnts (selflowline3dpoints, cross_section_pnts, "ice", icepoly) 
            else:
                arcpy.AddMessage("Deriving F factor from cross section...")
                #Try to use the curve fitting method to derive the F factor
                cross_section_width = AdjustFfactor_with_cross_section_pnts (selflowline3dpoints, cross_section_pnts, "ice", icepoly) 

            ##Recalculate the ice thickness for the flowpoints based on adjusted F and ss parameters
            Ice_Thickness_Calculation (selflowline3dpoints)

            ss0 = ss
            ss_ratio0 = ss_ratio
            nloop += 1
            
        if gid == 0: #The first time
            arcpy.CopyFeatures_management(selflowline3dpoints, outpoints)
            arcpy.CopyFeatures_management(icepoly, allicepolys)
            arcpy.CopyFeatures_management(ws, "in_memory\\watersheds")
            arcpy.CopyRaster_management(icesurface, icesurs)
        else:
            arcpy.Append_management(selflowline3dpoints, outpoints, "NO_TEST" )
            arcpy.Append_management(icepoly, allicepolys, "NO_TEST")
            arcpy.Append_management(ws, "in_memory\\watersheds")
            arcpy.Mosaic_management(icesurface, icesurs, "MEAN","","", "", "", "", "")

    ##Smooth the polgyons and extract the ice surface within the ice polygon
    arcpy.AddMessage("Generating final outputs...")

    arcpy.Dissolve_management(allicepolys, 'in_memory\\disolveicepoly', '#', '#', 'SINGLE_PART', '#')
    arcpy.cartography.SmoothPolygon("in_memory\\disolveicepoly", outIcePolys, "PAEK", cellsize_float * 5)

    '''
    ###Output method 1: directly output mosaic ice surface
    extSurface = ExtractByMask(icesurs, 'in_memory\\disolveicepoly')
    arcpy.CopyRaster_management (extSurface, outIceSurfaces)

    ##Interpret ice surface based on ice thickness
    icethickness = extSurface - BedDEM
    ##Should make sure that no negtive value for the thickness, is negative, set as zero?? Need to explore where are these negative values
    conIcethickness = Con(icethickness > 0, icethickness, 0) 
    #arcpy.CopyRaster_management (icethickness, outIceThickness)
    arcpy.CopyRaster_management (conIcethickness, outIceThickness)
    '''
    
    ##Output method 2: to be consistent and compariable with the Volta ice thickness program, first interpret the ice thickness and then reinterpret the ice surface
    ##if there is target features, using the target features as the zero contour lines?????
    singlepart_outlines = arcpy.MultipartToSinglepart_management(outIcePolys, "in_memory\\singlepart_outlines")
    outline_lines_in = arcpy.PolygonToLine_management(singlepart_outlines, "in_memory\\outlines_line_in")
    
    arcpy.AddField_management(outline_lines_in, "contour", "SHORT")
    arcpy.CalculateField_management(outline_lines_in,"contour",0)
    cellsize_interp = cellsize_float

    ##Select only the ice thick > 0 for the thickness interpretation to prevent the negative ice thickness interpretation
    outpoints_not_zero = "in_memory\\outpoints_not_zero"
    arcpy.Select_analysis (outpoints, outpoints_not_zero, "thick > 0")
    
    interpolated_ice_depth = TopoToRaster([TopoPointElevation([[outpoints_not_zero, 'thick']]), TopoContour([[outline_lines_in, 'contour']]), TopoBoundary ([singlepart_outlines])], cellsize_interp, "", '20', '0', '#', 'NO_ENFORCE', "SPOT", '1', '#', '1', '0', '0', '200') ##,"","","","",0.1)

    outThickness = Con(interpolated_ice_depth > 0, interpolated_ice_depth, 0)
    arcpy.CopyRaster_management (outThickness, outIceThickness)
    ##Get ice surface based on ice thickness
    icesurface_final = outThickness + BedDEM

    #arcpy.CopyRaster_management (interpolated_ice_depth, outIceThickness)

    ##Get ice surface based on ice thickness
    #icesurface_final = interpolated_ice_depth + BedDEM
    #icesurface_final = interpolated_ice_depth + fillDEM
    ##Should make sure that no negtive value for the thickness, is negative, set as zero?? Need to explore where are these negative values
    arcpy.CopyRaster_management (icesurface_final, outIceSurfaces)

    #icethickness = icesurface_final - BedDEM
    #arcpy.CopyRaster_management (icethickness, outIceThickness)
    
    ##Delete temp datasets
    arcpy.Delete_management(icesurs) 
    arcpy.Delete_management(oldsurface) 

    return outpoints, outIcePolys, outIceSurfaces, outIceThickness

####-------Start the main program-----------------------####
if __name__ == '__main__':

    #Define input data this is the core data
    BedDEM = arcpy.GetParameterAsText(0)
    inputflowline = arcpy.GetParameterAsText(1)
    Distance=int(arcpy.GetParameter(2))
    inwatershed = arcpy.GetParameterAsText(3)
    TargetFeatures = arcpy.GetParameterAsText(4)
    shearstress = float(arcpy.GetParameter(5))
    min_ss = float(arcpy.GetParameter(6))
    max_ss = float(arcpy.GetParameter(7))
    bFactorPolyfit = arcpy.GetParameter(8)
    outpoints=arcpy.GetParameterAsText(9)
    outIcePolys = arcpy.GetParameterAsText(10)
    outIceSurfaces = arcpy.GetParameterAsText(11)
    outIceThickness = arcpy.GetParameterAsText(12)

    arcpy.Delete_management("in_memory") ### Empty the in_memory


    PaleoIceReconstruction(BedDEM, inputflowline, Distance, inwatershed, TargetFeatures, shearstress, min_ss, max_ss, bFactorPolyfit, outpoints, outIcePolys, outIceSurfaces, outIceThickness)

    arcpy.Delete_management("in_memory") ### Empty the in_memory
