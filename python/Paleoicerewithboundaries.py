#-------------------------------------------------------------------------------
# Name: Paleo Ice Reconstruction with ice boundary
# Purpose: This tool calculates the ice thickness on flowlines based on the Excel flowline model introduced by the Benn and Houlton (2010)
# and the python model, GLaRe, developed by Pellitero et al.(2016). However, both of the above models require manually assign shear stress
# and F factors. It requires a lot of time and efforts to conduct a paleo ice reconstrcution work.
# 
# The new developed tool will automatically adjust shear stress and F factors based on the DEM and target features. For shear stress, the tool
# assumes one value for the whole glacier based on the recommendation of many previous studies. The shear stress is first derived based on a
# revised shearstress code from VOLTA based on the elevation distrbution of ice surface and then adjusted to reach the best fit to
# the ice boundary.
#
# This tool also automatically adjust the F factor (shape factor) based on the cross sections along the flowlines. There are two options to derive
# the F factor: one is based on the cross section area, ice-contact perimeter, and ice thickness; the other is based on the fit of the polynomial
# function introduced by Li et al (2012). This tool uses a maximum width to prevent the error long cross sections in some sections when
#  the direction is not paralell the valley direction and where the tributary valley joins the main valley. This tool also applied EU
# allocation method to make sure the cross section does not extend to the tributaries (May need to check if it is necessary because the width is
# already used to constrain the exent of the cross section).
# 
# Author:      Yingkui Li
# Created:     06/01/2021 to 06/02/2021
# Copyright:   (c) Yingkui Li 2020
# Department of Geography, University of Tennessee
# Knoxville, TN 37996, USA
#-------------------------------------------------------------------------------
#from __future__ import division
from SharedFunctions import *  ## 

from arcpy.sa import *


#------------------------------------------------------------------------------------------------------------
# This function calculates ice surface elevation for the points along the flowlines and assign the ice surface elevation
# of the these points to the closest cross section points. Then, the ice surface elevations of the the cross section
# points will be used to interpret the ice surface raster for the whole ice polygon boundary based on the Topo to Raster tool.
#------------------------------------------------------------------------------------------------------------
def IceSurfaceCalculation_with_crosssection_pnts_for_iceboundary (flowline_points, beddem, cellsize, icebnd, field, outicesurface):
    #arcpy.env.extent = beddem
    #arcpy.env.cellSize = beddem

    ##Convert the ice bondary to contour lines of zero depth
    dispoly = "in_memory\\dispoly"
    arcpy.Dissolve_management(icebnd, dispoly, '#', '#', 'SINGLE_PART', '#')

    singlepart_icebnd = arcpy.MultipartToSinglepart_management(dispoly, "in_memory\\singlepart_icebnd")
    icebnd_lines = arcpy.PolygonToLine_management(singlepart_icebnd, "in_memory\\icebnd_lines")
    arcpy.AddField_management(icebnd_lines, "contour", "SHORT")
    arcpy.CalculateField_management(icebnd_lines,"contour",0)

    #arcpy.CopyFeatures_management(flowline_points, "d:\\temp\\flowline_points.shp")
    #arcpy.CopyFeatures_management(icebnd_lines, "d:\\temp\\icebnd_lines.shp")
    #arcpy.CopyFeatures_management(singlepart_icebnd, "d:\\temp\\singlepart_icebnd.shp")
    

    palaeoThickness = TopoToRaster([TopoPointElevation([[flowline_points, field]]), TopoContour([[icebnd_lines, 'contour']]), TopoBoundary ([singlepart_icebnd])], cellsize, 
                       singlepart_icebnd, "#", "#", "#", "ENFORCE")

    icesurface = palaeoThickness + beddem

    arcpy.CopyRaster_management (icesurface, outicesurface)

    return outicesurface

#------------------------------------------------------------------------------------------------------------
# This fuction is the whole process to reconstruct paleoice based on DEM, input flowlines, ice boundary, and default shear stress
#------------------------------------------------------------------------------------------------------------
def PaleoIceReconstructionwithboundary(BedDEM, inputflowline, Distance, iceboundary, shearstress, min_ss, max_ss, bFactorPolyfit, outpoints, outIceSurfaces, outIceThickness):
    arcpy.env.extent = BedDEM
    arcpy.env.cellSize = BedDEM
    arcpy.env.snapRaster = BedDEM ##setup snap raster
    
    GlacierID = "GlacierID" ##This is an ID field in inputflowline to identify the flowline(s) for each glacier (maybe connected with multiple flowlines)
    spatialref=arcpy.Describe(inputflowline).spatialReference 
    arcpy.env.outputCoordinateSystem = spatialref

    flowlineType = arcpy.Describe(inputflowline).shapeType
    if not (flowlineType == "Polyline"): ##quit if not polyline features
        arcpy.AddError("You did not define a flowline polyline input")
        arcpy.GetMessage(0)
        exit()

    icepolyType =arcpy.Describe(iceboundary).shapeType
    if not (icepolyType == "Polygon" or icepolyType == "Polyline"): ##quit if not polygon features
        arcpy.AddError("You did not define the polygon or polyline ice boundary input")
        arcpy.GetMessage(0)
        exit()


    icebndpolys = "in_memory\\icebndpolys"
    if (icepolyType == "Polyline"):
        arcpy.FeatureToPolygon_management(iceboundary, icebndpolys)
        count_result = arcpy.GetCount_management(icebndpolys)
        if int(count_result.getOutput(0)) < 1:
            arcpy.AddError("The paleo ice boundary cannot form the polygons")
            arcpy.GetMessage(0)
            exit()
    else:
        arcpy.CopyFeatures_management(iceboundary, icebndpolys)

    ####Flow direction and accumulation analysis
    cellsize = arcpy.GetRasterProperties_management(BedDEM,"CELLSIZEX")
    cellsize_float = float(cellsize.getOutput(0)) # use float cell size


    ##remove potential spurious polygons from the icebndpolys
    count_result = arcpy.GetCount_management(icebndpolys)
    #arcpy.AddMessage("Number of Outlines: " + str(int(count_result.getOutput(0))))
    minArea = cellsize_float * cellsize_float * 5 ##5 pixels 
    with arcpy.da.UpdateCursor(icebndpolys, "SHAPE@AREA") as cursor:        
        for row in cursor:
            if row[0] < minArea:
                arcpy.AddMessage("Delete spurious polygon")
                cursor.deleteRow()  ##Just remove the max elevation and keep the lowest elevation
    del cursor, row

    count_result = arcpy.GetCount_management(icebndpolys)
    if int(count_result.getOutput(0))==0:
        arcpy.AddWarning("There is no outline")
        quit()
         

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

    #Make a copy of the flowline
    flowlines = "in_memory\\flowlines"
    flowline = "in_memory\\flowline"
    flowline3dpoints = "in_memory\\flowline3dpoints"
    selflowline3dpoints = "in_memory\\selflowline3dpoints"
    singepoint = "in_memory\\singepoint"
    ws = "in_memory\\ws"
    wsflowpoints = "in_memory\\wsflowpoints"
    icepolys = "in_memory\\icepolys"
    icepolyselect = "in_memory\\icepolyselect"
    mainflowline = "in_memory\\mainflowline"
    flowlinecopy = "in_memory\\flowlinecopy"
    allicepolys = "in_memory\\allicepolys"
    singeflowpnts = "in_memory\\singeflowpnts"
    icewatersheds = "in_memory\\icewatersheds"
    icesurs = arcpy.env.scratchFolder + "\\r" + "icesurs" ##the inmemory does not work for raster


    ##watershed delineations
    burninDEM = BedDEM - Power (cellsize_float / (cellsize_float + EucDistance(inputflowline) ), 2 ) * 10 ##Burn in the DEM to make sure the flow pass through the flowline start points
    ##Start to delineate the watershed
    #Hydro analysis
    fillDEM =Fill(burninDEM)  ##Fill the sink first
    fdir = FlowDirection(fillDEM,"NORMAL") ##Flow direction
        
    facc = FlowAccumulation(fdir) ##Flow accmulation    

    ##For some times this the following does not work
    arcpy.AddField_management(icebndpolys,GlacierID,"LONG", 10)
    arcpy.CalculateField_management(icebndpolys,GlacierID,str("!"+str(arcpy.Describe(icebndpolys).OIDFieldName)+"!"),"PYTHON_9.3")
    arcpy.SpatialJoin_analysis(inputflowline, icebndpolys, flowlines, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
    count_result = arcpy.GetCount_management(flowlines)
    #arcpy.AddMessage(int(count_result.getOutput(0)))

    if int(count_result.getOutput(0)) == 0: ##test the selection of flowlines by select layer by location
        outline_layer = arcpy.MakeFeatureLayer_management(icebndpolys, "in_memory\\outline_layer")
        flowline_layer = arcpy.MakeFeatureLayer_management(inputflowline, "in_memory\\flowline_layer")
        arcpy.SelectLayerByLocation_management(flowline_layer, "WITHIN", outline_layer)
        arcpy.CopyFeatures_management(flowline_layer, flowlines)
        count_result2 = arcpy.GetCount_management(flowlines)
        #arcpy.AddMessage(int(count_result2.getOutput(0)))
        if int(count_result2.getOutput(0)) == 0:
            arcpy.AddWarning("No flowlines inside glacier outlines!")
            quit()

    #arcpy.CopyFeatures_management(flowlines, "d:\\temp\\flowlines.shp")
    #exist_fields = [f.name for f in arcpy.ListFields(inputflowline)] #List of current field names in outline layer
    #if GlacierID not in exist_fields:
    #    arcpy.AddMessage("Assigning Glacier ID...")
    #    Add_GlacierID_by_Touches (inputflowline, GlacierID, flowlines)
    #else:
    #    arcpy.CopyFeatures_management(inputflowline, flowlines)

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

    del cursor, row

    ##Order the line geometries in the list
    arcpy.AddMessage("Ordering flowlines...")
    

    order = sorted(range(len(height)), key=lambda k: height[k])  ##order is the ID

    arcpy.AddField_management(flowlines,"ProcessID","LONG", 10)

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

    for gid in range(len(uniqueiceID)):
        query = GlacierID + " = " + str(uniqueiceID[gid])
        #arcpy.AddMessage(query)
        arcpy.AddMessage("Processing #" + str(gid+1) +"/" + str(len(uniqueiceID)) + " of reconstructed glaciers...")                                                                                       
        
        arcpy.Select_analysis (flowlines, flowline, query)
        #arcpy.CopyFeatures_management(flowline, "c:\\temp2\\flowline.shp")

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

        ##select the ice polygon corresponding to the flowline glacier ID
        ##Get the ice polygon with the considertation of the watershed
        startperppolygon = Startpoint_perpendiculars(mainflowline, int(cellsize_float), 300) ##Just create a short distance perpendicular polygon (need polygon because sometimes line does not cross the highest Facc) from the start point to get the highest Facc point, to avoid the snap pointpoints downstream
        outZonalStatistics = ZonalStatistics(startperppolygon, arcpy.Describe(startperppolygon).OIDFieldName, facc, "MAXIMUM") #Find the maximum flowaccumulation point on the moriane feature
        OutHighestFcc = Con(facc == outZonalStatistics,facc)  ##Determine the highest flowaccumuation part
        #Snap_Pour Points ## the purpose is to create the same extent raster with only the highest fcc point
        outSnapPour = SnapPourPoint(OutHighestFcc, facc, 0, "VALUE") ## Just create a pourpoint raster with the same extent of the input DEM
        outpntWs = Watershed(fdir, outSnapPour, "VALUE")

        ##make sure to have the watershed and iceboundary intersection that can cover all flowline points
        points_count = int(arcpy.GetCount_management(selflowline3dpoints).getOutput(0))
        arcpy.RasterToPolygon_conversion(outpntWs, ws)
        inside_point_count = 0
        i = 1
        while (points_count > inside_point_count): 
            arcpy.Clip_analysis (selflowline3dpoints, ws, "in_memory\\selectedflpoints")
            inside_point_count = int(arcpy.GetCount_management("in_memory\\selectedflpoints").getOutput(0))
            if inside_point_count < points_count:
                #arcpy.AddMessage("The watershed corresponding to the flowline does not include all flowline points! The flowline does not follow the streamline!")
                outSnapPour = SnapPourPoint(singepoint, facc, 100*i) ##each time increase 100 m downstream
                outpntWs = Watershed(fdir, outSnapPour)
                arcpy.RasterToPolygon_conversion(outpntWs, ws)
            ##add the maximum loop controls
            if i >= 10:
                #arcpy.AddMessage("Cannot find the watershed that include all flowline points. Use the input ice boundary instead!")
                arcpy.CopyFeatures_management(icebndpolys, ws)
                break
            i += 1

        arcpy.Clip_analysis (ws, icebndpolys, "in_memory\\icepolycliped")
        ##Convert to single parts
        arcpy.MultipartToSinglepart_management("in_memory\\icepolycliped", icepolyselect)
        #remove the potential spurious polygons
        with arcpy.da.UpdateCursor(icepolyselect, "SHAPE@AREA") as cursor:
            i = 0
            for row in cursor:
                i += 1
                if row[0] < minArea:
                    cursor.deleteRow()
                    #arcpy.AddMessage("delete one spurious polygon!")
        del cursor
        #arcpy.AddMessage("the number of icepolyselect is:" + str(i))
        if i> 0:
            del row


        cross_section_pnts = cross_section_points(selflowline3dpoints, flowline, icepolyselect, BedDEM, 1600, 60)

        ##Ice thickness calculation along the flowlines
        Ice_Thickness_Calculation (selflowline3dpoints)
        
        ss0 = shearstress ##each glacier only have one shear stress for all flowlines and flow points

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
            icesurface = IceSurfaceCalculation_with_crosssection_pnts_for_iceboundary (selflowline3dpoints, BedDEM, cellsize_float, icepolyselect, "thick", "in_memory\\icesurface")

            if bAdjustSS == False:
                #ss = shear_stress_calculation (mainflowline, icepoly, icesurface)
                ss = shear_stress_calculation(mainflowline, icepolyselect, icesurface, min_ss, max_ss)
            else:
                if abs(distance_list[-2] - distance_list[-1]) > 0: 
                    m = (ss_list[-2] - ss_list[-1])/(distance_list[-2] - distance_list[-1])
                    ss = ss_list[-1] - m * distance_list[-1]
                else:
                    ss = (ss_list[-2] - ss_list[-1]) / 2

                if ss < min_ss: ##If the shear stress is negative, using the default value
                    ss = min_ss
                if ss > max_ss:
                    ss = max_ss
                
            arcpy.AddMessage("The calculated shear stress is:" + str(ss) + " and the percentage difference with the previous value is " + str(abs(ss-ss0)/ss0)) ####+ " \%"
            ss_ratio = abs(ss - ss0)/ss0
            if abs(ss_ratio - ss_ratio0)< 0.005 or (ss_ratio < 0.01):  ##assuming ss always increase, to stop if decreasing; ##use 1.0% change as the threshold to stop the loop
                break
                
            ##update ss value
            arcpy.CalculateField_management(selflowline3dpoints,"SSTRESS",ss) ##update ss to the flowline points

            if bFactorPolyfit == True:
                arcpy.AddMessage("Deriving F factor based on the Polyfit of the cross section...")
                AdjustFfactor_Ployfit_with_cross_section_pnts (selflowline3dpoints, cross_section_pnts, "ice", icepolyselect)
                ##replace the ffactor == 0.8 as the average F factor of the whole section
                
            else:
                arcpy.AddMessage("Deriving F factor from cross section...")
                #Try to use the curve fitting method to derive the F factor
                AdjustFfactor_with_cross_section_pnts (selflowline3dpoints, cross_section_pnts, "ice", icepolyselect) 

            ##Recalculate the ice thickness for the flowpoints based on adjusted F and ss parameters
            Ice_Thickness_Calculation (selflowline3dpoints)

            ss0 = ss
            ss_ratio0 = ss_ratio
            nloop += 1
            
        if gid == 0: #The first time
            arcpy.CopyFeatures_management(selflowline3dpoints, outpoints)
            arcpy.CopyRaster_management(icesurface, icesurs)
        else:
            arcpy.Append_management(selflowline3dpoints, outpoints, "NO_TEST" )
            arcpy.Mosaic_management(icesurface, icesurs, "MEAN","","", "", "", "", "")

    ##Smooth the polgyons and extract the ice surface within the ice polygon
    arcpy.AddMessage("Generating final outputs...")

    ##Dissolve polyygons so that the shared boundaries are removed
    disolveicepoly = "in_memory\\disolveicepoly"
    arcpy.Dissolve_management(icebndpolys, disolveicepoly, '#', '#', 'SINGLE_PART', '#')

    
    ##Output method: to be consistent and compariable with the Volta ice thickness program, first interpret the ice thickness and then reinterpret the ice surface
    ##if there is target features, using the target features as the zero contour lines?????
    singlepart_outlines = arcpy.MultipartToSinglepart_management(disolveicepoly, "in_memory\\singlepart_outlines")
    outline_lines_in = arcpy.PolygonToLine_management(singlepart_outlines, "in_memory\\outlines_line_in")
    arcpy.AddField_management(outline_lines_in, "contour", "SHORT")
    arcpy.CalculateField_management(outline_lines_in,"contour",0)
    cellsize_interp = cellsize_float

    ##Select only the ice thick > 0 for the thickness interpretation to prevent the negative ice thickness interpretation
    outpoints_not_zero = "in_memory\\outpoints_not_zero"
    arcpy.Select_analysis (outpoints, outpoints_not_zero, "thick > 0")
    interpolated_ice_depth = TopoToRaster([TopoPointElevation([[outpoints_not_zero, 'thick']]), TopoContour([[outline_lines_in, 'contour']]), TopoBoundary ([singlepart_outlines])], cellsize_interp, singlepart_outlines, '#', '#', '#', 'ENFORCE')####, "SPOT", '1', '#', '1', '0', '0', '200') ##,"","","","",0.1)
    outThickness = Con(interpolated_ice_depth > 0, interpolated_ice_depth, 0)
    arcpy.CopyRaster_management (outThickness, outIceThickness)
    ##Get ice surface based on ice thickness
    icesurface_final = outThickness + BedDEM
    ##Should make sure that no negtive value for the thickness, is negative, set as zero?? Need to explore where are these negative values
    arcpy.CopyRaster_management (icesurface_final, outIceSurfaces)
    #except:
    #    arcpy.AddMessage("This is an error in TopoToRaster")
        
    ##Delete temp datasets
    arcpy.Delete_management(icesurs) 

    return outpoints, outIceSurfaces, outIceThickness

####-------Start the main program-----------------------####
if __name__ == '__main__':

    #Define input data this is the core data
    BedDEM = arcpy.GetParameterAsText(0)
    inputflowline = arcpy.GetParameterAsText(1)
    inputiceboundary = arcpy.GetParameterAsText(2)
    Distance=int(arcpy.GetParameter(3))
    shearstress = float(arcpy.GetParameter(4))
    min_ss = float(arcpy.GetParameter(5))
    max_ss = float(arcpy.GetParameter(6))
    bFactorPolyfit = arcpy.GetParameter(7)
    outpoints=arcpy.GetParameterAsText(8)
    outIceSurfaces = arcpy.GetParameterAsText(9)
    outIceThickness = arcpy.GetParameterAsText(10)

    arcpy.Delete_management("in_memory") ### Empty the in_memory

    PaleoIceReconstructionwithboundary(BedDEM, inputflowline, Distance, inputiceboundary, shearstress, min_ss, max_ss, bFactorPolyfit, outpoints, outIceSurfaces, outIceThickness)

    arcpy.Delete_management("in_memory") ### Empty the in_memory
