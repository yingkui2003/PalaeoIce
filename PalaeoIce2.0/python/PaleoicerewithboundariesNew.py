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
# Updated:     02/16/2023
# Copyright:   (c) Yingkui Li 2020
# Department of Geography, University of Tennessee
# Knoxville, TN 37996, USA
#-------------------------------------------------------------------------------
#from __future__ import division
from SharedFunctions import *  ## 

from arcpy.sa import *



#------------------------------------------------------------------------------------
# This function do the surface interpretation based on a set of points and a interpretation method
#------------------------------------------------------------------------------------
def surface_interpolation (inpoints, field, dem, boundary_polygon, method, outsurface):
    ####Flow direction and accumulation analysis
    #arcpy.env.cellSize = dem
    #arcpy.env.snapRaster = dem
    boundarypoints = "in_memory\\boundarypoints"
    arcpy.PolygonToLine_management(boundary_polygon, "in_memory\\boundary_lines")
    arcpy.FeatureVerticesToPoints_management("in_memory\\boundary_lines", boundarypoints, 'ALL')
    boundarypoints3D = ExtractValuesToPoints(boundarypoints, dem, "in_memory\\boundarypoints3D", "INTERPOLATE")
    arcpy.AddField_management(boundarypoints3D, field, "FLOAT",10,6)
    arcpy.CalculateField_management(boundarypoints3D, field, "!RASTERVALU!","PYTHON_9.3")
    arcpy.DeleteField_management(boundarypoints3D, 'RASTERVALU')

    ##Remove the inpoints close to the boundary line
    arcpy.Buffer_analysis("in_memory\\boundary_lines", "in_memory\\boundarydbuf", "90 Meter")
    arcpy.Erase_analysis(inpoints, "in_memory\\boundarydbuf", "in_memory\\inpoints_after_erase")

    #arcpy.CopyFeatures_management(inpoints, "c:\\test\\inpoints.shp")

    ##Remove the identical points
    arcpy.conversion.PointToRaster("in_memory\\inpoints_after_erase", field, "in_memory\\pntraster", "MEAN", "NONE", 30)###, "BUILD")
    arcpy.DeleteIdentical_management("in_memory\\inpoints_after_erase", "Shape", "30 Meter")
    ExtractValuesToPoints("in_memory\\inpoints_after_erase", "in_memory\\pntraster", "in_memory\\graded_points")
    arcpy.CalculateField_management("in_memory\\graded_points", field, "!RASTERVALU!","PYTHON_9.3")
    arcpy.DeleteField_management("in_memory\\graded_points", 'RASTERVALU')
    
    #arcpy.CopyFeatures_management("in_memory\\inpoints_after_erase", "c:\\test\\inpoints_after_erase.shp")
    #arcpy.CopyFeatures_management("in_memory\\graded_points", "c:\\test\\graded_points.shp")
    
    #arcpy.Append_management("in_memory\\inpoints_after_erase", boundarypoints3D, "NO_TEST")
    arcpy.Append_management("in_memory\\graded_points", boundarypoints3D, "NO_TEST")

    #arcpy.CopyFeatures_management(boundarypoints3D, "c:\\test\\boundarypoints3D.shp")

    surface3d = arcpy.env.scratchGDB + "\\surface3d"
    if method == "TopoToRaster":
        pointElevations = TopoPointElevation([[boundarypoints3D,field]])
        pointSurface = TopoToRaster([pointElevations])
        arcpy.CopyRaster_management (pointSurface, surface3d)
    elif method == "Kriging":
        arcpy.Kriging_3d(boundarypoints3D, field, surface3d,"Circular", "#", "Variable 5")
    elif method == "IDW":
        arcpy.Idw_3d(boundarypoints3D, field, surface3d, "#",2)
    elif method == "Spline":
        arcpy.Spline_3d(boundarypoints3D, field, surface3d, "", "REGULARIZED")###, 0.1)
    elif method == "Trend":
        arcpy.Trend_3d(boundarypoints3D, field, surface3d, "#", 2, "LINEAR")
    elif method == "NaturalNeighbor":
        arcpy.NaturalNeighbor_3d(boundarypoints3D, field, surface3d, "#")

    extSurface = ExtractByMask(surface3d, boundary_polygon)
    arcpy.CopyRaster_management (extSurface, outsurface)

    arcpy.Delete_management(surface3d)
    
    return outsurface   

#------------------------------------------------------------------------------------------------------------
# This function calculates ice surface elevation for the points along the flowlines and assign the ice surface elevation
# of the these points to the closest cross section points. Then, the ice surface elevations of the the cross section
# points will be used to interpret the ice surface raster for the whole ice polygon boundary based on the Topo to Raster tool.
#------------------------------------------------------------------------------------------------------------
def IceSurfaceCalculation_with_crosssection_pnts_for_iceboundary (flowline_points, beddem, icebnd, field, method, outicesurface):
    #arcpy.env.extent = beddem
    #arcpy.env.cellSize = beddem

    ##Convert the ice bondary to contour lines of zero depth
    #dispoly = "in_memory\\dispoly"
    #arcpy.Dissolve_management(icebnd, dispoly, '#', '#', 'SINGLE_PART', '#')

    #singlepart_icebnd = arcpy.MultipartToSinglepart_management(dispoly, "in_memory\\singlepart_icebnd")

    #icebnd_lines = arcpy.PolygonToLine_management(singlepart_icebnd, "in_memory\\icebnd_lines")
    #arcpy.AddField_management(icebnd_lines, "contour", "SHORT")
    #arcpy.CalculateField_management(icebnd_lines,"contour",0)
    #palaeoThickness = TopoToRaster([TopoPointElevation([[flowline_points, field]]), TopoContour([[icebnd_lines, 'contour']]), TopoBoundary ([singlepart_icebnd])])
    #icesurface = palaeoThickness + beddem
    #arcpy.CopyRaster_management (icesurface, outicesurface)


    ##Use the cross section points that already created to interpret the ice surfaces
    clip_flowline_points = "in_memory\\clip_flowline_points"
    arcpy.Clip_analysis (flowline_points, icebnd, clip_flowline_points)

    ##Check point1: should consider the target feature elevations here. If close to the target elevation, then, use the target elevation for the near
    #elev_control_points = "in_memory\\elev_control_points"
    #arcpy.CopyFeatures_management(flowline_points, elev_control_points)
    '''
    arcpy.Near_analysis (clip_CS_points, elev_control_points)#identify nearest flowline point
    
    #pointarray = arcpy.da.FeatureClassToNumPyArray(flowline_points,["OID@",field])#create array of ice values and populate a list with them
    pointarray = arcpy.da.FeatureClassToNumPyArray(elev_control_points,["OID@",field])#create array of ice values and populate a list with them
    fidval = np.array([item[0] for item in pointarray])
    iceval = np.array([item[1] for item in pointarray])

    ###Add field into the cross section points
    exist_fields = [f.name for f in arcpy.ListFields(clip_CS_points)] #List of current field names in outline layer
    if field not in exist_fields:
        arcpy.AddField_management(clip_CS_points, field, "FLOAT",10,6) #field for ice value

    with arcpy.da.UpdateCursor(clip_CS_points,("NEAR_FID", field, "PointZ")) as cursor:   #populate ice field with value from the nearest flowline point
        i = 0
        for row in cursor:
            fid = row[0]
            idx_result = np.where(fidval == fid)
            idx = idx_result[0][0]
            row[1]=iceval[idx]
            #row[2] = row[1] - row[3]                   
            cursor.updateRow(row)
            try:
                if row[1] < row[2]: ##if the ice surface elevation is lower than the bed elevation
                    cursor.deleteRow()
            except:
                #arcpy.AddMessage("PointZ is: " + str(row[3]))
                pass
            i += 1
    if i>0:   
        del row
    del cursor
    '''
    #modified on 2/16/2023
    #arcpy.CopyFeatures_management(clip_CS_points, "in_memory\\clip_CS_points_cp")    
    fieldlist=[field.name for field in arcpy.ListFields(clip_flowline_points)]
    if "RASTERVALU" in fieldlist:
        arcpy.DeleteField_management(clip_flowline_points, 'RASTERVALU')
   
    surface_interpolation (clip_flowline_points, field, beddem, icebnd, method, outicesurface)
    #arcpy.CopyFeatures_management(clip_CS_points, out_CS_points)

    return outicesurface

#------------------------------------------------------------------------------------------------------------
# This fuction creates perpendicular lines along the flowlines and then create a set of points along these
# perpendicular lines. This tool aslo extracts the elevation of each points. The max_width is used to limit
# the extent of these points. The division between the perp lines in different valleys are determined by the
# EU Allocation tool in ArcGIS.
#------------------------------------------------------------------------------------------------------------
def cross_section_points_with_TF(flowlinepoints, flowline, watershed, beddem, targetFC, cellsize_float, max_width, resolution): 

    #cellsize = arcpy.GetRasterProperties_management(beddem,"CELLSIZEX")
    #cellsize_float = float(cellsize.getOutput(0)) # use float cell size
    spatialref=arcpy.Describe(flowlinepoints).spatialReference
    
    width = max_width/2.0  ##only use the half the width for each side of the flowline

    oldextent = arcpy.env.extent
    arcpy.Buffer_analysis(watershed, "in_memory\\watershedbuf", (str(cellsize_float)+ " Meter"))
    arcpy.env.extent = "in_memory\\watershedbuf"

    #arcpy.CopyFeatures_management("in_memory\\watershedbuf", "c:\\test\\watershedbuf.shp")
    
    ##use the buffer so that the extract elevations of the points along the watershed boundary are not zero
    try:
        error = False
        extractDEM = ExtractByMask(beddem, "in_memory\\watershedbuf") 
    except:
        error = True

    ##Step 1: data preparation
    flowlinepointscp = arcpy.CopyFeatures_management(flowlinepoints, "in_memory\\flowlinepointscp")
    arcpy.AddField_management(flowlinepointscp, 'PntID', 'Long', 6) ##OFID is not FID, but it is related to FID
    arcpy.CalculateField_management(flowlinepointscp,"PntID",str("!"+str(arcpy.Describe(flowlinepointscp).OIDFieldName)+"!"),"PYTHON_9.3")

    arcpy.AddField_management(flowline, 'SegmentID', 'Long', 6) ##OFID is not FID, but it is related to FID
    arcpy.CalculateField_management(flowline,"SegmentID",str("!"+str(arcpy.Describe(flowline).OIDFieldName)+"!"),"PYTHON_9.3")
    flowlinepointscopy = "in_memory\\flowlinepointscopy"
    arcpy.SpatialJoin_analysis(flowlinepointscp, flowline, flowlinepointscopy, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")

    #arcpy.CopyFeatures_management(flowlinepointscopy, "c:\\test\\flowlinepointscopy.shp")

    arr=arcpy.da.FeatureClassToNumPyArray(flowlinepointscopy, ('SHAPE@X', 'SHAPE@Y', 'PntID', 'OFID', 'SegmentID'))
    segment_ids = np.array([item[4] for item in arr])
    unique_segment_ids = np.unique(segment_ids)

    ##Step 2: create perpendicular lines along the flowpoints
    distance = 10000
        
    perpendiculars = arcpy.CreateFeatureclass_management("in_memory", "perpendiculars","POLYLINE", "","","", spatialref)
    arcpy.AddField_management(perpendiculars, 'FlowlineID', 'Long', 6) ##OFID is not FID, but it is related to FID
    arcpy.AddField_management(perpendiculars, 'SegmentID', 'Long', 6) ##OFID is not FID, but it is related to FID
    arcpy.AddField_management(perpendiculars, 'FlowPntID', 'Long', 6) ##OFID is not FID, but it is related to FID
    
    exclude_PntID = []
    new_line_cursor = arcpy.da.InsertCursor(perpendiculars, ('SHAPE@', 'FlowPntID' , 'FlowlineID', 'SegmentID'))
    for row in range(len(arr)):
        if row > 0: ##not the first point
            startx = arr[row][0]
            starty = arr[row][1]
            creat_perp = 1
            if row <(len(arr)-1):
                if arr[row][3]==arr[row+1][3] and arr[row][3]==arr[row-1][3]: ##inside points
                    endx = arr[row+1][0]
                    endy = arr[row+1][1]
                elif arr[row][3]==arr[row+1][3]: ##start tributary points
                    exclude_PntID.append(arr[row][2])
                    creat_perp = 0
                else: ##end of the tributary
                    endx = arr[row-1][0]
                    endy = arr[row-1][1]
            else:
                endx = arr[row-1][0]
                endy = arr[row-1][1]

            if creat_perp > 0:
                if starty==endy or startx==endx:
                    if starty == endy:
                        y1 = starty + distance
                        y2 = starty - distance
                        x1 = startx
                        x2 = startx
                    if startx == endx:
                        y1 = starty
                        y2 = starty 
                        x1 = startx + distance
                        x2 = startx - distance     
                else:
                    m = ((starty - endy)/(startx - endx)) #get the slope of the line
                    negativereciprocal = -1*((startx - endx)/(starty - endy))    #get the negative reciprocal
                    if m > 0:
                        if m >= 1:
                            y1 = negativereciprocal*(distance)+ starty
                            y2 = negativereciprocal*(-distance) + starty
                            x1 = startx + distance
                            x2 = startx - distance
                        if m < 1:
                            y1 = starty + distance
                            y2 = starty - distance
                            x1 = (distance/negativereciprocal) + startx
                            x2 = (-distance/negativereciprocal)+ startx           
                    if m < 0:
                        if m >= -1:
                            y1 = starty + distance
                            y2 = starty - distance
                            x1 = (distance/negativereciprocal) + startx
                            x2 = (-distance/negativereciprocal)+ startx     
                        if m < -1:
                            y1 = negativereciprocal*(distance)+ starty
                            y2 = negativereciprocal*(-distance) + starty
                            x1 = startx + distance
                            x2 = startx - distance
                array = arcpy.Array([arcpy.Point(x1,y1),arcpy.Point(x2, y2)])
                polyline = arcpy.Polyline(array)
                lineID = arr[row][3]
                pntID = arr[row][2]
                segID = arr[row][4]
                new_line_cursor.insertRow([polyline, pntID, lineID, segID])
                
    del new_line_cursor

    ##Step 3: Use EU Alloation to cut perp lines
    ##divide watershed by different flowlines
    ##Delete the exclude_points
    if len(exclude_PntID) > 0:
        with arcpy.da.UpdateCursor(flowlinepointscopy, 'PntID') as cursor:
            for row in cursor:
                if row in exclude_PntID:
                    #arcpy.AddMessage("Delete one points!")
                    cursor.deleteRow()
        del row, cursor

    #arcpy.CopyFeatures_management(flowlinepointscopy, "c:\\test\\flowlinepointscopy2.shp")

    # Process: Euclidean Allocation
    eucAllocate = EucAllocation(flowlinepointscopy, "", "", cellsize_float,'SegmentID', "", "")
    
    arcpy.RasterToPolygon_conversion(eucAllocate, "in_memory\\eucAllocate", "SIMPLIFY", "VALUE")
    arcpy.Clip_analysis ("in_memory\\eucAllocate", watershed, "in_memory\\clipAllocation")

    selAllocation = "in_memory\\selAllocation"
    sel_perp = "in_memory\\sel_perp"
    clipedlines = "in_memory\\clipedlines"

    for segment_id in unique_segment_ids:
        query = "gridcode = " + str(segment_id) 
        arcpy.Select_analysis ("in_memory\\clipAllocation", selAllocation, query)
        query = "SegmentID = "+ str(segment_id) 
        arcpy.Select_analysis (perpendiculars, sel_perp, query)

        arcpy.Clip_analysis (sel_perp, selAllocation, "in_memory\\clip_perp")

        if segment_id == unique_segment_ids[0]:
            arcpy.CopyFeatures_management("in_memory\\clip_perp", clipedlines)
        else:
            arcpy.Append_management("in_memory\\clip_perp", clipedlines, "NO_TEST" )
        
    
    ##Multipart to singleparts
    singlepartlines = arcpy.MultipartToSinglepart_management(clipedlines, "in_memory\\singlepartlines")

    #arcpy.CopyFeatures_management("in_memory\\singlepartlines", "c:\\test\\singlepartlines.shp")

    ##Step 4: Clean perp lines, so that only the part intersect the flowpoints are preserved
    lines = []
    flowlineids = []
    flowpntids = []
    keep = []
    with arcpy.da.SearchCursor(singlepartlines, ['SHAPE@', 'FlowlineID', 'FlowPntID']) as perps:
        i = 0
        for perp in perps:
            lines.append(perp[0])
            flowlineids.append(perp[1]) 
            flowpntids.append(perp[2])
            keep.append(0)
            i += 1
    if i > 0:
        del perp
    del perps

    flowpntid_arr = np.array(flowpntids)

    with arcpy.da.SearchCursor(flowlinepointscopy, ['SHAPE@', 'OFID', 'OID@']) as pnts:
        for pnt in pnts:
            index = np.where(flowpntid_arr == pnt[2])
            indexlen = len(index[0])
            if indexlen > 0:
                if indexlen > 1: ##more than two lines for a point, only select the lines intersected with the point
                    for a in range (indexlen):
                        value = index[0][a]
                        within = pnt[0].within(lines[value])
                        if within==True: 
                            keep[value] = 1
                            break
                else:
                    value = index[0][0]
                    keep[value] = 1
    del pnts, pnt

    ##Delete the lines with keep = 0
    with arcpy.da.UpdateCursor(singlepartlines, 'SHAPE@') as perps:
        i = 0
        for perp in perps:
            if keep[i] == 0:
                perps.deleteRow()
            i += 1
    if i>0:
        del perp
    del perps

    #arcpy.CopyFeatures_management("in_memory\\singlepartlines", "c:\\test\\singlepartlines2.shp")


    ##Step 5: intersection of the perp lines
    ##Maybe this step is not necessary, instead, using a maximum width to limit the extremely long cross sections
    '''
    ##Solve the intersection issue of the perp lines
    splitted_perps = "in_memory\\splitted_perps"
    arcpy.PolygonToLine_management("in_memory\\clipAllocation", "in_memory\\clipAlloc_to_lines")
    arcpy.Append_management("in_memory\\clipAlloc_to_lines", singlepartlines, 'NO_TEST')
    arcpy.FeatureToPolygon_management(singlepartlines, "in_memory\\perppolys", "0.05 Meters")
    arcpy.PolygonToLine_management("in_memory\\perppolys", "in_memory\\perppolys_to_lines")
    arcpy.SplitLine_management("in_memory\\perppolys_to_lines", "in_memory\\splitline_perppolys")
    arcpy.SpatialJoin_analysis("in_memory\\splitline_perppolys", flowlinepointscopy, "in_memory\\final_perps", "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
    arcpy.SplitLineAtPoint_management("in_memory\\final_perps", flowlinepointscopy, splitted_perps, '1 Meters')

    #arcpy.CopyFeatures_management(splitted_perps, "c:\\test\\singlepartlines3.shp")
    '''
    ##Added on 2/13/2023
    ##Add the target elevation field that intersected with target features 
    #arcpy.CopyFeatures_management(singlepartlines, "c:\\test\\singlepartlines0.shp")
    #arcpy.AddField_management(singlepartlines, 'TargetElev', 'Long', 6) ##OFID is not FID, but it is related to FID
    ##intersect the singlepartline with the target features
    arcpy.Intersect_analysis([singlepartlines, targetFC], "in_memory\\intersect_Targetpoints", "ONLY_FID", "#", "POINT")
    arcpy.MultipartToSinglepart_management("in_memory\\intersect_Targetpoints", "in_memory\\singletargetpoints")    
    ExtractValuesToPoints("in_memory\\singletargetpoints", extractDEM, "in_memory\\points3D", "INTERPOLATE")
    ##Need to spatialjoin the point elevation back to the singlepartlines??
    ##How about if there are multiple intersected for a line?? Need to summarize the point elevation
    #arcpy.CopyFeatures_management("in_memory\\points3D", "c:\\test\\points3D.shp")
    #arcpy.CopyFeatures_management(singlepartlines, "c:\\test\\singlepartlines.shp")

    #RASTERVALU "RASTERVALU" true true false 13 Float 0 0,Median,#
    # Create a new fieldmappings and add the two input feature classes.
    fieldmappings = arcpy.FieldMappings()
    fieldmappings.addTable(singlepartlines)
    fieldmappings.addTable("in_memory\\points3D")
     
    RASTERVALU = fieldmappings.findFieldMapIndex("RASTERVALU")
    fieldmap = fieldmappings.getFieldMap(RASTERVALU)
     
    # Get the output field's properties as a field object
    field = fieldmap.outputField
     
    # Rename the field and pass the updated field object back into the field map
    field.name = "TargetEle"
    field.aliasName = "TargetEle"
    fieldmap.outputField = field
     
    # Set the merge rule to mean and then replace the old fieldmap in the mappings object
    # with the updated one
    fieldmap.mergeRule = "Median"
    fieldmappings.replaceFieldMap(RASTERVALU, fieldmap)
    
    arcpy.SpatialJoin_analysis(singlepartlines, "in_memory\\points3D", "in_memory\\lines_with_targetElev", "JOIN_ONE_TO_ONE", "KEEP_COMMON", fieldmappings, "INTERSECT", "#", "#")
    
    #arcpy.CalculateField_management("in_memory\\lines_with_targetElev","TargetElev", str("!"+"RASTERVALU"+"!"),"PYTHON_9.3")
    #arcpy.DeleteField_management("in_memory\\lines_with_targetElev", 'RASTERVALU')
    #arcpy.CopyFeatures_management("in_memory\\lines_with_targetElev", "c:\\test\\lines_with_targetElev.shp")

    #Spatialjoin target elevations to the flowline poitns as a targetElev field
    fieldmappings = arcpy.FieldMappings()
    fieldmappings.addTable(flowlinepoints)
    fieldmappings.addTable("in_memory\\lines_with_targetElev")
     
    TargetElev = fieldmappings.findFieldMapIndex("TargetEle")
    fieldmap = fieldmappings.getFieldMap(TargetElev)
     
    # Get the output field's properties as a field object
    field = fieldmap.outputField
     
    # Rename the field and pass the updated field object back into the field map
    field.name = "TargetE"
    field.aliasName = "TargetE"
    fieldmap.outputField = field
     
    # Set the merge rule to mean and then replace the old fieldmap in the mappings object
    # with the updated one
    fieldmap.mergeRule = "First"
    fieldmappings.replaceFieldMap(TargetElev, fieldmap)
    
    arcpy.SpatialJoin_analysis(flowlinepoints, "in_memory\\lines_with_targetElev", "in_memory\\flowlinepoints_with_targetElev", "JOIN_ONE_TO_ONE", "KEEP_All", fieldmappings, "INTERSECT", "#", "#")
    arcpy.CopyFeatures_management("in_memory\\flowlinepoints_with_targetElev", flowlinepoints)

    #arcpy.AddField_management(flowlinepoints, "TargetElev", "DOUBLE", 6, 2)

    ##Need to only use the top and low 10% flowline points as the target constrain
    pntarray = arcpy.da.FeatureClassToNumPyArray(flowlinepoints, ('TargetE', 'RASTERVALU'))
    target_elevs = np.array([item[0] for item in pntarray])
    #bed_elevs = np.array([item[1] for item in pntarray])
    #arcpy.AddMessage(target_elevs)
    #arcpy.AddMessage(bed_elevs)
    #threshold_ele = min(bed_elevs) + (max(bed_elevs) - min(bed_elevs))*0.25
    
    num_pnts = len(target_elevs)
    for i in range(num_pnts):
        if 'nan' in str(target_elevs[i]).lower():
            #arcpy.AddMessage("Target Elevation is Nan")
            target_elevs[i] = 0
        #if bed_elevs[i] > threshold_ele:
        #    target_elevs[i] = 0

    for i in range(1, num_pnts+1):
        #print ("loop:" +str(i))
        for j in range(num_pnts - i):
            if target_elevs[-i] > 0 and (target_elevs[j] - target_elevs[-i])> (resolution/30)*(num_pnts-i-j): ##less than 2 degree slope
                target_elevs[j] = 0

    #arcpy.AddMessage(target_elevs)

    with arcpy.da.UpdateCursor(flowlinepoints,("TargetElev")) as cursor:
        i = 0
        for row in cursor:
            row[0]=target_elevs[i]
            cursor.updateRow(row)
            i += 1
    del row, cursor
    
    #arcpy.CopyFeatures_management(flowlinepoints, "c:\\test\\flowlinepoints.shp")

    arcpy.DeleteField_management(flowlinepoints, 'TargetE')
    
    ##Step 6: Create a set of points along the perp lines on both side of the flowlinepointscopy use the steps of resolution; also add the maxmium width constrain
    ##the following three lines to replace step 5
    splitted_perps = "in_memory\\splitted_perps"
    #arcpy.SpatialJoin_analysis("in_memory\\lines_with_targetElev", flowlinepointscopy, "in_memory\\final_perps", "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
    arcpy.SpatialJoin_analysis(singlepartlines, flowlinepointscopy, "in_memory\\final_perps", "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
    arcpy.SplitLineAtPoint_management("in_memory\\final_perps", flowlinepointscopy, splitted_perps, '1 Meters')

    linearray = arcpy.da.FeatureClassToNumPyArray(splitted_perps, ('PntID', 'OFID'))  
    pnt_ids = np.array([item[0] for item in linearray])
    unique_pnt_ids = np.unique(pnt_ids)

    perp_points = arcpy.CreateFeatureclass_management("in_memory", "perp_points","POINT", "","","", spatialref)
    arcpy.AddField_management(perp_points, 'PntID', 'Long', 6) ##PntID related to the flowline point ID
    arcpy.AddField_management(perp_points, 'SortID', 'Long', 6) ## SortID: right from 0 to max, left from -min to -1, so that it can be used to derive the F factor
    arcpy.AddField_management(perp_points, 'PointZ', 'Double', 6, 2)
    
    sel_splitted_perps = "in_memory\\sel_splitted_perps"
    sel_flowline_pnt = "in_memory\\sel_flowline_pnt"
    perp_points_cursor = arcpy.da.InsertCursor(perp_points, ('SHAPE@', 'pntID', "SortID"))
    
    for pnt_id in unique_pnt_ids:
        query = "PntID = " + str(pnt_id)
        arcpy.Select_analysis (splitted_perps, sel_splitted_perps, query)
        arcpy.Select_analysis (flowlinepointscopy, sel_flowline_pnt, query)
        pntGeometry = arcpy.CopyFeatures_management(sel_flowline_pnt, arcpy.Geometry())
        pointx = pntGeometry[0].firstPoint.X
        pointy = pntGeometry[0].firstPoint.Y

        ##determine Left and right lines
        with arcpy.da.SearchCursor(sel_splitted_perps, ["SHAPE@","SHAPE@LENGTH"]) as cursor:
            for line in cursor:
                pntx = line[0].firstPoint.X
                pnty = line[0].firstPoint.Y
                length = line[1]
                effect_length = min(width,length)
                if (abs(pointx - pntx) < 0.01 and abs(pointy - pnty) < 0.01): ##start from the flowline point
                    start = 0
                    i = 0
                    while start < effect_length:
                        new_point = line[0].positionAlongLine(start)
                        perp_points_cursor.insertRow([new_point, pnt_id, i])
                        start = start + resolution
                        i += 1
                else:  #end at the flowline point
                    num_pnts = int(effect_length/resolution)
                    if length < width:
                        leftover = length - resolution * num_pnts
                    else:
                        leftover = length - width 
                    for i in range(num_pnts):
                        start = leftover + i * resolution
                        new_point = line[0].positionAlongLine(start)
                        perp_points_cursor.insertRow([new_point, pnt_id, (i-num_pnts)])
        del line, cursor            
    del perp_points_cursor

    ##Step 7: Extract Z values to points
    perp_points_with_Z = "in_memory\\perp_points_with_Z"
    if error == False:
        ExtractValuesToPoints(perp_points, extractDEM, perp_points_with_Z, "INTERPOLATE", "VALUE_ONLY")
    else:
        ExtractValuesToPoints(perp_points, beddem, perp_points_with_Z, "INTERPOLATE", "VALUE_ONLY")

    arcpy.CalculateField_management(perp_points_with_Z, 'PointZ', str("!"+"RASTERVALU"+"!"), 'PYTHON_9.3')
    arcpy.DeleteField_management(perp_points_with_Z, 'RASTERVALU')

    #arcpy.CopyFeatures_management(perp_points_with_Z, "c:\\test\\perp_points_with_Z.shp")
    arcpy.env.extent = oldextent

    ##Need to make sure that the poinZ is not none value
    with arcpy.da.UpdateCursor(perp_points_with_Z, "PointZ") as cursor:
        i = 0
        for row in cursor:
            #arcpy.AddMessage(str(row[0]))
            if str(row[0]) == "None":
                #arcpy.AddMessage(str(row[0]))
                cursor.deleteRow()  ##Just remove the max elevation and keep the lowest elevation
    if i > 0:
        del row
    del cursor
    
    return perp_points_with_Z


#------------------------------------------------------------------------------------------------------------
# This function calculates the ice thickness of each point along the flowline based on the Excel flowline model
# introduced by the Benn and Houlton (2010). It is revised from the codes by Pellitero et al.(2016) in GlaRe.
# Revised on 2/14/2023 to fit the target Elevations
#------------------------------------------------------------------------------------------------------------
def Ice_Thickness_Calculation(flowPnts, min_ss, max_ss):
    Array = arcpy.da.TableToNumPyArray(flowPnts,('RASTERVALU','OFID','NEAR_FID','POINT_X','POINT_Y',"SSTRESS", "ffactor","ice","thick", "TargetElev"))
    ##Find the middle elevation of the target
    #Elevations = np.array([item[9] for item in Array])
    #validElevs = Elevations[Elevations > 0]
    #mid_elev = (max(validElevs) + min(validElevs))/2
    #elev_bin = (max(validElevs) - min(validElevs))/10 ##divide the elevation range to 10 parts and each parts with a weight values 

    groundElevs = np.array([item[0] for item in Array])
    mid_elev = (max(groundElevs) + min(groundElevs))/2
    elev_bin = (max(groundElevs) - min(groundElevs))/10 ##divide the elevation range to 10 parts and each parts with a weight values 

    bLoop = True
    nLoop = 0
    while (bLoop):
        TerminusDist= 0
        ids = []
        OFIDs = []
        ratios = []
        weight2mids = []
        for i in range(len(Array)):
            if TerminusDist == 0:
                Array[i][7]=Array[i][0] ##ice surface elevation
                Array[i][8]=Array [i][7]-Array[i][0] ##thick value
                TerminusDist=1
            elif Array[i][1] == Array[i-1][1]: ##Along a same flowline
                #try:
                d=math.sqrt(math.pow((Array[i][3]-Array[i-1][3]),2)+math.pow((Array[i][4]-Array[i-1][4]),2))#euclidean distance
                b = -((Array[i-1][0]+Array[i][0]))
                c = Array[i-1][7]*(Array[i][0]-((Array[i-1][7]-Array[i-1][0])))-((2*d*(((Array[i][5]+Array[i-1][5])/2)/Array[i][6]))/(900*9.81))
                Array[i][7]=(-b+math.pow((math.pow(b,2))-(4*c),0.5))/2  ##ice surface
                #except:
                #    arcpy.AddMessage(str(b))
                #    arcpy.AddMessage(str(c))
                    
                Array[i][8]=Array [i][7]-Array[i][0]  ##ice thickness
                if Array[i][8] < 5: ##if the calculated thickness is less than 5
                    #arcpy.AddMessage("Negative thickness!")
                    Array[i][8] = Array[i-1][8] ##set the thickness as the previous thickness
                    Array[i][7] = Array[i][0] + Array[i][8]
                ##Calculate the difference between calcualted ice thickness and target elevation
                #arcpy.AddMessage(Array[i][9])
                #arcpy.AddMessage(str(len(Array[i][9])))
                if len(str(Array[i][9])) > 1 and Array[i][9] != "None": ##not an empty string
                    #arcpy.AddMessage(Array[i][9])
                    targetThickness = (float(Array[i][9]) - Array[i][0])
                    #arcpy.AddMessage("Target thickness is: " + str(targetThickness))
                    if targetThickness > 0:
                        ratio = targetThickness/(Array[i][8]+0.01) ##ratio of the target thickness to the derived thickness
                        #diff = max((float(Array[i][9]) - Array[i][7]),0) ##if 
                        diff = (float(Array[i][9]) - Array[i][7]) ##if
                        weight2mid = 1/int((abs(float(Array[i][9]) - mid_elev)/elev_bin) + 1) 
                        #weight2mid = 1/(abs(float(Array[i][9]) - mid_elev) + 1) 
                        ids.append(i)
                        OFIDs.append(Array[i][1])
                        ratios.append(ratio)
                        weight2mids.append(weight2mid)
            else:
                Array[i][7]=Array[(Array[i][2]-1)][7] ##ice surface 
                if 'nan' in str(Array[i][7]).lower():
                    Array[i][7] = Array[i][0]
                if Array[i][8] != 0:
                    Array[i][8]=Array [i][7]-Array[i][0]
                else:
                    Array[i][7] = Array[i][0]  ##reset the ice surface to bedDEM

        ratioArr = np.array(ratios)
        #diffsArr = np.array(diffs)
        OFIDsArr = np.array(OFIDs)
        weightsArr = np.array(weight2mids)
        
        uniqueOFIDsArr = np.unique(OFIDsArr)
        #print uniqueOFIDsArr
        mean_ratio = []
        #mean_diff = []
        #abs_diff = []
        for OFID in uniqueOFIDsArr:
            arr = ratioArr[OFIDsArr == OFID]
            weightArr = weightsArr[OFIDsArr == OFID]
            
            #arcpy.AddMessage("OFID is: " + str(OFID))
            #arcpy.AddMessage(arr)
            '''
            ##Need to calculate weighted mean
            num = len(arr)
            #weightArr = np.arange(num, 0, -1)

            middle_num = int(num/2+1)
            weights = []
            for i in range(1, num+1):
                weight = (i-middle_num)*(i-middle_num) - middle_num*middle_num
                weights.append(weight)
            weightArr = np.array(weights)
            '''
            
            total = np.sum(weightArr)
            weightArr = weightArr / total
            #arcpy.AddMessage(weightArr)

            weighted_ratio = weightArr * arr
            #arcpy.AddMessage(weighted_ratio)
            
            mean_ratio.append(np.sum(weighted_ratio))
            #arr2 = diffsArr[OFIDsArr == OFID]
            #mean_diff.append(np.mean(arr2))
            #abs_diff.append(np.mean(np.absolute(arr2)))    

        mean_ratioArr = np.array(mean_ratio)
        max_ratio_diff = np.max(np.abs((mean_ratioArr - 1 )))
        #print max_ratio_diff
        #max_abs_diff = max(abs_diff)
        #arcpy.AddMessage("Mean ratio is: ")
        #arcpy.AddMessage(mean_ratio)

        #arcpy.AddMessage("Ratio difference from 1.0 is: ")
        #arcpy.AddMessage(max_ratio_diff)
        
        #print max_abs_diff
        ##Adjust the shear stress based on the difference
        #err = 0
        if (max_ratio_diff > 0.01 and nLoop < 5):
            for i in range(len(Array)):
                OFID = Array[i][1]
                idx_result = np.where(uniqueOFIDsArr == OFID)
                try:
                    idx = idx_result[0][0]
                    #print idx
                    ##update Shearstress
                    ss = Array[i][5] * mean_ratio[idx]
                    ss = max(min_ss, ss)
                    ss = min(ss, max_ss)
                    Array[i][5] = ss
                    #print Array[i][5]
                except:
                    #print "error"
                    #err += 1
                    pass
            nLoop += 1
        else:
            bLoop = False
        #print err

           
    
    #Save and update the new calculations to flowPnts
    count=0
    with arcpy.da.UpdateCursor(flowPnts,("ice","thick", "SSTRESS")) as cursor:
        for row in cursor:
            row[0]=Array[count][7]
            row[1]=Array[count][8]
            row[2]=Array[count][5]
            cursor.updateRow(row)
            count+=1
    #delete cursors
    del row, cursor



#------------------------------------------------------------------------------------------------------------
# This fuction is the whole process to reconstruct paleoice based on DEM, input flowlines, ice boundary, and default shear stress
#------------------------------------------------------------------------------------------------------------
def PaleoIceReconstructionwithboundary(BedDEM, inputflowline, Distance, iceboundary, shearstress, min_ss, max_ss, bFactorPolyfit, method, outpoints, outIceSurfaces, outIceThickness):

    GlacierID = "GlacierID" ##This is an ID field in inputflowline to identify the flowline(s) for each glacier (maybe connected with multiple flowlines)


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
    minArea = cellsize_float * cellsize_float * 5 ##5 pixels 
    with arcpy.da.UpdateCursor(icebndpolys, "SHAPE@AREA") as cursor:        
        for row in cursor:
            if row[0] < minArea:
                cursor.deleteRow()  ##Just remove the max elevation and keep the lowest elevation
    del cursor, row
    

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
    #All_CS_ice_points = "in_memory\\All_CS_ice_points"
    icesurs = arcpy.env.scratchFolder + "\\r" + "icesurs.tif" ##the inmemory does not work for raster


    arcpy.env.extent = BedDEM
    arcpy.env.cellSize = BedDEM
    arcpy.env.snapRaster = BedDEM ##setup snap raster
    
    ##2/16/2023: Donot need the watershed because the ice polygons are already provided!!!!!
    '''
    ##watershed delineations
    burninDEM = BedDEM - Power (cellsize_float / (cellsize_float + EucDistance(inputflowline) ), 2 ) * 10 ##Burn in the DEM to make sure the flow pass through the flowline start points
    ##Start to delineate the watershed
    #Hydro analysis
    fillDEM =Fill(burninDEM)  ##Fill the sink first
    fdir = FlowDirection(fillDEM,"NORMAL") ##Flow direction
        
    facc = FlowAccumulation(fdir) ##Flow accmulation    
    '''
    
    #exist_fields = [f.name for f in arcpy.ListFields(inputflowline)] #List of current field names in outline layer
    #if GlacierID in exist_fields:
    #    arcpy.DeleteField_management(flowline3dpoints,GlacierID)

    #exist_fields = [f.name for f in arcpy.ListFields(inputflowline)] #List of current field names in outline layer
    #if GlacierID not in exist_fields:
    #    arcpy.AddMessage("Assigning Glacier ID...")
    #    Add_GlacierID_by_Touches (inputflowline, GlacierID, flowlines)
    #else:
    #    arcpy.CopyFeatures_management(inputflowline, flowlines)

    #arcpy.AddMessage("Assigning Glacier ID...")
    #Add_GlacierID_by_Touches (inputflowline, GlacierID, flowlines)

    '''
    if GlacierID not in exist_fields:
        arcpy.AddMessage("Assigning Glacier ID...")
        Add_GlacierID_by_Touches (inputflowline, GlacierID, flowlines)
    else:
        arcpy.CopyFeatures_management(inputflowline, flowlines)
    '''
    ##Added on 3/10/2023
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

    if not("TargetElev" in fieldlist):
        arcpy.AddField_management(flowline3dpoints, "TargetElev", "DOUBLE", 6, 2)

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
        try:
            query = GlacierID + " = " + str(uniqueiceID[gid])
            arcpy.AddMessage("Processing #" + str(gid+1) +"/" + str(len(uniqueiceID)) + " of reconstructed glaciers...")                                                                                       
            arcpy.Select_analysis (flowlines, flowline, query)

            ###Select the flowline points corresponding to the flowlines
            arcpy.Select_analysis (flowline3dpoints, selflowline3dpoints, query)
            

            ###Merge the branches and start to run the ice thickness calculation 
            arcpy.Near_analysis (selflowline3dpoints, selflowline3dpoints, Distance * 2) ##setup a search radius for near analysis (may be problematic to choose the upstream near point; should select the nearest two points and use the average value!! or keep the intersect points as part of the flowpoints)
            arcpy.AddXY_management(selflowline3dpoints)

            ##Get the ice polygon for the selflowline3dpoints by spatial join
            arcpy.SpatialJoin_analysis(iceboundary, selflowline3dpoints, icepolyselect, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "#", "#")
            

            '''
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

            arcpy.Clip_analysis (ws, icebndpolys, icepolyselect)
            #remove the potential spurious polygons
            with arcpy.da.UpdateCursor(icepolyselect, "SHAPE@AREA") as cursor:
                i = 0
                for row in cursor:
                    i += 1
                    if row[0] < minArea:
                        cursor.deleteRow()  
            del cursor
            #arcpy.AddMessage("the number of icepolyselect is:" + str(i))
            if i> 0:
                del row
            '''
            #cross_section_pnts = cross_section_points(selflowline3dpoints, flowline, icepolyselect, BedDEM, 1600, 60)
            arcpy.PolygonToLine_management(icepolyselect, "in_memory\\selTargetFC")
            cross_section_pnts = cross_section_points_with_TF(selflowline3dpoints, flowline, icepolyselect, BedDEM, "in_memory\\selTargetFC", cellsize_float, 10000, Distance)

            ##Ice thickness calculation along the flowlines
            #Ice_Thickness_Calculation (selflowline3dpoints)

            ##Get the shape factors first with the known ice polygon and target elevation for each flowline point
            if bFactorPolyfit == True:
                arcpy.AddMessage("Deriving F factor based on the Polyfit of the cross section...")
                AdjustFfactor_Ployfit_with_cross_section_pnts (selflowline3dpoints, cross_section_pnts, "TargetElev", icepolyselect)
                ##replace the ffactor == 0.8 as the average F factor of the whole section
                
            else:
                arcpy.AddMessage("Deriving F factor from cross section...")
                #Try to use the curve fitting method to derive the F factor
                AdjustFfactor_with_cross_section_pnts (selflowline3dpoints, cross_section_pnts, "TargetElev", icepolyselect) 
            
            arcpy.AddMessage("Ice thickness calculation along the flowlines...")
            Ice_Thickness_Calculation (selflowline3dpoints, min_ss, max_ss)

            ##No need to do surface interpolation for individual outlines becasue this will be done at the end 03/14/2013
            #arcpy.AddMessage("Interpolating ice thickness and surface rasters...")
            #icesurface = IceSurfaceCalculation_with_crosssection_pnts_for_iceboundary (selflowline3dpoints, BedDEM, icepolyselect, "ice", method, "in_memory\\icesurface")

            ###2/16/2023: Do not need to figure out the best shear stress because the revise ice thickness calculation already did that
            '''
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
                #icesurface = IceSurfaceCalculation_with_crosssection_pnts_for_iceboundary (selflowline3dpoints, BedDEM, icepolyselect, "thick", "in_memory\\icesurface")
                #2/16/2023
                icesurface = IceSurfaceCalculation_with_crosssection_pnts_for_iceboundary (selflowline3dpoints, BedDEM, icepolyselect, "ice", method, "in_memory\\icesurface")

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
                    
                arcpy.AddMessage("The calculated shear stress is:" + str(ss) + " and the difference with the previous value is " + str(abs(ss-ss0)/ss0))
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
            '''
            
            if gid == 0: #The first time
                arcpy.CopyFeatures_management(selflowline3dpoints, outpoints)
                ##arcpy.CopyFeatures_management(CS_ice_points, All_CS_ice_points)
                #arcpy.CopyRaster_management(icesurface, icesurs)
            else:
                arcpy.Append_management(selflowline3dpoints, outpoints, "NO_TEST" )
                ##arcpy.Append_management(CS_ice_points, All_CS_ice_points, "NO_TEST")
                #arcpy.Mosaic_management(icesurface, icesurs, "MEAN","","", "", "", "", "")
        except:
            arcpy.AddMessage("THere is an error inthe ice thickness calcualtion. Move to the next one!")
            pass

    ##Smooth the polgyons and extract the ice surface within the ice polygon
    arcpy.AddMessage("Generating final outputs...")

    ##Dissolve polyygons so that the shared boundaries are removed
    disolveicepoly = "in_memory\\disolveicepoly"
    ##Make sure to use the polygons that corresponding to the flowline points for the analysis
    arcpy.SpatialJoin_analysis(icebndpolys, outpoints, "in_memory\\outicepolys", "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
    arcpy.Dissolve_management("in_memory\\outicepolys", disolveicepoly, '#', '#', 'SINGLE_PART', '#')

    arcpy.CopyFeatures_management(outpoints, "in_memory\\outpoints_cp")    
    fieldlist=[field.name for field in arcpy.ListFields("in_memory\\outpoints_cp")]
    if "RASTERVALU" in fieldlist:
        arcpy.DeleteField_management("in_memory\\outpoints_cp", 'RASTERVALU')

    '''
    ##create a narrow buffer around the flowlines and concert it to points
    arcpy.Buffer_analysis(flowlines, "in_memory\\flowline_buf", "60 Meter", "", "", "ALL")
    ##convert the buffer to points
    arcpy.FeatureVerticesToPoints_management("in_memory\\flowline_buf", "in_memory\\flowline_buf_points", 'ALL')
    
    arcpy.Near_analysis ("in_memory\\flowline_buf_points", "in_memory\\outpoints_cp")#identify nearest flowline point
    
    pointarray = arcpy.da.FeatureClassToNumPyArray("in_memory\\outpoints_cp",["OID@","ice"])#create array of ice values and populate a list with them
    fidval = np.array([item[0] for item in pointarray])
    iceval = np.array([item[1] for item in pointarray])

    ###Add field into the cross section points
    arcpy.AddField_management("in_memory\\flowline_buf_points", "ice", "FLOAT",10,6) #field for ice value

    with arcpy.da.UpdateCursor("in_memory\\flowline_buf_points",("NEAR_FID", "ice")) as cursor:   #populate ice field with value from the nearest flowline point
        for row in cursor:
            fid = row[0]
            idx_result = np.where(fidval == fid)
            idx = idx_result[0][0]
            row[1]=iceval[idx]
            cursor.updateRow(row)
    del row, cursor

    ##append the buffer points
    arcpy.Append_management("in_memory\\flowline_buf_points", "in_memory\\outpoints_cp", "NO_TEST")
    '''
    
    surface_interpolation ("in_memory\\outpoints_cp", "ice", BedDEM, disolveicepoly, method, outIceSurfaces)

    #surface_interpolation ("in_memory\\outpoints_cp", "ice", BedDEM, disolveicepoly, method, outIceSurfaces)

    ##Interpret ice surface based on ice thickness
    icethickness = Raster(outIceSurfaces) - BedDEM
    ##Should make sure that no negtive value for the thickness, is negative, set as zero?? Need to explore where are these negative values
    conIcethickness = Con(icethickness > 0, icethickness, 0) 
    #arcpy.CopyRaster_management (icethickness, outIceThickness)
    arcpy.CopyRaster_management (conIcethickness, outIceThickness)

    '''
    ##Output method: to be consistent and compariable with the Volta ice thickness program, first interpret the ice thickness and then reinterpret the ice surface
    ##if there is target features, using the target features as the zero contour lines?????
    singlepart_outlines = arcpy.MultipartToSinglepart_management(disolveicepoly, "in_memory\\singlepart_outlines")
    outline_lines_in = arcpy.PolygonToLine_management(singlepart_outlines, "in_memory\\outlines_line_in")
    arcpy.AddField_management(outline_lines_in, "contour", "SHORT")
    arcpy.CalculateField_management(outline_lines_in,"contour",0)
    cellsize_interp = cellsize_float
    interpolated_ice_depth = TopoToRaster([TopoPointElevation([[outpoints, 'thick']]), TopoContour([[outline_lines_in, 'contour']]), TopoBoundary ([singlepart_outlines])], cellsize_interp, "", '20', '0', '#', 'NO_ENFORCE', "SPOT", '1', '#', '1', '0', '0', '200') ##,"","","","",0.1)
    arcpy.CopyRaster_management (interpolated_ice_depth, outIceThickness)

    ##Get ice surface based on ice thickness
    icesurface_final = interpolated_ice_depth + BedDEM
    ##Should make sure that no negtive value for the thickness, is negative, set as zero?? Need to explore where are these negative values
    arcpy.CopyRaster_management (icesurface_final, outIceSurfaces)
    '''
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
    Interpolation_Method = arcpy.GetParameterAsText(8)
    outpoints=arcpy.GetParameterAsText(9)
    outIceSurfaces = arcpy.GetParameterAsText(10)
    outIceThickness = arcpy.GetParameterAsText(11)

    arcpy.Delete_management("in_memory") ### Empty the in_memory

    PaleoIceReconstructionwithboundary(BedDEM, inputflowline, Distance, inputiceboundary, shearstress, min_ss, max_ss, bFactorPolyfit, Interpolation_Method, outpoints, outIceSurfaces, outIceThickness)

    arcpy.Delete_management("in_memory") ### Empty the in_memory
