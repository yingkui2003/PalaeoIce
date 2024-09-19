from __future__ import division
import locale
import os
#import operator
import arcpy
from arcpy import env
import numpy as np
from arcpy.sa import *
  
locale.setlocale(locale.LC_ALL,"")#sets local settings to decimals
arcpy.env.overwriteOutput = True
#arcpy.CheckOutExtension("Spatial")
#arcpy.CheckOutExtension("3D")

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


def ELA_AAR_MGE(dem, interval, ratio):
    #Maximum and minimum value of DEM
    upper=arcpy.GetRasterProperties_management(dem,"MAXIMUM")
    maximum = float(upper.getOutput(0))
    lower=arcpy.GetRasterProperties_management(dem,"MINIMUM")
    minimum = float(lower.getOutput(0))

    maxalt=int(maximum+interval)
    minalt=int(minimum-interval)

    volumetable = arcpy.env.scratchFolder + "\\xxELAtable.txt"

    # Create list of altitudes and populate primervalor
    Elelist = range(minalt, maxalt, interval)

    for plane in Elelist:
        arcpy.SurfaceVolume_3d(dem, volumetable, "ABOVE", plane)

    ##Step 2: Read the volume table for 3D area
    arr=arcpy.da.TableToNumPyArray(volumetable, ('AREA_3D'))
    Area3D_arr = np.array([item[0] for item in arr])

    ##Delete the volume table
    arcpy.Delete_management(volumetable)


    superf_total=max(Area3D_arr) # Get the total surface
    ELA=superf_total * ratio # Get the surface above the ELA
    kurowski= superf_total * 0.5

    # Create a list of the altitudes whose surface is less than ELA
    superf_en_ELA=[]
    superf_kurowski=[]
    for values in Area3D_arr:
        if values <= ELA and values<= kurowski:
            superf_en_ELA.append(values)
            superf_kurowski.append(values)
        elif values <= ELA and values> kurowski:
            superf_en_ELA.append(values)
        elif values > ELA and values<= kurowski:
            superf_kurowski.append(values)
        else:
            pass

    # Get the maximum surface value within the list
    ela=max(superf_en_ELA)
    kur=max(superf_kurowski)

    #arcpy.AddMessage("ela value is" + str(ela))
    #arcpy.AddMessage("kur value is" + str(kur))


    idx_result = np.where(Area3D_arr == ela)
    idx = idx_result[0][0]
    ELA_AAR=Elelist[idx]+(interval/2)
    #arcpy.AddMessage("The ELA AAR is:")
    #arcpy.AddMessage(ELA_AAR)

    idx_result = np.where(Area3D_arr == kur)
    idx = idx_result[0][0]
    ELA_MGE=Elelist[idx]+(interval/2)
    #arcpy.AddMessage("The ELA MGE is:")
    #arcpy.AddMessage(ELA_MGE)
    
    return ELA_AAR, ELA_MGE 

def ELA_AAR_MGE(EleArr, interval, ratio):

    #array = arcpy.RasterToNumPyArray(dem,"","","",0)
    #EleArr = array[array > 0].astype(int) ##Get the elevations greater than zero
    minimum = np.min(EleArr)
    maximum = np.max(EleArr)
      

    maxalt=int(maximum+interval)
    minalt=int(minimum-interval)

    # Create list of altitudes and populate primervalor
    Elelist = range(minalt, maxalt, interval)

    #H,X1 = np.histogram( EleArr, bins = Elelist, normed = True )
    H,X1 = np.histogram( EleArr, bins = Elelist, density = True )
    dx = X1[1] - X1[0]
    Area3D_arr = np.cumsum(H)*dx
    
    superf_total=max(Area3D_arr) # Get the total surface
    Area3D_arr = superf_total - Area3D_arr

    ELA=superf_total * ratio # Get the surface above the ELA
    kurowski= superf_total * 0.5

    # Create a list of the altitudes whose surface is less than ELA
    superf_en_ELA=[]
    superf_kurowski=[]
    for values in Area3D_arr:
        if values <= ELA and values<= kurowski:
            superf_en_ELA.append(values)
            superf_kurowski.append(values)
        elif values <= ELA and values> kurowski:
            superf_en_ELA.append(values)
        elif values > ELA and values<= kurowski:
            superf_kurowski.append(values)
        else:
            pass

    # Get the maximum surface value within the list
    ela=max(superf_en_ELA)
    kur=max(superf_kurowski)


    idx_result = np.where(Area3D_arr == ela)
    idx = idx_result[0][0]
    ELA_AAR=Elelist[idx]+(interval/2) + interval ##Add one interval to match the old AA value by Yingkui 10/08/2023

    idx_result = np.where(Area3D_arr == kur)
    idx = idx_result[0][0]
    ELA_MGE=Elelist[idx]+(interval/2) + interval ##Add one interval to match the old AA value
    
    return ELA_AAR, ELA_MGE 

def ELA_AA_AABR(EleArr, interval, ratio):
    #array = arcpy.RasterToNumPyArray(dem,"","","",0)
    #EleArr = array[array > 0].astype(int) ##Get the elevations greater than zero
    minimum = np.min(EleArr)
    maximum = np.max(EleArr)
   
    maxalt=int(maximum+interval)
    minalt=int(minimum-interval)

 
    # Create a list of altitudes
    list_altitudes=[]
    start_altitude=minalt+(interval/2)
    while start_altitude > minalt and start_altitude < maxalt:
        list_altitudes.append(start_altitude)
        start_altitude=start_altitude+interval

    Elelist = range(minalt, maxalt, interval)

    
    #H,X1 = np.histogram( EleArr, bins = Elelist, normed = True )
    H,X1 = np.histogram( EleArr, bins = Elelist, density = True )
    dx = X1[1] - X1[0]
    Area3D_arr = np.cumsum(H)*dx*100 ##times 100 to get the percentage
    #arcpy.AddMessage(Area3D_arr)

    # AA Calculation
    superf_total=max(Area3D_arr) # Get the total surface
    #arcpy.AddMessage(str(superf_total))

    resta=[int(x)-int(y) for (x,y) in zip(Area3D_arr[1:], Area3D_arr[0:])]

    multiplicacion=[int(x)*int (y) for (x,y) in zip (resta,list_altitudes)]

    finalmulti=sum(multiplicacion)

    ELA_AA=int(int(finalmulti)/int(superf_total)) + interval ##Add one interval to match the old AA value
    #arcpy.AddMessage("ELA_AA is" + str(ELA_AA))

    # AABR Calculation
    refinf=minalt
    valores_multi=[]
    valorAABR=[x*(y - refinf) for (x,y) in zip (resta, list_altitudes)]
    
    for valoracion in valorAABR:
        if valoracion<0:
            valores_multi.append(int (valoracion*ratio))
        else:
            valores_multi.append(int (valoracion))

    #arcpy.AddMessage(valoracion)

    valorAABRfinal=sum (valores_multi)

    while valorAABRfinal > 0:
        refinf = refinf + interval
        valores_multi=[]
        valorAABR=[x*(y - refinf) for (x,y) in zip (resta, list_altitudes)]
        #arcpy.AddMessage(valorAABR)

        for valoracion in valorAABR:
            if valoracion < 0:
                valores_multi.append(valoracion*ratio)
            else:
                valores_multi.append(valoracion)

        valorAABRfinal=sum (valores_multi)
        #arcpy.AddMessage(valorAABRfinal)

    ELA_AABR = refinf-(interval/2) + interval ##Add one interval to match the old AA value

    
    return ELA_AA, ELA_AABR 

# Define the input parameters
dem=arcpy.GetParameterAsText(0)
Glaciers = arcpy.GetParameterAsText(1)
interval= int(arcpy.GetParameterAsText(2))
AARsr=arcpy.GetParameterAsText(3)
AARratio=locale.atof(AARsr)
AABRsr=arcpy.GetParameterAsText(4)
AABRratio=locale.atof(AABRsr)

Flist = []
#outglaciers = "in_memory\\primary_flowline"
#arcpy.CopyFeatures_management(Glaciers, outglaciers)
#a "list" where the name of fields from the attributed table are copied in
ListFields=arcpy.ListFields(Glaciers)
for x in ListFields:
    Flist.append(x.baseName)
    
#arcpy.AddMessage(Flist)
#if fields with the same name are already in the attribute table, no need to create them again (which will crash the tool)

if "MGE" in Flist:
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    #arcpy.AddMessage("pass1")
    arcpy.AddField_management(Glaciers, "MGE", "LONG",6)


if "AAR" in Flist:
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(Glaciers, "AAR", "LONG",6)

if "AA" in Flist:
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(Glaciers, "AA", "LONG",6)

if "AABR" in Flist:
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(Glaciers, "AABR", "LONG",6)


fields = ("MGE","AAR","AA","AABR","SHAPE@")
#arcpy.AddMessage(fields)

i = 1
with arcpy.da.UpdateCursor(Glaciers, fields) as cursor:
    #arcpy.AddMessage("pass2")
    for row in cursor:
        arcpy.AddMessage("Processing glacier #" + str(i))
        galcierDEM = ExtractByMask(dem, row[4])

        array = arcpy.RasterToNumPyArray(galcierDEM,"","","",0)
        EleArr = array[array > 0].astype(int) ##Get the elevations greater than zero


        ela_aar, ela_mge = ELA_AAR_MGE(EleArr, interval, AARratio)
        row[0] = ela_mge
        row[1] = ela_aar

        ela_aa, ela_AABR = ELA_AA_AABR(EleArr, interval, AABRratio)
        row[2] = ela_aa
        row[3] = ela_AABR

        cursor.updateRow(row)
        i += 1
del row, cursor


