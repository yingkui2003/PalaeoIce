from __future__ import division
import locale
import os
#import operator
import arcpy
from arcpy import env
import numpy as np
from arcpy.sa import *

try:
    import numba
except:
    os.system("python -m pip install numba")
    #!pip install numba
    import numba

from numba import jit, prange
  
locale.setlocale(locale.LC_ALL,"")#sets local settings to decimals
arcpy.env.overwriteOutput = True

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

temp_workspace = "in_memory"  
if ArcGISPro:
    temp_workspace = "memory"

##Old codes  03/26/2025 no numba jit
def ELA_AAR_MGE(EleArr, interval, ratio):

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

##revised code with numba jit
#@jit(nopython=True)
@jit(nopython=True, parallel=True)
def ELA_AAR_MGE_jit(EleArr, interval, ratio):
    minimum = np.min(EleArr)
    maximum = np.max(EleArr)
    
    maxalt = int(maximum + interval)
    minalt = int(minimum - interval)

    # Create array of bin edges
    Elelist = np.arange(minalt, maxalt + interval, interval)
    
    # Calculate histogram
    H, X1 = np.histogram(EleArr, bins=Elelist)
    dx = X1[1] - X1[0]


    Area3D_arr = np.cumsum(H) * dx
    
    superf_total = np.max(Area3D_arr)  # Get the total surface
    Area3D_arr = superf_total - Area3D_arr

    ELA = superf_total * ratio  # Get the surface above the ELA
    kurowski = superf_total * 0.5

    # Find indices where values meet conditions
    ela_idx = -1
    kur_idx = -1
    min_ela_diff = np.inf
    min_kur_diff = np.inf
    
    for i in range(len(Area3D_arr)):
        val = Area3D_arr[i]
        
        # Check for ELA condition
        if val <= ELA:
            diff = ELA - val
            if diff < min_ela_diff:
                min_ela_diff = diff
                ela_idx = i
                
        # Check for Kurowski condition
        if val <= kurowski:
            diff = kurowski - val
            if diff < min_kur_diff:
                min_kur_diff = diff
                kur_idx = i

    # Calculate results
    ELA_AAR = Elelist[ela_idx] + (interval/2) + interval
    ELA_MGE = Elelist[kur_idx] + (interval/2) + interval
    
    return ELA_AAR, ELA_MGE

##Old codes  03/26/2025 no numba jit
def ELA_AA_AABR(EleArr, interval, ratio):
    minimum = np.min(EleArr)
    maximum = np.max(EleArr)
   
    maxalt=int(maximum+interval)
    minalt=int(minimum-interval)

    # Create a list of altitudes
    #num_bins = round((maxalt - minalt) / interval + 0.01) ##round up
    num_bins = int(np.ceil((maxalt - minalt) / interval))
    maxValue = minalt + interval * num_bins - interval/2
    
    list_altitudes = np.linspace(minalt + interval/2, maxValue, num_bins)

    #Elelist = range(minalt, maxalt, interval)
    Elelist = np.arange(minalt, maxalt, interval)
    # Create histogram bins
    #Elelist3 = np.linspace(minalt, minalt + interval * (num_bins-1), num_bins)

    H,X1 = np.histogram( EleArr, bins = Elelist, density = True )
    dx = X1[1] - X1[0]
    Area3D_arr = np.cumsum(H)*dx*100 ##times 100 to get the percentage

    arcpy.AddMessage(len(H))

    arcpy.AddMessage(len(Area3D_arr))    

    # AA Calculation
    superf_total=max(Area3D_arr) # Get the total surface

    resta=[int(x)-int(y) for (x,y) in zip(Area3D_arr[1:], Area3D_arr[0:])]
    arcpy.AddMessage(len(resta))    
    arcpy.AddMessage(len(list_altitudes))    
    arcpy.AddMessage(zip (resta,list_altitudes))
    multiplicacion=[int(x)*int (y) for (x,y) in zip (resta,list_altitudes)]

    finalmulti=sum(multiplicacion)

    ELA_AA=int(int(finalmulti)/int(superf_total)) + interval ##Add one interval to match the old AA value

    # AABR Calculation
    refinf=minalt
    valores_multi=[]
    valorAABR=[x*(y - refinf) for (x,y) in zip (resta, list_altitudes)]
    
    for valoracion in valorAABR:
        if valoracion<0:
            valores_multi.append(int (valoracion*ratio))
        else:
            valores_multi.append(int (valoracion))

    valorAABRfinal=sum (valores_multi)

    while valorAABRfinal > 0:
        refinf = refinf + interval
        valores_multi=[]
        valorAABR=[x*(y - refinf) for (x,y) in zip (resta, list_altitudes)]

        for valoracion in valorAABR:
            if valoracion < 0:
                valores_multi.append(valoracion*ratio)
            else:
                valores_multi.append(valoracion)

        valorAABRfinal=sum (valores_multi)

    ELA_AABR = refinf-(interval/2) + interval ##Add one interval to match the old AA value

    return ELA_AA, ELA_AABR 

##revised code with numba jit
@jit(nopython=True)
def ELA_AA_AABR_jit(EleArr, interval, AABRratio):
    minimum = np.min(EleArr)
    maximum = np.max(EleArr)
   
    maxalt = int(maximum + interval)
    minalt = int(minimum - interval)

    # Create array of altitudes using numpy instead of list
    num_bins = round((maxalt - minalt) / interval + 0.01) ##round up
    maxValue = minalt + interval * num_bins - interval/2
    list_altitudes = np.linspace(minalt + interval/2, maxValue, num_bins)
    
    # Create bin edges for histogram
    Elelist = np.arange(minalt, maxalt, interval)

    # Calculate histogram
    H, X1 = np.histogram(EleArr, bins=Elelist)
    dx = X1[1] - X1[0]
    
    # Calculate cumulative sum
    Area3D_arr = np.cumsum(H) * dx * 100  # times 100 to get percentage
    
    # AA Calculation
    superf_total = np.max(Area3D_arr)  # Get the total surface
    
    # Calculate differences between consecutive elements
    resta = np.diff(Area3D_arr)
    
    # Calculate weighted sum
    finalmulti = np.sum(resta * list_altitudes[1:])
    
    ELA_AA = int(finalmulti / superf_total) + interval
    
    # AABR Calculation
    refinf = minalt
    valorAABRfinal = np.inf  # Initialize with large value
    
    while valorAABRfinal > 0:
        # Calculate weighted differences
        valorAABR = resta * (list_altitudes[1:] - refinf)
        
        # Apply ratio to negative values
        valores_multi = np.where(valorAABR < 0, 
                                valorAABR * AABRratio, 
                                valorAABR)
        
        valorAABRfinal = np.sum(valores_multi)
        
        if valorAABRfinal > 0:
            refinf += interval

    ELA_AABR = refinf - (interval/2) + interval
    
    return ELA_AA, ELA_AABR

##revised code with numba jit
#@jit(nopython=True)
@jit(nopython=True, parallel=True)
def ELA_AA_AABR_jit2(EleArr, interval, AABRratio):
    # Calculate min/max with buffer
    minimum = np.min(EleArr)
    maximum = np.max(EleArr)
    maxalt = int(maximum + interval)
    minalt = int(minimum - interval)

    # Optimized bin calculation
    num_bins = int(np.ceil((maxalt - minalt) / interval))
    maxValue = minalt + interval * num_bins - interval/2
    list_altitudes = np.linspace(minalt + interval/2, maxValue, num_bins)
    
    # Create histogram bins
    Elelist = np.linspace(minalt, minalt + interval * (num_bins-1), num_bins)
    
    # Calculate histogram and cumulative area
    H, X1 = np.histogram(EleArr, bins=Elelist)
    dx = X1[1] - X1[0]
    Area3D_arr = np.cumsum(H) * dx * 100  # Convert to percentage

    # AA Calculation
    superf_total = np.max(Area3D_arr)
    resta = np.diff(Area3D_arr)
    finalmulti = np.sum(resta * list_altitudes[1:-1])
    ELA_AA = int(finalmulti / superf_total) ##+ interval
    
    # Optimized AABR Calculation
    refinf = minalt
    while True:
        # Vectorized calculation
        diff = list_altitudes[1:-1] - refinf
        weighted = resta * diff
        adjusted = np.where(weighted < 0, weighted * AABRratio, weighted)
        total = np.sum(adjusted)
        
        if total <= 0:
            break
        refinf += interval
    
    ELA_AABR = refinf - (interval/2) ##+ interval
    
    return ELA_AA, ELA_AABR


##Combined codes with numba jit
@jit(nopython=True, parallel=True)
def calculate_ELAs(EleArr, interval, AARratio, AABRratio):
    # Common calculations for all methods
    minimum = np.min(EleArr)
    maximum = np.max(EleArr)
    maxalt = int(maximum + interval)
    minalt = int(minimum - interval)

    # Optimized bin calculation
    num_bins = int(np.ceil((maxalt - minalt) / interval))
    maxValue = minalt + interval * num_bins - interval/2
    Elelist = np.linspace(minalt, minalt + interval * (num_bins-1), num_bins)

    # Create histogram bins
    H, X1 = np.histogram(EleArr, bins=Elelist)
    dx = X1[1] - X1[0]

    Area3D_arr = np.cumsum(H) * dx
    superf_total = np.max(Area3D_arr)
    Area3D_arr_rev = superf_total - Area3D_arr  # Reversed for AAR/MGE

    # --- AA/AABR Method ---
    list_altitudes = np.linspace(minalt + interval/2, maxValue, num_bins)

    resta = np.diff(Area3D_arr)
    
    # AA Calculation
    finalmulti = np.sum(resta * list_altitudes[1:-1])
    ELA_AA = int(finalmulti / superf_total) ##+ interval
    
    # AABR Calculation
    refinf = minalt
    valorAABRfinal = np.inf
    
    while valorAABRfinal > 0:
        valorAABR = resta * (list_altitudes[1:-1] - refinf)
        valores_multi = np.where(valorAABR < 0, 
                               valorAABR * AABRratio, 
                               valorAABR)
        valorAABRfinal = np.sum(valores_multi)
        if valorAABRfinal > 0:
            refinf += interval
    
    ELA_AABR = refinf - (interval/2) ##+ interval
    
    # --- AAR/MGE Method ---
    ELA = superf_total * AARratio
    kurowski = superf_total * 0.5
    
    # Find closest values in single pass
    ela_idx = aar_idx = mge_idx = 0
    min_ela_diff = min_aar_diff = min_mge_diff = np.inf
    
    for i in range(len(Area3D_arr_rev)):
        val = Area3D_arr_rev[i]
        
        # For AAR
        if val <= ELA:
            diff = ELA - val
            if diff < min_ela_diff:
                min_ela_diff = diff
                ela_idx = i
                
        # For MGE
        if val <= kurowski:
            diff = kurowski - val
            if diff < min_mge_diff:
                min_mge_diff = diff
                mge_idx = i
    
    ELA_AAR = Elelist[ela_idx] + (interval/2) + interval
    ELA_MGE = Elelist[mge_idx] + (interval/2) + interval
    
    return ELA_AA, ELA_AABR, ELA_AAR, ELA_MGE

# Define the input parameters
dem=arcpy.GetParameterAsText(0)
Glaciers = arcpy.GetParameterAsText(1)
interval= int(arcpy.GetParameterAsText(2))
AARsr=arcpy.GetParameterAsText(3)
AARratio=locale.atof(AARsr)
AABRsr=arcpy.GetParameterAsText(4)
AABRratio=locale.atof(AABRsr)
OutContours = arcpy.GetParameterAsText(5)

Flist = []
ListFields=arcpy.ListFields(Glaciers)
for x in ListFields:
    Flist.append(x.baseName)
    
if "MGE" in Flist:
    pass
else:
    arcpy.AddField_management(Glaciers, "MGE", "LONG",6)


if "AAR" in Flist:
    pass
else:
    arcpy.AddField_management(Glaciers, "AAR", "LONG",6)

if "AA" in Flist:
    pass
else:
    arcpy.AddField_management(Glaciers, "AA", "LONG",6)

if "AABR" in Flist:
    pass
else:
    arcpy.AddField_management(Glaciers, "AABR", "LONG",6)


ELAlines = arcpy.CreateFeatureclass_management(temp_workspace, "ELAlines", "POLYLINE", "","","",Glaciers)
arcpy.AddField_management(ELAlines, "LineID", 'Long', 6)
arcpy.AddField_management(ELAlines, "ELA", 'Long', 6)
arcpy.AddField_management(ELAlines, "Method", 'Text')
arcpy.AddField_management(ELAlines, "Ratio", 'DOUBLE',6, 3)

ELAline = temp_workspace + "\\ELAline"

ELAs = []
methods = []
ratios = []
fields = ("MGE","AAR","AA","AABR","SHAPE@")
with arcpy.da.UpdateCursor(Glaciers, fields) as cursor:
    i = 0
    j = 0
    for row in cursor:
        arcpy.AddMessage("Processing glacier #" + str(i+1))
        galcierDEM = ExtractByMask(dem, row[4])

        array = arcpy.RasterToNumPyArray(galcierDEM,"","","",0)
        EleArr = array[array > 0].astype(int) ##Get the elevations greater than zero

        
        ##MGE and AAR
        #ela_aar, ela_mge = ELA_AAR_MGE(EleArr, interval, AARratio)
        ela_aa, ela_AABR, ela_aar, ela_mge = calculate_ELAs(EleArr, interval, AARratio, AABRratio)
        #ela_aar, ela_mge = ELA_AAR_MGE_jit(EleArr, interval, AARratio)
        row[0] = ela_mge
        row[1] = ela_aar

        #Generate ELA contourlines
        ELAs.append(ela_mge)
        ELAs.append(ela_aar)
        methods.append("MGE")
        methods.append("AAR")
        ratios.append(0.5)
        ratios.append(AARratio)
        
        Contour(galcierDEM, ELAline, 10000, ela_mge)
        arcpy.AddField_management(ELAline, "LineID", 'Long', 6)
        arcpy.CalculateField_management(ELAline,"LineID",j,"PYTHON_9.3")
        arcpy.Append_management(ELAline, ELAlines, "NO_TEST")
        j += 1
        
        Contour(galcierDEM, ELAline, 10000, ela_aar)
        arcpy.AddField_management(ELAline, "LineID", 'Long', 6)
        arcpy.CalculateField_management(ELAline,"LineID",j,"PYTHON_9.3")
        arcpy.Append_management(ELAline, ELAlines, "NO_TEST")        
        j += 1
        
        ##AA and AABR
        #ela_aa, ela_AABR = ELA_AA_AABR(EleArr, interval, AABRratio)
        #ela_aa, ela_AABR = ELA_AA_AABR_jit2(EleArr, interval, AABRratio)
        row[2] = ela_aa
        row[3] = ela_AABR

        ELAs.append(ela_aa)
        ELAs.append(ela_AABR)
        methods.append("AA")
        methods.append("AABR")
        ratios.append(1.0)
        ratios.append(AABRratio)

        Contour(galcierDEM, ELAline, 10000, ela_aa)
        arcpy.AddField_management(ELAline, "LineID", 'Long', 6)
        arcpy.CalculateField_management(ELAline,"LineID",j,"PYTHON_9.3")
        arcpy.Append_management(ELAline, ELAlines, "NO_TEST")        
        j += 1

        Contour(galcierDEM, ELAline, 10000, ela_AABR)
        arcpy.AddField_management(ELAline, "LineID", 'Long', 6)
        arcpy.CalculateField_management(ELAline,"LineID",j,"PYTHON_9.3")
        arcpy.Append_management(ELAline, ELAlines, "NO_TEST")        
        j += 1

        cursor.updateRow(row)
        i += 1
del row, cursor

fields = ("LineID","ELA", "Method", "Ratio")
with arcpy.da.UpdateCursor(ELAlines, fields) as cursor:
    for row in cursor:
        lineid = row[0]
        row[1] = ELAs[lineid] 
        row[2] = methods[lineid] 
        row[3] = ratios[lineid]
        cursor.updateRow(row)
del row, cursor

arcpy.CopyFeatures_management(ELAlines, OutContours)
arcpy.Delete_management(temp_workspace) 
