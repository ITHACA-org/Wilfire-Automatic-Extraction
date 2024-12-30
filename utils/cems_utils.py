import arcpy
import os
import utm
import tempfile
from utils import globVar as glob
import geopandas as gpd

def getAprx(path='current'):
    aprx = arcpy.mp.ArcGISProject(path)
    
    return aprx

def getListMaps(path='current'):
    aprx = getAprx(path)
    listMaps = aprx.listMaps()
    
    return listMaps

def getMap(mapName='*', path='current'):
    aprx = getAprx(path)
    map = aprx.listMaps(mapName)[0]
    
    return map

def getMfLayers(path='current'):
    map = getMap(glob.MapDisp, path)[0]
    listLayers = map.listLayers()
    
    return listLayers

def getListLayers(mapName='*', path='current'):
    map = getMap(mapName, path)
    listLayers = map.listLayers()
    
    return listLayers

def getLayer(layerName, mapName='*', path='current'):
    map = getMap(mapName, path)
    layer = map.listLayers(layerName)[0]
    
    return layer

def getAoiLayer(path='current'):
    map = getMap(glob.MapDisp, path)
    layer = map.listLayers(glob.Aoi)[0]
    
    return layer

def addLayer(layer_path, mapName='*', path='current'):
    aprx = arcpy.mp.ArcGISProject(path)
    map = getMap(mapName)
    map.addDataFromPath(layer_path)

def removeLayer(layer_name, mapName='*', path='current'):
    aprx = arcpy.mp.ArcGISProject(path)
    map = getMap(mapName)
    layer = [layer for layer in map.listLayers() if layer.name == layer_name][0]
    map.removeLayer(layer)

def getLayerExt(layer):
    try:
        descDS = arcpy.Describe(layer.dataSource)
        extent = descDS.extent
        return extent
        
    except:
        arcpy.env.addOutputsToMap = True
        arcpy.management.MakeFeatureLayer(layer, 'layer')
        
        layer = getLayer('layer')
        descDS = arcpy.Describe(layer.dataSource)
        extent = descDS.extent
        removeLayer(layer)
        arcpy.env.addOutputsToMap = False
        return extent

def createTempFolder():
    out_folder = tempfile.mkdtemp()
    return out_folder
    
def createTempGdb():
    out_folder = createTempFolder()
    arcpy.CreateFileGDB_management(out_folder, "temp.gdb")
    int_gdb = os.path.join(out_folder, "temp.gdb")
    return int_gdb

def appendData(destFc, appFc, delFcs):
    arcpy.DeleteFeatures_management(destFc)
    arcpy.Append_management(appFc, destFc, schema_type="NO_TEST")
    arcpy.management.Delete(delFcs, 'FeatureClass')

def cutAlongBorder(self, fcA, mask, gdbPath):       
    fcClipped = os.path.join(gdbPath, 'Clipped')
    fcErased = os.path.join(gdbPath, 'Erased')
    arcpy.analysis.Clip(fcA, mask, fcClipped)
    arcpy.analysis.Erase(fcA, mask, fcErased)
    self.appendData(fcA, [fcClipped, fcErased], [fcClipped, fcErased])
        
def getUTMZone(layExt):
    """Getting the UTM zone from the AOI extent."""
    emisf = "327"
    x = (layExt.XMax + layExt.XMin) / 2
    y = (layExt.YMax + layExt.YMin) / 2
    if x < -180 or x > 180:
        arcpy.AddError("Error in getting coordinates of the center of the FC")
    if y >= 0:
        emisf = "326"
    utmZone = str(utm.from_latlon(y, x)[-2])
    utmCode = int("{}{}".format(emisf,utmZone))
    
    return utmCode 

def getUTMZoneGpd(featureClass):
    
    tempFolder = createTempFolder()
    Lyrjsonpath = os.path.join(tempFolder, "layer.geojson")
    arcpy.conversion.FeaturesToJSON(featureClass, Lyrjsonpath, geoJSON=True)
    lyrGdf = gpd.read_file(Lyrjsonpath)
    extent = lyrGdf.bounds
        
    emisf = "327"
    x = (extent['maxx']+ extent['minx']) / 2
    y = (extent['maxy']+ extent['miny']) / 2
    if x.item() < -180 or x.item() > 180:
            arcpy.AddError("Error in getting coordinates of the center of the FC")
    if y.item() >= 0:
        emisf = "326"
    utmZone = str(int(utm.from_latlon(y.item(), x.item())[-2]))
    utmCode = "{}{}".format(emisf,utmZone)
    
    return utmCode 