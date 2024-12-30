
for layer in QgsProject.instance().mapLayers().values():
    layerPath = layer.source()[:-4]
    sldPath = "{}.sld".format(layerPath)
    layer.loadSldStyle(sldPath)