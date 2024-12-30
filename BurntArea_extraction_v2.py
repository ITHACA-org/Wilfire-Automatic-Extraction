import rasterio
import geopandas as gpd
import json
from rasterio.enums import Resampling
from rasterio.mask import mask 
import arcpy
import os
import re
from utils import cems_utils as cems
from rasterio.features import shapes
from shapely.geometry import shape
from shapely.geometry import Polygon, MultiPolygon
import tkinter as tk
import time

class Classificator:

    def __init__(self):
        self.pre_image_arcgis = arcpy.GetParameterAsText(0)
        self.post_image_arcgis = arcpy.GetParameterAsText(1)
        self.clipping_feature_arcgis = arcpy.GetParameterAsText(2)
        self.geodatabase = arcpy.GetParameterAsText(11)
        self.extraction_method = arcpy.GetParameterAsText(3)
        self.pre_image_red = arcpy.GetParameter(4)
        self.pre_image_NIR = arcpy.GetParameter(5)
        self.pre_image_SWIR = arcpy.GetParameter(6)
        self.post_image_red = arcpy.GetParameter(7)
        self.post_image_NIR = arcpy.GetParameter(8)
        self.post_image_SWIR = arcpy.GetParameter(9)
        self.MMU = arcpy.GetParameter(10)
        self.ObjectDescription = arcpy.GetParameter(12)
        self.DamageSourceIdentifier = arcpy.GetParameter(13)
        self.ObjectDescriptionDict = {"Forest Fire": 81, "Land fire: Brush, bush, Pasture": 82, "Urban Fires": 83}
        self.pre_image = rasterio.open(self.pre_image_arcgis) # turn into rasterio object
        self.post_image = rasterio.open(self.post_image_arcgis)
        self.obs_event = os.path.join(self.geodatabase, "B1_observed_event_a")
        self.temp_folder = cems.createTempFolder() # create a temp folder where to store intermediate outputs
        self.aoi = os.path.join(self.geodatabase, "A1_area_of_interest_a")
        self.crs = int(cems.getUTMZoneGpd(self.aoi))
        self.clipping_feature = self.ReprojectClippingFeatureShp()
        self.clipping_feature = self.RefineClippingFeature() # only in case of MONITs
        arcpy.AddMessage("Clipping data to your custom extent...")
        self.pre_image_clip, self.pre_image_clip_transform = self.ClipImage(self.pre_image, "PreImageClip.tif")
        self.post_image_clip, self.post_image_clip_transform = self.ClipImage(self.post_image, "PostImageClip.tif")
        self.pre_image_clip_resampled, self.resampling_condition, self.transform = self.Resampling()
        self.index = self.CalculateIndex()
        self.classification_vector, self.threshold = self.ExtractClassificationVector()
        self.classification_vector_cleaned = self.CleanClassificationVector()
    
    def ReprojectClippingFeatureShp(self):
        aprx = arcpy.mp.ArcGISProject('CURRENT')
        map = aprx.listMaps('Map')[0]
        layer = map.listLayers(str(self.clipping_feature_arcgis))[0]
        path = layer.dataSource
        match = re.match(r'^(.*)\\([^\\]+)$', path)
        gdb_path = match.group(1)  
        layer_name = match.group(2)
        try:
            clipping_feature_shp = gpd.read_file(gdb_path, layer = layer_name)
            if 'DateTime' in clipping_feature_shp.columns:
                clipping_feature_shp = clipping_feature_shp.drop(columns = ['DateTime'])
            clipping_feature_shp = clipping_feature_shp.to_crs(self.crs)
            clipping_feature_shp_path = os.path.join(self.temp_folder, self.clipping_feature_arcgis + '.shp')
            clipping_feature_shp.to_file(clipping_feature_shp_path, driver= "ESRI Shapefile")
        except:
            arcpy.AddError("Only feature class objects are accepted as clipping feature!")
        return clipping_feature_shp
    
    def RefineClippingFeature(self):
        B1_obs_event_gdf = gpd.read_file(self.geodatabase, layer = "B1_observed_event_a")
        B1_obs_event_gdf = B1_obs_event_gdf.to_crs(self.crs)
        if B1_obs_event_gdf.shape[0] > 0:
            shapely_objects = B1_obs_event_gdf['geometry']
            erasing_feature = self.FillHoles(shapely_objects)
            erased_clipping_feature = self.clipping_feature.overlay(erasing_feature, how = 'difference')
            erased_buffered_clipping_feature = erased_clipping_feature.copy()
            erased_buffered_clipping_feature['geometry'] = erased_clipping_feature['geometry'].buffer(10)
            return erased_buffered_clipping_feature
        else: return self.clipping_feature

    def FillHoles(self, shapely_objects, filling_area = 100000):
        list_interiors = []
        for geom in shapely_objects:
            # Initialize a list to collect new polygons after filtering interiors
            filtered_polygons = []
            if isinstance(geom, Polygon):
                # Handle single Polygon
                filtered_interiors = []
                for interior in geom.interiors:
                    p = Polygon(interior)
                    if p.area > filling_area:
                        filtered_interiors.append(interior)
                new_polygon = Polygon(geom.exterior.coords, holes=filtered_interiors)
                filtered_polygons.append(new_polygon)
            elif isinstance(geom, MultiPolygon):
                # Handle MultiPolygon by processing each Polygon within it
                for poly in geom.geoms:
                    filtered_interiors = []
                    for interior in poly.interiors:
                        p = Polygon(interior)
                        if p.area > filling_area:
                            filtered_interiors.append(interior)
                    new_polygon = Polygon(poly.exterior.coords, holes=filtered_interiors)
                    filtered_polygons.append(new_polygon)
            # If at the end of the iteration one polygon is split in more parts save it as a Multipolygon, otherwise as a simple Polygon 
            if len(filtered_polygons) > 1:
                list_interiors.append(MultiPolygon(filtered_polygons))
            else:
                list_interiors.append(filtered_polygons[0])
        filled_vector  = gpd.GeoDataFrame({'geometry': list_interiors}, crs=self.crs)
        return filled_vector
   
    def ClipImage(self, image, clip_name: str):
        def getFeatures(gdf):
            return [json.loads(gdf.to_json())["features"][0]["geometry"]]
        coords = getFeatures(self.clipping_feature)
        image_clip, image_clip_transform = mask(dataset = image , shapes = coords, crop = True, nodata = 0)
        profile = image.meta.copy()
        profile.update({"driver": "GTiff", 'width' : image_clip.shape[2], 'height': image_clip.shape[1], 'transform': image_clip_transform, 'crs': self.crs, 'nodata': 0})
        image_clip_path = os.path.join(self.temp_folder, clip_name)
        with rasterio.open(image_clip_path, "w", **profile) as dst:
            dst.write(image_clip)
        return image_clip_path, image_clip_transform
    
    def Resampling(self):
        pre_image_clip = rasterio.open(self.pre_image_clip)
        post_image_clip = rasterio.open(self.post_image_clip)
        pre_image_res = pre_image_clip.res[0]
        post_image_res = post_image_clip.res[0]
        resampling_condition = True if pre_image_res != post_image_res else False
        if resampling_condition:
            arcpy.AddMessage("Pre and Post images have different resolutions! Rescaling Pre image to the resolution of Post image...")
            def ResampledBand(band_number: int):
                scaling_factor = pre_image_res / post_image_res
                resampled_band = pre_image_clip.read(band_number, resampling = Resampling.bilinear, out_shape = (int(pre_image_clip.height * scaling_factor), int(pre_image_clip.width * scaling_factor)))
                return resampled_band
            if self.extraction_method == "dNDVI":
                resampled_band1 = ResampledBand(self.pre_image_red)
                resampled_band2 = ResampledBand(self.pre_image_NIR)
            else:
                resampled_band1 = ResampledBand(self.pre_image_NIR)
                resampled_band2 = ResampledBand(self.pre_image_SWIR)
            transform = pre_image_clip.transform * pre_image_clip.transform.scale((pre_image_clip.width/resampled_band1.shape[-1]), (pre_image_clip.height / resampled_band1.shape[-2]))
            profile = pre_image_clip.profile
            profile.update({"count": 2, "height": resampled_band1.shape[0], "width": resampled_band1.shape[1], "transform": transform})
            resampled_pre_image_clip_path = os.path.join(self.temp_folder, "ResampledPreImageClip.tif")
            with rasterio.open(resampled_pre_image_clip_path, 'w', **profile) as dst:
                dst.write(resampled_band1, 1)
                dst.write(resampled_band2, 2)
            self.UploadData(resampled_pre_image_clip_path)
            arcpy.AddMessage("Resampling done!")
            return resampled_pre_image_clip_path, resampling_condition, transform
        else:
            arcpy.AddMessage("Pre and Post images have the same resolution, no resampling needed...")
            return self.pre_image_clip, resampling_condition, self.pre_image_clip_transform
    
    def CalculateIndex(self):
        arcpy.AddMessage(f"calculating {self.extraction_method} index...")
        if self.extraction_method == "dNDVI":
            if self.resampling_condition:
                NDVI_pre = arcpy.sa.NDVI(self.pre_image_clip_resampled, 2, 1)
                NDVI_post = arcpy.sa.NDVI(self.post_image_clip, self.post_image_NIR, self.post_image_red)
            else:
                NDVI_pre = arcpy.sa.NDVI(self.pre_image_clip_resampled, self.pre_image_NIR, self.pre_image_red)
                NDVI_post = arcpy.sa.NDVI(self.post_image_clip, self.post_image_NIR, self.post_image_red)
            index = NDVI_pre - NDVI_post
        elif self.extraction_method == "dNBR":
            if self.resampling_condition:
                NBR_pre = arcpy.sa.NBR(self.pre_image_clip_resampled, 2, 1)
                NBR_post = arcpy.sa.NBR(self.post_image_clip, self.post_image_SWIR, self.post_image_NIR)
            else:    
                NBR_pre = arcpy.sa.NBR(self.pre_image_clip_resampled, self.pre_image_SWIR, self.pre_image_NIR)
                NBR_post = arcpy.sa.NBR(self.post_image_clip, self.post_image_SWIR, self.post_image_NIR)
            index = NBR_pre - NBR_post
        index_path = os.path.join(self.temp_folder, self.extraction_method + ".tif")
        index.save(index_path)
        self.UploadData(index_path)
        return index 
    
    def ExtractClassificationVector(self, threshold = 0.2):
        start_time = time.time()
        binary_classification = self.index > threshold 
        binary_classification = arcpy.RasterToNumPyArray(binary_classification) # turn arcpy object into numpy object
        classification_mask = binary_classification == 1
        shapes_generator = shapes(binary_classification, mask = classification_mask, transform = self.post_image_clip_transform)
        geometries = [shape(geom) for geom, value in shapes_generator if value == 1]
        classification_vector = gpd.GeoDataFrame({'geometry': geometries}, crs = self.crs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        arcpy.AddMessage(f"Extracting classification vector took {elapsed_time} seconds")
        return classification_vector, threshold
    
    def CleanClassificationVector(self):
        # Erase burnt area from previous MONITs
        burnt_area_gdf = gpd.read_file(self.geodatabase, layer = "B1_observed_event_a")
        burnt_area_gdf = burnt_area_gdf.to_crs(self.crs)
        if burnt_area_gdf.shape[0] > 0:
            arcpy.AddMessage("erasing previously burnt areas...")
            self.classification_vector = self.classification_vector.overlay(burnt_area_gdf, how = "difference")
        # Erase Not-Analyzed
        footprint_gdf = gpd.read_file(self.geodatabase, layer = "A2_image_footprint_a")
        footprint_gdf = footprint_gdf.to_crs(self.crs) 
        clouds_gdf = footprint_gdf[(footprint_gdf["obj_type"] == 2) & (footprint_gdf["or_src_id"] == self.DamageSourceIdentifier)]
        if clouds_gdf.shape[0] > 0:
            arcpy.AddMessage("erasing clouds...")
            self.classification_vector = self.classification_vector.overlay(clouds_gdf, how = "difference")
        # Erase Hydrography_a 
        hydrography_gdf = gpd.read_file(self.geodatabase, layer = "D1_hydrography_a")
        hydrography_gdf = hydrography_gdf.to_crs(self.crs)
        if hydrography_gdf.shape[0] > 0:
            arcpy.AddMessage("erasing hydrography...")
            self.classification_vector = self.classification_vector.overlay(hydrography_gdf, how = "difference")
        # Remove small features 
        arcpy.AddMessage("removing feature smaller than MMU...")
        self.classification_vector['area'] = self.classification_vector.area
        classification_vector_filtered = self.classification_vector[self.classification_vector['area'] > self.MMU]
        # Fill small holes 
        arcpy.AddMessage("removing small holes from classification polygon...")
        shapely_objects = classification_vector_filtered['geometry']
        classification_vector_cleaned = self.FillHoles(shapely_objects, self.MMU)
        # Calculate attributes that will be placed in the feature class
        arcpy.AddMessage("calculating attributes of the classification vector...")
        classification_vector_cleaned["event_type"] = 8 
        classification_vector_cleaned["obj_desc"] = self.ObjectDescriptionDict[self.ObjectDescription]
        classification_vector_cleaned["det_method"] = 2
        classification_vector_cleaned["notation"] = 802
        classification_vector_cleaned["dmg_src_id"] = self.DamageSourceIdentifier
        classification_vector_cleaned['area'] = round((classification_vector_cleaned.area / 10000), 2)
        # Save classification vector
        extraction_name = f"extraction_thr_{self.threshold}.shp"
        extraction = os.path.join(self.temp_folder, extraction_name)
        classification_vector_cleaned.to_file(extraction, driver='ESRI Shapefile')
        return extraction
    
    def UploadData(self, data_path): 
        aprx = arcpy.mp.ArcGISProject('current')
        map = aprx.activeMap
        map.addDataFromPath(data_path)
    
    def RemoveTemporaryLayers(self):
        aprx = arcpy.mp.ArcGISProject('current')
        map = aprx.activeMap
        for layer in map.listLayers():
            if layer.name.startswith("extraction") or layer.name in [self.clipping_feature_arcgis, self.extraction_method + '.tif']:
                map.removeLayer(layer)

    def make_GUI(self):
        def press_yes_button():
            # This function is executed when "Yes" button is clicked
            new_label2.pack(pady=10)
            entry2.pack(pady=10)
            submit_button2.pack(pady=10)

        def press_no_button():
            # This function is executed when "No" button is clicked
            new_label.pack(pady=10)
            entry.pack(pady=10)
            submit_button.pack(pady=10)

        def recalculate():
            # This function is executed when the new value is submitted
            threshold_value = float(entry.get())
            root.destroy()
            arcpy.AddMessage("reclassifying...")  
            self.classification_vector, self.threshold = self.ExtractClassificationVector(threshold = threshold_value)
            self.classification_vector_cleaned = self.CleanClassificationVector()
            self.UploadData(self.classification_vector_cleaned)
            self.make_GUI()

        def accept_result():
            acceptance_value = entry2.get()
            root.destroy()
            extraction_name = f"extraction_thr_{acceptance_value}"
            arcpy.AddMessage("Appending data")
            append_data(extraction_name)
            arcpy.AddMessage("Append done")
        
        def append_data(extraction):
            arcpy.management.Append(
            inputs = extraction,
            target = self.obs_event,
            schema_type = "NO_TEST",
            field_mapping = r'event_type "Event Type" true true false 4 Long 0 0,First,#,extraction,event_type,-1,-1;obj_desc "Object Description" true true false 2 Short 0 0,First,#,extraction,obj_desc,-1,-1;det_method "Determination Method" true true false 2 Short 0 0,First,#,extraction,det_method,-1,-1;notation "Comment" true true false 255 Text 0 0,First,#,extraction,notation,0,254;dmg_src_id "Damage Source Identifier" true true false 4 Long 0 0,First,#,extraction,dmg_src_id,-1,-1;area "Affected Area (ha)" true true false 8 Double 0 0,First,#,extraction,area,-1,-1',
            subtype = "",
            expression = "",
            match_fields = None,
            update_geometry = "NOT_UPDATE_GEOMETRY"
            )

        root = tk.Tk()
        root.title("Sensitivity Analysis Dialog Box")
        window_width = 400
        window_height = 300
        screen_width = root.winfo_screenwidth()
        screen_height = root.winfo_screenheight()
        center_x = int(screen_width/2 - window_width/2)
        center_y = int(screen_height/2 - window_height/2)
        root.geometry(f'{window_width}x{window_height}+{center_x}+{center_y}')
        root.configure(bg="lightblue")
        label = tk.Label(root, text="Are you satisfied with the result?", font=("Helvetica", 14), bg="lightblue")
        label.pack(pady=20)
        button_frame = tk.Frame(root, bg="lightblue")
        button_frame.pack(pady=10)
        yes_button = tk.Button(button_frame, text="Yes", command=press_yes_button, font=("Helvetica", 12), bg="green", fg="white", padx=20, pady=10)
        yes_button.pack(side=tk.LEFT, padx=20)
        no_button = tk.Button(button_frame, text="No", command=press_no_button, font=("Helvetica", 12), bg="red", fg="white", padx=20, pady=10)
        no_button.pack(side=tk.RIGHT, padx=20)
        new_label = tk.Label(root, text="Please insert a new value:", font=("Helvetica", 12), bg="lightblue")
        entry = tk.Entry(root, font=("Helvetica", 12))
        submit_button = tk.Button(root, text="Submit", command=recalculate, font=("Helvetica", 12), bg="blue", fg="white", padx=20, pady=10)
        new_label2 = tk.Label(root, text="Choose the threshold you want to use:", font=("Helvetica", 12), bg="lightblue")
        entry2 = tk.Entry(root, font=("Helvetica", 12))
        submit_button2 = tk.Button(root, text="Submit", command=accept_result, font=("Helvetica", 12), bg="blue", fg="white", padx=20, pady=10)
        root.mainloop()
    
    def run(self):
        self.UploadData(self.classification_vector_cleaned)
        self.make_GUI()
        self.RemoveTemporaryLayers()
    
if __name__ == "__main__":
    Classificator().run()