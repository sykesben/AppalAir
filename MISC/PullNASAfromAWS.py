import geopandas as gpd
from shapely.geometry import Point

# Center point
lon = -81.54
lat = 35.91

# Create point in WGS84
gdf = gpd.GeoDataFrame(
    geometry=[Point(lon, lat)],
    crs="EPSG:4326"
)

# Project to a CRS with meter units
gdf_proj = gdf.to_crs("EPSG:5070")  # CONUS Albers

# Create 350 km buffer
buffer = gdf_proj.buffer(350000)

# Back to lat/lon
buffer_gdf = gpd.GeoDataFrame(
    geometry=buffer,
    crs="EPSG:5070"
).to_crs("EPSG:4326")

# Save geojson
buffer_gdf.to_file(
    r"C:\Users\bensy\Documents\Research\buffer_350km.geojson",
    driver="GeoJSON")