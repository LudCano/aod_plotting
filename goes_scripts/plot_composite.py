import pandas as pd
import os
import xarray as xr          # manejo de xarrays
import numpy as np
import datetime as dt
# LIBRERÍAS PARA MAPAS Y GEOGRAFÍA
from cartopy import crs as ccrs
import cartopy.feature as cfeature
# LIBRERÍAS PARA GRÁFICOS
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.offsetbox import AnchoredText
import matplotlib.ticker as ticker

dir_data = 'descarga_composite'
start = '1230'
end = '1430'
quality = 'top2'

lst_files = []
for n,i in enumerate(sorted(os.listdir(dir_data))):
    hr = int(i[34:38])
    if hr>=int(start) and hr<=int(end):
        #print(str(int(i[34:38])).zfill(4))
        lst_files.append(os.path.join(dir_data,i))
        #print(n)



# Calculate 2D latitude & longitude arrays from GOES ABI fixed grid projection
# If not using xarray, specify "dtype=np.int64" so squares of large integers don't overflow memory
def calculate_abi_lat_lon(ds):

    # ABI fixed grid projection constants
    abi_x = ds.x  # E/W scanning angle in radians
    abi_y = ds.y  # N/S elevation angle in radians
    lambda_0 = np.deg2rad(ds.goes_imager_projection.longitude_of_projection_origin)
    H = ds.goes_imager_projection.perspective_point_height+ds.goes_imager_projection.semi_major_axis
    r_eq = ds.goes_imager_projection.semi_major_axis
    r_pol = ds.goes_imager_projection.semi_minor_axis

    # Create 2D x,y arrays from 1D x,y data arrays
    abi_x_2d, abi_y_2d = np.meshgrid(abi_x, abi_y)

    # Geometry equations to calculate latitude and longitude
    a = np.square(np.sin(abi_x_2d))+np.square(np.cos(abi_x_2d))*(np.square(np.cos(abi_y_2d))+(np.square(r_eq)*np.square(np.sin(abi_y_2d))*np.reciprocal(np.square(r_pol))))
    b = np.negative(2*H)*np.cos(abi_x_2d)*np.cos(abi_y_2d)
    c = np.square(H)-np.square(r_eq)
    r_s = np.reciprocal(2*a)*(np.negative(b)-np.sqrt(np.square(b)-4*a*c))
    s_x = r_s*np.cos(abi_x_2d)*np.cos(abi_y_2d)
    s_y = np.negative(r_s)*np.sin(abi_x_2d)
    s_z = r_s*np.cos(abi_x_2d)*np.sin(abi_y_2d)

    lat = np.rad2deg(np.arctan(np.square(r_eq)*np.reciprocal(np.square(r_pol))*s_z*np.reciprocal(np.sqrt(np.square(H-s_x)+np.square(s_y)))))
    lon = np.rad2deg(lambda_0-np.arctan(s_y*np.reciprocal(H-s_x)))

    return lat, lon



ds = xr.open_dataset(lst_files[0])
latitude_array, longitude_array = calculate_abi_lat_lon(ds)

# Define constants for GOES-16 geostationary orbit & Earth reference system
satellite_height = ds.goes_imager_projection.perspective_point_height
semi_major_axis = ds.goes_imager_projection.semi_major_axis
semi_minor_axis = ds.goes_imager_projection.semi_minor_axis
central_longitude = ds.goes_imager_projection.longitude_of_projection_origin

# Define geostationary map projection using cartopy
globe = ccrs.Globe(semimajor_axis=semi_major_axis, semiminor_axis=semi_minor_axis)
geo_projection = ccrs.Geostationary(central_longitude=central_longitude,
                                    satellite_height=satellite_height,
                                    globe=globe, sweep_axis='x')

print('GOES-16 constants defined!')


lst_mtxs = []
for file in lst_files:
    with xr.open_dataset(file, engine='netcdf4') as ds:
        # Convert "DQF" DataArray to correct data type
        DQF = ds.DQF.astype('uint8')

        # Set AOD data quality for plot
        if quality == 'high':
            plot_quality = (DQF == 0)
        elif quality == 'top2':
            plot_quality = (DQF <= 1)
        elif quality == 'all':
            plot_quality = (DQF <= 2)
        lst_mtxs.append(ds.AOD.where(plot_quality))

aods = np.array(lst_mtxs)
#print(aods.shape)
aods = np.nanmean(aods, axis = 0)
print(aods.shape)

# Set up figure in Matplotlib
fig = plt.figure(figsize=(8, 6))

# Set map projection using cartopy
# Use geostationary projection set previously
ax = plt.axes(projection=geo_projection)

# Set geographic domain of map: [W_lon, E_lon, S_lat, N_lat]
# °E longitude > 0 > °W longitude, °N latitude > 0 > °S latitude
ax.set_extent([-82, -35, -29, 2], crs=ccrs.PlateCarree())   #small one

# Format lat/lon gridlines using cartopy
lon_ticks = np.arange(-90, -20, 5)
lat_ticks = np.arange(-35,0,5)
gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='silver')
gl.xlocator = ticker.FixedLocator(lon_ticks)
gl.ylocator = ticker.FixedLocator(lat_ticks)
gl.right_labels = None
gl.top_labels = None
gl.xlabel_style = {'size': 8}
gl.ylabel_style = {'size': 8}

# Add coastlines & borders, shade land & water polygons
# "zorder" argument sets order for plotting layers (larger zorder plots over smaller zorder)
ax.add_feature(cfeature.COASTLINE, linewidth=0.75, zorder=2)
ax.add_feature(cfeature.BORDERS, linewidth=0.75, zorder=2)
ax.add_feature(cfeature.LAKES, facecolor='lightgrey')
ax.add_feature(cfeature.LAND, facecolor='grey')
ax.add_feature(cfeature.OCEAN, facecolor='lightgrey')

# Set colormap and unique color for AOD > 1
cmap = plt.get_cmap('rainbow').with_extremes(over='darkred')

# Plot AOD data
# "transform=ccrs.PlateCarree()" argument is required b/c data are in geographic coordinates
plot = ax.pcolormesh(longitude_array, latitude_array, aods, cmap=cmap,
                    vmin=0, vmax=4, transform=ccrs.PlateCarree())

# Add colorbar
cb = fig.colorbar(plot, orientation='horizontal', fraction=0.2, pad=0.05, shrink=0.5,
                ticks=[0, 1, 2, 3, 4], extend='max')
cb.set_label(label='Aerosol Optical Depth at 550nm', size=8, weight='bold')
cb.ax.set_xticklabels(['0', '1','2','3','4'])
# Put extracted/reformated strings together to make title
def gen_hour(x):
    stri = x[0:2]  + ':' + x[2:]
    return stri
todaytitle0 = dt.datetime.now()
todaytitle = todaytitle0.date().strftime('%d %b %Y')
forfname = todaytitle0.date().strftime('%d%b')
image_title = f'GOES-16/ABI Aerosol Optical Depth  {todaytitle}  ' + gen_hour(start) + '-' + gen_hour(end) +' UTC Composite'
# Add plot title
plt.title(image_title, pad=8, size=8, weight='bold')


from matplotlib.patheffects import withStroke


cities = ['LP','CBBA','SCZ']
citiespos = [(-68.15,-16.5),(-66.15,-17.389),(-63.18,-17.786)]

# Definir el efecto de contorno
stroke_effect = withStroke(linewidth=1, foreground='black')
# Añadir una anotación con estilo

for c, p in zip(cities, citiespos):
    ax.annotate(
        c,
        xy=(p[0], p[1]-1),  # Coordenadas de la anotación
        xytext=(p[0], p[1]-1),  # Coordenadas del texto
        fontsize=7,
                color='white',  # Color de relleno
                path_effects=[stroke_effect],  # Efecto de borde
                ha='center', va='center', transform=ccrs.PlateCarree()
    )

    # Marcar un punto con estilo
    ax.scatter(p[0],p[1], color='white', edgecolor='black', s=20, linewidth=1, transform=ccrs.PlateCarree())


fig.savefig(f'composites/composite{forfname}.png', dpi = 500)
os.remove('current_plots/composite.png')
fig.savefig(f'current_plots/composite.png', dpi = 500)
for f in os.listdir(dir_data):
    pth = os.path.join(dir_data, f)
    os.remove(pth)
print('All AOD files were removed')
