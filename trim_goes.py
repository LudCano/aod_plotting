# ===================================================
# ---------------------------------------------------
#             DOWNLOADING GOES16 AOD DATA
# by: Ludving Cano, based on the course on aod and 
#     fire monitoring. 2024
# ---------------------------------------------------
# ===================================================

######################################################
###########   PARÁMETROS DE DESCARGA   ###############
######################################################
date0      = '2024-09-23'  #fecha inicial YYYY-MM-DD
datef      = '2024-09-23'  #fecha final YYYY-MM-DD
n_hour     = 6             #número de datos por hora, 1 (horario), 2(media hora), 3(20min), 6(todos)
outdir     = 'descarga_test'   #carpeta donde descargar los datos, se creará si no existe
trim_bf    = False         #Cortar en la nube? Puede tomar más tiempo (True/False)
flush_orig = True       #Eliminar los archivos originales (sin cortar)

######################################################
###########            DOMINIO         ###############
######################################################
# Subsetted domain settings
# Enter latitude/longitude values in degrees (integers or floats)
# °N latitude > 0 > °S latitude, °E longitude > 0 > °W longitude
# Set corners larger than desired map domain by ~5-10°
upper_left_latitude = 2  # Latitude of upper left corner
upper_left_longitude = -82  # Longitude of upper left corner
lower_right_latitude = -29  # Latitude of lower right corner
lower_right_longitude = -31  # Longitude of lower right corner

######################################################
##########            LIBRERIAS            ###########
######################################################
# Import modules and packages
import os                    # libreria para uso del sistema
import datetime as dt        # manejo de tiempos (timestamps)
from pathlib import Path     # manejo de rutas del sistema
import warnings              # alertas o errores
import s3fs                  # conexión al servidor
import xarray as xr          # manejo de xarrays
import numpy as np           # librería numérica
from PIL import Image        # manejo de imágenes
import pandas as pd          # tablas (y fechas)
from tqdm import tqdm        # barra de progreso
import h5netcdf

######################################################
#######   CONEXIÓN A AMAZON WEB SERVICES   ###########
######################################################
# Connect to AWS S3 anonymously
fs = s3fs.S3FileSystem(anon=True)
print('Connected to AWS...')
abi_path = Path(outdir)


bucket = 'noaa-goes16'  # string (INSTRUMENTO)
#product = 'ABI-L2-FDCF'  # string (PRODUCTO PARA FUEGOS - no disponible)
product = 'ABI-L2-AODF'  # string (PRODUCTO PARA AOD - F es FullDisk)

######################################################
#########   OBTENIENDO RUTAS A DESCARGAR   ###########
######################################################
# Obteniendo los índices de los datos a descargar (de los 6 disponibles)
if n_hour == 1:
   idxs = [0]
elif n_hour == 2:
   idxs = [0,3]
elif n_hour == 3:
   idxs = [0,2,4]
else:
   idxs = list(range(5))


if date0 == datef:
   lst_dates = [dt.datetime.strptime(date0, '%Y-%m-%d')]
else:
   d0 = dt.datetime.strptime(date0, '%Y-%m-%d')
   df = dt.datetime.strptime(datef, '%Y-%m-%d')
   lst_dates = list(pd.date_range(d0, df, freq='D'))
   

files_to_dwnld = []
for day in lst_dates:
    y = day.year
    m = day.month
    d = day.day
    # Obteniendo el dia del año 
    julian_day = dt.datetime(y, m, d).strftime('%j')
    for h in range(0,15):
        data_path = bucket + '/' + product + '/'  + str(y) + '/' + julian_day + '/' + str(h).zfill(2)
        fils = sorted(fs.ls(data_path))
        if len(fils) > 0:
            files_to_dwnld = files_to_dwnld + [fils[i] for i in idxs]



# Imprimiendo la cantidad total de datos a descargar
print(f'Total number of files to download: {len(files_to_dwnld)}')
# Imprimiendo mensaje de que los originales se van a borrar
if flush_orig:
    print('IMPORTANT: Original files will be removed after processed!')

######################################################
#######   FUNCIONES DE AMY - RECORTAR   ###########
######################################################
# Calculate ABI fixed grid N/S Elevation Angle (y) & E/W Scanning Angle (x) for given latitude & longitude
# xarray reads in all fixed grid constants as floasts, even H & r_eq (integers)
# If using another package, specify "dtype=np.int64" so squares of large integers don't overflow memory
# Also np.reciprocal can only be used with floats (not integers)
def calculate_abi_x_y(ds, latitude, longitude):

    # Convert entered latitude, longitude in degrees to radians
    phi = np.deg2rad(latitude)
    lamb = np.deg2rad(longitude)

    # ABI fixed grid projection constants
    lambda_0 = np.deg2rad(ds.goes_imager_projection.longitude_of_projection_origin)
    H = ds.goes_imager_projection.perspective_point_height+ds.goes_imager_projection.semi_major_axis
    r_eq = ds.goes_imager_projection.semi_major_axis
    r_pol = ds.goes_imager_projection.semi_minor_axis
    f = np.reciprocal(ds.goes_imager_projection.inverse_flattening)
    eccentricity = np.sqrt(f*(2-f))

    # Geometry equations to calculate y and x
    phi_c = np.arctan(np.square(r_pol)*np.reciprocal(np.square(r_eq))*np.tan(phi))
    r_c = r_pol*np.reciprocal(np.sqrt(1-(np.square(eccentricity)*np.square(np.cos(phi_c)))))
    s_x = H-r_c*np.cos(phi_c)*np.cos(lamb-lambda_0)
    s_y = np.negative(r_c)*np.cos(phi_c)*np.sin(lamb-lambda_0)
    s_z = r_c*np.sin(phi_c)

    abi_y = np.arctan(s_z*np.reciprocal(s_x))
    abi_x = np.arcsin(np.negative(s_y)*np.reciprocal(np.sqrt(np.square(s_x)+np.square(s_y)+np.square(s_z))))

    return abi_y, abi_x

# Subset ABI file xarray dataset to user-entered domain & save as new .nc file
def subset_abi_file(ds, upper_left_lat, upper_left_lon, lower_right_lat, lower_right_lon, fname):

    # Find y-index & x-index subsetted ranges in ABI Full Disk file
    # Get range of x and y corresponding to subsetted domain
    y_min, x_min = calculate_abi_x_y(ds, upper_left_lat, upper_left_lon)
    y_max, x_max = calculate_abi_x_y(ds, lower_right_lat, lower_right_lon)

    # Get indices of y_min, y_max & x_min, x_max
    # Finds 3-4 closest index values; select first one
    y_index_min = np.where(np.isclose(ds.y, y_min, atol=1e-04))[0][0]
    y_index_max = np.where(np.isclose(ds.y, y_max, atol=1e-04))[0][0]
    x_index_min = np.where(np.isclose(ds.x, x_min, atol=1e-04))[0][0]
    x_index_max = np.where(np.isclose(ds.x, x_max, atol=1e-04))[0][0]

    # Create new dataset: all variables w/ (y, x) dimensions subsetted to user-specified domain
    # Subsetted domain indices: [y_index_min:y_index_max, x_index_min:x_index_max]
    ds_sub = ds.isel(y=slice(y_index_min, y_index_max), x=slice(x_index_min, x_index_max))

    # Modify metadata for "geospatial_lat_lon_extent" variable to match subsetted domain
    ds_sub.geospatial_lat_lon_extent.attrs['geospatial_westbound_longitude'] = upper_left_longitude
    ds_sub.geospatial_lat_lon_extent.attrs['geospatial_eastbound_longitude'] = lower_right_longitude
    ds_sub.geospatial_lat_lon_extent.attrs['geospatial_northbound_latitude'] = upper_left_latitude
    ds_sub.geospatial_lat_lon_extent.attrs['geospatial_southbound_latitude'] = lower_right_latitude

    # Add info on subsetting indices to global file metadata
    ds_sub.attrs['File modified'] = 'Variables in this file with (y,x) dimensions were subsetted from ABI Full Disk using y[' + str(y_index_min) + ':' + str(y_index_max) + '], x[' + str(x_index_min) + ':' + str(x_index_max) + '] via code written by Dr. Amy Huff, IMSG at NOAA/NESDIS/STAR'

    # Save new dataset as an .nc file
    # Modify original file name to signify subsetted file ("-FSub")
    save_name = 'sub_' + fname
    # Silence warning from xarray about saving x, y with no fill value
    with warnings.catch_warnings():
      warnings.simplefilter('ignore')
      ds_sub.to_netcdf((abi_path / save_name).as_posix(), engine='h5netcdf')


if not os.path.exists(outdir):
    os.makedirs(outdir)


######################################################
#############   DESCARGANDO ARCHIVOS    ##############
######################################################
# Loop through list of ABI files on NODD
# Open each file remotely, subset variables & save as a new .nc file
orig_files_paths = []
lst_orig_files = os.listdir(outdir)
print(lst_orig_files)
for file in tqdm(lst_orig_files, desc='Trimming'):
    fpath = os.path.join(outdir,file)    
    print(fpath)
    with xr.open_dataset(fpath, engine = 'h5netcdf') as ds:
        subset_abi_file(ds, upper_left_latitude, upper_left_longitude,lower_right_latitude, lower_right_longitude, file)
    if os.path.exists(fpath) and flush_orig:
        os.remove(fpath)


print('Done!')
