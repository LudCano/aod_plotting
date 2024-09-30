# ===================================================
# ---------------------------------------------------
#          DOWNLOADING GOES16 AOD DATA (LIVE)
# by: Ludving Cano, based on the course on aod and 
#     fire monitoring. 2024
# Retrieves the AOD from important places for the LFA
# this should run from 7 (11) to 18 (22) local (UTC)
# ---------------------------------------------------
# ===================================================

# IMPORTANT!!!!
# This code will catch up the series, use it carefully since it can take a long time

######################################################
###########   PARÁMETROS DE DESCARGA   ###############
######################################################
outdir     = 'descarga_live_aod'   #carpeta donde descargar los datos, se creará si no existe
flush_orig = True       #Eliminar los archivos originales (sin cortar)
quality = 'top2'
date_start = '2024-09-24'
date_end   = '2024-09-27'
hrs_start = 7
hrs_end   = 18

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
# LIBRERÍAS PARA MAPAS Y GEOGRAFÍA
from cartopy import crs as ccrs
import cartopy.feature as cfeature
# LIBRERÍAS PARA GRÁFICOS
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.offsetbox import AnchoredText
import matplotlib.ticker as ticker

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
#########   FUNCION PARA OBTENER DATOS    ###########
######################################################
def get_aod_places(day):

    y = day.year
    m = day.month
    d = day.day
    h = day.hour + 4
    # Obteniendo el dia del año 
    julian_day = dt.datetime(y, m, d).strftime('%j')
    print("JULIAN DAY: ----",julian_day)
    data_path = bucket + '/' + product + '/'  + str(y) + '/' + julian_day + '/' + str(h).zfill(2)
    try:
        fils = sorted(fs.ls(data_path))
    except:
        print('No data, oh no')
        exit()
    if len(fils) == 0:
        exit()
    else:
        last_file = fils[-1]



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

    if not os.path.exists(outdir):
        os.makedirs(outdir)


    ######################################################
    #############   DESCARGANDO ARCHIVOS    ##############
    ######################################################
    # Loop through list of ABI files on NODD
    # Open each file remotely, subset variables & save as a new .nc file

    def getting_aodd(file):
        fname = file.split('/')[-1]
        subname = "sub_" + fname
        #print(file.split('/')[-1])  # Print the ABI file name
        fpath = (Path(abi_path) / file.split('/')[-1]).as_posix()
        trimmedpath = (Path(abi_path) / subname).as_posix()
        if not os.path.exists(trimmedpath):
            if not os.path.exists(fpath):
                fs.get(file, fpath)
                with xr.open_dataset(fpath, engine = 'h5netcdf') as ds:
                    subset_abi_file(ds, upper_left_latitude, upper_left_longitude,lower_right_latitude, lower_right_longitude, fname)
        if os.path.exists(fpath) and flush_orig:
            os.remove(fpath)
        else:
            print('file exists, yay!')

        print('Downloaded!')


        ds = xr.open_dataset(trimmedpath)
        latitude_array, longitude_array = calculate_abi_lat_lon(ds)
        print('GOES-16 constants defined!')


        with xr.open_dataset(trimmedpath, engine='netcdf4') as ds:
            # Convert "DQF" DataArray to correct data type
            DQF = ds.DQF.astype('uint8')

            # Set AOD data quality for plot
            if quality == 'high':
                plot_quality = (DQF == 0)
            elif quality == 'top2':
                plot_quality = (DQF <= 1)
            elif quality == 'all':
                plot_quality = (DQF <= 2)
            
            aod_arr = ds.AOD.where(plot_quality)


        ### GETTING THE INDEXES FOR THE PIXELS
        #aod data saved
        aod_data = open('aod_places.csv','+a')
        places = pd.read_csv('stations.csv')
        datestmp = file.split('/')[-1].split('_')[3][8:12]
        h_real = int(datestmp[:2])
        m_real = int(datestmp[2:])

        aods = []
        for _,s in places.iterrows():
            station = s.codename
            lat = s.lat
            lon = s.lon
            distances = np.sqrt((latitude_array - lat) ** 2 + (longitude_array - lon) ** 2)   
            idxs = np.argwhere(distances == np.min(distances))[0]
            aod_place = np.array(aod_arr)[idxs[0],idxs[1]]
            aods.append(aod_place)

        a = ','.join([str(i) for i in aods])
        a = dt.datetime.strftime(dt.datetime(y,m,d,h_real,m_real), '%Y-%m-%d %H:%M,') + a
        print(a, file = aod_data)

        os.remove(trimmedpath)
        
    getting_aodd(fils[-1])
    getting_aodd(fils[2])

######################################################
#########   OBTENIENDO DIAS A DESCARGAR    ###########
######################################################
d0 = dt.datetime.strptime(date_start, '%Y-%m-%d')
n_days = int(date_end.split('-')[-1]) - int(date_start.split('-')[-1]) + 1
for i in range(n_days):
    day0 = d0 + dt.timedelta(days = i)
    for h in range(hrs_start, hrs_end+1):
        day = day0 + dt.timedelta(hours = h)
        print(day)
        get_aod_places(day)

