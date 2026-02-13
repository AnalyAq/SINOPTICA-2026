import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import imageio.v2 as imageio
import os

# --- 1.CARGA DE DATOS ---
# Asegúrate de que estos nombres coincidan con tus archivos en la carpeta
file_accum = r"C:\Users\emmnu\OneDrive\Escritorio\Sinop_vera\final\data_stream-oper_stepType-accum.nc"

# Si tienes el viento y presión en otro archivo, cárgalo aquí:
ds_main = xr.open_dataset(r"C:\Users\emmnu\OneDrive\Escritorio\Sinop_vera\final\data_stream-oper_stepType-instant.nc") 
ds = xr.open_dataset(file_accum)

# --- 2. CONFIGURACIÓN DE VARIABLES Y REGIÓN ---
# Sudamérica y Pacífico Sur para ver los sistemas sinópticos
extent = [-100, -50, -55, -10] 

# Diccionario de configuración para cada GIF
# Formato: 'Variable': ['Nombre_en_NC', 'Unidades', 'Mapa_Color', 'Niveles']
config = {
    'Nevadas': ['sf', 'mm', 'Blues', np.linspace(0.1, 10, 21)],
    'Evaporacion': ['e', 'mm', 'PuOr', np.linspace(-1, 1, 21)],
    'Presion': ['sp', 'hPa', 'RdGy_r', np.arange(980, 1040, 4)],
    'Viento': ['ws', 'm/s', 'YlGnBu', np.linspace(0, 40, 21)],
    'Flujo_Humedad': ['q_flux', 'g/kg*m/s', 'Spectral_r', np.linspace(0, 150, 21)]
}

# --- 3. PROCESAMIENTO Y GENERACIÓN DE CUADROS ---
for nombre_gif, params in config.items():
    var_nc, unidad, cmap, levels = params
    print(f"Generando cuadros para: {nombre_gif}")
    
    frames = []
    
    # Verificamos si la variable existe en el dataset
    if var_nc not in ds and nombre_gif not in ['Viento', 'Flujo_Humedad', 'Presion']:
        print(f"Saltando {nombre_gif}: No se encontró {var_nc} en el archivo.")
        continue

    for t in range(len(ds.valid_time)):
        fig = plt.figure(figsize=(12, 8))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        
        # Capas geográficas
        ax.add_feature(cfeature.COASTLINE, linewidth=1.5)
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        
        # Etiquetas de Lat/Lon
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False; gl.right_labels = False
        gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER

        # --- LÓGICA POR VARIABLE ---
        data = None
        if var_nc in ds:
            data = ds[var_nc].isel(valid_time=t)
            if unidad == 'mm': data = data * 1000 # Conversión m a mm
            if unidad == 'hPa': data = data / 100 # Conversión Pa a hPa
        
        # Caso especial: Viento y Flujo (Si tienes u, v, q)
        if nombre_gif == 'Viento' and 'u' in ds:
            u = ds['u'].isel(valid_time=t); v = ds['v'].isel(valid_time=t)
            data = np.sqrt(u**2 + v**2)
            ax.quiver(ds.longitude[::10], ds.latitude[::10], u[::10,::10], v[::10,::10], transform=ccrs.PlateCarree())
        
        if data is not None:
            cf = ax.contourf(ds.longitude, ds.latitude, data, levels=levels, cmap=cmap, extend='both', transform=ccrs.PlateCarree())
            plt.colorbar(cf, label=f'{nombre_gif} ({unidad})', shrink=0.7)

        # Título con fecha y hora
        time_str = str(ds.valid_time.values[t])[:16]
        plt.title(f'Diagnóstico Sinóptico: {nombre_gif} a nivel de superficie \nFecha: {time_str} UTC', fontweight='bold')
        
        # Guardar frame temporal
        frame_name = f'temp_{nombre_gif}_{t:02d}.png'
        plt.savefig(frame_name, dpi=100, bbox_inches='tight')
        frames.append(frame_name)
        plt.close()

    # --- 4. CREAR EL GIF ---
    if frames:
        with imageio.get_writer(f'Atacama_{nombre_gif}.gif', mode='I', duration=300, loop=0) as writer:
            for f in frames:
                writer.append_data(imageio.imread(f))
                os.remove(f) # Limpieza
        print(f"¡GIF de {nombre_gif} creado con éxito!")

print("Análisis terminado.")
