import h5py
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import datetime

# --- 1. FUNCIÓN DE EXTRACCIÓN ---
def get_data(f, var, t_idx):
    data = f[var][t_idx].astype(float)
    if 'scale_factor' in f[var].attrs:
        data *= f[var].attrs['scale_factor']
    if 'add_offset' in f[var].attrs:
        data += f[var].attrs['add_offset']
    return data

# --- 2. CONFIGURACIÓN ---
archivo_accum = r"C:\Users\emmnu\OneDrive\Escritorio\Sinop_vera\final\data_stream-oper_stepType-accum.nc"
# Nota: Viento y Presión suelen estar en archivos de niveles de presión (instantáneos)
extent = [-120, 10, -60, 20] # Sudamérica Sur-Centro

# Índices para 12:00 UTC (cada 4 pasos si es 6-horario)
indices_12utc = [2, 6, 10, 14, 18, 22, 26] 
fechas = ["2025-06-24", "2025-06-25", "2025-06-26", "2025-06-27", "2025-06-28", "2025-06-29", "2025-06-30"]

with h5py.File(archivo_accum, 'r') as f:
    lats = f['latitude'][:]
    lons = f['longitude'][:]
    
    for idx, fecha in zip(indices_12utc, fechas):
        # Extraer variables del archivo 'accum'
        snow = get_data(f, 'sf', idx) * 1000 # a mm
        evap = get_data(f, 'e', idx) * 1000  # a mm
        
        # --- CREAR MAPA DE NEVADA ---
        fig = plt.figure(figsize=(10, 7))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(extent)
        ax.add_feature(cfeature.COASTLINE, linewidth=1)
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        
        # Etiquetas de Lat/Lon
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False; gl.right_labels = False
        gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER
        
        # Plot Nevada
        cf = ax.contourf(lons, lats, snow, levels=np.linspace(0.1, 20, 15), cmap='Blues', extend='max')
        plt.colorbar(cf, label='Snowfall (mm)')
        
        # Título Dinámico
        plt.title(f'Nevada en Atacama - {fecha} 12:00 UTC\nAnálisis de la Baja Segregada (DANA)', fontweight='bold')
        plt.savefig(f'mapa_nieve_{fecha}.png', dpi=150)
        plt.show()
        
        # --- MAPA DE EVAPORACIÓN ---
        plt.figure(figsize=(10, 7))
        ax2 = plt.axes(projection=ccrs.PlateCarree())
        ax2.set_extent(extent)
        ax2.add_feature(cfeature.COASTLINE)
        cf2 = ax2.contourf(lons, lats, evap, levels=20, cmap='RdBu_r')
        plt.colorbar(cf2, label='Evaporación (mm)')
        plt.title(f'Evaporación Superficial - {fecha} 12:00 UTC')
        plt.show()

print("Proceso completado. Revisa tus archivos PNG.")