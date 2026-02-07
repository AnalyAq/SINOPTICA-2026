import h5py
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# --- 1. FUNCIÓN PARA LEER DATOS ---
def get_var_data(f, var_name):
    var = f[var_name]
    data = var[:].astype(float)
    if 'scale_factor' in var.attrs:
        data *= var.attrs['scale_factor']
    if 'add_offset' in var.attrs:
        data += var.attrs['add_offset']
    return data

# --- 2. CARGA DE ARCHIVO ---
nombre_archivo = r'C:\Users\emmnu\OneDrive\Escritorio\Sinop_vera\guia.nc'

with h5py.File(nombre_archivo, 'r') as f:
    lats = f['latitude'][:]
    lons = f['longitude'][:]
    z = get_var_data(f, 'z')[0, 0] / 9.80665
    u = get_var_data(f, 'u')[0, 0]
    v = get_var_data(f, 'v')[0, 0]
    div = get_var_data(f, 'd')[0, 0] * 1e5
    omega = get_var_data(f, 'w')[0, 0]

extent = [260, 330, -60, 15] 

# --- FUNCIÓN PARA AGREGAR COORDENADAS AL MAPA ---
def agregar_coordenadas(ax):
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False   # Desactivar etiquetas arriba
    gl.right_labels = False # Desactivar etiquetas a la derecha
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}

# --- MAPA 1: VIENTO, GEOPOTENCIAL Y DIVERGENCIA ---
fig1 = plt.figure(figsize=(12, 8))
ax1 = plt.axes(projection=ccrs.PlateCarree())
ax1.set_extent(extent, crs=ccrs.PlateCarree())

# Añadir etiquetas de Lat/Lon
agregar_coordenadas(ax1)

ax1.add_feature(cfeature.COASTLINE, linewidth=1, edgecolor='black')
ax1.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='gray')

clevs_div = np.linspace(-6, 6, 25)
cf1 = ax1.contourf(lons, lats, div, levels=clevs_div, cmap='RdBu_r', extend='both', transform=ccrs.PlateCarree())
plt.colorbar(cf1, ax=ax1, label='Divergencia ($10^{-5} s^{-1}$)', shrink=0.7)

cs1 = ax1.contour(lons, lats, z, levels=np.arange(10000, 13000, 60), colors='black', linewidths=0.8, transform=ccrs.PlateCarree())
ax1.clabel(cs1, inline=True, fontsize=9, fmt='%1.0f')

skip = 12
ax1.quiver(lons[::skip], lats[::skip], u[::skip, ::skip], v[::skip, ::skip], 
           color='black', alpha=0.5, scale=800, transform=ccrs.PlateCarree())

ax1.set_title('Diagnóstico 200 hPa: Viento, Geopotencial y Divergencia\nSudamérica - 23/01/2026', fontsize=13)
plt.show()

# --- MAPA 2: OMEGA ---
fig2 = plt.figure(figsize=(12, 8))
ax2 = plt.axes(projection=ccrs.PlateCarree())
ax2.set_extent(extent, crs=ccrs.PlateCarree())

# Añadir etiquetas de Lat/Lon
agregar_coordenadas(ax2)

ax2.add_feature(cfeature.COASTLINE, linewidth=1, edgecolor='black')
ax2.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='gray')

clevs_w = np.linspace(-1.2, 1.2, 25)
cf2 = ax2.contourf(lons, lats, omega, levels=clevs_w, cmap='RdBu', extend='both', transform=ccrs.PlateCarree())
plt.colorbar(cf2, ax=ax2, label='Omega (Pa/s)', shrink=0.7)

ax2.set_title('Diagnóstico 200 hPa: Velocidad Vertical (Omega)\nAzul = Ascenso | Rojo = Subsidencia', fontsize=13)
plt.show()