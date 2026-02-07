import h5py
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import imageio.v2 as imageio
import os

# --- 1. FUNCIONES AUXILIARES ---
def get_var_data_subset(f, var_name, t_idx, lat_idx, lon_idx):
    """Extrae datos recortados compatible con las limitaciones de h5py"""
    var = f[var_name]
    # h5py no permite fancy indexing en dos ejes a la vez. 
    # Leemos el bloque completo que abarca nuestras latitudes y longitudes.
    l1, l2 = lat_idx[0], lat_idx[-1] + 1
    lo1, lo2 = lon_idx[0], lon_idx[-1] + 1
    
    # Extraemos el bloque (Slicing simple)
    data = var[t_idx, 0, l1:l2, lo1:lo2].astype(float)
    
    if 'scale_factor' in var.attrs:
        data *= var.attrs['scale_factor']
    if 'add_offset' in var.attrs:
        data += var.attrs['add_offset']
    return data

def agregar_coordenadas(ax):
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}

# --- 2. CONFIGURACIÓN ---
archivo = r"C:\Users\emmnu\OneDrive\Escritorio\Sinop_vera\23_todo.nc"
extent = [260, 330, -60, 15] # Sudamérica
lat_min, lat_max, lon_min, lon_max = -60, 15, 260, 330

with h5py.File(archivo, 'r') as f:
    lats_full = f['latitude'][:]
    lons_full = f['longitude'][:]
    tiempos = f['z'].shape[0]

    # Encontrar los rangos de índices
    lat_indices = np.where((lats_full >= lat_min) & (lats_full <= lat_max))[0]
    lon_indices = np.where((lons_full >= lon_min) & (lons_full <= lon_max))[0]
    
    # Coordenadas recortadas para el gráfico
    lats = lats_full[lat_indices]
    lons = lons_full[lon_indices]
    
    frames_div = []
    frames_omega = []

    print(f"Generando {tiempos} cuadros con etiquetas de coordenadas...")

    for t in range(tiempos):
        # Carga de datos optimizada para h5py
        z = get_var_data_subset(f, 'z', t, lat_indices, lon_indices) / 9.80665
        u = get_var_data_subset(f, 'u', t, lat_indices, lon_indices)
        v = get_var_data_subset(f, 'v', t, lat_indices, lon_indices)
        div = get_var_data_subset(f, 'd', t, lat_indices, lon_indices) * 1e5
        omega = get_var_data_subset(f, 'w', t, lat_indices, lon_indices)

        # --- GIF DIVERGENCIA ---
        fig1 = plt.figure(figsize=(10, 8))
        ax1 = plt.axes(projection=ccrs.PlateCarree())
        ax1.set_extent(extent, crs=ccrs.PlateCarree())
        agregar_coordenadas(ax1)
        ax1.add_feature(cfeature.COASTLINE, linewidth=1)
        ax1.add_feature(cfeature.BORDERS, linestyle=':', alpha=0.6)
        
        ax1.contourf(lons, lats, div, levels=np.linspace(-8, 8, 21), cmap='RdBu_r', extend='both', transform=ccrs.PlateCarree())
        cs1 = ax1.contour(lons, lats, z, levels=15, colors='black', linewidths=0.8, transform=ccrs.PlateCarree())
        ax1.clabel(cs1, inline=True, fontsize=8, fmt='%1.0f')
        
        # Ajustamos el salto de flechas para que se vean bien
        skip = 10
        ax1.quiver(lons[::skip], lats[::skip], u[::skip, ::skip], v[::skip, ::skip], transform=ccrs.PlateCarree(), alpha=0.5, scale=700)
        
        plt.title(f'Divergencia 200 hPa - {t:02d}:00 UTC')
        tmp_div = f'tmp_div_{t:02d}.png'
        plt.savefig(tmp_div, dpi=100, bbox_inches='tight')
        frames_div.append(tmp_div)
        plt.close()

        # --- GIF OMEGA ---
        fig2 = plt.figure(figsize=(10, 8))
        ax2 = plt.axes(projection=ccrs.PlateCarree())
        ax2.set_extent(extent, crs=ccrs.PlateCarree())
        agregar_coordenadas(ax2)
        ax2.add_feature(cfeature.COASTLINE, linewidth=1)
        ax2.add_feature(cfeature.BORDERS, linestyle=':', alpha=0.6)
        
        ax2.contourf(lons, lats, omega, levels=np.linspace(-1, 1, 21), cmap='RdBu', extend='both', transform=ccrs.PlateCarree())
        
        plt.title(f'Omega (Mov. Vertical) 200 hPa - {t:02d}:00 UTC')
        tmp_om = f'tmp_om_{t:02d}.png'
        plt.savefig(tmp_om, dpi=100, bbox_inches='tight')
        frames_omega.append(tmp_om)
        plt.close()

        if t % 6 == 0: print(f"Procesado cuadro {t}/{tiempos}")

# --- 3. COMPILACIÓN DE LOS GIFS ---
print("Creando GIFs finales...")
with imageio.get_writer('anim_divergencia_final.gif', mode='I', duration=250, loop=0) as writer:
    for f_name in frames_div:
        writer.append_data(imageio.imread(f_name))
        os.remove(f_name)

with imageio.get_writer('anim_omega_final.gif', mode='I', duration=250, loop=0) as writer:
    for f_name in frames_omega:
        writer.append_data(imageio.imread(f_name))
        os.remove(f_name)

print("¡Proceso terminado exitosamente!")


import os
print(os.getcwd())