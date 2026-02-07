# -*- coding: utf-8 -*-
"""
Created on Fri Feb  6 23:49:02 2026

@author: User
"""


#-------------------precipitacion-------------------------
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from shapely.geometry import Polygon  # Para crear polígonos geométricos

# 1. Cargar los datos ERA5 para un solo tiempo
ds_precip = xr.open_dataset(r"C:\Users\User\Desktop\mapitas\pp.nc")
ds_z200 = xr.open_dataset(r"C:\Users\User\Desktop\mapitas\vientolineas.nc")

# 2. Definir región de Sudamérica
lat_min, lat_max = -60, 20
lon_min_sa, lon_max_sa = -120, 10

# 3. Convertir longitudes de Sudamérica a rango 0-360 si es necesario
def convert_to_0_360(lon):
    """Convertir longitud de -180..180 a 0..360"""
    return lon % 360

lon_min_360 = convert_to_0_360(lon_min_sa)
lon_max_360 = convert_to_0_360(lon_max_sa)

# 4. Función para seleccionar región
def select_region_single_time(ds, lat_min, lat_max, lon_min_360, lon_max_360):
    """Seleccionar región para datos de un solo tiempo"""
    
    # Caso especial: región cruza el meridiano 0°
    if lon_min_360 > lon_max_360:
        print("  ⚠️ Región cruza el meridiano 0° - seleccionando dos partes")
        
        # Primera parte: lon_min_360 a 360
        region1 = ds.sel(
            latitude=slice(lat_max, lat_min),  # ERA5: lat de mayor a menor
            longitude=slice(lon_min_360, 360)
        )
        
        # Segunda parte: 0 a lon_max_360
        region2 = ds.sel(
            latitude=slice(lat_max, lat_min),
            longitude=slice(0, lon_max_360)
        )
        
        # Combinar
        ds_subset = xr.concat([region1, region2], dim='longitude')
        
    else:
        # Caso normal
        ds_subset = ds.sel(
            latitude=slice(lat_max, lat_min),
            longitude=slice(lon_min_360, lon_max_360)
        )
    
    # Convertir longitudes a -180..180 para graficar
    if 'longitude' in ds_subset.coords and ds_subset.longitude.values.max() > 180:
        lon_adj = ds_subset.longitude.values.copy()
        lon_adj[lon_adj > 180] = lon_adj[lon_adj > 180] - 360
        ds_subset = ds_subset.assign_coords(longitude=lon_adj)
        ds_subset = ds_subset.sortby('longitude')
    
    return ds_subset

# 5. Seleccionar región para ambos datasets
precip_subset = select_region_single_time(ds_precip, lat_min, lat_max, lon_min_360, lon_max_360)
z200_subset = select_region_single_time(ds_z200, lat_min, lat_max, lon_min_360, lon_max_360)


# 6. Extraer variables de precipitación
precip_var_name = 'tp'
precip_data_var = precip_subset[precip_var_name]

# Extraer datos (ya es un solo tiempo)
if 'valid_time' in precip_data_var.dims:
    precipitacion = precip_data_var.isel(valid_time=0)
else:
    precipitacion = precip_data_var

# 7. Extraer variables de altura geopotencial
z200_var_name = 'z'
z200_data_var = z200_subset[z200_var_name]

print(f"\n✓ Forma de Z200 original: {z200_data_var.shape}")

# IMPORTANTE: Seleccionar el nivel de presión (200 hPa)
# Como ya sabes que es 200 hPa (solo hay un nivel), puedes extraerlo directamente
if 'pressure_level' in z200_data_var.dims:
    print(f"  Niveles de presión disponibles: {z200_data_var.pressure_level.values}")
    # Seleccionar el primer (y único) nivel
    z200_at_level = z200_data_var.isel(pressure_level=0)
elif 'level' in z200_data_var.dims:
    print(f"  Niveles disponibles: {z200_data_var.level.values}")
    z200_at_level = z200_data_var.isel(level=0)
else:
    z200_at_level = z200_data_var

# Extraer el tiempo (primer y único tiempo)
if 'valid_time' in z200_at_level.dims:
    z200_single_time = z200_at_level.isel(valid_time=0)
else:
    z200_single_time = z200_at_level

print(f"  Forma de Z200 después de seleccionar nivel y tiempo: {z200_single_time.shape}")

# 8. Convertir unidades
# Convertir precipitación de m a mm
precipitacion_mm = precipitacion * 1000
precip_units = 'mm'

# Convertir altura geopotencial de m²/s² a metros geopotenciales
if 'units' in z200_single_time.attrs:
    print(f"  Unidades Z200: {z200_single_time.attrs['units']}")

g = 9.80665
z200_meters = z200_single_time / g
z200_units = 'mgp'


# 9. Preparar arrays para graficar
# Para precipitación
lon_precip = precipitacion_mm.longitude.values
lat_precip = precipitacion_mm.latitude.values
precip_data = precipitacion_mm.values

# Para Z200
lon_z200 = z200_meters.longitude.values
lat_z200 = z200_meters.latitude.values
z200_data = z200_meters.values  # Ahora debe ser 2D


# Crear mallas de coordenadas
if len(lon_precip.shape) == 1 and len(lat_precip.shape) == 1:
    lon_grid_precip, lat_grid_precip = np.meshgrid(lon_precip, lat_precip)
else:
    lon_grid_precip, lat_grid_precip = lon_precip, lat_precip

if len(lon_z200.shape) == 1 and len(lat_z200.shape) == 1:
    lon_grid_z200, lat_grid_z200 = np.meshgrid(lon_z200, lat_z200)
else:
    lon_grid_z200, lat_grid_z200 = lon_z200, lat_z200

# 10. Crear mapa
fig = plt.figure(figsize=(14, 10))
ax = plt.axes(projection=ccrs.PlateCarree())

# Configurar límites del mapa
ax.set_extent([lon_min_sa, lon_max_sa, lat_min, lat_max], crs=ccrs.PlateCarree())

# 11. Graficar precipitación
from matplotlib.colors import LogNorm
# ==============================================
# DELIMITACIÓN SIMPLE CON SHAPEFILE
# ==============================================
import geopandas as gpd

# 1. Cargar shapefile
amazonia_gdf = gpd.read_file(r"C:\Users\User\Desktop\mapitas\lim\Lim_Raisg.shp")

# 2. Asegurar proyección correcta
if amazonia_gdf.crs != 'EPSG:4326':
    amazonia_gdf = amazonia_gdf.to_crs('EPSG:4326')

# 3. Agregar al mapa
ax.add_geometries(amazonia_gdf.geometry,
                  crs=ccrs.PlateCarree(),
                  edgecolor='darkred',
                  facecolor='none',
                  linewidth=2.0,
                  linestyle='-',
                  alpha=0.9,
                  zorder=5)

# Configurar escala de colores para precipitación
vmin = 0.01
vmax = max(precip_data.max(), 10)

from matplotlib.colors import LinearSegmentedColormap

# Definir colores personalizados: desde un azul muy claro hasta un azul muy oscuro
colors = [
    '#F7FCF5',  # casi blanco
    '#E5F5E0',  # verde muy claro
    '#74C476',  # verde
    '#41B6C4',  # transición verde-azul
    '#2C7FB8',  # azul
    '#225EA8',  # azul más profundo
    '#081D58'   # azul muy oscuro
]


# Crear colormap personalizado
custom_cmap = LinearSegmentedColormap.from_list('azul_intenso', colors, N=256)

# Aplicar
img = ax.pcolormesh(lon_grid_precip, lat_grid_precip, precip_data,
                   transform=ccrs.PlateCarree(),
                   cmap=custom_cmap,  # Tu colormap personalizado
                   norm=LogNorm(vmin=0.01, vmax=vmax),
                   shading='auto',
                   alpha=0.8)

# 12. Graficar isohipsas (contornos de Z200)
# Determinar niveles automáticamente basado en los datos
z200_min = np.floor(z200_data.min() / 40) * 40
z200_max = np.ceil(z200_data.max() / 40) * 40
step = 39.05  # Paso en metros geopotenciales

levels_z200 = np.arange(z200_min, z200_max + step, step)



# Graficar contornos - AHORA con datos 2D
contours = ax.contour(lon_grid_z200, lat_grid_z200, z200_data,
                     levels=levels_z200,
                     colors= '#556B2F',  # Color rojo para mejor contraste
                     linewidths=1,
                     alpha=0.8,
                     transform=ccrs.PlateCarree())

# Etiquetar contornos
ax.clabel(contours, contours.levels[::2],  # Etiquetar cada segundo nivel
          inline=True, fontsize=8, fmt='%1.0f', colors='black')

# 13. Barra de colores para precipitación
cbar = plt.colorbar(img, ax=ax, orientation='horizontal',
                    pad=0.05, shrink=0.95, extend='both')
cbar.set_label(f'Precipitación Total ({precip_units})', fontsize=12, weight='bold')
cbar.ax.tick_params(labelsize=10)

# 14. Añadir características geográficas
ax.add_feature(cfeature.COASTLINE, linewidth=0.8, edgecolor='black')
ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='gray', alpha=0.7)

# Países
countries = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_0_countries',
    scale='50m',
    facecolor='none',
    edgecolor='gray',
    linewidth=0.3
)
ax.add_feature(countries, alpha=0.5)

# 15. Configurar grid
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray',
                  alpha=0.5, linestyle='--')

gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 10, 'color': 'black'}
gl.ylabel_style = {'size': 10, 'color': 'black'}

gl.xlocator = mticker.FixedLocator(np.arange(-120, 11, 10))
gl.ylocator = mticker.FixedLocator(np.arange(-60, 21, 10))

# 16. Título y etiquetas
fecha = "23-01-2026, 00:00 UTC"

titulo = f'Precipitación Total y Altura Geopotencial 200 hPa \n{fecha}'
plt.title(titulo, fontsize=14, weight='bold', pad=20)

# 17. Información adicional
info_text = (f"ERA5 | Precip: {precip_data.min():.2f}-{precip_data.max():.2f} mm")
ax.text(0.02, 0.02, info_text,
        transform=ax.transAxes, fontsize=9,
        verticalalignment='bottom',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# 18. Leyenda
from matplotlib.lines import Line2D

legend_elements = [
    Line2D([0], [0], color='#556B2F', lw=1.5, label=f'Isohipsas 200 hPa ({z200_units})'),
    Line2D([0], [0], color='blue', lw=4, alpha=0.8, label='Precipitación'),
    Line2D([0], [0], color='darkred', lw=4, alpha=0.8, label='Delimitacion de la Amazonia')
    ]

ax.legend(handles=legend_elements, loc='upper right', fontsize=9, framealpha=0.9)

# Distribución de precipitación
print(f"\n📈 Distribución de precipitación:")
rangos = [(0, 0.1), (0.1, 1), (1, 5), (5, 10), (10, 20), (20, 50), (50, 1000)]
for r_min, r_max in rangos:
    mask = (precip_data >= r_min) & (precip_data < r_max)
    porcentaje = (mask.sum() / precip_data.size) * 100
    if porcentaje > 0:
        print(f"  {r_min:5.1f}-{r_max:5.1f} mm: {porcentaje:5.1f}%")

# 20. Ajustar layout
plt.tight_layout()

# 21. Mostrar mapa
plt.show()

# 23. Cerrar datasets
ds_precip.close()
ds_z200.close()



#------------------altura geopotencial-------------

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import geopandas as gpd

# 1. Cargar los datos ERA5 para altura geopotencial a 200 hPa
ds_z200 = xr.open_dataset(r"C:\Users\User\Desktop\mapitas\vientolineas.nc")

# 2. Definir región de Sudamérica
lat_min, lat_max = -60, 20
lon_min_sa, lon_max_sa = -120, 10

# 3. Convertir longitudes de Sudamérica a rango 0-360 si es necesario
def convert_to_0_360(lon):
    """Convertir longitud de -180..180 a 0..360"""
    return lon % 360

lon_min_360 = convert_to_0_360(lon_min_sa)
lon_max_360 = convert_to_0_360(lon_max_sa)

# 4. Función para seleccionar región
def select_region_single_time(ds, lat_min, lat_max, lon_min_360, lon_max_360):
    """Seleccionar región para datos de un solo tiempo"""
    
    # Caso especial: región cruza el meridiano 0°
    if lon_min_360 > lon_max_360:
        print("  ⚠️ Región cruza el meridiano 0° - seleccionando dos partes")
        
        # Primera parte: lon_min_360 a 360
        region1 = ds.sel(
            latitude=slice(lat_max, lat_min),  # ERA5: lat de mayor a menor
            longitude=slice(lon_min_360, 360)
        )
        
        # Segunda parte: 0 a lon_max_360
        region2 = ds.sel(
            latitude=slice(lat_max, lat_min),
            longitude=slice(0, lon_max_360)
        )
        
        # Combinar
        ds_subset = xr.concat([region1, region2], dim='longitude')
        
    else:
        # Caso normal
        ds_subset = ds.sel(
            latitude=slice(lat_max, lat_min),
            longitude=slice(lon_min_360, lon_max_360)
        )
    
    # Convertir longitudes a -180..180 para graficar
    if 'longitude' in ds_subset.coords and ds_subset.longitude.values.max() > 180:
        lon_adj = ds_subset.longitude.values.copy()
        lon_adj[lon_adj > 180] = lon_adj[lon_adj > 180] - 360
        ds_subset = ds_subset.assign_coords(longitude=lon_adj)
        ds_subset = ds_subset.sortby('longitude')
    
    return ds_subset

# 5. Seleccionar región para altura geopotencial
z200_subset = select_region_single_time(ds_z200, lat_min, lat_max, lon_min_360, lon_max_360)

# 6. Extraer variables de altura geopotencial
z200_var_name = 'z'
z200_data_var = z200_subset[z200_var_name]

print(f"\n✓ Forma de Z200 original: {z200_data_var.shape}")

# Seleccionar el nivel de presión (200 hPa)
if 'pressure_level' in z200_data_var.dims:
    print(f"  Niveles de presión disponibles: {z200_data_var.pressure_level.values}")
    z200_at_level = z200_data_var.isel(pressure_level=0)
elif 'level' in z200_data_var.dims:
    print(f"  Niveles disponibles: {z200_data_var.level.values}")
    z200_at_level = z200_data_var.isel(level=0)
else:
    z200_at_level = z200_data_var

# Extraer el tiempo (primer y único tiempo)
if 'valid_time' in z200_at_level.dims:
    z200_single_time = z200_at_level.isel(valid_time=0)
else:
    z200_single_time = z200_at_level

print(f"  Forma de Z200 después de seleccionar nivel y tiempo: {z200_single_time.shape}")

# 7. Convertir unidades de altura geopotencial
if 'units' in z200_single_time.attrs:
    print(f"  Unidades Z200: {z200_single_time.attrs['units']}")

g = 9.80665
z200_meters = z200_single_time / g
z200_units = 'mgp'

# 8. Preparar arrays para graficar
lon_z200 = z200_meters.longitude.values
lat_z200 = z200_meters.latitude.values
z200_data = z200_meters.values

# Crear malla de coordenadas
if len(lon_z200.shape) == 1 and len(lat_z200.shape) == 1:
    lon_grid_z200, lat_grid_z200 = np.meshgrid(lon_z200, lat_z200)
else:
    lon_grid_z200, lat_grid_z200 = lon_z200, lat_z200

# 9. Crear mapa
fig = plt.figure(figsize=(14, 10))
ax = plt.axes(projection=ccrs.PlateCarree())

# Configurar límites del mapa
ax.set_extent([lon_min_sa, lon_max_sa, lat_min, lat_max], crs=ccrs.PlateCarree())

# 10. ================= Sombreado de altura geopotencial =================
# Crear niveles para el sombreado
z_min = np.floor(z200_data.min() / 10) * 10
z_max = np.ceil(z200_data.max() / 10) * 10
step_shading = 10  # Paso pequeño para sombreado suave

levels_shading = np.arange(z_min, z_max + step_shading, step_shading)

# Definir colormap para el sombreado
# Puedes usar: 'viridis', 'plasma', 'coolwarm', 'Spectral_r', 'terrain', 'gist_earth'
cmap_shading = plt.cm.Spectral_r  # O elige otro colormap

# Graficar sombreado (fondo)
shading = ax.contourf(lon_grid_z200, lat_grid_z200, z200_data,
                     levels=levels_shading,
                     cmap=cmap_shading,
                     extend='both',
                     transform=ccrs.PlateCarree(),
                     alpha=0.8)

# 11. ================= Isohipsas (líneas de contorno) =================
# Niveles más espaciados para contornos
step_contour = 40  # Paso mayor para contornos
levels_contour = np.arange(z_min, z_max + step_contour, step_contour)

# Graficar contornos
contours = ax.contour(lon_grid_z200, lat_grid_z200, z200_data,
                     levels=levels_contour,
                     colors='darkred',  # Color negro para contornos
                     linewidths=1.2,
                     alpha=0.9,
                     transform=ccrs.PlateCarree())

# Etiquetar contornos
ax.clabel(contours, contours.levels[::2],  # Etiquetar cada segundo nivel
          inline=True, fontsize=9, fmt='%1.0f', colors='darkred')


# 13. ================= Características geográficas =================
ax.add_feature(cfeature.COASTLINE, linewidth=0.8, edgecolor='black')
ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='darkgray', alpha=0.7)

# Añadir países
countries = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_0_countries',
    scale='50m',
    facecolor='none',
    edgecolor='gray',
    linewidth=0.3
)
ax.add_feature(countries, alpha=0.5)

# Añadir océanos y tierra para mejor contraste
ax.add_feature(cfeature.OCEAN, color='lightblue', alpha=0.2)
ax.add_feature(cfeature.LAND, color='lightgray', alpha=0.1)

# 14. ================= Configurar grid =================
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray',
                  alpha=0.5, linestyle='--')

gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 10, 'color': 'black'}
gl.ylabel_style = {'size': 10, 'color': 'black'}

gl.xlocator = mticker.FixedLocator(np.arange(-120, 11, 10))
gl.ylocator = mticker.FixedLocator(np.arange(-60, 21, 10))

# 15. ================= Barra de color para el sombreado =================
cbar = plt.colorbar(shading, ax=ax, orientation='horizontal',
                    pad=0.08, shrink=0.9, extend='both')
cbar.set_label(f'Altura Geopotencial 200 hPa ({z200_units})', fontsize=12, weight='bold')
cbar.ax.tick_params(labelsize=10)

# 16. ================= Título y etiquetas =================
# Obtener fecha del dataset
try:
    time_data = z200_single_time.time
    if hasattr(time_data, 'values'):
        fecha_str = str(time_data.values)[:19]
    else:
        fecha_str = "Fecha no disponible"
except:
    fecha_str = "23-01-2026, 00:00 UTC"  # Fecha por defecto

titulo = f'Altura Geopotencial 200 hPa\n{fecha_str}'
plt.title(titulo, fontsize=14, weight='bold', pad=20)



# Distribución de valores
print(f"\n📈 Distribución de valores:")
percentiles = [0, 10, 25, 50, 75, 90, 100]
for p in percentiles:
    value = np.percentile(z200_data, p)
    print(f"  Percentil {p:3d}%: {value:7.1f} {z200_units}")

# 20. ================= Ajustar layout y mostrar =================
plt.tight_layout()
plt.show()

# 21. ================= Cerrar dataset =================
ds_z200.close()

print("\n✓ Mapa generado exitosamente!")