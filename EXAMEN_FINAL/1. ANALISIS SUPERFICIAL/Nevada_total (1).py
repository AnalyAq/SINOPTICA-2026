import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# 1. Cargar datos del evento de Junio 2025
# El archivo subido contiene la variable 'sf' (Snowfall)
ds = xr.open_dataset(r"C:\Users\emmnu\OneDrive\Escritorio\Sinop_vera\final\data_stream-oper_stepType-accum.nc")
ds
# 2. Selección de la variable y conversión de unidades
# 'sf' viene en metros, lo multiplicamos por 1000 para tener milímetros (mm)
snow_mm = ds['sf'] * 1000

# 3. Cálculo de la Nieve Total Acumulada (24 al 30 de Junio)
# Sumamos todos los pasos de tiempo del archivo para ver el acumulado total
total_snow = snow_mm.sum(dim='valid_time')

# 4. Configuración del Mapa
fig = plt.figure(figsize=(12, 9))
ax = plt.axes(projection=ccrs.PlateCarree())

# Definir la extensión para ver Atacama y alrededores
# Longitud: 75°W a 65°W | Latitud: 35°S a 15°S
extent = [-76, -64, -35, -18]
ax.set_extent(extent, crs=ccrs.PlateCarree())

# Añadir características geográficas
ax.add_feature(cfeature.COASTLINE, linewidth=1.5)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=1)
ax.add_feature(cfeature.STATES, linestyle='--', alpha=0.5) # Para ver límites de regiones en Chile

# 5. Sombreado de la Nieve (Colormap específico para nieve: 'PuBu' o 'Blues')
levels = [0.1, 1, 2, 5, 10, 15, 20, 25, 30, 40, 50] # Niveles en mm
cf = ax.contourf(ds.longitude, ds.latitude, total_snow, 
                 levels=levels, 
                 cmap='PuBu', 
                 extend='max', 
                 transform=ccrs.PlateCarree())

# Barra de colores
cbar = plt.colorbar(cf, ax=ax, orientation='vertical', pad=0.03, shrink=0.8)
cbar.set_label('Nieve Total Acumulada (mm)', fontsize=12)

# 6. Títulos y etiquetas
plt.title('Nevada Extrema en el Desierto de Atacama\nAcumulado Total: 24 al 30 de Junio, 2025', 
          fontsize=14, fontweight='bold')

# Añadir etiquetas de Lat/Lon
gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False

plt.tight_layout()
plt.savefig('nevada_atacama_junio2025.png', dpi=300)
plt.show()

# Cerrar el dataset
ds.close()