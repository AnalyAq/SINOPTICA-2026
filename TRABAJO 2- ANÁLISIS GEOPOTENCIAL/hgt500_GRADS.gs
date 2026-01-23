'reinit'
'cd "C:\Users\ANALY AVALOS\Downloads"
'set display color white'


'sdfopen era5_2017.nc'


'set lon -90 30'
'set lat -60 15'
'set lev 500'


'set gxout shaded'
'set clevs 5200 5250 5300 5350 5400 5450 5500 5550 5600 5650 5700'
'set ccols 0 21 23 25 27 29 31 33 35 37 39'
'set mpdset hires'
'set grid off'
'set grads off'


t = 1
while (t <= 168)

  'set t 't
  'clear'


  'q time'
  tiempo = sublin(result,1)
  tiempo = substr(tiempo,8,11)

  

  'd z/9.80665'

  
  'run cbarn.gs'
  'draw string 1.2 0.9 Altura geopotencial (m)'

  
  'draw title Altura geopotencial a 500 hPa | ERA5 | 'tiempo

  
  fname = 'hgt500_'t'.png'
  'gxprint 'fname

  t = t + 1
endwhile

