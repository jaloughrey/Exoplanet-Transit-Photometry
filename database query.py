#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np

from astropy.visualization import astropy_mpl_style, quantity_support

plt.style.use(astropy_mpl_style)
quantity_support()


# In[2]:


import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time


# In[3]:


from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
NasaExoplanetArchive.query_object("wasp-52b") 


# In[4]:


NasaExoplanetArchive.query_criteria(table=("TD"), select="*")


# In[28]:


from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
NasaExoplanetArchive.TAP_TABLES


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[5]:


#INPUT SELECTED OBJECT 
obj = input('Enter object:')

#COORDINATES OF OBJECT (uses simbad)
obj_coord = SkyCoord.from_name(obj)
print(obj_coord)

#LOCATION AND DATE OF OBSERVATION
from datetime import datetime
from datetime import date

#cardiff location
Cardiff = EarthLocation(lat=51.5*u.deg, lon=-3.2*u.deg, height=40*u.m) 

# current date and time in form YYYY-mm-dd H:M:S
now = datetime.now()
current_datetime = now.strftime("%Y-%m-%d %H:%M:%S")
time = Time(current_datetime) 
today = date.today()


#GET OBJECT ALTITUDE 
objaltaz = obj_coord.transform_to(AltAz(obstime=time,location=Cardiff))
print(f"Object Altitude = {objaltaz.alt:.2}")


#CREATE TIME ARRAY AND FRAME FOR CURRENT DATE AND LOCATION 
midnight = Time(f'{today} 00:00:00') 
delta_midnight = np.linspace(-12, 12, 1000)*u.hour
frame_today = AltAz(obstime=midnight+delta_midnight,
                          location=Cardiff)
objaltazs_today = obj_coord.transform_to(frame_today)

#GET OBJECTS AIRMASS
objairmasss_today = objaltazs_today.secz


#GET SUN LOCATION
from astropy.coordinates import get_sun

#delta_midnight = np.linspace(-12, 12, 1000)*u.hour
times_today = midnight + delta_midnight
#frame_today = AltAz(obstime=times_today, location=Cardiff)
sunaltazs_today = get_sun(times_today).transform_to(frame_today)

#MOON LOCATION
from astropy.coordinates import get_body

moon_today = get_body("moon", times_today)
moonaltazs_today = moon_today.transform_to(frame_today)

#PLOT ALTITUDE OF OBJECT WITH TIME INCLUDING DARK AREAS REPRESENTING NIGHT. 
plt.plot(delta_midnight, sunaltazs_today.alt, color='r', label='Sun')
plt.plot(delta_midnight, moonaltazs_today.alt, color=[0.75]*3, ls='--', label='Moon')
plt.scatter(delta_midnight, objaltazs_today.alt,
            c=objaltazs_today.az, label=obj, lw=0, s=8,
            cmap='viridis')
plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                 sunaltazs_today.alt < -0*u.deg, color='0.5', zorder=0)
plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                 sunaltazs_today.alt < -18*u.deg, color='k', zorder=0)
plt.colorbar().set_label('Azimuth [deg]')
plt.legend(loc='upper left')
plt.xlim(-12*u.hour, 12*u.hour)
plt.xticks((np.arange(13)*2-12)*u.hour)
plt.ylim(0*u.deg, 90*u.deg)
plt.xlabel('Hours from GMT Midnight')
plt.ylabel('Altitude [deg]')
plt.show()


# In[6]:


#PLOT AIRMASS GRAPH 
plt.plot(delta_midnight, objairmasss_today)
plt.xlim(-12, 12)
plt.ylim(10, 1)
plt.xlabel('Hours from GMT Midnight')
plt.ylabel('Airmass [Sec(z)]')
plt.show()


# In[ ]:




