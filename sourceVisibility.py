"""Script from Vikram"""

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
import sys

from astropy.utils.iers import conf
conf.auto_max_age = None

def get_moon(time, location, pressure=None):
    """
    Position of the Earth's moon.

    Currently uses PyEphem to calculate the position of the moon by default
    (requires PyEphem to be installed). Set ``use_pyephem`` to `False` to
    calculate the moon position with jplephem (requires jplephem to be
    installed).

    Parameters
    ----------
    time : `~astropy.time.Time` or see below
        Time of observation. This will be passed in as the first argument to
        the `~astropy.time.Time` initializer, so it can be anything that
        `~astropy.time.Time` will accept (including a `~astropy.time.Time`
        object).

    location : `~astropy.coordinates.EarthLocation`
        Location of the observer on Earth

    pressure : `None` or `~astropy.units.Quantity` (optional)

    Returns
   -------

    moon_sc : `~astropy.coordinates.SkyCoord`
        Position of the moon at ``time``
    """
    if not isinstance(time, Time):
        time = Time(time)

    try:
        import ephem
    except ImportError:
        raise ImportError("The get_moon function currently requires "
                          "PyEphem to compute the position of the moon.")

    moon = ephem.Moon()
    obs = ephem.Observer()
    obs.lat = location.latitude.to_string(u.deg, sep=':')
    obs.lon = location.longitude.to_string(u.deg, sep=':')
    obs.elevation = location.height.to(u.m).value
    if pressure is not None:
        obs.pressure = pressure.to(u.bar).value*1000.0

    if time.isscalar:
        obs.date = time.datetime
        moon.compute(obs)
        moon_alt = float(moon.alt)
        moon_az = float(moon.az)
        moon_distance = moon.earth_distance
    else:
        moon_alt = []
        moon_az = []
        moon_distance = []
        for t in time:
            obs.date = t.datetime
            moon.compute(obs)
            moon_alt.append(float(moon.alt))
            moon_az.append(float(moon.az))
            moon_distance.append(moon.earth_distance)
    return SkyCoord(alt=moon_alt*u.rad, az=moon_az*u.rad,
                    distance=moon_distance*u.AU,
                    frame=AltAz(location=location, obstime=time))



# observatories and offsets
cerro_pachon = EarthLocation(lat=-30.240722*u.deg, lon=-70.736583*u.deg, height=2700.*u.m)
palomar = EarthLocation(lat=33.3558*u.deg, lon=-116.8639*u.deg, height=1700.*u.m)
mauna_kea = EarthLocation(lat=19.82636*u.deg, lon=-155.47501*u.deg, height=4145.*u.m)
ovro = EarthLocation(lat=37.2339*u.deg, lon=-118.282*u.deg, height=1222.*u.m)
vla = EarthLocation(lat=34.0784*u.deg, lon=-107.6184*u.deg, height=2124.*u.m)

observatory = {'cerro_pachon': cerro_pachon, 'palomar': palomar, 'mauna_kea': mauna_kea, 'ovro': ovro, 'vla': vla}
utcoffset = {'cerro_pachon': -3.*u.hour, 'palomar': -8.*u.hour, 'mauna_kea': -10.*u.hour, 'ovro': -8.*u.hour, 'vla': -7.*u.hour}
elevation_limit = 8.0*u.deg # degrees

# source, time and obs
obs = 'vla'
src = SkyCoord('13 09 00', '-23 23 00', unit=(u.hourangle, u.deg))
time = Time('2017-08-18 00:00:00')-utcoffset[obs]

# plotting
deltas = np.linspace(-12., 12., 100)*u.hour
times = time+deltas
altazframe = AltAz(obstime=times, location=observatory[obs])
sunaltazs = get_sun(times).transform_to(altazframe)
#moonaltazs = get_moon(times,observatory[obs]).transform_to(altazframe)
srcaltazs = src.transform_to(altazframe)  

# lst at ut date
#print (time+utcoffset[obs]).sidereal_time('mean',observatory[obs].longitude.degree)

srcaltazs.alt[srcaltazs.alt <= elevation_limit] = elevation_limit
plt.scatter(deltas, srcaltazs.alt, c=np.log10(srcaltazs.secz), label='Source', lw=0, s=8) 
plt.colorbar().set_label('log10 Airmass')  
plt.plot(deltas, sunaltazs.alt, 'y-', label='Sun')  
#plt.plot(deltas, moonaltazs.alt, 'g--', label='Moon')  
plt.fill_between(deltas.to('hr').value, 0, 90, sunaltazs.alt < -0*u.deg, color='0.6', zorder=0)  
plt.fill_between(deltas.to('hr').value, 0, 90, sunaltazs.alt < -18*u.deg, color='0.4', zorder=0)  

plt.legend(loc='upper left')  
plt.xlim(-12, 12)  
plt.xticks(np.arange(13)*2 -12)  
plt.ylim(0, 90)  
plt.xlabel('Hours from  local midnight')  
plt.ylabel('Altitude [deg]')  

plt.show()

