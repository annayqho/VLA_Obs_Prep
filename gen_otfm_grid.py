import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord


def get_coords():
    c1 = SkyCoord(11.659, 25.445, frame='icrs', unit='deg')
    c2 = SkyCoord(11.527, 24.975, frame='icrs', unit='deg') 
    c3 = SkyCoord(13.981, 24.140, frame='icrs', unit='deg')
    c4 = SkyCoord(14.121, 24.608, frame='icrs', unit='deg') 
    return c1,c2,c3,c4


def plot_rec():
    plt.scatter(c1.ra, c1.dec, c='blue', label="C1", lw=0, s=50)
    plt.scatter(c2.ra, c2.dec, c='magenta', label="C2", lw=0, s=50)
    plt.scatter(c3.ra, c3.dec, c='green', label="C3", lw=0, s=50)
    plt.scatter(c4.ra, c4.dec, c='black', label="C4", lw=0, s=50)
    plt.legend()


def fitline(x1, x2, y1, y2):
    m = (y2 - y1) / (x2 - x1)
    b = y2 - m*x2
    return m, b


def gen_left(c1,c2,c3,c4):
    # Starting points are the line from C1 to C2
    # + the line from C2 to C3
    m, b = fitline(c1.ra, c2.ra, c1.dec, c2.dec)
    x = np.linspace(c2.ra, c1.ra)
    y = m * x + b
    y1 = np.array(y)
    plt.plot(x, y1, lw=1, color='k')

    m, b = fitline(c2.ra, c3.ra, c2.dec, c3.dec)
    x = np.linspace(c2.ra, c3.ra)
    y = m * x + b
    y2 = np.array(y)
    plt.plot(x, y2, lw=1, color='k')
    return y1,y2


def gen_right(c1,c2,c3,c4):
    # End points are the line from C1 to C4
    # + the line from C3 to C4
    m, b = fitline(c1.ra, c4.ra, c1.dec, c4.dec)
    x = np.linspace(c1.ra, c4.ra)
    y = m * x + b
    y = np.array(y)
    plt.plot(x, y, lw=1, color='k')

    m, b = fitline(c3.ra, c4.ra, c3.dec, c4.dec)
    x = np.linspace(c3.ra, c4.ra)
    y = m * x + b
    y = np.array(y)
    plt.plot(x, y, lw=1, color='k')


def plot_fluxcal():
    # 3C48
    fluxcal = SkyCoord("01h37m41.1s", "+33d09m32s", frame='icrs')
    plt.scatter(fluxcal.ra, fluxcal.dec, c='black')
    plt.text(fluxcal.ra.value+0.1, fluxcal.dec.value+0.1, "3C48", fontsize=12)


def plot_cal():
    # J0122+2502
    cal = SkyCoord("01h22m38.82", "+25d02m31.8", frame='icrs')
    plt.scatter(cal.ra, cal.dec, c='black')
    plt.text(cal.ra.value+0.1, cal.dec.value+0.1, "J0122+2502", fontsize=12)


def run():
    c1 = SkyCoord(11.659, 25.445, frame='icrs', unit='deg')
    c2 = SkyCoord(11.527, 24.975, frame='icrs', unit='deg') 
    c3 = SkyCoord(13.981, 24.140, frame='icrs', unit='deg')
    c4 = SkyCoord(14.121, 24.608, frame='icrs', unit='deg') 
    #plot_rec()

    gen_left()
    gen_right()

    plt.xlabel("RA (deg)")
    plt.ylabel("Dec (deg)")

    ras = np.sort(np.array([c1.ra.value, c2.ra.value, c3.ra.value, c4.ra.value]))
    decs = np.sort(np.array([c1.dec.value, c2.dec.value, c3.dec.value, c4.dec.value]))

    dec_grid = np.arange(decs[2], decs[3]+0.1, 0.125)
    x = np.repeat(min(ras),len(dec_grid))
    plt.scatter(x, dec_grid, color='k', s=1)
    coord_grid = SkyCoord(x, dec_grid, unit='deg')
    for coord in coord_grid:
        print(coord.ra.hms)
        print(coord.dec.dms)

    dec_grid = np.arange(decs[0], decs[1]+0.1, 0.125)
    x = np.repeat(max(ras),len(dec_grid))
    plt.scatter(x, dec_grid, color='k', s=1)
    coord_grid = SkyCoord(x, dec_grid, unit='deg')
    for coord in coord_grid:
        print(coord.ra.hms)
        print(coord.dec.dms)


    plot_fluxcal()
    plot_cal()
    plt.show()
    plt.legend(loc="lower right")
    plt.ylim(24,34)

if __name__=="__main__":
    run()
