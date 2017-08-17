import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord


def get_coords():
    c1 = SkyCoord(157.035, -29.211, frame='icrs', unit='deg')
    c2 = SkyCoord(156.551, -28.839, frame='icrs', unit='deg') 
    c3 = SkyCoord(157.863, -26.870, frame='icrs', unit='deg')
    c4 = SkyCoord(158.343, -27.238, frame='icrs', unit='deg')
    c5 = SkyCoord(157.476, -29.111, frame='icrs', unit='deg')
    c6 = SkyCoord(157.014, -28.758, frame='icrs', unit='deg')
    c7 = SkyCoord(158.516, -26.503, frame='icrs', unit='deg')
    c8 = SkyCoord(158.974, -26.850, frame='icrs', unit='deg')
    c = [c1, c2, c3, c4, c5, c6, c7, c8]
    return c


def plot_rec(cs):
    colors = ['red', 'orange', 'yellow', 'green', 'cyan', 'blue', 'purple', 'black']
    for ii,c in enumerate(cs):
        plt.scatter(c.ra, c.dec, c=colors[ii], label=ii, lw=0, s=50)
    plt.legend()


def fitline(x1, x2, y1, y2):
    m = (y2 - y1) / (x2 - x1)
    b = y2 - m*x2
    return m, b


def gen_left(c):
    # Starting points are the line from C1 to C2
    # + the line from C2 to C3
    m, b = fitline(c[1].ra, c[2].ra, c[1].dec, c[2].dec)
    x = np.linspace(c[2].ra, c[1].ra)
    y = m * x + b
    y1 = np.array(y)
    plt.plot(x, y1, lw=1, color='k')

    m, b = fitline(c[1].ra, c[0].ra, c[1].dec, c[0].dec)
    x = np.linspace(c[1].ra, c[0].ra)
    y = m * x + b
    y2 = np.array(y)
    plt.plot(x, y2, lw=1, color='k')
    return y1,y2


def gen_right(c):
    # End points are the line from C1 to C4
    # + the line from C3 to C4
    m, b = fitline(c[4].ra, c[7].ra, c[4].dec, c[7].dec)
    x = np.linspace(c[4].ra, c[7].ra)
    y = m * x + b
    y = np.array(y)
    plt.plot(x, y, lw=1, color='k')

    m, b = fitline(c[6].ra, c[7].ra, c[6].dec, c[7].dec)
    x = np.linspace(c[6].ra, c[7].ra)
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
    c = get_coords()
    plot_rec(c)
    gen_left(c)
    gen_right(c)
# 
#     plt.xlabel("RA (deg)")
#     plt.ylabel("Dec (deg)")
# 
#     ras = np.sort(np.array([c1.ra.value, c2.ra.value, c3.ra.value, c4.ra.value]))
#     decs = np.sort(np.array([c1.dec.value, c2.dec.value, c3.dec.value, c4.dec.value]))
# 
#     dec_grid = np.arange(decs[2], decs[3]+0.1, 0.125)
#     x = np.repeat(min(ras),len(dec_grid))
#     plt.scatter(x, dec_grid, color='k', s=1)
#     coord_grid = SkyCoord(x, dec_grid, unit='deg')
#     for coord in coord_grid:
#         print(coord.ra.hms)
#         print(coord.dec.dms)
# 
#     dec_grid = np.arange(decs[0], decs[1]+0.1, 0.125)
#     x = np.repeat(max(ras),len(dec_grid))
#     plt.scatter(x, dec_grid, color='k', s=1)
#     coord_grid = SkyCoord(x, dec_grid, unit='deg')
#     for coord in coord_grid:
#         print(coord.ra.hms)
#         print(coord.dec.dms)
# 
# 
#     plot_fluxcal()
#     plot_cal()
    plt.show()
#     plt.legend(loc="lower right")
#     plt.ylim(24,34)
# 
if __name__=="__main__":
    run()
