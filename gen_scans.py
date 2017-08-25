import numpy as np
import itertools
import math
import matplotlib.pyplot as plt
from gen_otfm_grid import get_coords,gen_left,gen_right,plot_fluxcal,plot_cal,get_center
from astropy.coordinates import SkyCoord


def gen_rect_mosaic(theta_row, theta_hex, length, height):
    """ Generate the pointing centers for a rectangular mosaic """
    row_starts = gen_heights(theta_row, height)[::-1]
    nrows = len(row_starts)
    short_row = gen_row(theta_hex, length)
    long_row = np.append(short_row - theta_hex/2, short_row[-1] + theta_hex/2)
    nshort = math.ceil(nrows / 2)
    nlong = nrows - nshort
    npointings = nshort * len(short_row) + nlong * len(long_row)
    grid = np.zeros((npointings,2))
    
    # add pointings in order, starting at the top right, raster scan
    ii = 0
    jj = 0
    while ii < nrows:
        toadd = np.array(list(
            itertools.product(long_row[::-1], [row_starts[ii]])))
        grid[jj:jj+len(toadd)] = toadd
        jj += len(toadd)
        toadd = np.array(list(
            itertools.product(short_row, [row_starts[ii+1]])))
        grid[jj:jj+len(toadd)] = toadd
        jj += len(toadd)
        ii += 2

    return grid


def gen_heights(theta_row, height):
    """ Generate the vertical starting points for each row of the mosaic
    Center it on 0, 0 """
    dim = height/2
    return np.arange(-dim, dim, theta_row)


def gen_row(theta_hex, length):
    """ Generate pointing centers for an individual row 
    Center it on 0, 0 """
    dim = (length+theta_hex)/2
    return np.arange(-dim, dim, theta_hex) 


def rotate_mosaic(coords, theta):
    """ 
    Rotate a mosaic by an angle theta
    
    Parameters
    ----------
    coords: (ncoords, 2)
    theta: angle to rotate through 
    """
    R = np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]])
    return np.dot(coords,R)


if __name__=="__main__":
    # MAXI target coordinates
    c = get_coords()
    x_cen, y_cen = get_center(c)

    gen_left(c, 0, 1, 1, 2)
    gen_right(c, 2, 3, 3, 0)
    top_left = np.array([c[1].ra.value, c[1].dec.value])
    bottom_left = np.array([c[0].ra.value, c[0].dec.value])
    bottom_right = np.array([c[3].ra.value, c[3].dec.value])
    top_right = np.array([c[2].ra.value, c[2].dec.value])
    #plot_fluxcal()
    #plot_cal()

    # spacing between pointings in all rows
    theta_hex = 9.9 / 60 # degrees

    # spacing between rows
    theta_row = 8.57 / 60 # degrees

    # dimensions of mosaic
    length = 2.2 # degrees
    height = 0.45 # degrees
    coords = gen_rect_mosaic(theta_row, theta_hex, length, height)

    # put in the right order
    # in this case, start at the top right, so 

    # angle of mosaic
    delta_dec = (top_right[1]-top_left[1])
    delta_ra = (top_right[0]-top_left[0])
    theta_rot = np.arctan(delta_dec/delta_ra)

    # rotate
    rotated_coords = rotate_mosaic(coords,theta_rot)

    # translate
    x = rotated_coords[:,0] + x_cen
    y = rotated_coords[:,1] + y_cen

    # print scan file
    outputf = open("scanlist.txt", "w")
    coords_raw = SkyCoord(x,y,unit="deg")
    coords =coords_raw.to_string("hmsdms")
    for ii,coord in enumerate(coords):
        ra = coord.split()[0]
        dec = coord.split()[1][1:]
        outputf.write("ptng%s;;;;%s;%s;;;;;\n" %(ii,ra,dec))
        #outputf.write("%s;;;;%s;%s;;;;;+\n" (str(ii),ra,dec))
    outputf.close()

    # plot
    plt.scatter(x, y, c='k')
    plt.scatter(x_cen, y_cen, c='r')
    plt.xlabel("RA (deg)")
    plt.ylabel("Dec (deg)")
    plt.savefig("pointing_centers.png")
