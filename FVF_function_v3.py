from matplotlib.patches import Circle, Rectangle
import numpy as np

from random import *
from math import *
from sklearn.neighbors import KDTree
from scipy import spatial
from astropy.stats import RipleysKEstimator
from matplotlib import interactive, pyplot as plt


def intersection(r1, cent1, r2, cent2, tol):
    distance = ((cent1[0] - cent2[0]) ** 2 + (cent1[1] - cent2[1])
                ** 2) ** 0.5  # pythagoras thm
    min_distance = r1 + r2 + tol  # minimum distance between two centres of circle

    if distance >= min_distance:
        overlap = False
    else:
        overlap = True
    return overlap


def check_a_circle(cent1, centlist, tol):

    r1 = cent1[2]
    for cent2 in centlist:
        r2 = cent2[2]
        temp = intersection(r1, cent1, r2, cent2, tol)
        if temp == True:
            return temp
    return temp


def ran(x1, x2):
    # randomise
    # used to generate random points in the chosen region
    return ((x2 - x1) * np.random.rand() + x1)


def add_vff(radius, width, height):
    return ((pi * radius ** 2) / (width*height))


def plot_pdfgraph(x_range_int, y_range, y, line_label, x_label, y_label, line_type="--ro"):
    # x_range_int is a list
    # y_range is a list
    plt.plot(np.arange(x_range_int[0], x_range_int[1],
                       x_range_int[2]), y, line_type, label=line_label)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.ylim(y_range[0], y_range[1])
    plt.legend(loc='upper right')
    plt.show()


def modelprint(width, height, lw, vff, centlist):
    # generate model
    print(centlist)
    print('FVF:', vff, '\n')   # The final FVF

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rectan = Rectangle((0, 0), width, height, linewidth=lw, edgecolor='red', facecolor='w',
                       linestyle='solid', alpha=1.0)
    ax.add_patch(rectan)

    for circle in centlist:
        X = circle[0]
        Y = circle[1]
        R = circle[2]
        draw_circle = Circle(xy=(X, Y), radius=R, color='black', alpha=1.0)
        ax.add_patch(draw_circle)
        plt.axis('scaled')
        # change limits of x or y axis so that equal increemets of x and y have the same length
        plt.axis('equal')
    plt.show()


def fillcircle(centlist, vff):
    og = len(centlist)
    print("start", og)
    i = 0
    while i < len(centlist):
        # overlap
        # tilt one degree
        if centlist[i][1] > 0.8 or centlist[i][1] < 0 or centlist[i][0] > 0.8 or centlist[i][0] < 0:
            pass
        else:
            radius = 2*centlist[i][2] + tol
            angle = 0
            while angle < 360:
                x_diff = radius * cos(angle/180*pi)
                y_diff = radius * sin(angle/180*pi)
                imagcent = [centlist[i][0] + x_diff,
                            centlist[i][1] + y_diff, centlist[i][2]]
                # check imagcent

                if imagcent[1] > 0.8 or imagcent[1] < 0 or imagcent[0] > 0.8 or imagcent[0] < 0:
                    pass
                else:
                    indicator = check_a_circle(imagcent, centlist, tol)
                    if indicator == False:
                        # never overlap
                        centlist = centlist + [imagcent]
                        vff = vff + add_vff(imagcent[2], width, height)
                        print("increase")
                angle += 1
            print('done checking ', i, "/", len(centlist), og)
        i += 1
    return centlist, vff


"""START OF COMPUTATION"""
# Variable parameters
count = 0
rmax = 0.016  # the maxmum radius of the fibres
rmin = 0.016  # the minmum radius of the fibres
width = rmax*50  # 0.8   # the length of the RVE
height = rmax*50  # 0.8 # the width of the RVE
tol = 0.0005  # the minmum distance of two circles (except for the radius)
vf = 0.55    # the FVF in the RVE
# with every circle added to the model
# cap at 0.55
#

# Algorithm of generating random distributing fibres
# Ratio of each part to the whole RVE area
ratio_PartI = (width - 2 * rmax) * (height - 2 * rmax) / \
    (width + 2 * rmax) * (height + 2 * rmax)
# ratio for centre
ratio_PartII = (2 * 2 * rmax * (height - 2 * rmax)) / \
    (width + 2 * rmax) * (height + 2 * rmax)
# ratio for vertical sides
ratio_PartIII = (2 * 2 * rmax * (width - 2 * rmax)) / \
    (width + 2 * rmax) * (height + 2 * rmax)
# ratio for horizontal sides
ratio_PartIV = (4 * 2 * rmax * 2 * rmax) / \
    (width + 2 * rmax) * (height + 2 * rmax)
# ratio for 4 corners
#

centlist = [[ran(0, width), ran(0, height), ran(rmin, rmax)]]
vff = (pi * centlist[0][2] ** 2) / (width*height)
# vff = total area of circle / area of box

print("begin")
print("vff: ", vff, "\n Vf: ", vf)
print("Ratio 1: ", ratio_PartI)
print("Ratio 2: ", ratio_PartII)
print("Ratio 3: ", ratio_PartIII)
print("Ratio 4: ", ratio_PartIV)

# NO ISSUES UP TILL HERE

ratio_centre = ratio_PartI
ratio_heights = ratio_PartI + ratio_PartII
ratio_widths = ratio_PartI + ratio_PartII + ratio_PartIII
ratio_corners = ratio_PartI + ratio_PartII + ratio_PartIII + ratio_PartIV


while vff < vf:
    prob = random()
    if prob < ratio_centre:
        X = ran(0 + 1 * rmax, width - 1 * rmax)
        Y = ran(0 + 1 * rmax, height - 1 * rmax)
        R = ran(rmin, rmax)
        newcircle = [X, Y, R]
        check = check_a_circle(newcircle, centlist, tol)
        if check == False:
            centlist = centlist + [newcircle]
            vff = vff + add_vff(newcircle[2], width, height)

    elif ratio_centre <= prob < ratio_heights:
        X = ran(0 - 1 * rmax, 0 + 1 * rmax)
        Y = ran(0 + 1 * rmax, height - 1 * rmax)
        R = ran(rmin, rmax)
        newcircle_1 = [X, Y, R]
        newcircle_2 = [X + width, Y, R]  # second circle at the opposing end
        check_1 = check_a_circle(newcircle_1, centlist, tol)
        check_2 = check_a_circle(newcircle_2, centlist, tol)

        if check_1 == False and check_2 == False:
            # never intersect with other circles
            centlist = centlist + [newcircle_1] + [newcircle_2]
            vff = vff + \
                add_vff(newcircle_1[2], width, height) + \
                add_vff(newcircle_2[2], width, height)

    elif ratio_heights <= prob < ratio_widths:
        X = ran(0 + 1 * rmax, width - 1 * rmax)
        Y = ran(0 - 1 * rmax, 0 + 1 * rmax)
        R = ran(rmin, rmax)
        newcircle_3 = [X, Y, R]
        newcircle_4 = [X, Y + height, R]  # second circle at the opposing end
        check_3 = check_a_circle(newcircle_3, centlist, tol)
        check_4 = check_a_circle(newcircle_4, centlist, tol)

        if check_3 == False and check_4 == False:
            # never intersect with other circles
            centlist = centlist + [newcircle_3] + [newcircle_4]
            vff = vff + add_vff(newcircle_3[2], width, height) + \
                add_vff(newcircle_4[2], width, height)

    else:
        X = ran(0 - 1 * rmax, 0 + 1 * rmax)
        Y = ran(0 - 1 * rmax, 0 + 1 * rmax)
        R = ran(rmin, rmax)
        newcircle_5 = [X, Y, R]
        newcircle_6 = [X + width, Y, R]
        newcircle_7 = [X, Y + height, R]
        newcircle_8 = [X + width, Y + height, R]

        check_5 = check_a_circle(newcircle_5, centlist, tol)
        check_6 = check_a_circle(newcircle_6, centlist, tol)
        check_7 = check_a_circle(newcircle_7, centlist, tol)
        check_8 = check_a_circle(newcircle_8, centlist, tol)

        if check_5 == False and check_6 == False and check_7 == False and check_8 == False:
            centlist = centlist + [newcircle_5] + \
                [newcircle_6] + [newcircle_7] + [newcircle_8]

            vff = vff + add_vff(newcircle_5[2], width, height) + add_vff(newcircle_6[2], width, height) + add_vff(
                newcircle_7[2], width, height) + add_vff(newcircle_8[2], width, height)

        print(vff, ",", prob)  # FVF for each step
modelprint(width, height, 0.8, vff, centlist)
centlist, vff = fillcircle(centlist, vff)
modelprint(width, height, 0.8, vff, centlist)
