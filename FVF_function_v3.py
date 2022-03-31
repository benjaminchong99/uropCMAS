from astropy.stats import RipleysKEstimator
from scipy import spatial
from matplotlib.patches import Circle, Rectangle
import numpy as np
import time

from random import *
from math import *
from matplotlib import pyplot as plt

from sklearn.neighbors import KDTree


def intersection(r1, cent1, r2, cent2, tol):
    """inside check_a_circle function to determine whether 2 circles intersect"""
    distance = ((cent1[0] - cent2[0]) ** 2 + (cent1[1] -
                                              cent2[1]) ** 2) ** 0.5  # pythagoras thm
    min_distance = r1 + r2 + tol  # minimum distance between two centres of circle

    if distance >= min_distance:
        overlap = False
    else:
        overlap = True
    return overlap


def check_a_circle(cent1, centlist, tol):
    """to determine if a circle intersects with any circles in the list"""
    if len(centlist) == 0:
        return False
    r1 = cent1[2]
    for cent2 in centlist:
        r2 = cent2[2]
        check_overlap = intersection(r1, cent1, r2, cent2, tol)
        if check_overlap == True:
            return check_overlap

    return False


def add_vff(radius, width, height):
    """shorten the adding of vff per circle"""
    return ((pi * radius ** 2) / (width*height))


def model(width, height, lw, vff, centlist, dictionary):
    """ plot function, generate model"""

    print('FVF:', vff, '\n')   # The current FVF

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rectan = Rectangle((0, 0), width, height, linewidth=lw, edgecolor='red', facecolor='w',
                       linestyle='solid', alpha=1.0)
    ax.add_patch(rectan)

    for element in centlist:
        X = element[0]
        Y = element[1]
        R = element[2]
        draw_circle = Circle(xy=(X, Y), radius=R, color='black', alpha=1.0)
        ax.add_patch(draw_circle)
        plt.axis('scaled')
        # change limits of x or y axis so that equal increemets of x and y have the same length
        plt.axis('equal')

    plt.show()


def fillcircle(centlist, vff, dictionary, rmax):
    """fill possible additional circles if possible"""
    og = len(centlist)
    print("start", og)
    i = 0
    while i < len(centlist):
        # tilt one degree
        if centlist[i][1] > 0.8-rmax or centlist[i][1] < 0+rmax or centlist[i][0] > 0.8-rmax or centlist[i][0] < 0+rmax:
            pass
        else:
            radius = 2*centlist[i][2] + tol
            angle = 0
            while angle < 360:
                x_diff = radius * cos(angle/180*pi)
                y_diff = radius * sin(angle/180*pi)
                imaginary_cent = [centlist[i][0] + x_diff,
                                  centlist[i][1] + y_diff, centlist[i][2]]
                # check imaginary_cent

                if imaginary_cent[1] > 0.8-rmax or imaginary_cent[1] < 0+rmax or imaginary_cent[0] > 0.8-rmax or imaginary_cent[0] < 0+rmax:
                    pass
                else:
                    indicator = check_a_circle(imaginary_cent, centlist, tol)
                    if indicator == False:
                        # never overlap
                        centlist = centlist + [imaginary_cent]
                        dictionary[len(dictionary.items())] = imaginary_cent
                        vff = vff + add_vff(imaginary_cent[2], width, height)
                        print("increase")
                angle += 1  # go all 360 degrees
        print('done checking ', i, "/", len(centlist), og)
        if i == round(0.8*og):
            # marker, ratio to be determined
            if len(centlist)-og <= 2:
                i = len(centlist)  # stop and shake again to save time
        i += 1
    return centlist, vff


def random_movement(centlist, tol, vff, rmax, dictionary):
    """
    automated random movement of circles in centlist
    Supposedly generated until the vff hit
    """

    # KDTree
    arr = np.array(centlist)
    kdt = KDTree(arr, leaf_size=30, metric="euclidean")
    # ind there to seperate out distanec and indexes
    dist, ind = kdt.query(arr, k=7)  # 7 columns including original
    avail_range = dist[:, 1:]  # distance away from the closest neighbour
    avail_range_ls = np.ndarray.tolist(avail_range)
    avail_idx = ind[:, 1:]  # index
    avail_idx_ls = np.ndarray.tolist(avail_idx)
    print('before ', len(centlist), vff, ' before')
    # print(dictionary)
    # shake it thee times before finding space to insert circles
    for n in range(5):
        for i in range(len(centlist)):
            if centlist[i][1] > 0.8-centlist[i][2] or centlist[i][1] < 0+centlist[i][2] or centlist[i][0] > 0.8-centlist[i][2] or centlist[i][0] < 0+centlist[i][2]:
                pass
            else:
                # x coord of nn - xcoord of current
                ref_idx = avail_idx_ls[i]  # find the list of neighbours index
                neighcount = 0
                # print(str(dictionary[ref_idx]))
                while neighcount <= 5:
                    # if dictionary[ref_idx[neighcount]][1] > 0.8-centlist[i][2] or dictionary[ref_idx[neighcount]][1] < 0+centlist[i][2] or dictionary[ref_idx[neighcount]][0] > 0.8-centlist[i][2] or dictionary[ref_idx[neighcount]][0] < 0+centlist[i][2]:
                    #     neighcount += 1
                    # else:
                    # print(ref_idx)
                    # print(dictionary[ref_idx])
                    avail_x = abs(
                        dictionary[ref_idx[neighcount]][0] - centlist[i][0])
                    avail_y = abs(
                        dictionary[ref_idx[neighcount]][1] - centlist[i][1])
                    temp_x = centlist[i][0] + \
                        round(uniform(-avail_x, avail_x), 5)
                    temp_y = centlist[i][1] + \
                        round(uniform(-avail_y, avail_y), 5)
                    tempcircle = [temp_x, temp_y, centlist[i][2]]
                    if ((rmax < tempcircle[0] < 0.8-rmax) and (rmax < tempcircle[1] < 0.8-rmax)):
                        # check intercepting; revert back to normal checking
                        comparelist = centlist[:]
                        # comparelist.append(tempcircle)
                        comparelist.remove(comparelist[i])
                        check_overlap = check_a_circle(
                            tempcircle, comparelist, tol)
                        if check_overlap == False:
                            centlist[i] = tempcircle
                            dictionary[i] = tempcircle
                            neighcount = 10
                        else:
                            neighcount += 1
                    else:
                        neighcount += 1
    centlist, vff = fillcircle(centlist, vff, dictionary, rmax)
    print(len(centlist), vff)
    # model(width, height, 0.8, vff, centlist)
    # print_checkpoint = input("returning centlist vff")
    return centlist, vff


def maincheck(X, Y, R, centlist, vff):
    """MOST PARTS OF THE ALGORITHM HAPPENS HERE, SHOULD NOT HAVE MUCH CHANGES"""
    if (rmax <= X <= (width-rmax)) and (rmax <= Y <= (height-rmax)):
        newcircle = [X, Y, R]
        check = check_a_circle(newcircle, centlist, tol)
        if check == False:
            centlist = centlist + [newcircle]
            dictionary[len(dictionary.items())] = newcircle
            # center area
            vff = vff + add_vff(newcircle[2], width, height)

    # left height
    elif (-rmax <= X <= rmax) and (rmax <= Y <= (height-rmax)):
        # left side and right side
        newcircle_1 = [X, Y, R]
        newcircle_2 = [X + width, Y, R]  # second circle at the opposing end
        check_1 = check_a_circle(newcircle_1, centlist, tol)
        check_2 = check_a_circle(newcircle_2, centlist, tol)

        if check_1 == False and check_2 == False:
            # never intersect with other circles
            centlist = centlist + [newcircle_1] + [newcircle_2]
            dictionary[len(dictionary.items())] = newcircle_1
            dictionary[len(dictionary.items())] = newcircle_2

            # left side and right side
            vff = vff + add_vff(newcircle_1[2], width, height)

    # right height
    elif ((width-rmax) <= X <= (width+rmax)) and (rmax <= Y <= (height-rmax)):
        # left side and right side
        newcircle_1 = [X, Y, R]
        newcircle_2 = [X - width, Y, R]  # second circle at the opposing end
        check_1 = check_a_circle(newcircle_1, centlist, tol)
        check_2 = check_a_circle(newcircle_2, centlist, tol)

        if check_1 == False and check_2 == False:
            # never intersect with other circles
            centlist = centlist + [newcircle_1] + [newcircle_2]
            dictionary[len(dictionary.items())] = newcircle_1
            dictionary[len(dictionary.items())] = newcircle_2
            # left side and right side
            vff = vff + add_vff(newcircle_1[2], width, height)

    # bottom width
    elif (rmax <= X <= (width-rmax) and (-rmax <= Y <= rmax)):
        newcircle_3 = [X, Y, R]
        newcircle_4 = [X, Y + height, R]  # second circle at the opposing end
        check_3 = check_a_circle(newcircle_3, centlist, tol)
        check_4 = check_a_circle(newcircle_4, centlist, tol)

        if check_3 == False and check_4 == False:
            # never intersect with other circles
            centlist = centlist + [newcircle_3] + [newcircle_4]
            dictionary[len(dictionary.items())] = newcircle_3
            dictionary[len(dictionary.items())] = newcircle_4

            # topside and bottomside
            vff = vff + add_vff(newcircle_3[2], width, height)
    # top width
    elif ((rmax <= X <= (width-rmax)) and (height-rmax) <= Y <= (height+rmax)):
        newcircle_3 = [X, Y, R]
        newcircle_4 = [X, Y - height, R]  # second circle at the opposing end
        check_3 = check_a_circle(newcircle_3, centlist, tol)
        check_4 = check_a_circle(newcircle_4, centlist, tol)

        if check_3 == False and check_4 == False:
            # never intersect with other circles
            centlist = centlist + [newcircle_3] + [newcircle_4]
            dictionary[len(dictionary.items())] = newcircle_3
            dictionary[len(dictionary.items())] = newcircle_4
            # topside and bottomside
            vff = vff + add_vff(newcircle_3[2], width, height)

    # 4 corners
    # bottom left corner
    elif (-rmax <= X <= rmax) and (-rmax <= Y <= rmax):
        newcircle_5 = [X, Y, R]  # bottom left corner
        newcircle_6 = [X + width, Y, R]  # bottom right corner
        newcircle_7 = [X, Y + height, R]  # top left corner
        newcircle_8 = [X + width, Y + height, R]  # top right corner

        check_5 = check_a_circle(newcircle_5, centlist, tol)
        check_6 = check_a_circle(newcircle_6, centlist, tol)
        check_7 = check_a_circle(newcircle_7, centlist, tol)
        check_8 = check_a_circle(newcircle_8, centlist, tol)

        if check_5 == False and check_6 == False and check_7 == False and check_8 == False:
            centlist = centlist + [newcircle_5] + \
                [newcircle_6] + [newcircle_7] + [newcircle_8]
            dictionary[len(dictionary.items())] = newcircle_5
            dictionary[len(dictionary.items())] = newcircle_6
            dictionary[len(dictionary.items())] = newcircle_7
            dictionary[len(dictionary.items())] = newcircle_8
            # 4 corners
            vff = vff + add_vff(newcircle_5[2], width, height)

    # bottom right corner
    elif ((width-rmax) <= X <= (width+rmax)) and (-rmax <= Y <= rmax):
        newcircle_5 = [X, Y, R]  # bottom right corner
        newcircle_6 = [X - width, Y, R]  # bottom left corner
        newcircle_7 = [X, Y + height, R]  # top right corner
        newcircle_8 = [X - width, Y + height, R]  # top left corner

        check_5 = check_a_circle(newcircle_5, centlist, tol)
        check_6 = check_a_circle(newcircle_6, centlist, tol)
        check_7 = check_a_circle(newcircle_7, centlist, tol)
        check_8 = check_a_circle(newcircle_8, centlist, tol)

        if check_5 == False and check_6 == False and check_7 == False and check_8 == False:
            centlist = centlist + [newcircle_5] + \
                [newcircle_6] + [newcircle_7] + [newcircle_8]
            dictionary[len(dictionary.items())] = newcircle_5
            dictionary[len(dictionary.items())] = newcircle_6
            dictionary[len(dictionary.items())] = newcircle_7
            dictionary[len(dictionary.items())] = newcircle_8
            # 4 corners
            vff = vff + add_vff(newcircle_5[2], width, height)

    # top left corner
    elif (-rmax <= X <= rmax) and ((height-rmax) <= Y <= (height+rmax)):
        newcircle_5 = [X, Y, R]  # top left corner
        newcircle_6 = [X + width, Y, R]  # top right corner
        newcircle_7 = [X, Y - height, R]  # bottom left corner
        newcircle_8 = [X + width, Y - height, R]  # bottom right corner

        check_5 = check_a_circle(newcircle_5, centlist, tol)
        check_6 = check_a_circle(newcircle_6, centlist, tol)
        check_7 = check_a_circle(newcircle_7, centlist, tol)
        check_8 = check_a_circle(newcircle_8, centlist, tol)

        if check_5 == False and check_6 == False and check_7 == False and check_8 == False:
            centlist = centlist + [newcircle_5] + \
                [newcircle_6] + [newcircle_7] + [newcircle_8]
            dictionary[len(dictionary.items())] = newcircle_5
            dictionary[len(dictionary.items())] = newcircle_6
            dictionary[len(dictionary.items())] = newcircle_7
            dictionary[len(dictionary.items())] = newcircle_8
            # 4 corners
            vff = vff + add_vff(newcircle_5[2], width, height)

    # top right corner
    elif ((width-rmax) <= X <= (width+rmax)) and ((height-rmax) <= Y <= (height+rmax)):
        newcircle_5 = [X, Y, R]  # top right corner
        newcircle_6 = [X - width, Y, R]  # top left corner
        newcircle_7 = [X, Y - height, R]  # bottom right corner
        newcircle_8 = [X - width, Y - height, R]  # bottom left corner

        check_5 = check_a_circle(newcircle_5, centlist, tol)
        check_6 = check_a_circle(newcircle_6, centlist, tol)
        check_7 = check_a_circle(newcircle_7, centlist, tol)
        check_8 = check_a_circle(newcircle_8, centlist, tol)

        if check_5 == False and check_6 == False and check_7 == False and check_8 == False:

            centlist = centlist + [newcircle_5] + \
                [newcircle_6] + [newcircle_7] + [newcircle_8]
            dictionary[len(dictionary.items())] = newcircle_5
            dictionary[len(dictionary.items())] = newcircle_6
            dictionary[len(dictionary.items())] = newcircle_7
            dictionary[len(dictionary.items())] = newcircle_8
            # 4 corners
            vff = vff + add_vff(newcircle_5[2], width, height)

    print(vff)  # FVF for each step
    return centlist, vff

################################################################################
    """START OF COMPUTATION"""


################################################################################
start = time.time()


# Variable parameters
rmax = 0.016  # the maxmum radius of the fibres
rmin = 0.016  # the minmum radius of the fibres
width = rmax*50  # 0.8   # the length of the RVE
height = rmax*50  # 0.8 # the width of the RVE
tol = 0.001  # the minmum distance of two circles (except for the radius)
target_vf = 0.65  # target FVF to hit
# preset_vf = 0.45    # the FVF in the RVE; to hit initially
centlist = []  # initial empty space
vff = 0
# vff = total area of circle / area of box
# with every circle added to the model
# cap at 0.55 if don't want to wait

dictionary = {
}
# put into the dictionary accordingly,
# shift based on corresponding indexing in the list
# not gg for time complexity but for functionality , making use of dictionary for indexing


# Algorithm of generating random distributing fibres
print("begin")
print("vff: ", vff, "\n Vf: ", target_vf)
counting = 0
# while vff < preset_vf:
while counting < 50:
    # change to stright compare x and y
    X = uniform(0 - 1 * rmax, width + 1 * rmax)
    Y = uniform(0 - 1 * rmax, height + 1 * rmax)
    R = uniform(rmin, rmax)
    centlist, vff = maincheck(X, Y, R, centlist, vff)
    counting += 1
centlist, vff = fillcircle(centlist, vff, dictionary, rmax)
model(width, height, 0.8, vff, centlist, dictionary)
print(centlist)
while vff < target_vf:
    centlist, vff = random_movement(centlist, tol, vff, rmax, dictionary)
# namespaced(centlist, tol, vff, width, height)
print('check space!')
print(centlist)
print(vff)
model(width, height, 0.8, vff, centlist, dictionary)

# 470 for 0.55
# 513 for 0.6
# records the time you take to view the model photo also
end = time.time()
print("Time: ", end - start)
################################################################################
""" END OF COMPUTATION """
################################################################################

"""# KDTree
# One searchig methods of KNN
arr = np.array(centlist)
print(arr)
kdt = KDTree(arr, leaf_size=30, metric='euclidean')

dist, ind = kdt.query(arr, k=2)
print("kdtree: ", dist, ind)  # distance, then index of the nearest circle


# the 1st nearest neighbor distance, B0 reprents the point itself
B1 = dist[:, 1]
print(B1, '\n')
ls = np.ndarray.tolist(B1-0.016)
print(len(ls), len(centlist))"""


"""
Pseudocode
# cases if circle is at the outskirts

# zone left side right side
if left side:
    second_circle = first_circle [x+width], first_circle [y]
    
if right side:
    second_circle = first_circle [x-width], first_circle [y]
# zone top bottom
if top:
    second_circle = first_circle [x], first_circle [y-height]
if bottom:
    second_circle = first_circle [x-width], first_circle [y+height]
check circle

# zone four corners




"""

#
###############################################################################
#  2 - Nearest neighbor orientation                                           #
###############################################################################
#
#
O = []
A = centlist
# Variable parameters
start_l = 0
end_l = 2 * pi
n_step = 10
interval_step = (end_l - start_l) / n_step
# print(interval_step)
#
M = np.zeros(n_step)
# for each fibre a, search for the nearest neighbor fibre b
for i in range(len(A)):
    a = A[i]  # list a   =====================================================#
    B = A[:]
    B.remove(a)
    array_a = np.array(a)
    array_B = np.array(B)
    print("ArrayA: ", array_a)
    print("ArrayB: ", array_B)
    array_minP = array_B[spatial.KDTree(array_B).query(array_a)[1]]
    # list b   ========================================#
    b = list(array_minP)
    # print(b)
    minD = ((list(array_a)[0] - list(array_minP)[0]) ** 2.0 +
            (list(array_a)[1] - list(array_minP)[1]) ** 2.0) ** 0.5
    # print(minD)
    x_dis = abs(a[0] - b[0])  # cosine edge
    cos_alpha = x_dis / minD
    # print(cos_alpha)    # the cosine value
    alpha_o = acos(cos_alpha)  # min_alpha for each quadrant
    #
    # alpha for different quadrant using the alpha_o
    if b[0] > a[0] and b[1] < a[1]:
        alpha = alpha_o
    if b[0] < a[0] and b[1] < a[1]:
        alpha = pi - alpha_o
    if b[0] < a[0] and b[1] > a[1]:
        alpha = pi + alpha_o
    if b[0] > a[0] and b[1] > a[1]:
        alpha = 2 * pi - alpha_o
    # print(i, alpha)   # the angle of the two fibres
#=============================================================================#
    i = 0
    for k in np.arange(start_l, end_l, interval_step):
        if alpha >= k and alpha < k + interval_step:
            M[i] += 1
        i += 1
# print(M, '\n')   # counts of alpha in each interval of angles
#
O_bank = [0]
for i in range(len(M)):
    x = sum(M[: i + 1])
    # print(x)
    prob_o = x / len(A)
    O_bank.append(prob_o)
# cumulative probability of alpha with increasing of intervals
print('2: cumulative probability of alpha:\n', O_bank, '\n')
O.append(O_bank)
print(O)
#============================================================================#
#


def curv_plt(x_axis, y_array, l_color, l_label):
    l_mean = np.mean(y_array, axis=0)
    l_max = np.max(y_array, axis=0)
    l_min = np.min(y_array, axis=0)
    plt.plot(x_axis, l_mean, color=l_color, label=l_label, marker='.')
    plt.fill_between(x_axis, l_min, l_max, color=l_color, alpha=0.2)
    plt.legend(fontsize=16)


#
# plot the nearest neighbor orientation
plt.figure()
curv_plt(np.arange(0, n_step + 1, 1), O, 'r', 'Nearest neighbor orientation')
plt.xlabel('Orientation [degree]')
plt.ylabel('Cumulative Distribution Function')
plt.xlim(0, 10)
plt.ylim(0, 1.0)
plt.legend(loc='upper left')
plt.show()
#
###############################################################################
#  End of function                                                            #
###############################################################################


#
###############################################################################
#  3 - Ripley's K function                                                    #
###############################################################################
#
#
R = []
x = np.array(centlist)   # array of list of CentList
z = x[:, [0, 1]]         # extract the first two columns
# print(z, '\n')
#
Kest = RipleysKEstimator(area=0.8**2, x_max=0.8 + rmax,
                         y_max=0.8 + rmax, x_min=-rmax, y_min=- rmax)
# start from 0, end with 0.24, 100 floats, but how to determine the r value? 0.3 * x_max
h = np.linspace(0, 0.3, 30)
y_r = Kest(data=z, radii=h, mode='ripley')
R.append(y_r)
#=============================================================================#
#


def curv_plt(x_axis, y_array, l_color, l_label):
    l_mean = np.mean(y_array, axis=0)
    l_max = np.max(y_array, axis=0)
    l_min = np.min(y_array, axis=0)
    plt.plot(x_axis, l_mean, color=l_color, label=l_label, marker='.')
    plt.fill_between(x_axis, l_min, l_max, color=l_color, alpha=0.2)
    plt.legend(fontsize=16)


#
# plot the Poisson and Ripley's K function
plt.figure()
plt.plot(h/rmax, Kest.poisson(h), color='green',
         ls=':', label=r'$K_{pois}$')  # poisson curve
# print('3-1: y coordinate values on poisson curve:\n', Kest.poisson(h), '\n')
curv_plt(h/rmax, np.array(R), 'r', '$K_{ripley}$')  # ripley curve
# print('3-2: y coordinate values on Ripley curve:\n', Kest(data=z, radii=h, mode='ripley'), '\n')  # y coordinate value
plt.xlabel('h/R')
plt.ylabel('K(h)')
plt.legend()
plt.show()
#
###############################################################################
#  End of function                                                            #
###############################################################################


# for randomness analysis
R_mean = np.average(R, axis=0)   # average value of D1 for 20 calculations
print(list(R_mean))
