from matplotlib.patches import Circle, Rectangle
import numpy as np

from random import *
from math import *
from matplotlib import pyplot as plt

from sklearn.neighbors import KDTree


def intersection(r1, cent1, r2, cent2, tol):
    """inside check_a_circle function to determine whether 2 circles intersect"""
    distance = ((cent1[0] - cent2[0]) ** 2 + (cent1[1] - cent2[1])
                ** 2) ** 0.5  # pythagoras thm
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
        temp = intersection(r1, cent1, r2, cent2, tol)
        if temp == True:
            return temp
    return temp


def add_vff(radius, width, height):
    """shorten the adding of vff per circle"""
    return ((pi * radius ** 2) / (width*height))


def model(width, height, lw, vff, centlist, dictionary):
    """ plot function, generate model"""

    print('FVF:', vff, '\n')   # The final FVF

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rectan = Rectangle((0, 0), width, height, linewidth=lw, edgecolor='red', facecolor='w',
                       linestyle='solid', alpha=1.0)
    ax.add_patch(rectan)

    # for values in dictionary.values():
    #     for v in values:
    #         X = v[0]
    #         Y = v[1]
    #         R = v[2]
    #         draw_circle = Circle(xy=(X, Y), radius=R, color='black', alpha=1.0)
    #         ax.add_patch(draw_circle)
    #         plt.axis('scaled')
    #         # change limits of x or y axis so that equal increemets of x and y have the same length
    #         plt.axis('equal')
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


def fillcircle(centlist, vff):
    """fill possible additional circles if possible"""
    og = len(centlist)
    print("start", og)
    i = 0
    while i < len(centlist):
        # tilt one degree
        if centlist[i][1] > 0.8 or centlist[i][1] < 0 or centlist[i][0] > 0.8 or centlist[i][0] < 0:
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

                if imaginary_cent[1] > 0.8 or imaginary_cent[1] < 0 or imaginary_cent[0] > 0.8 or imaginary_cent[0] < 0:
                    pass
                else:
                    indicator = check_a_circle(imaginary_cent, centlist, tol)
                    if indicator == False:
                        # never overlap
                        centlist = centlist + [imaginary_cent]
                        vff = vff + add_vff(imaginary_cent[2], width, height)
                        print("increase")
                angle += 1
        print('done checking ', i, "/", len(centlist), og)
        if i == round(0.6*og):
            # marker, ratio to be determined
            if len(centlist)-og <= 2:
                i = len(centlist)  # stop and shake again to save time
        i += 1
    return centlist, vff


def namespaced(centlist, tol, vff, width, height):
    """plot function, name the circles"""

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rectan = Rectangle((0, 0), width, height, linewidth=0.8,
                       edgecolor='red', facecolor='w', linestyle='solid', alpha=1.0)
    ax.add_patch(rectan)
    for i in range(len(centlist)):
        X = centlist[i][0]
        Y = centlist[i][1]
        R = centlist[i][2]
        draw_circle = Circle(xy=(X, Y), radius=R, color='black', alpha=1.0)
        ax.add_patch(draw_circle)
        plt.plot(centlist[i][0], centlist[i][1], c='red', marker="o")
        plt.text(centlist[i][0], centlist[i][1],
                 str(i), color='blue', fontsize=12)
        plt.axis('scaled')
        # change limits of x or y axis so that equal increemets of x and y have the same length
        plt.axis('equal')
    print(f"vff: {vff}, circles: {len(centlist)}")
    plt.show()


def random_movement(centlist, tol, vff, rmax):
    """
    automated random movement of circles in centlist
    Supposedly generated until the vff hit
    """

    # KDTree
    arr = np.array(centlist)
    kdt = KDTree(arr, leaf_size=30, metric="euclidean")

    # ind there to seperate out distanec and indexes
    dist, ind = kdt.query(arr, k=2)
    avail_range = dist[:, 1]  # distance away from the closest neighbour
    avail_range_ls = np.ndarray.tolist(avail_range)
    print('before ', len(centlist), vff, ' before')

    for i in range(len(centlist)):
        if centlist[i][1] > 0.8-centlist[i][2] or centlist[i][1] < 0+centlist[i][2] or centlist[i][0] > 0.8-centlist[i][2] or centlist[i][0] < 0+centlist[i][2]:
            pass
        else:
            temp_x = centlist[i][0] + \
                round(uniform(-avail_range_ls[i] -
                              rmax, avail_range_ls[i]-rmax), 4)
            temp_y = centlist[i][1] + \
                round(uniform(-avail_range_ls[i] -
                              rmax, avail_range_ls[i]-rmax), 4)
            tempcircle = [temp_x, temp_y, centlist[i][2]]
            if (-rmax < temp_x < 0.8+rmax) and (-rmax < temp_y < 0.8+rmax):
                comparelist = centlist[:]
                comparelist.remove(comparelist[i])
                result = check_a_circle(tempcircle, comparelist, tol)
                if result == False:
                    centlist[i] = tempcircle
            print(f"Done: {i}/{len(centlist)}")
    centlist, vff = fillcircle(centlist, vff)
    print(len(centlist), vff)
    # model(width, height, 0.8, vff, centlist)
    # print_checkpoint = input("returning centlist vff")
    return centlist, vff


def maincheck(X, Y, R, centlist, vff):
    if (rmax <= X <= (width-rmax)) and (rmax <= Y <= (height-rmax)):
        newcircle = [X, Y, R]
        check = check_a_circle(newcircle, centlist, tol)
        if check == False:
            centlist = centlist + [newcircle]
            dictionary["cc"].append(newcircle)
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
            dictionary["lh"].append(newcircle_1)
            dictionary["rh"].append(newcircle_2)
            # left side and right side
            vff = vff + \
                add_vff(newcircle_1[2], width, height)

    # right height
    elif ((width-rmax) <= X <= (width+rmax)) and (rmax <= Y <= (height-rmax)):
        # left side and right side
        newcircle_1 = [X, Y, R]
        newcircle_2 = [X - width, Y, R]  # second circle at the opposing end
        check_1 = check_a_circle(newcircle_1, centlist, tol)
        check_2 = check_a_circle(newcircle_2, centlist, tol)

        if check_1 == False and check_2 == False:
            # never intersect with other circles
            dictionary["rh"].append(newcircle_1)
            dictionary["lh"].append(newcircle_2)
            centlist = centlist + [newcircle_1] + [newcircle_2]
            # left side and right side
            vff = vff + \
                add_vff(newcircle_1[2], width, height)

    # bottom width
    elif (rmax <= X <= (width-rmax) and (-rmax <= Y <= rmax)):
        newcircle_3 = [X, Y, R]
        newcircle_4 = [X, Y + height, R]  # second circle at the opposing end
        check_3 = check_a_circle(newcircle_3, centlist, tol)
        check_4 = check_a_circle(newcircle_4, centlist, tol)

        if check_3 == False and check_4 == False:
            # never intersect with other circles
            dictionary["bw"].append(newcircle_3)
            dictionary["tw"].append(newcircle_4)
            centlist = centlist + [newcircle_3] + [newcircle_4]
            # topside and bottomside
            vff = vff + add_vff(newcircle_3[2], width, height)
    # top width
    elif ((rmax <= X <= (width-rmax)) and (height-rmax) <= Y <= (height+rmax)):
        newcircle_3 = [X, Y, R]
        newcircle_4 = [X, Y - height, R]  # second circle at the opposing end
        check_3 = check_a_circle(newcircle_3, centlist, tol)
        check_4 = check_a_circle(newcircle_4, centlist, tol)

        if check_3 == False and check_4 == False:
            dictionary["tw"].append(newcircle_3)
            dictionary["bw"].append(newcircle_4)
            # never intersect with other circles
            centlist = centlist + [newcircle_3] + [newcircle_4]
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
            dictionary["ltc"].append(newcircle_7)
            dictionary["rtc"].append(newcircle_8)
            dictionary["lbc"].append(newcircle_5)
            dictionary["rbc"].append(newcircle_6)
            centlist = centlist + [newcircle_5] + \
                [newcircle_6] + [newcircle_7] + [newcircle_8]
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
            dictionary["ltc"].append(newcircle_8)
            dictionary["rtc"].append(newcircle_7)
            dictionary["lbc"].append(newcircle_6)
            dictionary["rbc"].append(newcircle_5)
            centlist = centlist + [newcircle_5] + \
                [newcircle_6] + [newcircle_7] + [newcircle_8]
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
            dictionary["ltc"].append(newcircle_5)
            dictionary["rtc"].append(newcircle_6)
            dictionary["lbc"].append(newcircle_7)
            dictionary["rbc"].append(newcircle_8)
            centlist = centlist + [newcircle_5] + \
                [newcircle_6] + [newcircle_7] + [newcircle_8]
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
            dictionary["ltc"].append(newcircle_6)
            dictionary["rtc"].append(newcircle_5)
            dictionary["lbc"].append(newcircle_8)
            dictionary["rbc"].append(newcircle_7)

            centlist = centlist + [newcircle_5] + \
                [newcircle_6] + [newcircle_7] + [newcircle_8]
            # 4 corners
            vff = vff + add_vff(newcircle_5[2], width, height)

    print(vff)  # FVF for each step
    return centlist, vff


def specialmovements(circle, radius, tol, vff):
    """I need to store the locations properly first such that I can search it later"""
    """FUNCTION NOT IN USED YET BECAUSE OF ^"""
    # for movement in part 2 3 and 4
    # if y > 0.8-radius or y < 0 + rad or x > 0.8-radius or x < 0+rad

    if (-radius <= circle[0] <= radius and radius < circle[1] < 0.8-radius) or (circle[0] > 0.8-radius and radius < circle[1] < 0.8-radius):
        # part 2
        # if it is part 2 (vertical sides), part 3 horizontal side, part 4 corners
        pass

    if (radius < circle[0] < 0.8-radius and circle[1] <= radius) or (radius < circle[0] < 0.8-radius and circle[1] >= 0.8-radius):
        # part 3
        pass

    if (circle[0] <= radius and circle[1] <= radius) or (circle[0] <= radius and circle[1] >= 0.8-radius) or (circle[0] >= 0.8-radius and circle[1] <= radius) or (circle[0] <= 0.8-radius and circle[1] >= 0.8-radius):
        # if you shift one corner, you shift adjacent one as well
        # check all to see if got overlap or not
        # add inside
        pass

################################################################################
        """START OF COMPUTATION"""
################################################################################


# Variable parameters
rmax = 0.016  # the maxmum radius of the fibres
rmin = 0.016  # the minmum radius of the fibres
width = rmax*50  # 0.8   # the length of the RVE
height = rmax*50  # 0.8 # the width of the RVE
tol = 0.001  # the minmum distance of two circles (except for the radius)
target_vf = 0.55  # target FVF to hit
preset_vf = 0.45    # the FVF in the RVE; to hit initially
centlist = []  # initial empty space
vff = 0
# vff = total area of circle / area of box
# with every circle added to the model
# cap at 0.55 if don't want to wait

dictionary = {
    "cc": [],
    "lh": [],
    "rh": [],
    "tw": [],
    "bw": [],
    "ltc": [],
    "rtc": [],
    "lbc": [],
    "rbc": []
}
# put into the dictionary accordingly,
# shift based on corresponding indexing in the list
# not gg for time complexity but for functionality , making use of dictionary for indexing


# Algorithm of generating random distributing fibres
print("begin")
print("vff: ", vff, "\n Vf: ", target_vf)

while vff < preset_vf:
    # change to stright compare x and y
    X = uniform(0 - 1 * rmax, width + 1 * rmax)
    Y = uniform(0 - 1 * rmax, height + 1 * rmax)
    R = uniform(rmin, rmax)
    centlist, vff = maincheck(X, Y, R, centlist, vff)
# print(dictionary)
centlist, vff = fillcircle(centlist, vff)
model(width, height, 0.8, vff, centlist, dictionary)
print(centlist)

while vff < target_vf:
    centlist, vff = random_movement(centlist, tol, vff, rmax)
namespaced(centlist, tol, vff, width, height)
print('check space!')
print(centlist)
model(width, height, 0.8, vff, centlist, dictionary)

# 470 for 0.55
# 513 for 0.6
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
