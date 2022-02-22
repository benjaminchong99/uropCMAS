from matplotlib.patches import Circle, Rectangle
import numpy as np

from random import *
from math import *
from matplotlib import pyplot as plt


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


def model(width, height, lw, vff, centlist):
    """ plot function, generate model"""

    print(len(centlist))
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

    # centlist = sorted(centlist)
    # model(width, height, 0.8, vff, centlist)
    # namespaced(centlist, 0.0005, vff)

    plt.show()


def fillcircle(centlist, vff):
    """fill possible additional circles if possible"""
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
#    movespaced(centlist, vff, tol, width, height)


def random_movement(centlist, tol, vff):
    """automated random movement of circles in centlist"""
    print('before ', len(centlist), vff, ' before')
    for i in range(len(centlist)):
        if centlist[i][1] > 0.8-centlist[i][2] or centlist[i][1] < 0+centlist[i][2] or centlist[i][0] > 0.8-centlist[i][2] or centlist[i][0] < 0+centlist[i][2]:
            pass
        else:

            temp_x = centlist[i][0] + round(uniform(-0.01, 0.01), 4)
            temp_y = centlist[i][1] + round(uniform(-0.01, 0.01), 4)
            tempcircle = [temp_x, temp_y, centlist[i][2]]
            comparelist = centlist[:]
            comparelist.remove(comparelist[i])
            result = check_a_circle(tempcircle, comparelist, tol)
            if result == False:
                centlist[i] = tempcircle
            print(f"Done: {i}/{len(centlist)}")
    centlist, vff = fillcircle(centlist, vff)
    print(len(centlist), vff)
    # model(width, height, 0.8, vff, centlist)
    # useless_checkpoint = input("returning centlist vff")
    return centlist, vff


"""manual moving of individual selected circles"""
# def movespaced(centlist, vff, tol, width, height):

#     answer = input(
#         f"type in the numbered circle (current max = {len(centlist)}) : ")
#     while not answer.isnumeric():
#         answer = input(
#             "invalid input, try again. Type in the numbered circle: ")

#     if int(answer) == len(centlist):
#         # fill circle
#         centlist, vff = fillcircle(centlist, vff)
#         model(width, height, 0.8, vff, centlist)
#         namespaced(centlist, tol, vff, width, height)

#     elif int(answer) > len(centlist):
#         print("out of range, end")
#         model(width, height, 0.8, vff, centlist)
#         plt.show()

#     else:
#         index = int(answer)
#         # store og position first
#         storage = centlist[index]
#         movement = [width, height]
#         while abs(movement[0]) >= width or abs(movement[1]) >= height:
#             answer2 = input(
#                 f"Chosen circle {answer}, coordinates: {centlist[index][0], centlist[index][1]}. Enter valid movement in x and y direction: ")
#             movement = answer2.split(" ")

#             if len(movement) < 2:
#                 movement = [width, height]
#             else:
#                 for i in range(len(movement)):
#                     movement[i] = float(movement[i])
#                 # need to add distance to the point
#                 centlist[index] = [centlist[index][0] +
#                                    float(movement[0]), centlist[index][1] + float(movement[1]), centlist[index][2]]
#                 # check valid movement
#                 if 0 < centlist[index][0] < width and 0 < centlist[index][1] < height:
#                     pass
#                 else:
#                     centlist[index] = storage
#                     movement = [width, height]
#         print(centlist[index])
#         print(movement)
#         # check the validity of the movement before you move forward with moving the point
#         comparelist = centlist[:]
#         comparelist.remove(centlist[index])
#         indicator = check_a_circle(centlist[index], comparelist, tol)
#         if indicator == False:
#             print('enter when no overlap')
#             nonsense = input('enter to confirm')
#             centlist, vff = single_fillcircle(centlist[index], centlist, vff)

#         # else reject and ask to move the point again
#         else:
#             print("invalid movement, try again")
#             centlist[index] = storage
#             nonsense = input('enter to confirm')
#         model(width, height, 0.8, vff, centlist)
#         namespaced(centlist, tol, vff, width, height)

#         plt.show()

# def single_fillcircle(circle, centlist, vff):
#     """
#     fill possible circles to an individual selected circle
#     under movespaced function only
#     y then x
#     """
#     if circle[1] > 0.8-circle[2] or circle[1] < 0+circle[2] or circle[0] > 0.8-circle[2] or circle[0] < 0+circle[2]:
#         pass
#     else:
#         radius = 2*circle[2] + tol
#         angle = 0
#         while angle < 360:
#             x_diff = radius * cos(angle/180*pi)
#             y_diff = radius * sin(angle/180*pi)
#             imagcent = [circle[0] + x_diff,
#                         circle[1] + y_diff, circle[2]]
#             # check imagcent
#             if imagcent[1] > 0.8 or imagcent[1] < 0 or imagcent[0] > 0.8 or imagcent[0] < 0:
#                 pass
#             else:
#                 indicator = check_a_circle(imagcent, centlist, tol)
#                 if indicator == False:
#                     # never overlap
#                     centlist = centlist + [imagcent]
#                     vff = vff + add_vff(imagcent[2], width, height)
#                     print("increase")
#             angle += 1
#     return centlist, vff
""""""

"""START OF COMPUTATION"""
# Variable parameters
rmax = 0.016  # the maxmum radius of the fibres
rmin = 0.016  # the minmum radius of the fibres
width = rmax*50  # 0.8   # the length of the RVE
height = rmax*50  # 0.8 # the width of the RVE
tol = 0.001  # the minmum distance of two circles (except for the radius)
vf = 0.45    # the FVF in the RVE; to hit initially
centlist = []  # initial empty space
vff = 0
# vff = total area of circle / area of box
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
'''^change to use x and y and then compare straight away'''
"""vff calculation only add one after every case"""

center_list = []
heightsLeft_list = []
heightsRight_list = []
widthsLeft_list = []
widthsRight_list = []
corners_list = []

allCircles = {
    "centre": [],
    "leftside": [],
    "rightside": [],
    "topside": [],
    "bottomside": [],
    "bottomleft": [],
    "bottomright": [],
    "topleft": [],
    "topright": []
}


while vff < vf:
    prob = random()
    if prob < ratio_centre:
        X = uniform(0 + 1 * rmax, width - 1 * rmax)
        Y = uniform(0 + 1 * rmax, height - 1 * rmax)
        R = uniform(rmin, rmax)
        newcircle = [X, Y, R]
        check = check_a_circle(newcircle, centlist, tol)
        if check == False:
            centlist = centlist + [newcircle]
            # center area
            allCircles["centre"].append(newcircle)
            vff = vff + add_vff(newcircle[2], width, height)

    elif ratio_centre <= prob < ratio_heights:
        # left side and right side
        X = uniform(0 - 1 * rmax, 0 + 1 * rmax)
        Y = uniform(0 + 1 * rmax, height - 1 * rmax)
        R = uniform(rmin, rmax)
        newcircle_1 = [X, Y, R]
        newcircle_2 = [X + width, Y, R]  # second circle at the opposing end
        check_1 = check_a_circle(newcircle_1, centlist, tol)
        check_2 = check_a_circle(newcircle_2, centlist, tol)

        if check_1 == False and check_2 == False:
            # never intersect with other circles
            centlist = centlist + [newcircle_1] + [newcircle_2]
            # left side and right side
            allCircles["leftside"].append(newcircle_1)
            allCircles["rightside"].append(newcircle_2)
            vff = vff + \
                add_vff(newcircle_1[2], width, height)

    elif ratio_heights <= prob < ratio_widths:
        # top and bottom
        X = uniform(0 + 1 * rmax, width - 1 * rmax)
        Y = uniform(0 - 1 * rmax, 0 + 1 * rmax)
        R = uniform(rmin, rmax)
        newcircle_3 = [X, Y, R]
        newcircle_4 = [X, Y + height, R]  # second circle at the opposing end
        check_3 = check_a_circle(newcircle_3, centlist, tol)
        check_4 = check_a_circle(newcircle_4, centlist, tol)

        if check_3 == False and check_4 == False:
            # never intersect with other circles
            centlist = centlist + [newcircle_3] + [newcircle_4]
            # topside and bottomside
            allCircles["bottomside"].append(newcircle_3)
            allCircles["topside"].append(newcircle_4)
            vff = vff + add_vff(newcircle_3[2], width, height)

    else:
        # 4 corners
        X = uniform(0 - 1 * rmax, 0 + 1 * rmax)
        Y = uniform(0 - 1 * rmax, 0 + 1 * rmax)
        R = uniform(rmin, rmax)
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
            # 4 corners
            allCircles["bottomleft"].append(newcircle_5)
            allCircles["bottomright"].append(newcircle_6)
            allCircles["topleft"].append(newcircle_7)
            allCircles["topright"].append(newcircle_8)
            vff = vff + add_vff(newcircle_5[2], width, height)

        print(vff, ",", prob)  # FVF for each step
for keys in allCircles:
    print(len(allCircles[keys]))

print("ALL CIRCLES :", len(allCircles))
useless = input("ENDS here")
centlist, vff = fillcircle(centlist, vff)
model(width, height, 0.8, vff, centlist)
# namespaced(centlist, tol, vff, width, height)
# model(width, height, 0.8, vff, centlist)


while vff < 0.6:
    centlist, vff = random_movement(centlist, tol, vff)
    # useless = input('shaken once, enter for result')
namespaced(centlist, tol, vff, width, height)
print('check space!')
print(centlist)
model(width, height, 0.8, vff, centlist)
