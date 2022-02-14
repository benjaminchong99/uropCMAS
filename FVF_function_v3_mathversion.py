from matplotlib.patches import Circle, Rectangle
import numpy as np

from random import *
from math import *
from matplotlib import pyplot as plt


def intersection(r1, cent1, r2, cent2, tol):
    # inside check_a_circle function to determine whether 2 circles intersect
    distance = ((cent1[0] - cent2[0]) ** 2 + (cent1[1] - cent2[1])
                ** 2) ** 0.5  # pythagoras thm
    min_distance = r1 + r2 + tol  # minimum distance between two centres of circle

    if distance >= min_distance:
        overlap = False
    else:
        overlap = True
    return overlap


def check_a_circle(cent1, centlist, tol):
    # to determine if a circle intersects with any circles in the list
    r1 = cent1[2]
    for cent2 in centlist:
        r2 = cent2[2]
        temp = intersection(r1, cent1, r2, cent2, tol)
        if temp == True:
            return temp
    return temp


'''
def ran(x1, x2):
    # randomise
    # replaced with uniform() (random.uniform)
    # used to generate random centlist[i]s in the chosen region
    return ((x2 - x1) * np.random.rand() + x1)
'''


def add_vff(radius, width, height):
    # shorten the adding of vff per circle
    return ((pi * radius ** 2) / (width*height))


def modelprint(width, height, lw, vff, centlist):
    # plot function
    # generate model
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
    # modelprint(width, height, 0.8, vff, centlist)
    # namespaced(centlist, 0.0005, vff)

    plt.show()


def fillcircle(centlist, vff):
    # fill possible additional circles if possible
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
    # plot function
    # name the circles
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


def movespaced(centlist, vff, tol, width, height):
    # manual moving of individual selected circles
    answer = input(
        f"type in the numbered circle (current max = {len(centlist)}) : ")
    while not answer.isnumeric():
        answer = input(
            "invalid input, try again. Type in the numbered circle: ")

    if int(answer) == len(centlist):
        # fill circle
        centlist, vff = fillcircle(centlist, vff)
        modelprint(width, height, 0.8, vff, centlist)
        namespaced(centlist, tol, vff, width, height)

    elif int(answer) > len(centlist):
        print("out of range, end")
        modelprint(width, height, 0.8, vff, centlist)
        plt.show()

    else:
        index = int(answer)
        # store og position first
        storage = centlist[index]
        movement = [width, height]
        while abs(movement[0]) >= width or abs(movement[1]) >= height:
            answer2 = input(
                f"Chosen circle {answer}, coordinates: {centlist[index][0], centlist[index][1]}. Enter valid movement in x and y direction: ")
            movement = answer2.split(" ")

            if len(movement) < 2:
                movement = [width, height]
            else:
                for i in range(len(movement)):
                    movement[i] = float(movement[i])
                # need to add distance to the point
                centlist[index] = [centlist[index][0] +
                                   float(movement[0]), centlist[index][1] + float(movement[1]), centlist[index][2]]
                # check valid movement
                if 0 < centlist[index][0] < width and 0 < centlist[index][1] < height:
                    pass
                else:
                    centlist[index] = storage
                    movement = [width, height]
        print(centlist[index])
        print(movement)
        # check the validity of the movement before you move forward with moving the point
        comparelist = centlist[:]
        comparelist.remove(centlist[index])
        indicator = check_a_circle(centlist[index], comparelist, tol)
        if indicator == False:
            print('enter when no overlap')
            nonsense = input('enter to confirm')
            centlist, vff = single_fillcircle(centlist[index], centlist, vff)

        # else reject and ask to move the point again
        else:
            print("invalid movement, try again")
            centlist[index] = storage
            nonsense = input('enter to confirm')
        modelprint(width, height, 0.8, vff, centlist)
        namespaced(centlist, tol, vff, width, height)

        plt.show()


def single_fillcircle(circle, centlist, vff):
    # fill possible circles to an individual selected circle
    # under movespaced function only
    if circle[1] > 0.8 or circle[1] < 0 or circle[0] > 0.8 or circle[0] < 0:
        pass
    else:
        radius = 2*circle[2] + tol
        angle = 0
        while angle < 360:
            x_diff = radius * cos(angle/180*pi)
            y_diff = radius * sin(angle/180*pi)
            imagcent = [circle[0] + x_diff,
                        circle[1] + y_diff, circle[2]]
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
    return centlist, vff


def setbounds(circle_i, centlist, rmax):
    limit = 4*rmax
    xcentre = circle_i
    comparelist = centlist[:]
    comparelist.remove(xcentre)
    possiblecircles = []
    for cir in range(len(comparelist)):
        # hypotenuse the radius, compare it such that it must be less than 4r, add to new list
        distance = ((xcentre[0] - comparelist[cir][0])**2 +
                    (xcentre[1] - comparelist[cir][1])**2)**0.5
        if distance < limit:
            possiblecircles.append([comparelist[cir], cir])
        else:
            pass
    # cap at a certian number of neighbours
    return possiblecircles


def closein(circle, idx, centlist, req_dis, tol):
    # create linear equation
    xcentre = circle
    possiblecircles = setbounds(xcentre, centlist, xcentre[2])
    if len(possiblecircles) == 0:
        return centlist
    else:
        # [[circle,index],[circle,index],[circle,index]...]
        """cannot random.uniform, must find space instead of randomise, if not its useless"""
        theta = (pi)*uniform(-1, 1)
        gradient = tan(theta)
        constant_c = xcentre[0]*gradient - xcentre[1]
        # linear_eqn = gradient*x + constant
        print("theta:", theta)
        print(gradient)

        x_1, y_1 = xcentre[0], xcentre[1]
        h_list = []
        for neighbours in possiblecircles:
            x_2, y_2 = neighbours[0][0], neighbours[0][1]

            # abs(tan(theta)*x_2 - y_2 + (y_1 - tan(theta)*x_1))/sqrt(tan(theta)**2 + (-1)**2)
            upperequation = abs(tan(theta)*x_2 - y_2 + (y_1 - tan(theta)*x_1))
            lowerequation = (tan(theta)**2 + (-1)**2)**0.5
            h = upperequation/lowerequation
            h_list = h_list + [h]
        print("h_list: ", h_list)
        min_h = min(h_list)
        index_nearest = h_list.index(min(h_list))
        x_3, y_3 = possiblecircles[index_nearest][0][0], possiblecircles[index_nearest][0][1]

        # pythagoras
        # part1eqn = (req_dis)**2 - (gradient*x_1 + y_1-gradient*x_1)
        # circumference eqn with varied angle alpha
        alpha = 2*pi - acos((min_h/(req_dis)) % 1)
        print('alpha: ', alpha)

        circum_x = x_3 + req_dis*cos(alpha)
        circum_y = y_3 + req_dis*sin(alpha)
        comparelist = centlist[:]
        comparelist.remove(comparelist[idx])
        current_cir = [circum_x, circum_y,
                       possiblecircles[index_nearest][0][2]]

        indicator = check_a_circle(
            current_cir, comparelist, tol)
        if indicator == True:
            pass
        else:
            centlist[idx] = current_cir
        return centlist
        # modelprint(width, height, 0.8, centlist, rmax, possiblecircles)


def recursion(index, centlist, possiblecircles, vff):
    if index == len(centlist):
        modelprint(width, height, 0.8, vff, centlist)
    else:
        centlist = closein(centlist[index], index, centlist, rmax*2+tol, tol)
        index = index + 1
        recursion(index, centlist, possiblecircles, vff)


"""START OF COMPUTATION"""
# Variable parameters
count = 0
possiblecircles = []
rmax = 0.016  # the maxmum radius of the fibres
rmin = 0.016  # the minmum radius of the fibres
width = rmax*50  # 0.8   # the length of the RVE
height = rmax*50  # 0.8 # the width of the RVE
tol = 0.001  # the minmum distance of two circles (except for the radius)
vf = 0.45    # the FVF in the RVE
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

centlist = []
vff = 0
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
        X = uniform(0 + 1 * rmax, width - 1 * rmax)
        Y = uniform(0 + 1 * rmax, height - 1 * rmax)
        R = uniform(rmin, rmax)
        newcircle = [X, Y, R]
        check = check_a_circle(newcircle, centlist, tol)
        if check == False:
            centlist = centlist + [newcircle]
            vff = vff + add_vff(newcircle[2], width, height)

    elif ratio_centre <= prob < ratio_heights:
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
            vff = vff + \
                add_vff(newcircle_1[2], width, height)

    elif ratio_heights <= prob < ratio_widths:
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
            vff = vff + add_vff(newcircle_3[2], width, height)

    else:
        X = uniform(0 - 1 * rmax, 0 + 1 * rmax)
        Y = uniform(0 - 1 * rmax, 0 + 1 * rmax)
        R = uniform(rmin, rmax)
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

            vff = vff + add_vff(newcircle_5[2], width, height)

        print(vff, ",", prob)  # FVF for each step
# centlist, vff = fillcircle(centlist, vff)
modelprint(width, height, 0.8, vff, centlist)
# namespaced(centlist, tol, vff, width, height)
# modelprint(width, height, 0.8, vff, centlist)

index = 0
recursion(index, centlist, possiblecircles, vff)

# while vff < 0.6:
# centlist, vff = random_movement(centlist, tol, vff)
# useless = input('shaken once, enter for result')
# namespaced(centlist, tol, vff, width, height)
# print('check space!')
modelprint(width, height, 0.8, vff, centlist)
