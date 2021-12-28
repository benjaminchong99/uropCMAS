from matplotlib.patches import Circle, Rectangle
import numpy as np

from random import *
from math import *
from sklearn.neighbors import KDTree
from scipy import spatial
from astropy.stats import RipleysKEstimator
from matplotlib import pyplot as plt

'''
Problems:
1. place all codes outside into functions. OOP maybe?
2. import part need to standardise, got repeated imports
3. Have a function to plot the graphs
4. MOST TIME SPENT ON FINDING THE FVF FOR SOME REASON, NEED TO CHECK ON THAT

'''


# Function to check the Intersection of TWO circles


def Intersection(R1, Cen1, R2, Cen2, Tol):
    Distan = ((Cen1[0] - Cen2[0]) ** 2.0 + (Cen1[1] - Cen2[1]) ** 2.0) ** 0.5
    Min_Dis = R1 + R2 + Tol

    if Distan > Min_Dis:
        overlap = 'No'
    else:
        overlap = 'Yes'
    return overlap
#
# Function to check the Intersection of a circle with a list of circles


# !!! WARNING !!! VARIABLE NAME SAME AS FUNCTION
def Intersec(Cen1, CentList, Tol):
    Intersect = 'No'
    R1 = Cen1[2]
    for Cen2 in CentList:
        R2 = Cen2[2]
        Check = Intersection(R1, Cen1, R2, Cen2, Tol)
        if Check == 'Yes':
            Intersect = 'Yes'
    return Intersect
#
# Function of ran

# put a seed value to standardise


def ran(X1, X2):
    # randomise
    # not sure if this is needed, might have built in function alrdy
    value = (X2 - X1) * np.random.rand() + X1
    return value
#


# new function
# np_arrange --> format your x values and plot points
# should be done, just need test out
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

def recount_limit(CentList, Vff, record, section, count):
    if len(record) > count:
        if section == 1:
            CentList = CentList + [Newcircle]
            Vff = Vff + pi * Newcircle[2] ** 2.0 / (Wcomp * Hcomp)
            record = [Vff]
        elif section == 2:
            CentList = CentList + [Newcircle_1]
            Vff = Vff + pi * Newcircle_1[2] ** 2.0 / (Wcomp * Hcomp)
            CentList = CentList + [Newcircle_2]
            # this Vff will overwrite the previous one, add recursion?
            Vff = Vff + pi * Newcircle_2[2] ** 2.0 / (Wcomp * Hcomp)
            record = [Vff]
        elif section == 3:
            CentList = CentList + [Newcircle_3]
            Vff = Vff + pi * Newcircle_3[2] ** 2.0 / (Wcomp * Hcomp)
            CentList = CentList + [Newcircle_4]
            # this Vff will overwrite the previous one, add recursion?
            Vff = Vff + pi * Newcircle_4[2] ** 2.0 / (Wcomp * Hcomp)
            record = [Vff]
        elif section == 4:
            CentList = CentList + [Newcircle_5]
            Vff = Vff + pi * Newcircle_5[2] ** 2.0 / (Wcomp * Hcomp)
            CentList = CentList + [Newcircle_6]
            # this Vff will overwrite the previous one, add recursion?
            Vff = Vff + pi * Newcircle_6[2] ** 2.0 / (Wcomp * Hcomp)
            CentList = CentList + [Newcircle_7]
            # this Vff will overwrite the previous one, add recursion?
            Vff = Vff + pi * Newcircle_7[2] ** 2.0 / (Wcomp * Hcomp)
            CentList = CentList + [Newcircle_8]
            # this Vff will overwrite the previous one, add recursion?
            Vff = Vff + pi * Newcircle_8[2] ** 2.0 / (Wcomp * Hcomp)
            record = [Vff]
    else:
        pass
    return record, Vff

class Queue:
    def __init__(self):
        self._items = []

    def enqueue(self, item):
        self._items.insert(0, item)

    def dequeue(self):
        return self._items.pop()

    def is_empty(self):
        return self._items == []

    def peek(self):
        return self._items[-1]



"""START OF COMPUTATION"""
# Variable parameters
Wcomp = 0.8   # the length of the RVE
Hcomp = 0.8   # the width of the RVE
Rmax = 0.016  # the maxmum radius of the fibres
Rmin = 0.016  # the minmum radius of the fibres
Tol = 0.001  # the minmum distance of two circles (except for the radius)
Vf = 0.65    # the FVF in the RVE
# cap at 0.55
#
# Algorithm of generating random distributing fibres
# Ratio of each part to the whole RVE area
ratio_PartI = (Wcomp - 2 * Rmax) * (Hcomp - 2 * Rmax) / \
    (Wcomp + 2 * Rmax) * (Hcomp + 2 * Rmax)
ratio_PartII = (2 * 2 * Rmax * (Hcomp - 2 * Rmax)) / \
    (Wcomp + 2 * Rmax) * (Hcomp + 2 * Rmax)
ratio_PartIII = (2 * 2 * Rmax * (Wcomp - 2 * Rmax)) / \
    (Wcomp + 2 * Rmax) * (Hcomp + 2 * Rmax)
ratio_PartIV = (4 * 2 * Rmax * 2 * Rmax) / \
    (Wcomp + 2 * Rmax) * (Hcomp + 2 * Rmax)
#
#
CentList = [[ran(0, Wcomp), ran(0, Hcomp), ran(Rmin, Rmax)]]
Vff = pi * CentList[0][2] ** 2.0 / (Wcomp * Hcomp)
# an indication for the user to know whether the codes is correct or not
print("begin")
print("Vff: ", Vff, "\n Vf: ", Vf)
print("Ratio 1: ", ratio_PartI)
print("Ratio 2: ", ratio_PartII)
print("Ratio 3: ", ratio_PartIII)
print("Ratio 4: ", ratio_PartIV)

# ^ up till here still fast
#count = 0
#randlist = np.random.random(size=20000)
# print(randlist)
#listofrandom = randlist.tolist()
# for i in listofrandom


firstRatio = ratio_PartI
secondRatio = ratio_PartI + ratio_PartII
thirdRatio = ratio_PartI + ratio_PartII + ratio_PartIII
fourthRatio = ratio_PartI + ratio_PartII + ratio_PartIII + ratio_PartIV


# !!! WARNIING BIG WHILE LOOP HERE !!!
count = 0
record_count = 0
record = [0]

# Initialize queue
Vffqueue = Queue()
Vffqueue.enqueue(Vff)
current_Vff = Vffqueue.dequeue()
print(current_Vff)

while current_Vff < Vf:
    count = count + 1
    prob = random()
    # print(prob)  # prob for each step when Vff < Vf

    # for Part I, Rmax < X < Wcomp - Rmax, Rmax < Y < Hcomp - Rmax:
    if prob < ratio_PartI:
        X = ran(0 + 1 * Rmax, Wcomp - 1 * Rmax)
        Y = ran(0 + 1 * Rmax, Hcomp - 1 * Rmax)
        R = ran(Rmin, Rmax)
        Newcircle = [X, Y, R]
        Check = Intersec(Newcircle, CentList, Tol)
        if Check == 'No':
            CentList = CentList + [Newcircle]
            Vff = Vff + pi * Newcircle[2] ** 2.0 / (Wcomp * Hcomp)
            # for chck the code sentence by sentence

        record.append(Vff)

        if record[-1] != record[0]:
            record = [Vff]
        else:
            record += [Vff]
            record, Vff = recount_limit(CentList, Vff, record, 1, 200)
            # if len(record) > 2000:
            #    CentList = CentList + [Newcircle]
            #    Vff = Vff + pi * Newcircle[2] ** 2.0 / (Wcomp * Hcomp)
            #    record = [Vff]
        Vffqueue.enqueue(Vff)

    # for Part II,  -Rmax <= X <= Rmax, Rmax <= Y <= Hcomp - Rmax:
    elif ratio_PartI <= prob < (ratio_PartI + ratio_PartII):
        X = ran(0 - 1 * Rmax, 0 + 1 * Rmax)
        Y = ran(0 + 1 * Rmax, Hcomp - 1 * Rmax)
        R = ran(Rmin, Rmax)
        Newcircle_1 = [X, Y, R]
        Newcircle_2 = [X + Wcomp, Y, R]
        Check_1 = Intersec(Newcircle_1, CentList, Tol)
        Check_2 = Intersec(Newcircle_2, CentList, Tol)

        if Check_1 == 'No' and Check_2 == 'No':
            CentList = CentList + [Newcircle_1]
            Vff = Vff + pi * Newcircle_1[2] ** 2.0 / (Wcomp * Hcomp)
            CentList = CentList + [Newcircle_2]
            Vffqueue.enqueue(Vff)
            # this Vff will overwrite the previous one, add recursion?
            Vff = Vff + pi * Newcircle_2[2] ** 2.0 / (Wcomp * Hcomp)
            Vffqueue.enqueue(Vff)

        record.append(Vff)

        if record[-1] != record[0]:
            record = [Vff]
        else:
            record += [Vff]
            record, Vff = recount_limit(CentList, Vff, record, 2, 200)
            # if len(record) > 2000:
            #     CentList = CentList + [Newcircle_1]
            #     Vff = Vff + pi * Newcircle_1[2] ** 2.0 / (Wcomp * Hcomp)
            #     CentList = CentList + [Newcircle_2]
            #     # this Vff will overwrite the previous one, add recursion?
            #     Vff = Vff + pi * Newcircle_2[2] ** 2.0 / (Wcomp * Hcomp)
            #     record = [Vff]
        Vffqueue.enqueue(Vff)

    # for Part III, Rmax <= X <= Wcomp -Rmax, -Rmax <= Y <= Rmax:
    elif (ratio_PartI + ratio_PartII) <= prob < (ratio_PartI + ratio_PartII + ratio_PartIII):
        X = ran(0 + 1 * Rmax, Wcomp - 1 * Rmax)
        Y = ran(0 - 1 * Rmax, 0 + 1 * Rmax)
        R = ran(Rmin, Rmax)
        Newcircle_3 = [X, Y, R]
        Newcircle_4 = [X, Y + Hcomp, R]
        Check_3 = Intersec(Newcircle_3, CentList, Tol)
        Check_4 = Intersec(Newcircle_4, CentList, Tol)

        if Check_3 == 'No' and Check_4 == 'No':
            CentList = CentList + [Newcircle_3]
            Vff = Vff + pi * Newcircle_3[2] ** 2.0 / (Wcomp * Hcomp)
            CentList = CentList + [Newcircle_4]
            Vffqueue.enqueue(Vff)
            # this Vff will overwrite the previous one, add recursion?
            Vff = Vff + pi * Newcircle_4[2] ** 2.0 / (Wcomp * Hcomp)
            Vffqueue.enqueue(Vff)

        record.append(Vff)

        if record[-1] != record[0]:
            record = [Vff]
        else:
            record += [Vff]
            record, Vff = recount_limit(CentList, Vff, record, 3, 200)
            # if len(record) > 2000:
            #     CentList = CentList + [Newcircle_3]
            #     Vff = Vff + pi * Newcircle_3[2] ** 2.0 / (Wcomp * Hcomp)
            #     CentList = CentList + [Newcircle_4]
            #     # this Vff will overwrite the previous one, add recursion?
            #     Vff = Vff + pi * Newcircle_4[2] ** 2.0 / (Wcomp * Hcomp)
            #     record = [Vff]
        Vffqueue.enqueue(Vff)

    # for Part IV, Rmax <= X <= -Rmax, -Rmax <= Y <= Rmax:
    else:
        X = ran(0 - 1 * Rmax, 0 + 1 * Rmax)
        Y = ran(0 - 1 * Rmax, 0 + 1 * Rmax)
        R = ran(Rmin, Rmax)
        Newcircle_5 = [X, Y, R]
        Newcircle_6 = [X + Wcomp, Y, R]
        Newcircle_7 = [X, Y + Hcomp, R]
        Newcircle_8 = [X + Wcomp, Y + Hcomp, R]

        Check_5 = Intersec(Newcircle_5, CentList, Tol)
        Check_6 = Intersec(Newcircle_6, CentList, Tol)
        Check_7 = Intersec(Newcircle_7, CentList, Tol)
        Check_8 = Intersec(Newcircle_8, CentList, Tol)

        if Check_5 == 'No' and Check_6 == 'No' and Check_7 == 'No' and Check_8 == 'No':
            CentList = CentList + [Newcircle_5]
            # this Vff will overwrite the previous one, add recursion?
            Vff = Vff + pi * Newcircle_5[2] ** 2.0 / (Wcomp * Hcomp)
            CentList = CentList + [Newcircle_6]
            Vffqueue.enqueue(Vff)
            # this Vff will overwrite the previous one, add recursion?
            Vff = Vff + pi * Newcircle_6[2] ** 2.0 / (Wcomp * Hcomp)
            CentList = CentList + [Newcircle_7]
            Vffqueue.enqueue(Vff)
            # this Vff will overwrite the previous one, add recursion?
            Vff = Vff + pi * Newcircle_7[2] ** 2.0 / (Wcomp * Hcomp)
            CentList = CentList + [Newcircle_8]
            Vffqueue.enqueue(Vff)
            # this Vff will overwrite the previous one, add recursion?
            Vff = Vff + pi * Newcircle_8[2] ** 2.0 / (Wcomp * Hcomp)

            record.append(Vff)

            if record[-1] != record[0]:
                record = [Vff]
            else:
                record += [Vff]
                record, Vff = recount_limit(CentList, Vff, record, 4, 200)
                # if len(record) > 2000:
                #     CentList = CentList + [Newcircle_5]
                #     Vff = Vff + pi * Newcircle_5[2] ** 2.0 / (Wcomp * Hcomp)
                #     CentList = CentList + [Newcircle_6]
                #     # this Vff will overwrite the previous one, add recursion?
                #     Vff = Vff + pi * Newcircle_6[2] ** 2.0 / (Wcomp * Hcomp)
                #     CentList = CentList + [Newcircle_7]
                #     # this Vff will overwrite the previous one, add recursion?
                #     Vff = Vff + pi * Newcircle_7[2] ** 2.0 / (Wcomp * Hcomp)
                #     CentList = CentList + [Newcircle_8]
                #     # this Vff will overwrite the previous one, add recursion?
                #     Vff = Vff + pi * Newcircle_8[2] ** 2.0 / (Wcomp * Hcomp)
                #     record = [Vff]
                Vffqueue.enqueue(Vff)

        Vffqueue.enqueue(Vff)

    current_Vff = Vffqueue.dequeue()
    print(Vff)   # FVF for each step
# print(CentList)  # The coordinate values of each fibre
print('FVF:', Vff, '\n')   # The final FVF
print(CentList)
print(count)


fig = plt.figure()
ax = fig.add_subplot(111)
rectan = Rectangle((0, 0), Wcomp, Hcomp, linewidth=0.8, edgecolor='red', facecolor='w',
                   linestyle='solid', alpha=1.0)
ax.add_patch(rectan)

for Cir in CentList:
    X = Cir[0]
    Y = Cir[1]
    R = Cir[2]
    Cir = Circle(xy=(X, Y), radius=R, color='black', alpha=1.0)
    ax.add_patch(Cir)
    plt.axis('scaled')
    # change limits of x or y axis so that equal increemets of x and y have the same length
    plt.axis('equal')
plt.show()
#
