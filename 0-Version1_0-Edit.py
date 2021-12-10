###############################################################################
#                                                                             #
#                Function to Quantitatively Characterize                      #
#                                                                             #
#                    Randomness of Fibre Distribution                         #
#                                                                             #
###############################################################################
#                                                                             #
#                                V1.0                                         #
#                                                                             #
#                Yongfeng Ding - yongfeng_ding@mymail.sutd.edu.sg             #
#                              November 2021                                  #
#                                                                             #
###############################################################################
#                                                                             #
#  Modification log:                                                          #
#                                                                             #
#       2021/11/23 - Add Pair distribution function, V1.0 release             #
#       2021/11/17 - Add Ripley's K function                                  #
#       2021/11/16 - Add Nearest neighbor orientation                         #
#       2021/11/09 - Improvement of Nearest neighbor distance                 # 
#       2021/10/29 - Add Nearest neighbor distance - 1st nearest              #
#       2021/10/16 - plot RVE image in Python                                 #
#       2021/09/22 - Add periodicity on RVE algorithm in ABAQUS               #
#       2021/09/14 - Add RVE algorithm without periodicity in ABAQUS          #                          
#                                                                             #                                            
###############################################################################
#  Pre-processing (Algorithm of generating RVE microstructure)                #
###############################################################################
#
from matplotlib.patches import Ellipse, Circle, Rectangle
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

from random import *
from math import *

# Variable parameters
Wcomp = 0.8   # the length of the RVE 
Hcomp = 0.8   # the width of the RVE
Rmax = 0.016  # the maxmum radius of the fibres
Rmin = 0.016  # the minmum radius of the fibres
Tol = 0.0005  # the minmum distance of two circles (except for the radius)
Vf = 0.55      # the FVF in the RVE
#
# Algorithm of generating random distributing fibres
# Ratio of each part to the whole RVE area
ratio_PartI = (Wcomp - 2 * Rmax) * (Hcomp - 2 * Rmax) / (Wcomp + 2 * Rmax) * (Hcomp + 2 * Rmax)
ratio_PartII = (2 * 2 * Rmax * (Hcomp - 2 * Rmax)) / (Wcomp + 2 * Rmax) * (Hcomp + 2 * Rmax)
ratio_PartIII = (2 * 2 * Rmax * (Wcomp - 2 * Rmax)) / (Wcomp + 2 * Rmax) * (Hcomp + 2 * Rmax)
ratio_PartIV = (4 * 2 * Rmax * 2 * Rmax) / (Wcomp + 2 * Rmax) * (Hcomp + 2 * Rmax)
#
# Function to check the Intersection of TWO circles
def Intersection (R1, Cen1, R2, Cen2, Tol):
    Distan = ((Cen1[0] - Cen2[0]) ** 2.0 + (Cen1[1] - Cen2[1]) ** 2.0) ** 0.5
    Min_Dis = R1 + R2 + Tol
  
    if Distan > Min_Dis:
        Intersection = 'No'
    else:
        Intersection = 'Yes'
    return Intersection
#
# Function to check the Intersection of a circle with a list of circles  
def Intersec(Cen1, CentList, Tol):
    Intersec = 'No'
    R1 = Cen1[2]
    for Cen2 in CentList:
        R2 = Cen2[2]
        Check = Intersection(R1, Cen1, R2, Cen2, Tol)
        if Check == 'Yes':
            Intersec = 'Yes'
    return Intersec
#
# Function of ran
def ran(X1, X2):
    value = (X2 - X1) * np.random.rand() + X1
    return value
#
CentList = [[ran(0, Wcomp), ran(0, Hcomp), ran(Rmin, Rmax)]]  
Vff = pi * CentList[0][2] ** 2.0 / (Wcomp * Hcomp)
print("begin")    # an indication for the user to know whether the codes is correct or not 
while Vff < Vf:
    prob = random()
    # print(prob)  # prob for each step when Vff < Vf
#    
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
        continue        # for chck the code sentence by sentence
        
# for Part II,  -Rmax <= X <= Rmax, Rmax <= Y <= Hcomp - Rmax:     
    if prob < (ratio_PartI + ratio_PartII) and prob >= ratio_PartI:
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
            Vff = Vff + pi * Newcircle_2[2] ** 2.0 / (Wcomp * Hcomp)
        continue
     
# for Part III, Rmax <= X <= Wcomp -Rmax, -Rmax <= Y <= Rmax:
    if prob < (ratio_PartI + ratio_PartII + ratio_PartIII) and prob >= (ratio_PartI + ratio_PartII):   
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
            Vff = Vff + pi * Newcircle_4[2] ** 2.0 / (Wcomp * Hcomp)
        continue
         
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
           Vff = Vff + pi * Newcircle_5[2] ** 2.0 / (Wcomp * Hcomp)
           CentList = CentList + [Newcircle_6]
           Vff = Vff + pi * Newcircle_6[2] ** 2.0 / (Wcomp * Hcomp)
           CentList = CentList + [Newcircle_7]
           Vff = Vff + pi * Newcircle_7[2] ** 2.0 / (Wcomp * Hcomp)
           CentList = CentList + [Newcircle_8]
           Vff = Vff + pi * Newcircle_8[2] ** 2.0 / (Wcomp * Hcomp)
     
    # print(Vff)   # FVF for each step 
# print(CentList)  # The coordinate values of each fibre
print('FVF:', Vff, '\n')   # The final FVF
#
###############################################################################
#  0 - RVE image                                                              #
############################################################################### 
#
fig = plt.figure()
ax = fig.add_subplot(111)
rectan = Rectangle((0, 0), Wcomp, Hcomp, linewidth = 0.8, edgecolor = 'red', facecolor = 'w',
                   linestyle = 'solid', alpha = 1.0)
ax.add_patch(rectan)

for Cir in CentList:
    X = Cir[0]
    Y = Cir[1]
    R = Cir[2]
    Cir = Circle(xy = (X, Y), radius = R, color = 'black', alpha = 1.0)   
    ax.add_patch(Cir)
    plt.axis('scaled')
    plt.axis('equal')  # change limits of x or y axis so that equal increemets of x and y have the same length
plt.show()
#
############################################################################### 
#  1 - Nearest neighbor distance (including 1st, 2nd, 3rd nearest dist.)      #
############################################################################### 
#
# search the 1,2,3 nearest neighbor distance
from sklearn.neighbors import KDTree

A = np.array(CentList)   # convert list to array
kdt = KDTree(A, leaf_size=30, metric='euclidean')   # One searchig methods of KNN
#
dist, ind = kdt.query(A, k=4)                
# print(ind)  # indices of 4 closest neighbors
# print(dist)  # distances to 4 closest neighbors
#
#=============================================================================#
# sort the 1st, 2nd, 3rd nearest neighbor point and put in separate arrays
B1 = dist[:, 1]   # the 1st nearest neighbor distance, B0 reprents the point itself
C1 = B1 / R       # h / R for 1st
# print(C1, '\n')
#
B2 = dist[:, 2]   # the 2nd nearest neighbor distance
C2 = B2 / R       # h / R for 2nd
# print(C2, '\n')
#
B3 = dist[:, 3]   # the 3rd nearest neighbor distance
C3 = B3 / R       # h / R for 3rd
# print(C3, '\n')
#=============================================================================#
#
max_C1 = np.amax(C1)
min_C1 = np.amin(C1)   # the minmum value
print('1-1: max. and min. distance of 1st nearest neighbor distance:\n', max_C1, min_C1)
range_C1 = np.ptp(C1)   # be equivalent to maxY - minY, the ptp sentence can be used for list and array
# print(range_C1)
interval = range_C1 / 10
# print(interval)
#
# Generating probability of each interval
D1 = []
for i in range(10):
    num_Inter = ((C1 >= min_C1 + i * interval) & (C1 <= min_C1 + (i + 1) * interval)).sum()
    # print(num_Inter)
    prob_Inter = num_Inter / len(A)
    D1.append(prob_Inter)
print('Probability of 1st nearest neighbor distance in each interval:\n', D1, '\n')
#=============================================================================#
# search for range of B2 and calculate pdf of B2 (C2)
#
max_C2 = np.amax(C2)
min_C2 = np.amin(C2) 
print('1-2: max. and min. distance of 2nd nearest neighbor distance:\n', max_C2, min_C2)  
range_C2 = np.ptp(C2)   # be equivalent to maxY - minY
# print(range_C2)
interval = range_C2 / 10
# print(interval)
#
# Generating probability of each interval
D2 = []
for i in range(10):
    num_Inter = ((C2 >= min_C2 + i * interval) & (C2 <= min_C2 + (i + 1) * interval)).sum()
    # print(num_Inter)
    prob_Inter = num_Inter / len(A)
    D2.append(prob_Inter)
print('Probability of 2nd nearest neighbor distance in each interval:\n', D2, '\n')
#=============================================================================#
# search for range of B3 and calculate pdf of B3 (C3)
#
max_C3 = np.amax(C3)
min_C3 = np.amin(C3)   
print('1-3: max. and min. distance of 3rd nearest neighbor distance:\n', max_C3, min_C3)
range_C3 = np.ptp(C3)   # be equivalent to maxY - minY
# print(range_C3)
#
interval = range_C3 / 10
# print(interval)
#
# Generating probability of each interval
D3 = []
for i in range(10):
    num_Inter = ((C3 >= min_C3 + i * interval) & (C3 <= min_C3 + (i + 1) * interval)).sum()
    # print(num_Inter)
    prob_Inter = num_Inter / len(A)
    D3.append(prob_Inter)
print('Probability of 3rd nearest neighbor distance in each interval:\n', D3, '\n')
#=============================================================================#
# plot the 1st pdf
plt.plot(np.arange(2, 3.5, 0.15), D1, '--ro', label = 'Real structure-1st')
plt.xlabel('h/R')
plt.ylabel('Probability Density Function')
plt.ylim(0, 0.8)
plt.legend(loc = 'upper right')
plt.show()
#
# plot the 2nd pdf
plt.plot(np.arange(2, 4.5, 0.25), D2, '--ro', label = 'Real structure-2nd')
plt.xlabel('h/R')
plt.ylabel('Probability Density Function')
plt.ylim(0, 0.8)
plt.legend(loc = 'upper right')
plt.show()
#
# plot the 3rd pdf
plt.plot(np.arange(2, 5, 0.3), D3, '--ro', label = 'Real structure-3rd')
plt.xlabel('h/R')
plt.ylabel('Probability Density Function')
plt.ylim(0, 0.8)
plt.legend(loc = 'upper right')
plt.show()
#
############################################################################### 
#  2 - Nearest neighbor orientation                                           #
############################################################################### 
#
from scipy import spatial

A = CentList  # for simplity
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
     
    array_minP = array_B[spatial.KDTree(array_B).query(array_a)[1]]
    b = list(array_minP)   # list b   ========================================# 
    # print(b)
    minD = ((list(array_a)[0] - list(array_minP)[0]) ** 2.0 +
            (list(array_a)[1] - list(array_minP)[1]) ** 2.0) ** 0.5
    # print(minD)
    
    x_dis = abs(a[0] - b[0])  #   cosine edge
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
        alpha = 2 * pi -alpha_o
    # print(i, alpha)   # the angle of the two fibres 
#=============================================================================#
    i = 0
    for k in np.arange(start_l, end_l, interval_step):
        if alpha >= k and alpha < k + interval_step:
            M[i] += 1
        i += 1     
# print(M, '\n')   # counts of alpha in each interval of angles
#
o = [0]
for i in range(len(M)):
    x = sum(M[: i + 1])
    # print(x)
    prob_o = x / len(A)
    o.append(prob_o)
print('2: cumulative probability of alpha:\n', o, '\n')   # cumulative probability of alpha with increasing of intervals 
#=============================================================================#
# plot the nearest neighbor orientation
plt.plot(np.arange(0, n_step + 1, 1), o,'--ro', label = 'Nearest neighbor orientation')
plt.xlabel('Orientation [degree]')
plt.ylabel('Cumulative Distribution Function')
plt.xlim(0, 10)
plt.ylim(0, 1.0)
plt.legend(loc = 'upper left')
plt.show()
#
############################################################################### 
#  3 - Ripley's K function                                                    #
############################################################################### 
#
from astropy.stats import RipleysKEstimator
from matplotlib import pyplot as plt

x = np.array(CentList)   # array of list of CentList
z = x[:, [0, 1]]         # extract the first two columns 
# print(z, '\n')
#
Kest = RipleysKEstimator(area = Wcomp * Hcomp, x_max = Wcomp + Rmax, 
                         y_max = Hcomp + Rmax, x_min = -Rmax, y_min = - Rmax)
h = np.linspace(0, 0.3, 30)   # start from 0, end with 0.24, 100 floats, but how to determine the r value? 0.3 * x_max
plt.plot(h/R, Kest.poisson(h), color='green', ls=':', label=r'$K_{pois}$')  # poisson curve
print('3-1: y coordinate values on poisson curve:\n', Kest.poisson(h), '\n')
plt.plot(h/R, Kest(data=z, radii=h, mode='ripley'), '--ro', label=r'$K_{ripley}$')  # ripley curve
print('3-2: y coordinate values on Ripley curve:\n', Kest(data=z, radii=h, mode='ripley'), '\n')  # y coordinate value
plt.xlabel('h/R')
plt.ylabel('K(h)')
plt.legend()
plt.show()
#
############################################################################### 
#  4 - Pair distributiong function                                            #
############################################################################### 
#
y = Kest(data=z, radii=h, mode='ripley')  # y coordinates value of ripley curve
p = np.polyfit(h, y, 10)   # polynomial coefficients for fitting the ripley k function curve, order=10  
dp = [10 * p[0], 9 * p[1], 8 * p[2], 7 * p[3], 6 * p[4], 5 * p[5], 4 * p[6], 3 * p[7], 2 * p[8], 1 * p[9]]
y_deriv = np.polyval(dp, h)
g_r = y_deriv / (2 * pi * h)
print('4: g(r) values:\n', g_r)
#=============================================================================#
plt.plot(np.arange(0.0/R, 0.3/R, 0.01/R), g_r, '--r.', label = 'Pair distribution function')
plt.xlabel('h/R')
plt.ylabel('g(r)')
plt.ylim(0.0, 1.5)
plt.legend(loc = 'upper right')
plt.show()
#
###############################################################################
#  End of function                                                            #
###############################################################################
