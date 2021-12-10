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
    value = (X2 - X1) * np.random.rand() + X1
    return value
#


# Variable parameters
Wcomp = 0.8   # the length of the RVE
Hcomp = 0.8   # the width of the RVE
Rmax = 0.016  # the maxmum radius of the fibres
Rmin = 0.016  # the minmum radius of the fibres
Tol = 0.0005  # the minmum distance of two circles (except for the radius)
Vf = 0.55      # the FVF in the RVE
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

print(ratio_PartI, ratio_PartII, ratio_PartIII, ratio_PartIV)

CentList = [[ran(0, Wcomp), ran(0, Hcomp), ran(Rmin, Rmax)]]
Vff = pi * CentList[0][2] ** 2.0 / (Wcomp * Hcomp)
# an indication for the user to know whether the codes is correct or not
print("begin")

print(CentList)

print(np.random.uniform(0, Wcomp))