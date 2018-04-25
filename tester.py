#!/usr/bin/env python3

import csv
import os.path
import sys
import numpy as np
import matplotlib as plt
sys.path.insert(0,os.path.abspath(sys.path[0]))

from vpm.vpmfuncs import *
from vpm.examples import *

def tester1():
    X = []
    Y = []
    cl = []

    m = 0/100
    p = 0/100
    t = 12/100
    c = 1
    N = 100
    alpha = 1
    Cp = []

    relpath = sys.path[0] + '/tests/test_data/NACA0012.csv'
    with open(relpath,'r') as csvfile:
        spamreader = csv.reader(csvfile,delimiter=',')
        for row in spamreader:
            if (X == []):
                X = row
            elif (Y == []):
                Y = row
            elif (cl == []):
                cl = row
            elif (Cp == []):
                Cp = row

    # Conversion to numpy arrays
    npX = np.array(X)
    npX = npX.astype(np.float)
    npY = np.array(Y)
    npY = npY.astype(np.float)
    npcl = float(cl[0])
    npCp = np.array(Cp)
    npCp = npCp.astype(np.float)

    XTest, YTest = Get_AirfoilCoordinates(m,p,t,c,N,False)
    clTest, CpTest = Get_PanelCoefficients(XTest,YTest,N,alpha,'',False)
    XTest = np.around(XTest,decimals=5)
    npX = np.around(npX,decimals=5)
    print("ClTest: " + str(clTest) + " CpTest: " + str(CpTest))
    for i in range(len(npX)):
        if (XTest[i] != npX[i]):
            print("Index: " + str(i) + " Not Equal.")
            print("XTest: " + str(XTest[i]) + " npX: " + str(npX[i]))
            input()

def StringExtractor():
    testString = 'NACA0012'
    Data = []
    for i in range(4,len(testString)-1):
        if (i <= 5):
            Data.append(float(testString[i]))
        else:
            Data.append(float(testString[i] + testString[i+1]))

    print(Data)

def Plotting_Tester():
    name = 'NACA0012'
    plotColor = '#00ff08'
    coeffFIG = plt.figure()
    c = 1
    N = 100
    alpha = 5
    plotColor = '#00ff08'
    NACA0012 = Airfoil(name,c,N,alpha,coeffFIG,plotColor)
    NACA0012.Get_ParsedData()
    NACA0012.Get_AirfoilCoordinates()
    NACA0012.Get_PanelCoefficients(True)
    newFig = plt.figure()
    alpha = 2
    plotColor = '#ff1f00'
    NACA0012.Get_ParsedData()
    NACA0012.Set_AngleOfAttack(alpha)
    NACA0012.Set_PlotColor(plotColor)
    NACA0012.Set_Figure(newFig)
    NACA0012.Get_PanelCoefficients(True)
    alpha = 10
    plotColor = '#0074ff'
    NACA0012.Get_ParsedData()
    NACA0012.Set_AngleOfAttack(alpha)
    NACA0012.Set_PlotColor(plotColor)
    NACA0012.Get_PanelCoefficients(True)
    plt.show()

def exampleTester():
    example_NACA0012()



def main():
    exampleTester()


if __name__ == "__main__":
    main()
