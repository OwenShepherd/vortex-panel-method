#!/usr/bin/env python3

import csv
import os.path
import sys
import numpy as np
import matplotlib as plt

from vpm import vpmfuncs

def test_NACA0012():
    X = []
    Y = []
    cl = []

    NACA0012 = vpmfuncs.Airfoil('NACA0012')

    file_directory = os.path.abspath(os.path.join(__file__,'..'))
    sys.path.insert(0,file_directory)
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

    Airfoil.get_panel_coefficients()
    Airfoil.get_airfoil_coordinates()
    Airfoil.get_panel_coefficients()

    XTest = Airfoil.x_boundary_points
    YTest = Airfoil.y_boundary_points
    clTest = self.full_coefficient_lift
    CpTest = self.pressure_coefficient

    XTest = np.around(XTest,decimals=5)
    npX = np.around(npX,decimals=5)
    assert XTest == npX

def tester1():
    X = []
    Y = []
    cl = []

    file_directory = os.path.abspath(os.path.join(__file__,'..'))
    sys.path.insert(0,file_directory)
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

    XTest, YTest = vpmfuncs.Get_AirfoilCoordinates(m,p,t,c,N,False)
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

def exampleTester():
    StringExtractor()
    Plotting_Tester()



def main():
    exampleTester()


if __name__ == "__main__":
    main()
