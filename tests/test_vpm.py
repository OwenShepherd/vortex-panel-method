import os.path
import sys
import csv
import numpy as np

scriptDir = sys.path[0]
filepath = os.path.abspath(os.path.join(scriptDir,".."))
sys.path.insert(0,filepath)

from vpm.vpmfuncs import *

def test_NACA0012():
    testName = 'NACA0012'
    XPy, YPy, XMat, YMat, clPy, CpPy, clMat, cPMat = Data_Comparison(testName)
    for i in range(len(XPy)):
        assert(abs(XPy[i]-XMat[i])<=0.001)
        assert(abs(YPy[i]-YMat[i])<=0.001)
        if (i<=(len(XPy)-2)):
            assert(abs(cPMat[i]-CpPy[i])<=0.001)
    assert(abs(clPy-clMat)<=0.001)

def test_NACA1408():
    testName = 'NACA1408'
    XPy, YPy, XMat, YMat, clPy, CpPy, clMat, cPMat = Data_Comparison(testName)
    for i in range(len(XPy)):
        assert(abs(XPy[i]-XMat[i])<=0.001)
        assert(abs(YPy[i]-YMat[i])<=0.001)
        if (i<=(len(XPy)-2)):
            assert(abs(cPMat[i]-CpPy[i])<=0.001)
    assert(abs(clPy-clMat)<=0.001)

def test_NACA2410():
    testName = 'NACA2410'
    XPy, YPy, XMat, YMat, clPy, CpPy, clMat, cPMat = Data_Comparison(testName)
    for i in range(len(XPy)):
        assert(abs(XPy[i]-XMat[i])<=0.001)
        assert(abs(YPy[i]-YMat[i])<=0.001)
        if (i<=(len(XPy)-2)):
            assert(abs(cPMat[i]-CpPy[i])<=0.001)
    assert(abs(clPy-clMat)<=0.001)

def test_NACA4415():
    testName = 'NACA4415'
    XPy, YPy, XMat, YMat, clPy, CpPy, clMat, cPMat = Data_Comparison(testName)
    for i in range(len(XPy)):
        assert(abs(XPy[i]-XMat[i])<=0.001)
        assert(abs(YPy[i]-YMat[i])<=0.001)
        if (i<=(len(XPy)-2)):
            assert(abs(cPMat[i]-CpPy[i])<=0.001)
    assert(abs(clPy-clMat)<=0.001)

def test_NACA6412():
    testName = 'NACA6412'
    XPy, YPy, XMat, YMat, clPy, CpPy, clMat, cPMat = Data_Comparison(testName)
    for i in range(len(XPy)):
        assert(abs(XPy[i]-XMat[i])<=0.001)
        assert(abs(YPy[i]-YMat[i])<=0.001)
        if (i<=(len(XPy)-2)):
            assert(abs(cPMat[i]-CpPy[i])<=0.001)
    assert(abs(clPy-clMat)<=0.001)

def Data_Comparison(NACA_Name):
    X = []
    Y = []
    cl = []
    Cp = []
    airfoil_parameters = []
    alpha = 1

    for i in range(4,len(NACA_Name)-1):
        if (i <= 5):
            airfoil_parameters.append(float(NACA_Name[i]))
        else:
            airfoil_parameters.append(float(NACA_Name[i] + NACA_Name[i+1]))

    m = airfoil_parameters[0]/100
    p = airfoil_parameters[1]/100
    t = airfoil_parameters[2]/100
    c = 1
    N = 100

    relpath = sys.path[0] + '/tests/test_data/' + NACA_Name + '.csv'

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
    npcl = np.array(cl)
    npcl = npcl.astype(np.float)
    npCp = np.array(Cp)
    npCp = npCp.astype(np.float)

    XTest, YTest = Get_AirfoilCoordinates(m,p,t,c,N,False)
    newcl, newCp = Get_PanelCoefficients(XTest,YTest,N,alpha,'',False)
    XTest = np.around(XTest,decimals=4)
    npX = np.around(npX,decimals=4)
    YTest = np.around(YTest,decimals=4)
    npY = np.around(npY,decimals=4)
    newcl = np.around(newcl,decimals=4)
    npcl = np.around(npcl,decimals=4)
    npCp = np.around(npCp,decimals=4)
    newCp = np.around(newCp,decimals=4)

    return XTest, YTest, npX, npY, newcl, newCp, npcl, npCp
