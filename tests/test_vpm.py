import os.path
import sys
import csv
import numpy as np

scriptDir = sys.path[0]
filepath = os.path.abspath(os.path.join(scriptDir,".."))
sys.path.insert(0,filepath)

from vpm.vpmfuncs import *

def test_NACA0012():
    X = []
    Y = []
    cl = []
    Cp = []

    m = 0/100
    p = 0/100
    t = 12/100
    c = 1
    N = 100


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
    XTest = np.around(XTest,decimals=5)
    npX = np.around(npX,decimals=5)
    for i in range(len(npX)):
        assert(XTest[i] == npX[i])



class TestClass(object):
    def test_one(self):
        x = 'this'
        assert 'h' in x
