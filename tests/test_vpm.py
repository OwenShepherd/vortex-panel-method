import os.path
import sys

def func(x):
    return (x+1)


def test_answer():
    assert func(3) == 5


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
    npX = npX.astype(npX)
    npY = np.array(Y)
    npY = npY.astype(npY)
    npcl = float(cl[0])
    npCp = np.array(Cp)
    npCp = npCp.astype(np.float)

    XTest, YTest = Get_AirfoilCoordinates(m,p,t,c,N,False)

    print(X)
    input()
    print(XTest)




class TestClass(object):
    def test_one(self):
        x = 'this'
        assert 'h' in x

    def test_two(self):
        x = 'hello'
        assert hasattr(x, 'check')
