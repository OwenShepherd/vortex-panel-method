import math

def NACAAirfoil(m,p,t,c,N=100,PLOT=False):
    """ Computes the x-y locations of a 2-D NACA Airfoil"""


    # Splitting the airfoil into equal-length x-locations ALA K&C
    dt = 2*math.pi/N

    # Creating some control values and storage variables
    NOTDONE = True
    th = []
    th.append(0)
    i = 1

    # Want angles to go from zero to pi and split the boundary x-y locations
    # into upper and lower locations
    while (NOTDONE):
        value = th[i-1] + dt

        if (abs(value-math.pi)<=0.0000001):
            th.append(value)
            NOTDONE = False
        else:
            th.append(value)
            i = i+1

    # Now here we determine the x-lcations based off of the thetas
    R = 0.5*c
    xb = [R+R*math.cos(k) for k in th]

    # Now we shall determine the points of interest on the airfoil (the y-
    # locations)
    yt = [t*c/0.2*(0.2969*math.sqrt(k/c)-0.1260*(k/c)-0.3516*(k/c)*(k/c)+0.2843*pow(k/c,3)-0.1036*pow(k/c,4)) for k in xb]


    # Here we can determine the mean camber line
    yc = []
    dyc = []

    for x in xb:
        if (x<=(p*c) and p!=0):
            yc.append(m*x/pow(p,2)*(2*p-x/c))
            dyc.append(m*x/pow(p,2)*(-1/float(c))+m/pow(p,2)*(2*p-x/c))
        else:
            yc.append(m*(c-x)/(1-p*p)*(1+x/c-2*p))
            dyc.append(m*(c-x)/(1-p*p)*1/float(c)+-1*m/(1-p*p)*(1+x)/(c-2*p))

    # Now we can define parameter zeta
    zeta = [math.atan(k) for k in dyc]



    XU = []
    YU = []
    XL = []
    YL = []

    # Creating the upper and lower x,y locations
    for x in range(0,len(xb)):
        XU.append(xb[x])
        YU.append(yc[x]+yt[x]*math.cos(zeta[x]))

        XL.append(xb[x])
        YL.append(yc[x]-yt[x]*math.cos(zeta[x]))



    # This ensures that the boundary points go from the trailing edge, clockwise
    # back to the trailing edge with two boundary points at the trailing edge
    # and one at the leading edge
    XU.pop()
    YU.pop()


    X = XL+list(reversed(XU))
    Y = YL+list(reversed(YU))


    return X,Y
