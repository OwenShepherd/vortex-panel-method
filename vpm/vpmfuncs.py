import math
import numpy as np
import matplotlib.pyplot as plt




def Get_LiftCoefficients(V,S,M):
    """
    Function: Calculate_LiftCoefficients

    Purpose: This function computes the coefficient of lift for an airfoil as to
    be used in the vortex panel method "Calculate_PanelCoefficients"

    Parameters:
        V - The dimensionlesss velocity at each control point
        S - The dimensionless length of each of the control points
        M - The number of panels

    Returns:
        cl - The coefficient of lift
    """

    gamma = 0
    for j in range(M):
        gamma = gamma + V[j]*S[j]

    cl = 2*gamma

    return cl



def Get_PanelCoefficients(XB, YB, M, alpha, NACA, PLOT, color):
    """
    Function: CalculatePanelCoefficients
    Purpose: Formulates the system of equations for the vortex paneling method.

    Parameters:
        XB - The dimensionless boundary x locations on the airfoil
        YB - The dimensionless boundary y locations on the airfoil
        M - The total number of panels
        alpha - The angle of attack
        NACA - The string that specifies the 4-digit NACA airfoil, ex: 0012
        PLOT - Whether or not to plot Cp vs. x/c; if the plot is desired,
            Plot == 1, if the plot is not desired then PLOT == 0

    Returns:
        cl - The lift coefficient
    """

    alpha = math.radians(alpha)

    MP1 = M+1 # The trailing edge requries an extra point

    """ TODO: Find better variable names for these arrays """
    X = np.zeros((1,M))
    Y = np.zeros((1,M))
    RHS = np.zeros((M,1))
    theta = np.zeros((1,M))
    S = np.zeros((1,M))

    for i in range(M):
        IP1 = i + 1

        X[i] = 0.5 * (XB[i] + XB[IP1])
        Y[i] = 0.5 * (YB[i] + YB[IP1])

        S[i] = math.sqrt(math.pow(XB[IP1] - XB[i],2) + math.pow(YB[IP1] - YB[i],2))

        theta[i] = math.atan2(YB[IP1] - YB[i], XB[IP1] - XB[i])

        RHS[i,1] = math.sin(theta[i] - alpha)


    CN1 = np.zeros((M,M))
    CN2 = np.zeros((M,M))
    CT1 = np.zeros((M,M))
    CT2 = np.zeros((M,M))
    An = np.zeros((MP1))
    At = np.zeros((MP1))


    for i in range(M):
        for j in range(M):

            if (i == j):
                CN1[i,j] = -1
                CN2[i,j] = 1
                CT1[i,j] = 0.5 * math.pi
                CT2[i,j] = 0.5 * math.pi
            else:
                A = -1*(X[i]-XB[j])*math.cos(theta[j])-(Y[i]-YB[j])*math.sin(theta[j])
                B = (X[i]-XB[j])**2+(Y[i]-YB[j])**2
                C = math.sin(theta[i]-theta[j])
                D = math.cos(theta[i]-theta[j])
                E = (X[i]-XB[j])*math.sin(theta[j])-(Y[i]-YB[j])*math.cos(theta[j])
                F = math.log(1+S[j]*(S[j]+2*A)/B)
                G = math.atan2((E*S[j]),(B+A*S[j]))
                P = (X[i]-XB[j])*math.sin(theta[i]-2*theta[j])+(Y[i]-YB[j])*math.cos(theta[i]-2*theta[j])
                Q = (X[i]-XB[j])*math.cos(theta[i]-2*theta[j])-(Y[i]-YB[j])*math.sin(theta[i]-2*theta[j])

                CN2[i,j] = D + 0.5 * Q * F/S[j] - (A*C + D*E) * G/S[j]
                CN1[i,j] = 0.5*D*F + C*G - CN2[i,j]
                CT2[i,j] = C + 0.5*P*F/S[j] + (A*D - C*E) * G/S[j]
                CT1[i,j] = 0.5*C*F - D*G - CT2[i,j]

    for i in range(M):
        An[i,1] = CN1[i,1]
        An[i,MP1] = CN1[i,1]
        At[i,1] = CT1[i,1]
        At[i,MP1] = CT2[i,M]

        for j in range(1,M):
            An[i,j] = CN1[i,j] + CN2[i,(j-1)]
            At[i,j] = CT1[i,j] + CT2[i,(j-1)]

    # Trailing edge conditions
    An[MP1,1] = 1
    An[MP1,MP1] = 1

    for j in range(1,M):
        An[MP1,j] = 0

    RHS[MP1] = 0

    gamma = np.linalg.solve(An,RHS)

    V = np.zeros((1,M))
    Cp = zeros((1,M))


    for i in range(M):
        V[i] = math.cos(theta[i] - alpha)

        for j in range(MP1):
            V[i] = V[i] + At[i,j]*gamma[j]
            C[i] = 1 - (V[i])**2

    cl = Calculate_LiftCoefficients(V,S,M)

    CpLower = Cp[:M/2]
    CpUpper = Cp[M/2:]

    if (PLOT):
        """ TODO: Plotting """


    return cl,Cp





def Get_AirfoilCoordinates(m,p,t,c,N,PLOT):
    """
    Function: Get_AirfoilCoordinates

    Purpose: Calculates the points on a NACA 4-digit series airfoil.

    Parameters:
        m - The maximum camber
        p - The location of the maximum camber
        t - The thickness
        c - The chord length
        N - The number of employed panels
        PLOT - Whether or not to plot Cp vs. x/c; if the plot is desired,
            PLOT == 1

    Returns:
        X - The x-locations on the airfoil
        Y - The y-locations on the airfoil
    """

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
    yt = [t*c/0.2*(0.2969*math.sqrt(k/c)-0.1260*(k/c)-0.3516*(k/c)*(k/c)+0.2843*(k/c)**3-0.1036**4) for k in xb]


    # Here we can determine the mean camber line
    yc = np.zeros((1,len(xb)))
    dyc = np.zeros((1,len(xb)))

    if ((p!=0) or (m!=0)):
        for i in range(len(xb)):
            if (xb[i]<=(p*c)):
                yc[0,i] = (m*xb[i]/(p**2)*(2*p-xb[i]/c))
                dyc[0,i] = (m*xb[i]/(p**2)*(-1/float(c))+m/(p**2)*(2*p-xb[i]/c))
            else:
                yc[0,i] = (m*(c-xb[i])/(1-p*p)*(1+xb[i]/c-2*p))
                dyc[0,i] = (m*(c-xb[i])/(1-p*p)*1/float(c)+-1*m/(1-p*p)*(1+xb[i])/(c-2*p))

    # Now we can define parameter zeta
    zeta = [math.atan(k) for k in dyc[0]]

    XU = np.zeros((1,len(xb)))
    YU = np.zeros((1,len(xb)))
    XL = np.zeros((1,len(xb)))
    YL = np.zeros((1,len(xb)))

    # Creating the upper and lower x,y locations
    for i in range(len(xb)):
        XU[0,i] = (xb[i])
        YU[0,i] = (yc[0,i]+yt[i]*math.cos(zeta[i]))

        XL[0,i] = (xb[i])
        YL[0,i] = (yc[0,i]-yt[i]*math.cos(zeta[i]))



    # This ensures that the boundary points go from the trailing edge, clockwise
    # back to the trailing edge with two boundary points at the trailing edge
    # and one at the leading edge
    XU = np.delete(XU,-1)
    YU = np.delete(YU,-1)
    XU = XU[::-1]
    YU = YU[::-1]
    XL = XL[0,:]
    YL = YL[0,:]

    X = np.concatenate([XL,XU])
    Y = np.concatenate([YL,YU])

    return X,Y
