from __future__ import division
import math
import numpy as np
import matplotlib.pyplot as plt
import pdb
from vpm.aeroFuncs import *


class Airfoil:
    """ Airfoil handles calculations made on 4-digit series NACA airfoil.

    The Airfoil class calculates the x-y coordinates of all boundary points on
    an NACA 4-digit series airfoil.  In addition, plotting of the pressure
    coefficients is also an option.

    Parameters:
        NACA_ID: The NACA 4-digit series airfoil name i.e. "NACA0012"
        chord: The chord length of the airfoil
        NUM_SAMPLES: The number of samples / panels considered
        angle_of_attack: The angle of attack of the airfoil
        Figure_To_Plot: The figure handle from matplotlib to plot on
        PlotColor: The current color handle for matplotlib to plot with
        HASPLOT: Keeps track of whether each fig has been plotted yet
        max_camber: The maximum camber of the airfoil i.e "0/100, 1/100"
        position_maxCamber: The position of the max. camber i.e. "4/100"
        thickness: The maximum thickness of the airfoil i.e. "12/100", "08/100"
        BoundaryPoints_X: The x-locations of each boundary point on the airfoil
        BoundaryPoints_Y: The y-locations of each boundary point on the airfoil
        full_coefficientLift: The coefficient of lift of the airfoil
        pressure_coefficient: The pressure coefficient at each (x,y) point

    """


    def __init__(self,NACA_Name,Chord_Length,NUM_SAMPLES,Angle_Of_Attack,
                 FigID=None,plotColor=None):
        """ The Airfoil class constructor.

        This function will convert the standard name of a naca airfoil
        (ex: "NACA0012" or "NACA1408", etc) into its corresponding parameters.

        """
        self.NACA_ID = NACA_Name
        self.chord = Chord_Length
        self.NUM_SAMPLES = NUM_SAMPLES
        self.angle_of_attack = Angle_Of_Attack
        self.Figure_To_Plot = FigID
        self.PlotColor = plotColor
        self.HASPLOT = {}
        if (FigID != None):
            self.HASPLOT[str(self.Figure_To_Plot.number)] = False
        self.max_camber = None
        self.position_maxCamber = None
        self.thickness = None
        self.BoundaryPoints_X = None
        self.BoundaryPoints_Y = None
        self.full_coefficientLift = None
        self.pressure_coefficient = None


    def Set_Figure(self,FigID):
        """ Sets the figure to plot on.

        Changes the figure identifier for the class.

        Args:
            FigID: The matplotlib figure identifier

        """
        self.Figure_To_Plot = FigID
        if (self.Figure_To_Plot.number not in self.HASPLOT):
            self.HASPLOT[str(self.Figure_To_Plot.number)] = False

    def Set_PlotColor(self,plotColor):
        """ Sets the color of the pressure coefficient plot to use.

        Args:
            plotColor: The new matplotlib plot color

        """

        self.PlotColor = plotColor

    def Set_AngleOfAttack(self,alpha):
        """ Sets the angle of attack to use for next calculations.

        Args:
            alpha: The new angle of attack

        """

        self.angle_of_attack = alpha

    def Get_ParsedData(self):
        """ Converts a 'NACAXXX' string into its parameters.

        This function will convert the standard name of a naca airfoil
        (ex: "NACA0012" or "NACA1408", etc) into its corresponding parameters.

        Args:
            NACA_ID: The standard name of a NACA airfoil

        """
        name = self.NACA_ID
        airfoil_parameters = []

        for i in range(4,len(name)-1):
            if (i<=5):
                airfoil_parameters.append(float(name[i]))
            else:
                airfoil_parameters.append(float(name[i] + name[i+1]))


        self.max_camber = airfoil_parameters[0]/100
        self.position_maxCamber = airfoil_parameters[1]/100
        self.thickness = airfoil_parameters[2]/100

    def Get_AirfoilCoordinates(self):
        """ Calculates the points on a NACA 4-digit series airfoil.

        Args:
            max_camber: The maximum camber
            position_maxCamber: The location of the maximum camber
            thickness: The thickness
            chord: The chord length
            NUM_SAMPLES: The number of employed panels

        """

        m = self.max_camber
        p = self.position_maxCamber
        t = self.thickness
        c = self.chord
        N = self.NUM_SAMPLES

        # Splitting the airfoil into equal-length x-locations ALA K&C
        dt = 2*math.pi/N

        # Creating some control values and storage variables
        NOTDONE = True
        th = [0] * int(N/2+1)
        th[0] = 0
        i = 1

        # Want angles to go from zero to pi and split the boundary x-y locations
        # into upper and lower locations
        while (NOTDONE):
            value = th[i-1] + dt

            if (abs(value-math.pi)<=0.0000001):
                th[i] = value
                NOTDONE = False
            else:
                th[i] = value
                i = i+1

        # Now here we determine the x-lcations based off of the thetas
        R = 0.5*c
        xb = [R+R*math.cos(k) for k in th]

        # Now we shall determine the points of interest on the airfoil (the y-
        # locations)
        yt = [t*c/0.2*
              (0.2969*pow((k/c),0.5)-
              0.1260*(k/c)-
              0.3516*pow((k/c),2)+
              0.2843*pow((k/c),3)-
              0.1036*pow(k/c,4)) for k in xb]

        # Here we can determine the mean camber line
        yc = np.zeros((1,len(xb)))
        dyc = np.zeros((1,len(xb)))

        if ((p!=0) or (m!=0)):
            for i in range(len(xb)):
                if (xb[i]<=(p*c)):
                    yc[0,i] = (m*xb[i]/(pow(p,2))
                               * (2*p-xb[i]/c))

                    dyc[0,i] = (m*xb[i]/(pow(p,2))
                                * (-1/float(c))
                                + m/(pow(p,2))
                                * (2*p-xb[i]/c))

                else:
                    yc[0,i] = (m*(c-xb[i])
                               / ((1-p)**2)
                               * (1+xb[i]/c-2*p))

                    dyc[0,i] = ((m*(c-xb[i])/(1-p*p))
                                * 1/float(c)
                                + (-1*m)/(1-p*p)
                                * (1+xb[i]/c-2*p))


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



        self.BoundaryPoints_X = X
        self.BoundaryPoints_Y = Y


    def Get_PanelCoefficients(self,PLOT=False):
        """ Forumlates the system of equations for the vortex panel method.

        Args:
            BoundaryPoints_X: The dimensionless boundary x locations on the
                airfoil
            BoundaryPoints_Y: The dimensionless boundary y locations on the
                airfoil
            NUM_SAMPLES: The total number of panels
            angle_of_attack: The angle of attack
            NACA_ID: The string that specifies the 4-digit NACA airfoil, for
                example: 0012
            PLOT: Whether or not to plot Cp vs. x/c; if the plot is desired,
                PLOT == True
            PlotColor: The matplotlib color handle

        """
        XB = self.BoundaryPoints_X
        YB = self.BoundaryPoints_Y
        M = self.NUM_SAMPLES
        alpha = self.angle_of_attack
        NACA = self.NACA_ID[4:]
        FigID = self.Figure_To_Plot
        plotColor = self.PlotColor



        alphad = alpha
        alpha = math.radians(alpha)

        MP1 = M+1 # The trailing edge requries an extra point

        X = np.zeros(M)
        Y = np.zeros(M)
        RHS = np.zeros(MP1)
        theta = np.zeros(M)
        S = np.zeros(M)

        for i in range(M):
            IP1 = i + 1

            X[i] = 0.5 * (XB[i] + XB[IP1])
            Y[i] = 0.5 * (YB[i] + YB[IP1])

            S[i] = math.sqrt(pow(XB[IP1] - XB[i],2) + pow(YB[IP1] - YB[i],2))

            theta[i] = math.atan2(YB[IP1] - YB[i], XB[IP1] - XB[i])

            RHS[i] = math.sin(theta[i] - alpha)

        CN1 = np.zeros((M,M))
        CN2 = np.zeros((M,M))
        CT1 = np.zeros((M,M))
        CT2 = np.zeros((M,M))
        An = np.zeros((MP1,MP1))
        At = np.zeros((MP1,MP1))


        for i in range(M):
            for j in range(M):

                if (i == j):
                    CN1[i,j] = -1
                    CN2[i,j] = 1
                    CT1[i,j] = 0.5 * math.pi
                    CT2[i,j] = 0.5 * math.pi
                else:
                    A = (-1*(X[i]-XB[j])
                         * math.cos(theta[j])
                         - (Y[i]-YB[j])
                         * math.sin(theta[j]))

                    B = ((X[i]-XB[j])**2
                         + (Y[i]-YB[j])**2)

                    C = math.sin(theta[i]-theta[j])
                    D = math.cos(theta[i]-theta[j])

                    E = ((X[i]-XB[j])
                         * math.sin(theta[j])
                         - (Y[i]-YB[j])
                         * math.cos(theta[j]))

                    F = math.log(1+S[j]*(S[j]+2*A)/B)
                    G = math.atan2((E*S[j]),(B+A*S[j]))

                    P = ((X[i]-XB[j])
                         * math.sin(theta[i]
                                    - 2*theta[j])
                         + (Y[i]-YB[j])
                         * math.cos(theta[i]
                                    - 2*theta[j]))

                    Q = ((X[i]-XB[j])
                         * math.cos(theta[i]
                                    - 2*theta[j])
                         - (Y[i]-YB[j])
                         * math.sin(theta[i]
                                    - 2*theta[j]))

                    CN2[i,j] = D + 0.5 * Q * F/S[j] - (A*C + D*E) * G/S[j]
                    CN1[i,j] = 0.5*D*F + C*G - CN2[i,j]
                    CT2[i,j] = C + 0.5*P*F/S[j] + (A*D - C*E) * G/S[j]
                    CT1[i,j] = 0.5*C*F - D*G - CT2[i,j]

        for i in range(M):
            An[i,0] = CN1[i,0]
            An[i,MP1-1] = CN2[i,M-1]
            At[i,0] = CT1[i,0]
            At[i,MP1-1] = CT2[i,M-1]

            for j in range(1,M):
                An[i,j] = CN1[i,j] + CN2[i,(j-1)]
                At[i,j] = CT1[i,j] + CT2[i,(j-1)]

        # Trailing edge conditions
        An[MP1-1,0] = 1
        An[MP1-1,MP1-1] = 1

        for j in range(1,M):
            An[MP1-1,j] = 0

        RHS[MP1-1] = 0

        gamma = np.linalg.solve(An,RHS)
        V = np.zeros(M)
        Cp = np.zeros(M)


        for i in range(M):
            V[i] = math.cos(theta[i] - alpha)

            for j in range(MP1):
                V[i] = V[i] + At[i,j]*gamma[j]
                Cp[i] = 1 - (V[i])**2

        cl = Get_LiftCoefficients(V,S,M)

        CpLower = Cp[0:int(M/2)]
        CpUpper = Cp[int(M/2):]
        Cp = Cp.astype(np.float)
        newcl = np.array(cl)
        newcl = newcl.astype(np.float)

        if (PLOT):
            self.HASPLOT = Plot_PressureCoefficients(X,Y,alphad,CpUpper,CpLower,
                                                     M,NACA,FigID,plotColor,
                                                     self.HASPLOT)

        self.full_coefficientLift = newcl
        self.pressure_coefficient = Cp
