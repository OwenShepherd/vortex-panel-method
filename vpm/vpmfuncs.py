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
        figure_to_plot: The figure handle from matplotlib to plot on
        plot_color: The current color handle for matplotlib to plot with
        HASPLOT: Keeps track of whether each fig has been plotted yet
        max_camber: The maximum camber of the airfoil i.e "0/100, 1/100"
        position_max_camber: The position of the max. camber i.e. "4/100"
        thickness: The maximum thickness of the airfoil i.e. "12/100", "08/100"
        x_boundary_points: The x-locations of each boundary point on the airfoil
        y_boundary_points: The y-locations of each boundary point on the airfoil
        full_coefficient_lift: The coefficient of lift of the airfoil
        pressure_coefficient: The pressure coefficient at each (x,y) point

    """


    def __init__(self,NACA_Name,Chord_Length,NUM_SAMPLES,Angle_Of_Attack,
                 figure_ID=None,plot_color=None):
        """ The Airfoil class constructor.

        This function will convert the standard name of a naca airfoil
        (ex: "NACA0012" or "NACA1408", etc) into its corresponding parameters.

        """
        self.x_boundary_points = None
        self.y_boundary_points = None
        self.NUM_SAMPLES = NUM_SAMPLES
        self.NACA_ID = NACA_Name
        self.chord = Chord_Length
        self.angle_of_attack = Angle_Of_Attack
        self.figure_to_plot = figure_ID
        self.plot_color = plot_color
        self.HASPLOT = {}
        if (figure_ID != None):
            self.HASPLOT[str(self.figure_to_plot.number)] = False
        self.max_camber = None
        self.position_max_camber = None
        self.thickness = None
        self.full_coefficient_lift = None
        self.pressure_coefficient = None


    def set_figure(self,figure_ID):
        """ Sets the figure to plot on.

        Changes the figure identifier for the class.

        Args:
            figure_ID: The matplotlib figure identifier

        """
        self.figure_to_plot = figure_ID
        if (self.figure_to_plot.number not in self.HASPLOT):
            self.HASPLOT[str(self.figure_to_plot.number)] = False

    def set_plot_color(self,plot_color):
        """ Sets the color of the pressure coefficient plot to use.

        Args:
            plot_color: The new matplotlib plot color

        """

        self.plot_color = plot_color

    def set_angle_of_attack(self,angle):
        """ Sets the angle of attack to use for next calculations.

        Args:
            self.angle_of_attack: The new angle of attack

        """

        self.angle_of_attack = angle

    def get_parsed_data(self):
        """ Converts a 'NACAXXX' string into its parameters.

        This function will convert the standard name of a naca airfoil
        (ex: "NACA0012" or "NACA1408", etc) into its corresponding parameters.

        Args:
            NACA_ID: The standard name of a NACA airfoil

        """
        airfoil_parameters = []

        for i in range(4,len(self.NACA_ID)-1):
            if (i<=5):
                airfoil_parameters.append(float(self.NACA_ID[i]))
            else:
                airfoil_parameters.append(float(self.NACA_ID[i] + self.NACA_ID[i+1]))


        self.max_camber = airfoil_parameters[0]/100
        self.position_max_camber = airfoil_parameters[1]/100
        self.thickness = airfoil_parameters[2]/100

    def get_airfoil_coordinates(self):
        """ Calculates the points on a NACA 4-digit series airfoil.

        Args:
            max_camber: The maximum camber
            position_max_camber: The location of the maximum camber
            thickness: The thickness
            chord: The chord length
            NUM_SAMPLES: The number of employed panels

        """

        # Splitting the airfoil into equal-length x-locations ALA K&C
        dt = 2*math.pi/self.NUM_SAMPLES

        # Creating some control values and storage variables
        NOTDONE = True
        th = [0] * int(self.NUM_SAMPLES/2+1)
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
        R = 0.5*self.chord
        x_coordinates = [R+R*math.cos(k) for k in th]

        # Now we shall determine the points of interest on the airfoil (the y-
        # locations)
        yt = [self.thickness*self.chord/0.2*
              (0.2969*pow((k/self.chord),0.5)-
              0.1260*(k/self.chord)-
              0.3516*pow((k/self.chord),2)+
              0.2843*pow((k/self.chord),3)-
              0.1036*pow(k/self.chord,4)) for k in x_coordinates]

        # Here we can determine the mean camber line
        yc = np.zeros((1,len(x_coordinates)))
        dyc = np.zeros((1,len(x_coordinates)))

        if ((self.position_max_camber!=0) or (self.max_camber!=0)):
            for i in range(len(x_coordinates)):
                if (x_coordinates[i]<=(self.position_max_camber*self.chord)):
                    yc[0,i] = (self.max_camber*x_coordinates[i]/(pow(self.position_max_camber,2))
                               * (2*self.position_max_camber-x_coordinates[i]/self.chord))

                    dyc[0,i] = (self.max_camber*x_coordinates[i]/(pow(self.position_max_camber,2))
                                * (-1/float(self.chord))
                                + self.max_camber/(pow(self.position_max_camber,2))
                                * (2*self.position_max_camber-x_coordinates[i]/self.chord))

                else:
                    yc[0,i] = (self.max_camber*(self.chord-x_coordinates[i])
                               / ((1-self.position_max_camber)**2)
                               * (1+x_coordinates[i]/self.chord-2*self.position_max_camber))

                    dyc[0,i] = ((self.max_camber*(self.chord-x_coordinates[i])/(1-self.position_max_camber*self.position_max_camber))
                                * 1/float(self.chord)
                                + (-1*self.max_camber)/(1-self.position_max_camber*self.position_max_camber)
                                * (1+x_coordinates[i]/self.chord-2*self.position_max_camber))


        # Now we can define parameter zeta
        zeta = [math.atan(k) for k in dyc[0]]



        XU = np.zeros((1,len(x_coordinates)))
        YU = np.zeros((1,len(x_coordinates)))
        XL = np.zeros((1,len(x_coordinates)))
        YL = np.zeros((1,len(x_coordinates)))

        # Creating the upper and lower x,y locations
        for i in range(len(x_coordinates)):
            XU[0,i] = (x_coordinates[i])
            YU[0,i] = (yc[0,i]+yt[i]*math.cos(zeta[i]))

            XL[0,i] = (x_coordinates[i])
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



        self.x_boundary_points = X
        self.y_boundary_points = Y


    def get_panel_coefficients(self,PLOT=False):
        """ Forumlates the system of equations for the vortex panel method.

        Args:
            x_boundary_points: The dimensionless boundary x locations on the
                airfoil
            y_boundary_points: The dimensionless boundary y locations on the
                airfoil
            NUM_SAMPLES: The total number of panels
            angle_of_attack: The angle of attack
            NACA_ID: The string that specifies the 4-digit NACA airfoil, for
                example: 0012
            PLOT: Whether or not to plot Cp vs. x/c; if the plot is desired,
                PLOT == True
            plot_color: The matplotlib color handle

        """
        x_coordinates = self.x_boundary_points
        y_coordinates = self.y_boundary_points
        number_panels = self.NUM_SAMPLES
        NACA_4digit = self.NACA_ID[4:]
        figure_ID = self.figure_to_plot
        plot_color = self.plot_color



        alphad = self.angle_of_attack
        self.angle_of_attack = math.radians(self.angle_of_attack)

        MP1 = number_panels + 1 # The trailing edge requries an extra point

        X = np.zeros(number_panels)
        Y = np.zeros(number_panels)
        RHS = np.zeros(MP1)
        theta = np.zeros(number_panels)
        S = np.zeros(number_panels)

        for i in range(number_panels):
            IP1 = i + 1

            X[i] = 0.5 * (x_coordinates[i] + x_coordinates[IP1])
            Y[i] = 0.5 * (y_coordinates[i] + y_coordinates[IP1])

            S[i] = math.sqrt(pow(x_coordinates[IP1] - x_coordinates[i],2) + pow(y_coordinates[IP1] - y_coordinates[i],2))

            theta[i] = math.atan2(y_coordinates[IP1] - y_coordinates[i], x_coordinates[IP1] - x_coordinates[i])

            RHS[i] = math.sin(theta[i] - self.angle_of_attack)

        CN1 = np.zeros((number_panels,number_panels))
        CN2 = np.zeros((number_panels,number_panels))
        CT1 = np.zeros((number_panels,number_panels))
        CT2 = np.zeros((number_panels,number_panels))
        An = np.zeros((MP1,MP1))
        At = np.zeros((MP1,MP1))


        for i in range(number_panels):
            for j in range(number_panels):

                if (i == j):
                    CN1[i,j] = -1
                    CN2[i,j] = 1
                    CT1[i,j] = 0.5 * math.pi
                    CT2[i,j] = 0.5 * math.pi
                else:
                    A = (-1*(X[i]-x_coordinates[j])
                         * math.cos(theta[j])
                         - (Y[i]-y_coordinates[j])
                         * math.sin(theta[j]))

                    B = ((X[i]-x_coordinates[j])**2
                         + (Y[i]-y_coordinates[j])**2)

                    C = math.sin(theta[i]-theta[j])
                    D = math.cos(theta[i]-theta[j])

                    E = ((X[i]-x_coordinates[j])
                         * math.sin(theta[j])
                         - (Y[i]-y_coordinates[j])
                         * math.cos(theta[j]))

                    F = math.log(1+S[j]*(S[j]+2*A)/B)
                    G = math.atan2((E*S[j]),(B+A*S[j]))

                    P = ((X[i]-x_coordinates[j])
                         * math.sin(theta[i]
                                    - 2*theta[j])
                         + (Y[i]-y_coordinates[j])
                         * math.cos(theta[i]
                                    - 2*theta[j]))

                    Q = ((X[i]-x_coordinates[j])
                         * math.cos(theta[i]
                                    - 2*theta[j])
                         - (Y[i]-y_coordinates[j])
                         * math.sin(theta[i]
                                    - 2*theta[j]))

                    CN2[i,j] = D + 0.5 * Q * F/S[j] - (A*C + D*E) * G/S[j]
                    CN1[i,j] = 0.5*D*F + C*G - CN2[i,j]
                    CT2[i,j] = C + 0.5*P*F/S[j] + (A*D - C*E) * G/S[j]
                    CT1[i,j] = 0.5*C*F - D*G - CT2[i,j]

        for i in range(number_panels):
            An[i,0] = CN1[i,0]
            An[i,MP1-1] = CN2[i,self.NUM_SAMPLES-1]
            At[i,0] = CT1[i,0]
            At[i,MP1-1] = CT2[i,self.NUM_SAMPLES-1]

            for j in range(1,self.NUM_SAMPLES):
                An[i,j] = CN1[i,j] + CN2[i,(j-1)]
                At[i,j] = CT1[i,j] + CT2[i,(j-1)]

        # Trailing edge conditions
        An[MP1-1,0] = 1
        An[MP1-1,MP1-1] = 1

        for j in range(1,self.NUM_SAMPLES):
            An[MP1-1,j] = 0

        RHS[MP1-1] = 0

        gamma = np.linalg.solve(An,RHS)
        V = np.zeros(number_panels)
        Cp = np.zeros(number_panels)


        for i in range(number_panels):
            V[i] = math.cos(theta[i] - self.angle_of_attack)

            for j in range(MP1):
                V[i] = V[i] + At[i,j]*gamma[j]
                Cp[i] = 1 - (V[i])**2

        cl = Get_LiftCoefficients(V,S,self.NUM_SAMPLES)

        CpLower = Cp[0:int(self.NUM_SAMPLES/2)]
        CpUpper = Cp[int(self.NUM_SAMPLES/2):]
        Cp = Cp.astype(np.float)
        newcl = np.array(cl)
        newcl = newcl.astype(np.float)

        if (PLOT):
            self.HASPLOT = Plot_PressureCoefficients(X,Y,alphad,CpUpper,CpLower,
                                                     self.NUM_SAMPLES,NACA_4digit,figure_ID,plot_color,
                                                     self.HASPLOT)

        self.full_coefficient_lift = newcl
        self.pressure_coefficient = Cp
