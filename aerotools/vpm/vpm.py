import math
import numpy as np

class Airfoil:
    """ Airfoil handles calculations made on 4-digit series NACA airfoil.

    The Airfoil class calculates the x-y coordinates of all boundary points on
    an NACA 4-digit series airfoil.

    Parameters:
        NACA_ID: The NACA 4-digit series airfoil name i.e. "NACA0012"
        chord: The chord length of the airfoil
        NUM_SAMPLES: The number of samples / panels considered
        angle_of_attack: The angle of attack of the airfoil
        max_camber: The maximum camber of the airfoil i.e "0/100, 1/100"
        position_max_camber: The position of the max. camber i.e. "4/100"
        thickness: The maximum thickness of the airfoil i.e. "12/100", "08/100"
        x_boundary_points: The x-locations of each boundary point on the airfoil
        y_boundary_points: The y-locations of each boundary point on the airfoil
        full_coefficient_lift: The coefficient of lift of the airfoil
        pressure_coefficient: The pressure coefficient at each (x,y) point

    """

############################## Constructor #####################################

    def __init__(self,NACA_Name,Chord_Length=1,NUM_SAMPLES=100,Angle_Of_Attack=0
                ):
        """ The Airfoil class constructor.

        This function will convert the standard name of a naca airfoil
        (ex: "NACA0012" or "NACA1408", etc) into its corresponding parameters.

        """
        self.x_boundary_points = None
        self.y_boundary_points = None
        self.x_panel_points = None
        self.y_panel_points = None
        self.NUM_SAMPLES = NUM_SAMPLES
        self.NACA_ID = NACA_Name
        self.chord = Chord_Length
        self.angle_of_attack = Angle_Of_Attack
        self.max_camber = None
        self.position_max_camber = None
        self.thickness = None
        self.full_coefficient_lift = None
        self.pressure_coefficient = None

        self._calculate_parameters()

####################### Private Functions ##############################
    def _mean_camber_and_gradient(self, x):
        # Here we can determine the mean camber line
        yc = []
        dyc = []
        yc_value = 0
        dyc_value = 0

        if ((self.position_max_camber!=0) or (self.max_camber!=0)):
            for i in range(len(x)):
                if (x[i]<=(self.position_max_camber*self.chord)):
                    yc_value = (self.max_camber*x[i]/(pow(self.position_max_camber,2))
                               * (2*self.position_max_camber-x[i]/self.chord))

                    dyc_value = (self.max_camber*x[i]/(pow(self.position_max_camber,2))
                                * (-1/float(self.chord))
                                + self.max_camber/(pow(self.position_max_camber,2))
                                * (2*self.position_max_camber-x[i]/self.chord))

                else:
                    yc_value = (self.max_camber*(self.chord-x[i])
                               / ((1-self.position_max_camber)**2)
                               * (1+x[i]/self.chord-2*self.position_max_camber))

                    dyc_value = ((self.max_camber*(self.chord-x[i])/(1-self.position_max_camber*self.position_max_camber))
                                * 1/float(self.chord)
                                + (-1*self.max_camber)/(1-self.position_max_camber*self.position_max_camber)
                                * (1+x[i]/self.chord-2*self.position_max_camber))
                yc.append(yc_value)
                dyc.append(dyc_value)

        else:
            yc = [0] * len(x)
            dyc = yc

        return yc, dyc

    def _x_coordinates(self):
        """ Determines the x-locations by splitting an airfoil into even arcs.

        Returns:
            float[]: the x-coordinates of the points of interest.
        """
        dt = 2*math.pi/self.NUM_SAMPLES
        th = []
        th.append(0)

        for i in range(self.NUM_SAMPLES//2):
            value = th[i] + dt
            th.append(value)

        R = 0.5*self.chord
        x_coordinates = [R+R*math.cos(k) for k in th]

        return x_coordinates

    def _thickness_distribution(self, x_coordinates):
        """ Determines the thickness distribution of the airfoil.

        Returns:
            float[]: y-location of the thickness distribution
        """
        # Now we shall determine the points of interest on the airfoil (the y-
        # locations)
        yt = [self.thickness*self.chord/0.2*
              (0.2969*pow((k/self.chord),0.5)-
              0.1260*(k/self.chord)-
              0.3516*pow((k/self.chord),2)+
              0.2843*pow((k/self.chord),3)-
              0.1036*pow(k/self.chord,4)) for k in x_coordinates]

        return yt

    def _calculate_parameters(self):
        """ Sets the parsed data, calculates the geometric shape, and then
        produces the coefficients of pressure and the coefficient of lift for
        the airfoil.
        """

        self._parsed_data()
        self._airfoil_coordinates()
        self._panel_coefficients()

    def _parsed_data(self):
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

    def _airfoil_coordinates(self):
        """ Calculates the points on a NACA 4-digit series airfoil.

        Args:
            max_camber: The maximum camber
            position_max_camber: The location of the maximum camber
            thickness: The thickness
            chord: The chord length
            NUM_SAMPLES: The number of employed panels

        """
        x_coordinates = self._x_coordinates()

        yt = self._thickness_distribution(x_coordinates)

        yc, dyc = self._mean_camber_and_gradient(x_coordinates)

        zeta = [math.atan(k) for k in dyc]

        XU = []
        YU = []
        XL = []
        YL = []

        # Creating the upper and lower x,y locations
        for i in range(len(x_coordinates)):
            XU.append(x_coordinates[i])
            YU.append(yc[i]+yt[i]*math.cos(zeta[i]))

            XL.append(x_coordinates[i])
            YL.append(yc[i]-yt[i]*math.cos(zeta[i]))

        # This ensures that the boundary points go from the trailing edge, clockwise
        # back to the trailing edge with two boundary points at the trailing edge
        # and one at the leading edge
        del XU[-1]
        del YU[-1]
        XU.reverse()
        YU.reverse()

        X = XL+XU
        Y = YL+YU

        self.x_boundary_points = X
        self.y_boundary_points = Y



    def _panel_coefficients(self):
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
        """
        x_coordinates = self.x_boundary_points
        y_coordinates = self.y_boundary_points
        number_panels = self.NUM_SAMPLES
        NACA_4digit = self.NACA_ID[4:]
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

        cl = self._lift_coefficients(V,S,self.NUM_SAMPLES)
        CpLower = Cp[0:int(self.NUM_SAMPLES/2)]
        CpUpper = Cp[int(self.NUM_SAMPLES/2):]
        Cp = Cp.astype(np.float)
        newcl = np.array(cl)
        newcl = newcl.astype(np.float)

        self.x_panel_points = X
        self.y_panel_points = Y

        self.full_coefficient_lift = newcl
        self.pressure_coefficient = Cp


    def _lift_coefficients(self,V,S,M):
        """
        Function: Calculate_LiftCoefficients

        Purpose: This function computes the coefficient of lift for an airfoil as to
        be used in the vortex panel method "Calculate_PanelCoefficients"

        Args:
            V: The dimensionlesss velocity at each control point
            S: The dimensionless length of each of the control points
            M: The number of panels

        Returns:
            float: The coefficient of lift, Cl.
        """

        gamma = 0
        for j in range(M):
            gamma = gamma + V[j]*S[j]
            cl = 2*gamma

        return cl

######################## Public Functions ######################################

    def set_angle_of_attack(self,angle):
        """ Sets the angle of attack to use for next calculations.

        Args:
            angle: The new angle of attack
        """

        self.angle_of_attack = angle
        self._calculate_parameters()

    def set_chord_length(self, length):
        """ Sets the chord length used in the calculations.

        Args:
            length: The chord length of the airfoil.
        """
        self.chord = length
        self._calculate_parameters()

    def set_num_samples(self, samples):
        """ Sets the number of panels used for sampling the airfoil.

        Args:
            samples: The number of samples / panels used for airfoil calculations.
        """
        self.NUM_SAMPLES = samples
        self._calculate_parameters()

    def set_airfoil(self, NACA_ID):
        """Sets the airfoil type used.

        Args:
            NACA_ID: The 4-digit series airfoil name.  Example: 'NACA0012'
        """

        self.NACA_ID = NACA_ID
        self._calculate_parameters()

    def get_airfoil_coordinates(self):
        """ Returns the coordinates of the airfoil.

        Returns:
            tuple: An array of x-coordinates and y-coordinates of boundary
            points.  [X,Y]
        """
        return self.x_boundary_points, self.y_boundary_points

    def get_panel_coordinates(self):
        """ Returns the coordinates of the midpoints of the panels.

        Returns:
            tuple: An array of x-coordinates and y-coordiantes of the panel mid-
            points.  [X, Y]
        """

        return self.x_panel_points, self.y_panel_points

    def get_coefficient_lift(self):
        """ Returns the coefficient of lift.

        Returns:
            float: The airfoil's coefficient of lift (per meter span), Cl.
        """
        return self.full_coefficient_lift

    def get_pressure_coefficients(self):
        """ Returns the pressure coefficient at the midpoint of each panel.

        Returns:
            float[]: Pressure coefficient at each boundary point, Cp.
        """
        return self.pressure_coefficient
