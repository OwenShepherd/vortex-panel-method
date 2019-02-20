import numpy as np
import matplotlib.pyplot as plt
from context import aerotools
from aerotools import vpm

def example_NACA0012():
    """ Example plotting a NACA0012 at a few angles of attack.

    Will plot the pressure coefficients of a NACA0012 at three angles of attack
    (-5, 0, 5, 10) and each at a different color.

    """

    # Initialization
    name = 'NACA0012'
    c = 1
    N = 100
    alphas = [-5, 0, 5, 10]
    Cp = []
    NACA0012 = vpm.Airfoil(name,c,N,0)
    colors = ['b', 'g', 'm', 'c', 'r']

    # Setting up plots
    plt.figure()
    plt.plot(1.2,0,'|',color='#7d7d7d',visible=True, label='Upper Surface',markersize=7)
    plt.plot(1.2,0,'2',color='#7d7d7d',visible=True, label='Lower Surface',markersize=7)
    plt.grid()
    temp_gca = plt.gca()
    temp_gca.invert_yaxis()

    [X, Y] = NACA0012.get_panel_coordinates()
    xUpper = X[int(N/2):]
    xLower = X[0:int(N/2)]

    for angle in alphas:
        NACA0012.set_angle_of_attack(angle)
        pressure = NACA0012.get_pressure_coefficients()
        pressureUpper = pressure[int(N/2):]
        pressureLower = pressure[0:int(N/2)]
        plot_Label = name + " AOA: " + str(angle) + "$^\circ$"
        angle_color = colors.pop()

        plt.plot(xUpper,pressureUpper,'-|',color=angle_color, label=plot_Label,markersize=7)
        plt.plot(xLower,pressureLower,'-2',color=angle_color, markersize=9)
        plt.xlim(-0.1,1.1)
        plt.legend()
        plt.xlabel("Dimensionless Chord Location [X/C]")
        plt.ylabel("Pressure Coefficient, Cp")

    plt.show()

def example_noplot():
    """ Demonstrates how to use the airfoil class without plotting.

    This function demonstrates how to use the airfoil class without
    plotting anything.

    Returns:
        numpy.ndarray: Numpy array containing the pressure coefficients at each
            calculated boundary point on the airfoil.

    """

    name = 'NACA2412'
    c = 1
    N = 100
    alpha = 1
    NACA0012 = vpm.Airfoil(name,c,N,alpha)

if __name__=="__main__":
    example_NACA0012()
