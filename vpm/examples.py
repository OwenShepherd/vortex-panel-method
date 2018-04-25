import os.path
import sys
import numpy as np
import matplotlib as plt
sys.path.insert(0,os.path.abspath(sys.path[0]))

from vpm.vpmfuncs import *

def example_NACA0012():
    """
    Function: example_NACA0012

    Purpose:
        This function demonstrates how to plot a NACA0012 and its pressure
        coefficients at 3 different angles of attack, ensuring that the plots
        have a different color each time.

    Parameters:
        none

    Returns:
        none
    """

    name = 'NACA0012'
    plotColor = '#00ff08'
    coeffFIG = plt.figure()
    c = 1
    N = 100
    alpha = 5
    plotColor = '#00ff08'
    NACA0012 = Airfoil(name,c,N,alpha,coeffFIG,plotColor)
    NACA0012.Get_ParsedData()
    NACA0012.Get_AirfoilCoordinates()
    NACA0012.Get_PanelCoefficients(True)
    alpha = 2
    plotColor = '#ff1f00'
    NACA0012.Get_ParsedData()
    NACA0012.Set_AngleOfAttack(alpha)
    NACA0012.Set_PlotColor(plotColor)
    NACA0012.Get_PanelCoefficients(True)
    alpha = 10
    plotColor = '#0074ff'
    NACA0012.Get_ParsedData()
    NACA0012.Set_AngleOfAttack(alpha)
    NACA0012.Set_PlotColor(plotColor)
    NACA0012.Get_PanelCoefficients(True)
    plt.show()

def example_NoPlot():
    """
    Function: example_NoPlot

    Purpose:
        This function demonstrates how to use the airfoil class without
        plotting anything.

    Parameters:
        none

    Returns:
        none
    """

    name = 'NACA0012'
    c = 1
    N = 100
    alpha = 5
    NACA0012 = Airfoil(name,c,N,alpha)
    NACA0012.Get_ParsedData()
    NACA0012.Get_AirfoilCoordinates()
    NACA0012.Get_PanelCoefficients()
