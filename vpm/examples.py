import os.path
import sys
import numpy as np
import matplotlib as plt
import pdb
RootPath = os.path.abspath(os.path.join(sys.path[0],".."))
sys.path.insert(0,os.path.abspath(RootPath))

from vpm.vpmfuncs import *


def example_NACA0012():
    """ Example plotting a NACA0012 at a few angles of attack.

    Will plot the pressure coefficients of a NACA0012 at three angles of attack
    (5, 2, 10) and each at a different color.

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
    """ Demonstrates how to use the airfoil class without plotting.

    This function demonstrates how to use the airfoil class without
    plotting anything.

    Returns:
        numpy.ndarray: Numpy array containing the pressure coefficients at each
            calculated boundary point on the airfoil.

    """

    name = 'NACA0012'
    c = 1
    N = 100
    alpha = 5
    NACA0012 = Airfoil(name,c,N,alpha)
    NACA0012.Get_ParsedData()
    NACA0012.Get_AirfoilCoordinates()
    NACA0012.Get_PanelCoefficients()
    return NACA0012.pressure_coefficient


if __name__=="__main__":
    FunctionCall = sys.argv[1]
    possibles = globals().copy()
    possibles.update(locals())
    method = possibles.get(FunctionCall)
    if not method:
        raise NotImplementedError("Method %s not implemented." %FunctionCall)
    method()
