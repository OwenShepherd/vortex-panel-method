from __future__ import division
import math
import numpy as np
import matplotlib.pyplot as plt
import pdb

def Get_LiftCoefficients(V,S,M):
    """ Calculates the lift coefficient of the wing

    This function computes the coefficient of lift for an airfoil as to
    be used in the vortex panel method "Get_PanelCoefficients"

    Args:
        V: The dimensionlesss velocity at each control point
        S: The dimensionless length of each of the control points
        M: The number of panels

    Returns:
        float: The coefficient of lift

    """

    gamma = 0
    for j in range(M):
        gamma = gamma + V[j]*S[j]

    cl = (2*gamma)

    return cl

def Plot_PressureCoefficients(X,Y,alpha,CpUpper,CpLower,M,NACA,FigID,plotColor,HASPLOT):
    """ Calcultes and may plot the pressure coefficients of an airfoil.

    This function will calculate the pressure coefficient at each point on
    the chosen airfoil, and will also plot if chosen.

    Args:
        X: The x-locations of each boundary point on the airfoil
        Y: The y-locations of each boundary point on the airfoil
        alpha: The angle of attack of the airfoil
        CpUpper: The numpy array containing the upper pressure coefficient
            values
        CpLower: The numpy array containing the lower pressure coefficient
            values
        M: The number of samples / panels considered
        NACA: The NACA 4-digit series airfoil name i.e. "NACA0012"
        FigID: The figure handle from matplotlib to plot on
        PlotColor: The current color handle for matplotlib to plot with
        HASPLOT: Dictionary that keeps track of which figs have an existing plot.
            I.e: HASPLOT["1"] = True, if the fig with fig.number == 1 already has
            an existing plot

    Returns:
        dict: Updated version of the HASPLOT dictionary

    """

    plt.figure(FigID.number)

    if not HASPLOT[str(FigID.number)]:
        plt.plot(1.2,0,'|',color='#7d7d7d',visible=True,
                 label='Upper Surface',markersize=7)

        plt.plot(1.2,0,'2',color='#7d7d7d',visible=True,
                 label='Lower Surface',markersize=7)

        plt.grid(FigID)
        temp_gca = plt.gca()
        temp_gca.invert_yaxis()
        HASPLOT[str(FigID.number)] = True

    x_upper = X[int(M/2):]
    x_lower = X[0:int(M/2)]

    plot_Label = NACA + " AOA: " + str(alpha) + "$^\circ$"


    plt.plot(x_upper,CpUpper,'-|',
             color=plotColor, label=plot_Label,markersize=7)

    plt.plot(x_lower,CpLower,'-2',
             color=plotColor,markersize=9)

    plt.xlim(-0.1,1.1)
    plt.legend()
    plt.xlabel("Dimensionless Chord Location [X/C]")
    plt.ylabel("Pressure Coefficient, Cp")

    return HASPLOT
