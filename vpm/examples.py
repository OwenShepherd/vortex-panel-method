import numpy as np

def example_NACA0012():
    """ Example plotting a NACA0012 at a few angles of attack.

    Will plot the pressure coefficients of a NACA0012 at three angles of attack
    (5, 2, 10) and each at a different color.

    """

    name = 'NACA0012'
    c = 1
    N = 10
    alpha = 5
    NACA0012 = Airfoil(name,c,N,alpha)
    NACA0012.set_num_samples(16)
    NACA0012.set_num_samples(18)

def example_noplot():
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

if __name__=="__main__":
    example_NACA0012()
