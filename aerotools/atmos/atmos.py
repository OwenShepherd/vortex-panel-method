import pdb

def get_pressure():
    pass

def get_temperature(geopotential_height):
    """
    """
    # Temperature height profile information from 0 to 86 geometric km
    Lmb = [-6.5, 0, 1.0, 2.8, 0, -2.8, -2.0] # Molecular-scale temperature gradient; Lm,b [K/km']
    Hb = [0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0, 84852.0] # Geopotential height [km']
    Hb = Hb*1000
    Tb = [288.15] # Sea-level kinetic temperature [K]
    Tm = 0

    Lmb = [i/1000 for i in Lmb]



    for index, layer in enumerate(Hb):
        if geopotential_height <= Hb[index+1]:
            Tm = temperature(layer, Tb[-1], Lmb[index], geopotential_height)
            break

        Tb.append(temperature(Hb[index], Tb[-1], Lmb[index], Hb[index+1]))

    return Tm

def temperature(layer_height, layer_temp, layer_gradient, geopotential_height):
    return layer_temp + layer_gradient * (geopotential_height - layer_height)
