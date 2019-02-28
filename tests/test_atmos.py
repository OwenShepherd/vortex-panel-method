from .context import aerotools
from aerotools import atmos

def test_temperature():
    temperature = atmos.get_temperature(11000)
    assert round(temperature) == round(216.65)
