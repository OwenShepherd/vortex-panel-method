from .context import aerotools
from aerotools import atmos

def test_temperature():
    temperature = atmos.get_temperature(81000)
    assert temperature == 287.51
