from Modules.CTMinus import *
from Modules.Module import *
from Modules.ETangleStrands import *
import timeit

from Tangles import TangleRenderer


# - 130/10
def test_da_speed():
    cup = ETangle(ETangle.Type.OVER, (-1, 1), 1)
    for _ in range(5):
        da = type_da(cup)
    assert True
