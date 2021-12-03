from Modules.CTMinus import type_da
from Tangles.Tangle import ETangle


# - 109/10 on high performance mode
# - 106/10 on high performance mode, plugged in
def test_da_speed():
    # cup = ETangle(ETangle.Type.OVER, (-1, 1, -1), 1)
    cup = ETangle(ETangle.Type.OVER, (-1, 1), 1)
    for _ in range(1):
        da = type_da(cup)
    assert True
