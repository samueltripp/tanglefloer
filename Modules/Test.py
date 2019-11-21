from Tangles.Tangle import *
from Modules.CTMinus import *

import unittest


class TestCTMinus(unittest.TestCase):

    def test_dplus(self):
        over1 = ETangle(ETangle.Type.OVER, (1, 1, 1, 1), 2)
        sd_over1_1 = StrandDiagram(over1, {1: 4, 2: 3}, {0: 0, 1: 1, 2: 2})
        sd_over1_1_out = Bimodule.Element(
            {StrandDiagram(over1, {1: 3, 2: 4}, {0: 0, 1: 1, 2: 2}): over1.polyring.one()})
        sd_over1_2 = StrandDiagram(over1, {1: 3, 4: 0}, {1: 1, 2: 2, 4: 4})
        sd_over1_2_out = Bimodule.Element(
            {StrandDiagram(over1, {1: 0, 4: 3}, {1: 1, 2: 2, 4: 4}): over1.polyring['U2'] * over1.polyring['U3']})
        over2 = ETangle(ETangle.Type.OVER, (1, -1, 1, 1), 2)
        sd_over2_1 = StrandDiagram(over2, {1: 4, 2: 3}, {0: 0, 1: 1, 2: 2})
        sd_over2_1_out = Bimodule.Element(
            {StrandDiagram(over2, {1: 3, 2: 4}, {0: 0, 1: 1, 2: 2}): over2.polyring.one()})
        sd_over2_2 = StrandDiagram(over2, {1: 3, 4: 0}, {1: 1, 2: 2, 4: 4})
        sd_over2_2_out = Bimodule.Element(
            {StrandDiagram(over2, {1: 0, 4: 3}, {1: 1, 2: 2, 4: 4}): over2.polyring.zero()})

        under1 = ETangle(ETangle.Type.UNDER, (1, 1, 1, 1), 2)
        sd_under1_1 = StrandDiagram(under1, {2: 2, 0: 4}, {0: 0, 1: 1, 3: 3})
        sd_under1_1_out = Bimodule.Element(
            {StrandDiagram(under1, {0: 2, 2: 4}, {0: 0, 1: 1, 3: 3}): under1.polyring['U3']}
        )
        sd_under1_2 = StrandDiagram(under1, {2: 1, 0: 4}, {0: 0, 2: 2, 3: 3})
        sd_under1_2_out = Bimodule.Element(
            {StrandDiagram(under1, {0: 1, 2: 4}, {0: 0, 2: 2, 3: 3}): under1.polyring['U3']}
        )
        under2 = ETangle(ETangle.Type.UNDER, (1, -1, 1, 1), 2)
        sd_under2_1 = StrandDiagram(under2, {3: 1, 0: 4}, {0: 0, 2: 2, 3: 3})
        sd_under2_1_out = Bimodule.Element(
            {StrandDiagram(under2, {0: 1, 3: 4}, {0: 0, 2: 2, 3: 3}): under2.polyring.zero()}
        )

        cup1 = ETangle(ETangle.Type.CUP, (1, 1, -1, 1), 2)
        sd_cup1_1 = StrandDiagram(cup1, {1: 4, 2: 3}, {0: 0, 1: 1})
        sd_cup1_1_out = Bimodule.Element(
            {StrandDiagram(cup1, {1: 3, 2: 4}, {0: 0, 1: 1}): cup1.polyring['U4']}
        )
        sd_cup1_2 = StrandDiagram(cup1, {0: 4, 2: 0}, {3: 3, 1: 1})
        sd_cup1_2_out = Bimodule.Element(
            {StrandDiagram(cup1, {0: 0, 2: 4}, {3: 3, 1: 1}): cup1.polyring['U1'] * cup1.polyring['U4']}
        )

        self.assertEqual(dplus(sd_over1_1), sd_over1_1_out)
        self.assertEqual(dplus(sd_over1_2), sd_over1_2_out)
        self.assertEqual(dplus(sd_over2_1), sd_over2_1_out)
        self.assertEqual(dplus(sd_over2_2), sd_over2_2_out)

        self.assertEqual(dplus(sd_under1_1), sd_under1_1_out)
        self.assertEqual(dplus(sd_under1_2), sd_under1_2_out)
        self.assertEqual(dplus(sd_under2_1), sd_under2_1_out)

        self.assertEqual(dplus(sd_cup1_1), sd_cup1_1_out)
        self.assertEqual(dplus(sd_cup1_2), sd_cup1_2_out)

    def test_dminus(self):
        under1 = ETangle(ETangle.Type.UNDER, (-1, -1, -1, -1), 2)
        sd_under1_1 = StrandDiagram(under1, {0: 0, 1: 1, 3: 3}, {2: 2, 4: 4})
        sd_under1_1_out = Bimodule.Element(
            {StrandDiagram(under1, {0: 0, 1: 1, 3: 3}, {2: 4, 4: 2}): under1.polyring['U3'] * under1.polyring['U4']}
        )

        self.assertEqual(dminus(sd_under1_1), sd_under1_1_out)


if __name__ == '__main__':
    unittest.main()
