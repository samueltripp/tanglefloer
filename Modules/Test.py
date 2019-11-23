from Tangles.Tangle import *
from Modules.CTMinus import *

import unittest


class TestCTMinus(unittest.TestCase):

    def test_dplus(self):
        under1 = ETangle(ETangle.Type.UNDER, (1, 1, 1, 1), 2)
        sd_under1_1 = StrandDiagram(under1, {0: 0, 3: 3, 4: 4}, {1: 4, 2: 3})
        sd_under1_1_out = Bimodule.Element(
            {StrandDiagram(under1, {0: 0, 3: 3, 4: 4}, {1: 3, 2: 4}): under1.polyring.one()})
        sd_under1_2 = StrandDiagram(under1, {0: 0, 2: 2, 3: 3}, {1: 3, 4: 0})
        sd_under1_2_out = Bimodule.Element(
            {StrandDiagram(under1, {0: 0, 2: 2, 3: 3}, {1: 0, 4: 3}): under1.polyring['U2'] * under1.polyring['U3']})
        under2 = ETangle(ETangle.Type.UNDER, (1, -1, 1, 1), 2)
        sd_under2_1 = StrandDiagram(under2, {0: 0, 3: 3, 4: 4}, {1: 4, 2: 3})
        sd_under2_1_out = Bimodule.Element(
            {StrandDiagram(under2, {0: 0, 3: 3, 4: 4}, {1: 3, 2: 4}): under2.polyring.one()})
        sd_under2_2 = StrandDiagram(under2, {0: 0, 2: 2, 3: 3}, {1: 3, 4: 0})
        sd_under2_2_out = Bimodule.Element(
            {StrandDiagram(under2, {0: 0, 2: 2, 3: 3}, {1: 0, 4: 3}): under2.polyring.zero()})

        over1 = ETangle(ETangle.Type.OVER, (1, 1, 1, 1), 2)
        sd_over1_1 = StrandDiagram(over1, {4: 4, 1: 1, 3: 3}, {2: 2, 0: 4})
        sd_over1_1_out = Bimodule.Element(
            {StrandDiagram(over1, {4: 4, 1: 1, 3: 3}, {0: 2, 2: 4}): over1.polyring['U2']}
        )
        sd_over1_2 = StrandDiagram(over1, {1: 1, 4: 4, 3: 3}, {2: 1, 0: 4})
        sd_over1_2_out = Bimodule.Element(
            {StrandDiagram(over1, {1: 1, 4: 4, 3: 3}, {0: 1, 2: 4}): over1.polyring['U2']}
        )
        over2 = ETangle(ETangle.Type.OVER, (1, 1, -1, 1), 2)
        sd_over2_1 = StrandDiagram(over2, {1: 1, 0: 0, 4: 4}, {3: 1, 2: 4})
        sd_over2_1_out = Bimodule.Element(
            {StrandDiagram(over2, {1: 1, 0: 0, 4: 4}, {3: 4, 2: 1}): over2.polyring.zero()}
        )
        # Figure 9 from "An introduction..."
        over3 = ETangle(ETangle.Type.OVER, (1, 1, -1, -1), 2)
        sd_over3_1 = StrandDiagram(over3, {1: 2, 2: 1, 3: 4}, {0: 1, 3: 2})
        sd_over3_1_out = Bimodule.Element()

        cap1 = ETangle(ETangle.Type.CAP, (1, 1, -1, 1), 2)
        sd_cap1_1 = StrandDiagram(cap1, {0: 0, 1: 1}, {4: 1, 3: 2})
        sd_cap1_1_out = Bimodule.Element(
            {StrandDiagram(cap1, {0: 0, 1: 1}, {4: 2, 3: 1}): cap1.polyring['U4']}
        )
        sd_cap1_2 = StrandDiagram(cap1, {3: 3, 1: 1}, {4: 0, 0: 2})
        sd_cap1_2_out = Bimodule.Element(
            {StrandDiagram(cap1, {3: 3, 1: 1}, {4: 2, 0: 0}): cap1.polyring['U1'] * cap1.polyring['U4']}
        )

        self.assertEqual(sd_under1_1_out, dplus(sd_under1_1))
        self.assertEqual(sd_under1_2_out, dplus(sd_under1_2))
        self.assertEqual(sd_under2_1_out, dplus(sd_under2_1))
        self.assertEqual(sd_under2_2_out, dplus(sd_under2_2))

        self.assertEqual(sd_over1_1_out, dplus(sd_over1_1))
        self.assertEqual(sd_over1_2_out, dplus(sd_over1_2))
        self.assertEqual(sd_over2_1_out, dplus(sd_over2_1))
        self.assertEqual(sd_over3_1_out, dplus(sd_over3_1))

        self.assertEqual(sd_cap1_1_out, dplus(sd_cap1_1))
        self.assertEqual(sd_cap1_2_out, dplus(sd_cap1_2))

    def test_dminus(self):
        under1 = ETangle(ETangle.Type.OVER, (-1, -1, -1, -1), 2)
        sd_under1_1 = StrandDiagram(under1, {2: 2, 4: 4}, {0: 0, 1: 1, 3: 3})
        sd_under1_1_out = Bimodule.Element(
            {StrandDiagram(under1, {2: 4, 4: 2}, {0: 0, 1: 1, 3: 3}): under1.polyring['U3'] * under1.polyring['U4']}
        )

        # Figure 9 from "An introduction..."
        over3 = ETangle(ETangle.Type.OVER, (1, 1, -1, -1), 2)
        sd_over3_1 = StrandDiagram(over3, {1: 2, 2: 1, 3: 4}, {0: 1, 3: 2})
        sd_over3_1_out = Bimodule.Element(
            {StrandDiagram(over3, {1: 2, 2: 4, 3: 1}, {0: 1, 3: 2}): over3.polyring['U3'],
             StrandDiagram(over3, {1: 4, 2: 1, 3: 2}, {0: 1, 3: 2}): over3.polyring['U3']}
        )

        self.assertEqual(sd_under1_1_out, dminus(sd_under1_1))
        self.assertEqual(sd_over3_1_out, dminus(sd_over3_1))

    def test_dmixed(self):
        # Figure 9 from "An introduction..."
        over3 = ETangle(ETangle.Type.OVER, (1, 1, -1, -1), 2)
        sd_over3_1 = StrandDiagram(over3, {1: 2, 2: 1, 3: 4}, {0: 1, 3: 2})
        sd_over3_1_out = Bimodule.Element(
            {StrandDiagram(over3, {1: 1, 2: 2, 3: 4}, {0: 1, 3: 2}): over3.polyring['U2'],
             StrandDiagram(over3, {1: 2, 2: 3, 3: 4}, {0: 1, 1: 2}): over3.polyring['U2'] * over3.polyring['U3'],
             StrandDiagram(over3, {1: 3, 2: 1, 3: 4}, {0: 1, 2: 2}): over3.polyring['U3'],
             StrandDiagram(over3, {1: 2, 2: 1, 3: 3}, {0: 1, 4: 2}): over3.polyring.one()}
        )

        self.assertEqual(sd_over3_1_out, dmixed(sd_over3_1))


if __name__ == '__main__':
    unittest.main()
