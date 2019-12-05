from Modules.CTMinus import *
import unittest


class TestCTMinus(unittest.TestCase):

    def test_d_plus(self):
        under1 = ETangle(ETangle.Type.UNDER, (1, 1, 1, 1), 2)
        sd_under1_1 = ETangleStrands(under1, {0: 0, 3: 3, 4: 4}, {1: 4, 2: 3})
        sd_under1_1_out = Bimodule.Element(
            {ETangleStrands(under1, {0: 0, 3: 3, 4: 4}, {1: 3, 2: 4}): under1.polyring.one()})
        sd_under1_2 = ETangleStrands(under1, {0: 0, 2: 2, 3: 3}, {1: 3, 4: 0})
        sd_under1_2_out = Bimodule.Element(
            {ETangleStrands(under1, {0: 0, 2: 2, 3: 3}, {1: 0, 4: 3}): under1.polyring['U2'] * under1.polyring['U3']})
        under2 = ETangle(ETangle.Type.UNDER, (1, -1, 1, 1), 2)
        sd_under2_1 = ETangleStrands(under2, {0: 0, 3: 3, 4: 4}, {1: 4, 2: 3})
        sd_under2_1_out = Bimodule.Element(
            {ETangleStrands(under2, {0: 0, 3: 3, 4: 4}, {1: 3, 2: 4}): under2.polyring.one()})
        sd_under2_2 = ETangleStrands(under2, {0: 0, 2: 2, 3: 3}, {1: 3, 4: 0})
        sd_under2_2_out = Bimodule.Element(
            {ETangleStrands(under2, {0: 0, 2: 2, 3: 3}, {1: 0, 4: 3}): under2.polyring.zero()})

        over1 = ETangle(ETangle.Type.OVER, (1, 1, 1, 1), 2)
        sd_over1_1 = ETangleStrands(over1, {4: 4, 1: 1, 3: 3}, {2: 2, 0: 4})
        sd_over1_1_out = Bimodule.Element(
            {ETangleStrands(over1, {4: 4, 1: 1, 3: 3}, {0: 2, 2: 4}): over1.polyring['U2']}
        )
        sd_over1_2 = ETangleStrands(over1, {1: 1, 4: 4, 3: 3}, {2: 1, 0: 4})
        sd_over1_2_out = Bimodule.Element(
            {ETangleStrands(over1, {1: 1, 4: 4, 3: 3}, {0: 1, 2: 4}): over1.polyring['U2']}
        )
        over2 = ETangle(ETangle.Type.OVER, (1, 1, -1, 1), 2)
        sd_over2_1 = ETangleStrands(over2, {1: 1, 0: 0, 4: 4}, {3: 1, 2: 4})
        sd_over2_1_out = Bimodule.Element(
            {ETangleStrands(over2, {1: 1, 0: 0, 4: 4}, {3: 4, 2: 1}): over2.polyring.zero()}
        )
        # Figure 9 from "An introduction..."
        over3 = ETangle(ETangle.Type.OVER, (1, 1, -1, -1), 2)
        sd_over3_1 = ETangleStrands(over3, {1: 2, 2: 1, 3: 4}, {0: 1, 3: 2})
        sd_over3_1_out = Bimodule.Element()

        cap1 = ETangle(ETangle.Type.CAP, (1, 1, -1, 1), 2)
        sd_cap1_1 = ETangleStrands(cap1, {0: 0, 1: 1}, {4: 1, 3: 2})
        sd_cap1_1_out = Bimodule.Element(
            {ETangleStrands(cap1, {0: 0, 1: 1}, {4: 2, 3: 1}): cap1.polyring['U3']}
        )
        sd_cap1_2 = ETangleStrands(cap1, {3: 3, 1: 1}, {4: 0, 0: 2})
        sd_cap1_2_out = Bimodule.Element(
            {ETangleStrands(cap1, {3: 3, 1: 1}, {4: 2, 0: 0}): cap1.polyring['U1'] * cap1.polyring['U3']}
        )

        self.assertEqual(sd_under1_1_out, d_plus(sd_under1_1))
        self.assertEqual(sd_under1_2_out, d_plus(sd_under1_2))
        self.assertEqual(sd_under2_1_out, d_plus(sd_under2_1))
        self.assertEqual(sd_under2_2_out, d_plus(sd_under2_2))

        self.assertEqual(sd_over1_1_out, d_plus(sd_over1_1))
        self.assertEqual(sd_over1_2_out, d_plus(sd_over1_2))
        self.assertEqual(sd_over2_1_out, d_plus(sd_over2_1))
        self.assertEqual(sd_over3_1_out, d_plus(sd_over3_1))

        self.assertEqual(sd_cap1_1_out, d_plus(sd_cap1_1))
        self.assertEqual(sd_cap1_2_out, d_plus(sd_cap1_2))

    def test_d_minus(self):
        under1 = ETangle(ETangle.Type.OVER, (-1, -1, -1, -1), 2)
        sd_under1_1 = ETangleStrands(under1, {2: 2, 4: 4}, {0: 0, 1: 1, 3: 3})
        sd_under1_1_out = Bimodule.Element(
            {ETangleStrands(under1, {2: 4, 4: 2}, {0: 0, 1: 1, 3: 3}): under1.polyring['U3'] * under1.polyring['U4']}
        )

        # Figure 9 from "An introduction..."
        over3 = ETangle(ETangle.Type.OVER, (1, 1, -1, -1), 2)
        sd_over3_1 = ETangleStrands(over3, {1: 2, 2: 1, 3: 4}, {0: 1, 3: 2})
        sd_over3_1_out = Bimodule.Element(
            {ETangleStrands(over3, {1: 2, 2: 4, 3: 1}, {0: 1, 3: 2}): over3.polyring['U3'],
             ETangleStrands(over3, {1: 4, 2: 1, 3: 2}, {0: 1, 3: 2}): over3.polyring['U3']}
        )

        self.assertEqual(sd_under1_1_out, d_minus(sd_under1_1))
        self.assertEqual(sd_over3_1_out, d_minus(sd_over3_1))

    def test_d_mixed(self):
        # Figure 9 from "An introduction..."
        over3 = ETangle(ETangle.Type.OVER, (1, 1, -1, -1), 2)
        sd_over3_1 = ETangleStrands(over3, {1: 2, 2: 1, 3: 4}, {0: 1, 3: 2})
        sd_over3_1_out = Bimodule.Element(
            {ETangleStrands(over3, {1: 1, 2: 2, 3: 4}, {0: 1, 3: 2}): over3.polyring['U2'],
             ETangleStrands(over3, {1: 2, 2: 3, 3: 4}, {0: 1, 1: 2}): over3.polyring['U2'] * over3.polyring['U3'],
             ETangleStrands(over3, {1: 3, 2: 1, 3: 4}, {0: 1, 2: 2}): over3.polyring['U3'],
             ETangleStrands(over3, {1: 2, 2: 1, 3: 3}, {0: 1, 4: 2}): over3.polyring.one()}
        )

        self.assertEqual(sd_over3_1_out, d_mixed(sd_over3_1))

    def test_m2(self):
        # Figure 10 from "An introduction..."
        cap2 = ETangle(ETangle.Type.CAP, (1, -1, -1), 1)
        sd_cap2_1 = ETangleStrands(cap2, {0: 0, 1: 3}, {2: 0})
        sd_cap2_1_out = Bimodule.Element(
            {ETangleStrands(cap2, {0: 0, 1: 3}, {2: 1}): cap2.polyring.one()}
        )
        algebra1 = AMinus((-1,))
        elt1 = algebra1.gen({0: 1})
        idem = algebra1.gen({0: 0})

        self.assertEqual(sd_cap2_1_out, m2(sd_cap2_1, elt1))
        self.assertEqual(Bimodule.Element({sd_cap2_1: cap2.polyring.one()}), m2(sd_cap2_1, idem))

        cup1 = ETangle(ETangle.Type.CUP, (1, -1), 1)
        sd_cup1_1 = ETangleStrands(cup1, {}, {2: 1, 0: 0})
        algebra2 = AMinus((1, -1))
        idem2 = algebra2.gen({0: 0, 1: 1})
        sd_cup1_2 = ETangleStrands(cup1, {}, {2: 0, 0: 1})

        self.assertEqual(Bimodule.Element({sd_cup1_1: cup1.polyring.one()}), m2(sd_cup1_1, idem2))
        self.assertEqual(Bimodule.Element({sd_cup1_2: cup1.polyring.one()}), m2(sd_cup1_2, idem2))

    def test_delta_ell_case_1(self):
        over1 = ETangle(ETangle.Type.OVER, (-1, -1), 1)
        x = ETangleStrands(over1, {0: 2}, {0: 0, 1: 1})
        a1 = 1
        a2 = 2
        y = x
        out = Bimodule.Edge(
            x, y, over1.polyring['U2'], (over1.left_algebra.gen({1: 2, 2: 1}),), tuple()
        )
        self.assertEqual(out, delta_ell_case_1(x, a1, a2))

    def test_delta_ell_case_2(self):
        over1 = ETangle(ETangle.Type.OVER, (1, 1), 1)
        x = ETangleStrands(over1, {0: 1, 1: 0}, {2: 2})
        a1 = 0
        a2 = 1
        y = ETangleStrands(over1, {0: 0, 1: 1}, {2: 2})
        out = Bimodule.Edge(
            x, y, over1.polyring['U1'], (over1.left_algebra.gen({2: 2}),), tuple()
        )
        self.assertEqual(out, delta_ell_case_2(x, a1, a2))

    def test_delta_ell_case_3(self):
        over1 = ETangle(ETangle.Type.OVER, (-1, 1), 1)
        x = ETangleStrands(over1, {1: 1, 2: 2}, {0: 0})
        a1 = 0
        a2 = 1
        y = ETangleStrands(over1, {0: 1, 2: 2}, {0: 0})
        out = Bimodule.Edge(
            x, y, over1.polyring['U1'], (over1.left_algebra.gen({0: 1}),), tuple()
        )
        self.assertEqual(out, delta_ell_case_3(x, a1, a2))

    def test_delta_ell_case_4(self):
        over1 = ETangle(ETangle.Type.OVER, (-1, -1), 1)
        x = ETangleStrands(over1, {0: 1, 1: 0}, {2: 2})
        a1 = 0
        a2 = 2
        y = ETangleStrands(over1, {2: 1, 1: 0}, {2: 2})
        out = Bimodule.Edge(
            x, y, over1.polyring['U2'], (over1.left_algebra.gen({2: 0}),), tuple()
        )
        self.assertEqual(out, delta_ell_case_4(x, a1, a2))

    def test_delta_ell(self):
        over1 = ETangle(ETangle.Type.OVER, (-1, -1, 1), 2)
        x = ETangleStrands(over1, {1: 0, 2: 1}, {3: 0, 2: 2})
        y = ETangleStrands(over1, {0: 0, 2: 1}, {3: 0, 2: 2})
        c = over1.polyring.one()
        elt = over1.left_algebra.gen({0: 1, 3: 3})
        out = [Bimodule.Edge(x, y, c, (elt,), tuple())]
        self.assertEqual(out, delta_ell(x))

    def test_type_da(self):
        et = ETangle(ETangle.Type.CUP, (1, -1), 1)
        da = type_da(et)
        gv = da.to_agraph(idempotents=False)
        gv.draw('output/test_type_da.png')


if __name__ == '__main__':
    unittest.main()
