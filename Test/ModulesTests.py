from Modules.CTMinus import *
from Modules.Module import *
from Modules.ETangleStrands import *
import unittest
import timeit


class TestCTMinus(unittest.TestCase):

    @staticmethod
    def empty_da_module(etangle: ETangle):
        return TypeDA(etangle.ring, etangle.left_algebra, etangle.right_algebra,
                      etangle.right_scalar_action)

    def test_d_plus(self):
        under1 = ETangle(ETangle.Type.UNDER, (1, 1, 1, 1), 2)
        under1_module = TestCTMinus.empty_da_module(under1)
        sd_under1_1 = ETangleStrands(under1, {0: 0, 3: 3, 4: 4}, {1: 4, 2: 3})
        sd_under1_1_out = ETangleStrands(under1, {0: 0, 3: 3, 4: 4}, {1: 3, 2: 4}).to_generator(under1_module)
        sd_under1_2 = ETangleStrands(under1, {0: 0, 2: 2, 3: 3}, {1: 3, 4: 0})
        sd_under1_2_out = under1.ring['U2'] * under1.ring['U3'] * \
                          ETangleStrands(under1, {0: 0, 2: 2, 3: 3}, {1: 0, 4: 3}).to_generator(under1_module)
        under2 = ETangle(ETangle.Type.UNDER, (1, -1, 1, 1), 2)
        under2_module = TestCTMinus.empty_da_module(under2)
        sd_under2_1 = ETangleStrands(under2, {0: 0, 3: 3, 4: 4}, {1: 4, 2: 3})
        sd_under2_1_out = ETangleStrands(under2, {0: 0, 3: 3, 4: 4}, {1: 3, 2: 4}).to_generator(under2_module)
        sd_under2_2 = ETangleStrands(under2, {0: 0, 2: 2, 3: 3}, {1: 3, 4: 0})
        sd_under2_2_out = under2_module.zero()

        over1 = ETangle(ETangle.Type.OVER, (1, 1, 1, 1), 2)
        over1_module = TestCTMinus.empty_da_module(over1)
        sd_over1_1 = ETangleStrands(over1, {4: 4, 1: 1, 3: 3}, {2: 2, 0: 4})
        sd_over1_1_out = over1.ring['U2'] * \
                         ETangleStrands(over1, {4: 4, 1: 1, 3: 3}, {0: 2, 2: 4}).to_generator(over1_module)
        sd_over1_2 = ETangleStrands(over1, {1: 1, 4: 4, 3: 3}, {2: 1, 0: 4})
        sd_over1_2_out = over1.ring['U2'] * \
                         ETangleStrands(over1, {1: 1, 4: 4, 3: 3}, {0: 1, 2: 4}).to_generator(over1_module)
        over2 = ETangle(ETangle.Type.OVER, (1, 1, -1, 1), 2)
        over2_module = TestCTMinus.empty_da_module(over2)
        sd_over2_1 = ETangleStrands(over2, {1: 1, 0: 0, 4: 4}, {3: 1, 2: 4})
        sd_over2_1_out = over2_module.zero()
        # Figure 9 from "An introduction..."
        over3 = ETangle(ETangle.Type.OVER, (1, 1, -1, -1), 2)
        over3_module = TestCTMinus.empty_da_module(over3)
        sd_over3_1 = ETangleStrands(over3, {1: 2, 2: 1, 3: 4}, {0: 1, 3: 2})
        sd_over3_1_out = over3_module.zero()

        cap1 = ETangle(ETangle.Type.CAP, (1, 1, -1, 1), 2)
        cap1_module = TestCTMinus.empty_da_module(cap1)
        sd_cap1_1 = ETangleStrands(cap1, {0: 0, 1: 1}, {4: 1, 3: 2})
        sd_cap1_1_out = cap1.ring['U3'] * ETangleStrands(cap1, {0: 0, 1: 1}, {4: 2, 3: 1}).to_generator(cap1_module)
        sd_cap1_2 = ETangleStrands(cap1, {3: 3, 1: 1}, {4: 0, 0: 2})
        sd_cap1_2_out = cap1.ring['U1'] * cap1.ring['U3'] * \
                        ETangleStrands(cap1, {3: 3, 1: 1}, {4: 2, 0: 0}).to_generator(cap1_module)

        self.assertEqual(sd_under1_1_out, d_plus(under1_module, sd_under1_1))
        self.assertEqual(sd_under1_2_out, d_plus(under1_module, sd_under1_2))
        self.assertEqual(sd_under2_1_out, d_plus(under2_module, sd_under2_1))
        self.assertEqual(sd_under2_2_out, d_plus(under2_module, sd_under2_2))

        self.assertEqual(sd_over1_1_out, d_plus(over1_module, sd_over1_1))
        self.assertEqual(sd_over1_2_out, d_plus(over1_module, sd_over1_2))
        self.assertEqual(sd_over2_1_out, d_plus(over2_module, sd_over2_1))
        self.assertEqual(sd_over3_1_out, d_plus(over3_module, sd_over3_1))

        self.assertEqual(sd_cap1_1_out, d_plus(cap1_module, sd_cap1_1))
        self.assertEqual(sd_cap1_2_out, d_plus(cap1_module, sd_cap1_2))

    def test_d_minus(self):
        under1 = ETangle(ETangle.Type.OVER, (-1, -1, -1, -1), 2)
        under1_module = TestCTMinus.empty_da_module(under1)
        sd_under1_1 = ETangleStrands(under1, {2: 2, 4: 4}, {0: 0, 1: 1, 3: 3})
        sd_under1_1_out = under1.ring['U3'] * under1.ring['U4'] * \
                          ETangleStrands(under1, {2: 4, 4: 2}, {0: 0, 1: 1, 3: 3}).to_generator(under1_module)

        # Figure 9 from "An introduction..."
        over3 = ETangle(ETangle.Type.OVER, (1, 1, -1, -1), 2)
        over3_module = TestCTMinus.empty_da_module(over3)
        sd_over3_1 = ETangleStrands(over3, {1: 2, 2: 1, 3: 4}, {0: 1, 3: 2})
        sd_over3_1_out = over3.ring['U3'] * \
                         ETangleStrands(over3, {1: 2, 2: 4, 3: 1}, {0: 1, 3: 2}).to_generator(over3_module) + \
                         over3.ring['U3'] * \
                         ETangleStrands(over3, {1: 4, 2: 1, 3: 2}, {0: 1, 3: 2}).to_generator(over3_module)

        self.assertEqual(sd_under1_1_out, d_minus(under1_module, sd_under1_1))
        self.assertEqual(sd_over3_1_out, d_minus(over3_module, sd_over3_1))

    def test_d_mixed(self):
        # Figure 9 from "An introduction..."
        over3 = ETangle(ETangle.Type.OVER, (1, 1, -1, -1), 2)
        over3_module = TestCTMinus.empty_da_module(over3)
        sd_over3_1 = ETangleStrands(over3, {1: 2, 2: 1, 3: 4}, {0: 1, 3: 2})

        sd_over3_1_out = over3.ring['U2'] * \
                         ETangleStrands(over3, {1: 1, 2: 2, 3: 4}, {0: 1, 3: 2}).to_generator(over3_module) + \
                         over3.ring['U2'] * over3.ring['U3'] * \
                         ETangleStrands(over3, {1: 2, 2: 3, 3: 4}, {0: 1, 1: 2}).to_generator(over3_module) + \
                         over3.ring['U3'] * \
                         ETangleStrands(over3, {1: 3, 2: 1, 3: 4}, {0: 1, 2: 2}).to_generator(over3_module) + \
                         over3.ring.one() * \
                         ETangleStrands(over3, {1: 2, 2: 1, 3: 3}, {0: 1, 4: 2}).to_generator(over3_module)

        self.assertEqual(sd_over3_1_out, d_mixed(over3_module, sd_over3_1))

    def test_m2(self):
        # Figure 10 from "An introduction..."
        cap2 = ETangle(ETangle.Type.CAP, (1, -1, -1), 1)
        cap2_module = TestCTMinus.empty_da_module(cap2)
        sd_cap2_1 = ETangleStrands(cap2, {0: 0, 1: 3}, {2: 0})
        sd_cap2_1_out = ETangleStrands(cap2, {0: 0, 1: 3}, {2: 1}).to_generator(cap2_module)
        algebra1 = AMinus((-1,))
        elt1 = algebra1.generator({0: 1})
        idem = algebra1.generator({0: 0})

        self.assertEqual(sd_cap2_1_out, m2(cap2_module, sd_cap2_1, elt1))
        self.assertEqual(sd_cap2_1.to_generator(cap2_module), m2(cap2_module, sd_cap2_1, idem))

        cup1 = ETangle(ETangle.Type.CUP, (1, -1), 1)
        cup1_module = TestCTMinus.empty_da_module(cup1)
        sd_cup1_1 = ETangleStrands(cup1, {}, {2: 1, 0: 0})
        algebra2 = AMinus((1, -1))
        idem2 = algebra2.generator({0: 0, 1: 1})
        elt = algebra2.generator({1: 0, 2: 2})
        sd_cup1_2 = ETangleStrands(cup1, {}, {2: 2, 0: 1})
        sd_cup1_2_out = cup1_module.ring['U1'] * ETangleStrands(cup1, {}, {2: 2, 0: 0}).to_generator(cup1_module)

        self.assertEqual(sd_cup1_1.to_generator(cup1_module), m2(cup1_module, sd_cup1_1, idem2))
        self.assertEqual(sd_cup1_2_out, m2(cup1_module, sd_cup1_2, elt))

        cup2 = ETangle(ETangle.Type.CUP, (-1, 1), 1)
        cup2_module = TestCTMinus.empty_da_module(cup2)
        sd_cup2_1 = ETangleStrands(cup2, {}, {0: 0, 2: 1})
        elt2 = cup2.right_algebra.generator({0: 0, 1: 2})
        sd_cup2_1_out = cup2_module.ring['U1'] * ETangleStrands(cup2, {}, {0: 0, 2: 2}).to_generator(cup2_module)

        self.assertEqual(sd_cup2_1_out, m2(cup2_module, sd_cup2_1, elt2))

    def test_delta_ell_case_1(self):
        over1 = ETangle(ETangle.Type.OVER, (-1, -1), 1)
        over1_module = TestCTMinus.empty_da_module(over1)
        x = ETangleStrands(over1, {0: 2}, {0: 0, 1: 1})
        a1 = 1
        a2 = 2
        y = x
        out = over1.left_algebra.generator({1: 2, 2: 1}).to_element() ** \
              (over1.ring['U2'] * y.to_generator(over1_module))
        self.assertEqual(out, delta_ell_case_1(over1_module, x, a1, a2))

    def test_delta_ell_case_2(self):
        over1 = ETangle(ETangle.Type.OVER, (1, 1), 1)
        over1_module = TestCTMinus.empty_da_module(over1)
        x = ETangleStrands(over1, {0: 1, 1: 0}, {2: 2})
        a1 = 0
        a2 = 1
        y = ETangleStrands(over1, {0: 0, 1: 1}, {2: 2})
        out = (over1.left_algebra.ring['U1'] * over1.left_algebra.generator({2: 2}).to_element()) ** \
              y.to_generator(over1_module)
        self.assertEqual(out, delta_ell_case_2(over1_module, x, a1, a2))

    def test_delta_ell_case_3(self):
        over1 = ETangle(ETangle.Type.OVER, (-1, 1), 1)
        over1_module = TestCTMinus.empty_da_module(over1)
        x = ETangleStrands(over1, {1: 1, 2: 2}, {0: 0})
        a1 = 0
        a2 = 1
        y = ETangleStrands(over1, {0: 1, 2: 2}, {0: 0})
        out = over1.left_algebra.generator({0: 1}).to_element() ** (over1.ring['U1'] * y.to_generator(over1_module))
        self.assertEqual(out, delta_ell_case_3(over1_module, x, a1, a2))

    def test_delta_ell_case_4(self):
        over1 = ETangle(ETangle.Type.OVER, (-1, -1), 1)
        over1_module = TestCTMinus.empty_da_module(over1)
        x = ETangleStrands(over1, {0: 1, 1: 0}, {2: 2})
        a1 = 0
        a2 = 2
        y = ETangleStrands(over1, {2: 1, 1: 0}, {2: 2})
        out = over1.left_algebra.generator({2: 0}).to_element() ** (over1.ring['U2'] * y.to_generator(over1_module))
        self.assertEqual(out, delta_ell_case_4(over1_module, x, a1, a2))

    def test_delta_ell(self):
        over1 = ETangle(ETangle.Type.OVER, (-1, -1, 1), 2)
        over1_module = TestCTMinus.empty_da_module(over1)
        x = ETangleStrands(over1, {1: 0, 2: 1}, {3: 0, 2: 2})
        y = ETangleStrands(over1, {0: 0, 2: 1}, {3: 0, 2: 2})
        c = over1.ring.one()
        elt = over1.left_algebra.generator({0: 1, 3: 3}).to_element()
        out = elt ** (c * y.to_generator(over1_module))
        self.assertEqual(out, delta_ell(over1_module, x))

    def test_type_da(self):
        idempotents = False
        cup = ETangle(ETangle.Type.CUP, (-1, 1), 1)
        cup_da = type_da(cup)
        da_list = [cup_da]
        for da in da_list:
            self.assertEqual(12, len(da.graph.nodes))
            self.assertEqual(2, len(da.decomposed()))

    # def test_tensor(self):
    #     idempotents = False
    #     cup = ETangle(ETangle.Type.CUP, (-1, 1), 1)
    #     cup_da = type_da(cup)
    #     cap = ETangle(ETangle.Type.CAP, (-1, 1), 1)
    #     cap_da = type_da(cap)
    #     cap_da.to_agraph(idempotents=idempotents).draw('output/test_cap.svg')
    #     unknot_da = cup_da ** cap_da
    #     unknot_da.to_agraph(idempotents=idempotents).draw('output/test_unknot_t.svg')

    # def test_type_da_reduced(self):
        # idempotents = False
        #
        # # some simple examples
        # cup = ETangle(ETangle.Type.CUP, (1, -1), 1)
        # over = ETangle(ETangle.Type.OVER, (1, -1), 1)
        # cap = ETangle(ETangle.Type.CAP, (-1, 1), 1)
        #
        # cup_da = type_da(cup)
        # cup_da_reduced = cup_da.reduce()
        # over_da = type_da(over)
        # over_da_reduced = over_da.reduce()
        # cap_da = type_da(cap)
        # cap_da_reduced = cap_da.reduce()
        # cap_da.to_agraph(idempotents=False).draw('output/cap_da.svg')
        # cap_da_reduced.to_agraph(idempotents=True).draw('output/cap_da_reduced.svg')
        #
        # (cup_da_reduced ** (over_da_reduced ** cap_da_reduced)) \
        #     .to_agraph(idempotents=idempotents).draw('output/test0.svg')
        # (cup_da_reduced ** (over_da_reduced ** cap_da_reduced)).to_chain_complex().write_m2_def('output/test0.m2')
        # ((cup_da_reduced ** over_da_reduced) ** cap_da_reduced) \
        #     .to_agraph(idempotents=idempotents).draw('output/test1.svg')
        # ((cup_da_reduced ** over_da_reduced) ** cap_da_reduced).to_chain_complex().write_m2_def('output/test1.m2')

    def test_trefoil(self):
        t1 = ETangle(ETangle.Type.CUP, (-1, 1), 1)
        t2 = ETangle(ETangle.Type.CUP, (-1, 1, -1, 1), 3)
        t3 = ETangle(ETangle.Type.OVER, (-1, 1, -1, 1), 2)
        t4 = ETangle(ETangle.Type.UNDER, (-1, -1, 1, 1), 1)
        t5 = ETangle(ETangle.Type.OVER, (-1, -1, 1, 1), 2)
        t6 = ETangle(ETangle.Type.CAP, (-1, 1, -1, 1), 1)
        t7 = ETangle(ETangle.Type.CAP, (-1, 1), 1)
        t = t1+t2+t3+t4+t5+t6+t7

        da = reduced_type_da(t)
        da.to_chain_complex().write_m2_def('output/trefoil.m2')

    # def test_priority(self):
    #     t = ETangle(ETangle.Type.CAP, (1, -1), 1)
    #     da = type_da(t)
    #     da_reduced = da.priority_reduced()
    #     da_reduced.to_agraph(idempotents=False).draw('output/test_pool.svg')

    # def test_timer(self):
    #     print(timeit.timeit('TestCTMinus.trial1()', setup='from Test.ModulesTests import TestCTMinus', number=1))
    #     print(timeit.timeit('TestCTMinus.trial2()', setup='from Test.ModulesTests import TestCTMinus', number=1))
    #
    # @staticmethod
    # def trial1():
    #     cup = ETangle(ETangle.Type.CUP, (1, -1), 1)
    #     over = ETangle(ETangle.Type.OVER, (1, -1), 1)
    #     cap = ETangle(ETangle.Type.CAP, (-1, 1), 1)
    #     t = cup + over + cap
    #     da = reduced_type_da(t)
    #
    # @staticmethod
    # def trial2():
    #     cup = ETangle(ETangle.Type.CUP, (1, -1), 1)
    #     over = ETangle(ETangle.Type.OVER, (1, -1), 1)
    #     cap = ETangle(ETangle.Type.CAP, (-1, 1), 1)
    #     t = cup + over + cap
    #     da = reduced_type_da2(t)

    # def test_test(self):
    #     t1 = ETangle(ETangle.Type.CUP, (1, -1), 1)
    #     t2 = ETangle(ETangle.Type.OVER, (1, -1), 1)
    #     t3 = ETangle(ETangle.Type.UNDER, (1, -1), 1)
    #     t4 = ETangle(ETangle.Type.CAP, (1, -1), 1)
    #     t = t1 + t2 + t3 + t4
    #     da = reduced_type_da(t)
    #     da.to_agraph().draw('output/test.svg')

    # def test_macaulay2(self):
    #     cup = ETangle(ETangle.Type.CUP, (1, -1), 1)
    #     cup_da = type_da(cup)
    #     over = ETangle(ETangle.Type.OVER, (1, -1), 1)
    #     over_da = type_da(over)
    #     cap = ETangle(ETangle.Type.CAP, (-1, 1), 1)
    #     cap_da = type_da(cap)
    #     unknot_da_tt = cup_da ** over_da ** cap_da
    #     unknot_da_tt.to_agraph(idempotents=False).draw('output/test_unknot_tt.svg')
    #     unknot_da_ttr = unknot_da_tt.reduce()
    #     unknot_da_ttr.to_agraph(idempotents=False).draw('output/test_unknot_ttr.svg')
    #     unknot_da_ttr.write_m2_def('output/unknot_da_ttr.m2')


if __name__ == '__main__':
    unittest.main()
