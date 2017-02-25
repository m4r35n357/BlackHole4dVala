from decimal import Decimal
from unittest import TestCase, skip, main

import KdS

class KdSTest(TestCase):

    # kds_sym_b1 = KdS.BhSymp(0.0, 0.9, 1.0, 1.0, 4.0, 0.0, 4.0, 0.0)
    # kds_sym_b2 = KdS.BhSymp(0.0, 0.9, 1.0, 1.0, 4.0, 0.0, 4.0, 0.0)
    # kds_sym_b4 = KdS.BhSymp(0.0, 0.9, 1.0, 1.0, 4.0, 0.0, 4.0, 0.0)
    #
    # def setUp(self):
    #     pass
    #
    # @skip("temporarily skipped")
    # def test_good_integrator_types(self):
    #     self.assertIsInstance(self.kds_sym_b2, KdS.BhSymp)
    #     self.assertIsInstance(self.kds_sym_b4, KdS.BhSymp)
    #
    # @skip("temporarily skipped")
    # def test_bad_integrator_types(self):
    #     self.assertRaises(Exception, KdS.BhRk4, 0.0, 0.9, 1.0, 1.0, 4.0, 0.0, 4.0, 0.0, 0.0, 100.0, 0.001, 500, "xxx")
    #     self.assertRaises(Exception, KdS.BhSymp, 0.0, 0.9, 1.0, 1.0, 4.0, 0.0, 4.0, 0.0, 0.0, 100.0, 0.001, 500, "xxx")
    #
    # @skip("temporarily skipped")
    # def test_solve_symp_polar(self):
    #     start = 0.0
    #     end = 100.0
    #     step = 0.00001
    #     interval = 10
    #     counts = KdS.BhSymp(0.0, 1.0, 1.0, 0.96210432940242041, 5.6843449527674236e-13, 15.914691393798241, 12.0, 0.0).solve(KdS.Symplectic(step, 'sb2'), step, start, end, interval)
    #     self.assertEqual(69174, counts[0])
    #     self.assertEqual(6918, counts[1])

    # @skip("temporarily skipped")
    def test_solve_symp_light(self):
        start = Decimal('0.0')
        end = Decimal('10.0')
        step = Decimal('0.001')
        interval = 1
        counts = KdS.BhSymp(Decimal('0.0'),
                            Decimal('1.0'),
                            Decimal('0.0'),
                            Decimal('1.0'),
                            Decimal('-2.0'),
                            Decimal('27.0'),
                            Decimal('3.0'),
                            Decimal('0.0')).solve(KdS.Symplectic(step, 'sb1'), step, start, end, interval)
        self.assertEqual(1058, counts[0])
        self.assertEqual(1058, counts[1])

    # @skip("temporarily skipped")
    def test_solve_symp_0(self):
        start = Decimal('0.0')
        end = Decimal('10.0')
        step = Decimal('0.001')
        interval = 1
        counts = KdS.BhSymp(Decimal('0.0'),
                            Decimal('0.8'),
                            Decimal('1.0'),
                            Decimal('0.94550509567490792'),
                            Decimal('1.4343745095317371'),
                            Decimal('7.9787599589278697'),
                            Decimal('7.5'),
                            Decimal('0.0')).solve(KdS.Symplectic(step, 'sb2'), step, start, end, interval)
        self.assertEqual(220, counts[0])
        self.assertEqual(220, counts[1])

    # @skip("temporarily skipped")
    def test_solve_symp_non_0(self):
        start = Decimal('5.0')
        end = Decimal('10.0')
        step = Decimal('0.001')
        interval = 1
        counts = KdS.BhSymp(Decimal('0.0'),
                            Decimal('0.8'),
                            Decimal('1.0'),
                            Decimal('0.94550509567490792'),
                            Decimal('1.4343745095317371'),
                            Decimal('7.9787599589278697'),
                            Decimal('7.5'),
                            Decimal('0.0')).solve(KdS.Symplectic(step, 'sb4'), step, start, end, interval)
        self.assertEqual(220, counts[0])
        self.assertEqual(121, counts[1])

if __name__ == '__main__':
    main()
