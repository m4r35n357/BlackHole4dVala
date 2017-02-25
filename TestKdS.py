from decimal import Decimal
from unittest import TestCase, skip, main

import KdS

class KdSTest(TestCase):

    kds_sym_b1 = KdS.Symplectic(Decimal('1.0'), 'sb1')
    kds_sym_b2 = KdS.Symplectic(Decimal('1.0'), 'sb2')
    kds_sym_b4 = KdS.Symplectic(Decimal('1.0'), 'sb4')

    def setUp(self):
        pass

    # @skip("temporarily skipped")
    def test_good_integrator_types(self):
        self.assertIsInstance(self.kds_sym_b1, KdS.Symplectic)
        self.assertIsInstance(self.kds_sym_b2, KdS.Symplectic)
        self.assertIsInstance(self.kds_sym_b4, KdS.Symplectic)

    @skip("temporarily skipped")
    def test_bad_integrator_types(self):
        self.assertRaises(Exception, KdS.Symplectic(Decimal('1.0'), 'xxx'))

    # @skip("temporarily skipped")
    def test_solve_symp_polar(self):
        start = Decimal('0.0')
        end = Decimal('10.0')
        step = Decimal('0.0001')
        interval = 1
        counts = KdS.BhSymp(Decimal('0.0'),
                            Decimal('1.0'),
                            Decimal('1.0'),
                            Decimal('0.96210432940242041'),
                            Decimal('5.6843449527674236e-13'),
                            Decimal('15.914691393798241'),
                            Decimal('12.0'),
                            Decimal('0.0')).solve(KdS.Symplectic(step, 'sb2'), step, start, end, interval)
        self.assertEqual(695, counts[0])
        self.assertEqual(695, counts[1])

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
