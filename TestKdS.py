from decimal import Decimal
from unittest import TestCase, main

import KdS


class KdSTest(TestCase):

    def test_good_integrator_types(self):
        self.assertIsInstance(KdS.Symplectic(Decimal('1.0'), 'sb1'), KdS.Symplectic)
        self.assertIsInstance(KdS.Symplectic(Decimal('1.0'), 'sb2'), KdS.Symplectic)
        self.assertIsInstance(KdS.Symplectic(Decimal('1.0'), 'sb4'), KdS.Symplectic)

    def test_bad_integrator_types(self):
        self.assertRaises(Exception, KdS.Symplectic, Decimal('1.0'), 'xxx')

    def test_solve_symp_polar_sb2(self):
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
                            Decimal('0.0'),
                            True).solve(KdS.Symplectic(step, 'sb2'), start, end, interval)
        self.assertEqual(695, counts[0])
        self.assertEqual(695, counts[1])

    def test_solve_symp_light_sb1(self):
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
                            Decimal('0.0'),
                            True).solve(KdS.Symplectic(step, 'sb1'), start, end, interval)
        self.assertEqual(1058, counts[0])
        self.assertEqual(1058, counts[1])

    def test_solve_symp_0_sb2(self):
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
                            Decimal('0.0'),
                            True).solve(KdS.Symplectic(step, 'sb2'), start, end, interval)
        self.assertEqual(220, counts[0])
        self.assertEqual(220, counts[1])

    def test_solve_symp_non_0_sb1(self):
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
                            Decimal('0.0'),
                            True).solve(KdS.Symplectic(step, 'sb1'), start, end, interval)
        self.assertEqual(220, counts[0])
        self.assertEqual(121, counts[1])

    def test_solve_symp_non_0_sb2(self):
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
                            Decimal('0.0'),
                            True).solve(KdS.Symplectic(step, 'sb2'), start, end, interval)
        self.assertEqual(220, counts[0])
        self.assertEqual(121, counts[1])


    def test_solve_symp_non_0_sb4(self):
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
                            Decimal('0.0'),
                            True).solve(KdS.Symplectic(step, 'sb4'), start, end, interval)
        self.assertEqual(220, counts[0])
        self.assertEqual(121, counts[1])


if __name__ == '__main__':
    main()
