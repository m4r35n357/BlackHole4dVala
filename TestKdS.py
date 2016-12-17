from unittest import TestCase, skip, main

import KdS

class KdSTest(TestCase):

    kds_rk4 = KdS.BhRk4(0.0, 0.9, 1.0, 1.0, 4.0, 0.0, 4.0, 0.0, 0.0, 100.0, 0.001, 500, "rk4")
    kds_rk438 = KdS.BhRk4(0.0, 0.9, 1.0, 1.0, 4.0, 0.0, 4.0, 0.0, 0.0, 100.0, 0.001, 500, "rk438")
    kds_sym_b2 = KdS.BhSymp(0.0, 0.9, 1.0, 1.0, 4.0, 0.0, 4.0, 0.0, 0.0, 100.0, 0.001, 500, "sb2")
    kds_sym_b4 = KdS.BhSymp(0.0, 0.9, 1.0, 1.0, 4.0, 0.0, 4.0, 0.0, 0.0, 100.0, 0.001, 500, "sb4")
    kds_sym_c4 = KdS.BhSymp(0.0, 0.9, 1.0, 1.0, 4.0, 0.0, 4.0, 0.0, 0.0, 100.0, 0.001, 500, "sc4")
    kds_sym_c6 = KdS.BhSymp(0.0, 0.9, 1.0, 1.0, 4.0, 0.0, 4.0, 0.0, 0.0, 100.0, 0.001, 500, "sc6")

    def setUp(self):
        pass

    def test_good_integrator_types(self):
        self.assertIsInstance(self.kds_rk4, KdS.BhRk4)
        self.assertIsInstance(self.kds_rk438, KdS.BhRk4)
        self.assertIsInstance(self.kds_sym_b2, KdS.BhSymp)
        self.assertIsInstance(self.kds_sym_b4, KdS.BhSymp)
        self.assertIsInstance(self.kds_sym_c4, KdS.BhSymp)
        self.assertIsInstance(self.kds_sym_c6, KdS.BhSymp)

    def test_bad_integrator_types(self):
        self.assertRaises(Exception, KdS.BhRk4, 0.0, 0.9, 1.0, 1.0, 4.0, 0.0, 4.0, 0.0, 0.0, 100.0, 0.001, 500, "xxx")
        self.assertRaises(Exception, KdS.BhSymp, 0.0, 0.9, 1.0, 1.0, 4.0, 0.0, 4.0, 0.0, 0.0, 100.0, 0.001, 500, "xxx")

    def test_solve_rk4_0(self):
        start = 0.0
        end = 100.0
        step = 0.001
        interval = 10
        counts = KdS.BhRk4(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, start, end, step, interval, "rk4").solve()
        self.assertEqual(100000, counts[0])
        self.assertEqual(10000, counts[1])

    def test_solve_rk4_non_0(self):
        start = 50.0
        end = 100.0
        step = 0.001
        interval = 10
        counts = KdS.BhRk4(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, start, end, step, interval, "rk4").solve()
        self.assertEqual(100000, counts[0])
        self.assertEqual(5000, counts[1])

    def test_solve_symp_0(self):
        start = 0.0
        end = 100
        step = 0.00005
        interval = 10
        counts = KdS.BhSymp(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, start, end, step, interval, "sb4").solve()
        self.assertEqual(85532, counts[0])
        self.assertEqual(8554, counts[1])

    #@skip("temporarily skipped")
    def test_solve_symp_non_0(self):
        start = 50.0
        end = 100
        step = 0.00005
        interval = 10
        counts = KdS.BhSymp(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, start, end, step, interval, "sb4").solve()
        self.assertEqual(85532, counts[0])
        self.assertEqual(3215, counts[1])

if __name__ == '__main__':
    main()
