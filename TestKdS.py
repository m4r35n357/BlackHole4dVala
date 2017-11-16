from decimal import Decimal
from unittest import TestCase, main

from Bh3d import BhSymp
from Symplectic import Symplectic


class KdSTest(TestCase):
    def test_debug_output_order_2(self):
        Symplectic(None, 1.0, 'b2', stages = 5, debug=True).method()

    def test_debug_output_order_4(self):
        Symplectic(None, 1.0, 'b4', stages = 5, debug=True).method()

    def test_debug_output_order_6(self):
        Symplectic(None, 1.0, 'b6', stages = 5, debug=True).method()

    def test_debug_output_order_8(self):
        Symplectic(None, 1.0, 'b8', stages = 5, debug=True).method()

    def test_good_integrator_types(self):
        self.assertIsInstance(Symplectic(None, 1.0, 'b1', stages=5), Symplectic)
        self.assertIsInstance(Symplectic(None, 1.0, 'b2', stages=5), Symplectic)
        self.assertIsInstance(Symplectic(None, 1.0, 'b4', stages=5), Symplectic)

    def test_bad_integrator_types(self):
        self.assertRaises(Exception, Symplectic, 1.0, 'xxx')

    def test_solve_symp_polar_sb2(self):
        start = 0.0
        end = 10.0
        step = 0.0001
        interval = 1
        model = BhSymp(0.0,
                        1.0,
                        1.0,
                        0.96210432940242041,
                        5.6843449527674236e-13,
                        15.914691393798241,
                        12.0,
                        0.0,
                        True)
        counts = model.solve(Symplectic(model, step, 'b2', stages=5).method, step, start, end, interval)
        self.assertEqual(695, counts[0])
        self.assertEqual(695, counts[1])

    def test_solve_symp_light_sb1(self):
        start = 0.0
        end = 10.0
        step = 0.001
        interval = 1
        model = BhSymp(0.0,
                        1.0,
                        0.0,
                        1.0,
                        -2.0,
                        27.0,
                        3.0,
                        0.0,
                        True)
        counts = model.solve(Symplectic(model, step, 'b1', stages=5).method, step, start, end, interval)
        self.assertEqual(1058, counts[0])
        self.assertEqual(1058, counts[1])

    def test_solve_symp_0_sb2(self):
        start = 0.0
        end = 10.0
        step = 0.001
        interval = 1
        model = BhSymp(0.0,
                        0.8,
                        1.0,
                        0.94550509567490792,
                        1.4343745095317371,
                        7.9787599589278697,
                        7.5,
                        0.0,
                        True)
        counts = model.solve(Symplectic(model, step, 'b2', stages=5).method, step, start, end, interval)
        self.assertEqual(220, counts[0])
        self.assertEqual(220, counts[1])

    def test_solve_symp_non_0_sb1(self):
        start = 5.0
        end = 10.0
        step = 0.001
        interval = 1
        model = BhSymp(0.0,
                        0.8,
                        1.0,
                        0.94550509567490792,
                        1.4343745095317371,
                        7.9787599589278697,
                        7.5,
                        0.0,
                        True)
        counts = model.solve(Symplectic(model, step, 'b1', stages=5).method, step, start, end, interval)
        self.assertEqual(220, counts[0])
        self.assertEqual(121, counts[1])

    def test_solve_symp_non_0_sb2(self):
        start = 5.0
        end = 10.0
        step = 0.001
        interval = 1
        model = BhSymp(0.0,
                        0.8,
                        1.0,
                        0.94550509567490792,
                        1.4343745095317371,
                        7.9787599589278697,
                        7.5,
                        0.0,
                        True)
        counts = model.solve(Symplectic(model, step, 'b2', stages=5).method, step, start, end, interval)
        self.assertEqual(220, counts[0])
        self.assertEqual(121, counts[1])

    def test_solve_symp_non_0_sb4(self):
        start = 5.0
        end = 10.0
        step = 0.001
        interval = 1
        model = BhSymp(0.0,
                        0.8,
                        1.0,
                        0.94550509567490792,
                        1.4343745095317371,
                        7.9787599589278697,
                        7.5,
                        0.0,
                        True)
        counts = model.solve(Symplectic(model, step, 'b4', stages=5).method, step, start, end, interval)
        self.assertEqual(220, counts[0])
        self.assertEqual(121, counts[1])


if __name__ == '__main__':
    main()
