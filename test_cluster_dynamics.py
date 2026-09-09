"""cluster_dynamics 的解析校验用例（纯 unittest，无第三方依赖）。

关键交叉验证（仿真 vs 闭式解）：
  * GCO 同环两星间距 |d|^2 = 8 rho^2 (1 - cos(d_alpha))，与时间无关
  * PCO 单链路距离摆幅恒为 sqrt(5/4) - 1 = 11.8034%
  * 高斯光束在 z = z_R 处 w = w0 * sqrt(2)
  * 远场极限退化回平方反比
"""

import contextlib
import io
import math
import unittest

import cluster_dynamics as cd
import orbit_solver as osv


class TestRelativeOrbit(unittest.TestCase):
    def test_bounded_solution_shape(self):
        sat = cd.RelativeOrbit("s", 100.0, 0.0, 100.0, 0.0)
        x, y, z = sat.position(1.0, 0.0)
        self.assertAlmostEqual(x, 0.0, places=12)
        self.assertAlmostEqual(y, 200.0, places=12)
        self.assertAlmostEqual(z, 0.0, places=12)

    def test_along_track_amplitude_is_twice_radial(self):
        sat = cd.RelativeOrbit("s", 100.0, 0.0, 0.0, 0.0)
        n = 1e-3
        xs = [abs(sat.position(n, t)[0]) for t in range(0, 7000, 7)]
        ys = [abs(sat.position(n, t)[1]) for t in range(0, 7000, 7)]
        self.assertAlmostEqual(max(ys) / max(xs), 2.0, places=3)

    def test_rejects_negative_amplitude(self):
        with self.assertRaises(ValueError):
            cd.RelativeOrbit("s", -1.0, 0.0, 0.0, 0.0)


class TestGcoRigidity(unittest.TestCase):
    """GCO 构型的核心性质：整个集群像刚体一样旋转。"""

    def test_range_from_centre_is_constant(self):
        sat = cd.RelativeOrbit("s", 500.0, 30.0, math.sqrt(3) * 500.0, 30.0)
        n = osv.mean_motion(osv.R_E + 650)
        period = osv.period_from_a(osv.R_E + 650)
        ranges = [sat.range_from_centre_m(n, period * k / 32) for k in range(32)]
        self.assertAlmostEqual(max(ranges) - min(ranges), 0.0, places=9)
        self.assertAlmostEqual(ranges[0], 1000.0, places=9)

    def test_pair_distance_matches_closed_form(self):
        """|d|^2 = 8 rho^2 (1 - cos(d_alpha))，且不含 t。"""
        rho, d_alpha = 500.0, 11.25
        expected = math.sqrt(8 * rho ** 2 * (1 - math.cos(math.radians(d_alpha))))
        a, b = (cd.RelativeOrbit("a", rho, 0.0, math.sqrt(3) * rho, 0.0),
                cd.RelativeOrbit("b", rho, d_alpha, math.sqrt(3) * rho, d_alpha))
        n = osv.mean_motion(osv.R_E + 650)
        for k in range(16):
            t = osv.period_from_a(osv.R_E + 650) * k / 16
            pa, pb = a.position(n, t), b.position(n, t)
            d = math.dist(pa, pb)
            self.assertAlmostEqual(d, expected, places=9)
        self.assertAlmostEqual(expected, 196.0, delta=0.1)

    def test_whole_cluster_is_rigid(self):
        r = cd.propagate_period(cd.build_cluster(1000.0, mode="gco"), 650.0,
                                steps=48)
        self.assertTrue(r["is_rigid"])
        self.assertAlmostEqual(r["per_link"]["worst_swing_pct"], 0.0, places=9)
        self.assertAlmostEqual(r["per_link"]["worst_coupling_swing_db"], 0.0,
                               places=9)

    def test_cluster_span_is_twice_the_radius(self):
        r = cd.propagate_period(cd.build_cluster(1000.0, mode="gco"), 650.0,
                                steps=8)
        self.assertAlmostEqual(r["span_spread"]["max"], 2000.0, places=6)


class TestPcoBreathing(unittest.TestCase):
    def test_worst_link_swing_matches_closed_form(self):
        """PCO 同环链路 |d|^2 = 4 sin^2(phi) (4 + cos^2(psi)) → 摆幅 sqrt(5)/2。"""
        r = cd.propagate_period(cd.build_cluster(1000.0, mode="pco"), 650.0,
                                steps=180)
        expected_pct = (math.sqrt(5.0 / 4.0) - 1.0) * 100.0
        self.assertAlmostEqual(r["per_link"]["worst_swing_pct"], expected_pct,
                               places=2)
        self.assertAlmostEqual(expected_pct, 11.8034, places=3)

    def test_pco_is_not_rigid(self):
        r = cd.propagate_period(cd.build_cluster(1000.0, mode="pco"), 650.0,
                                steps=48)
        self.assertFalse(r["is_rigid"])

    def test_ensemble_mean_hides_the_per_link_swing(self):
        """回归测试：集合统计量会把逐链路波动平均掉，两者必须分开报。"""
        r = cd.propagate_period(cd.build_cluster(1000.0, mode="pco"), 650.0,
                                steps=48)
        self.assertLess(r["nearest_spread"]["swing_pct"], 1.0)
        self.assertGreater(r["per_link"]["worst_swing_pct"], 10.0)


class TestClusterBuilder(unittest.TestCase):
    def test_default_rings_give_81_satellites(self):
        self.assertEqual(len(cd.build_cluster(1000.0)), 81)

    def test_ring_radii_are_evenly_spaced(self):
        sats = cd.build_cluster(1000.0, [1, 8, 16], mode="gco")
        radii = sorted({round(2 * s.rho_x_m, 6) for s in sats})
        self.assertEqual(radii, [0.0, 500.0, 1000.0])

    def test_gco_and_pco_differ_only_in_cross_track_shape(self):
        g = cd.build_cluster(1000.0, [1, 8], "gco")[-1]
        p = cd.build_cluster(1000.0, [1, 8], "pco")[-1]
        self.assertAlmostEqual(g.rho_x_m, p.rho_x_m, places=12)
        self.assertAlmostEqual(g.rho_z_m / g.rho_x_m, math.sqrt(3), places=12)
        self.assertAlmostEqual(p.rho_z_m / p.rho_x_m, 2.0, places=12)

    def test_rejects_bad_input(self):
        with self.assertRaises(ValueError):
            cd.build_cluster(0.0)
        with self.assertRaises(ValueError):
            cd.build_cluster(1000.0, mode="square")
        with self.assertRaises(ValueError):
            cd.build_cluster(1000.0, [1, 0, 8])


class TestGaussianBeam(unittest.TestCase):
    def test_rayleigh_range_closed_form(self):
        w0, lam = 0.025, 1.55e-6
        self.assertAlmostEqual(cd.rayleigh_range_m(w0, lam),
                               math.pi * w0 ** 2 / lam, places=9)
        self.assertAlmostEqual(cd.rayleigh_range_m(w0, lam), 1266.8, delta=1.0)

    def test_beam_radius_at_waist_and_rayleigh_range(self):
        w0, lam = 0.025, 1.55e-6
        z_r = cd.rayleigh_range_m(w0, lam)
        self.assertAlmostEqual(cd.beam_radius_m(0.0, w0, lam), w0, places=12)
        self.assertAlmostEqual(cd.beam_radius_m(z_r, w0, lam),
                               w0 * math.sqrt(2), places=12)

    def test_near_field_coupling_is_nearly_flat(self):
        """2.5 cm 口径下 100 m -> 200 m 的耦合损失可忽略。"""
        w = 0.025
        ratio = cd.coupling_efficiency(100.0, w) / cd.coupling_efficiency(200.0, w)
        self.assertLess(10 * math.log10(ratio), 0.1)

    def test_far_field_recovers_inverse_square(self):
        """z >> z_R 时距离翻倍 → 耦合降为 1/4（6 dB）。"""
        w = 0.001
        z_r = cd.rayleigh_range_m(w, cd.DEFAULT_WAVELENGTH_M)
        far = 200 * z_r
        ratio = cd.coupling_efficiency(far, w) / cd.coupling_efficiency(2 * far, w)
        self.assertAlmostEqual(ratio, 4.0, delta=0.05)

    def test_coupling_is_monotonic_and_bounded(self):
        w = 0.025
        prev = 1.1
        for d in (0.0, 100.0, 500.0, 2000.0, 10000.0):
            eta = cd.coupling_efficiency(d, w)
            self.assertTrue(0.0 < eta <= 1.0)
            self.assertLess(eta, prev)
            prev = eta

    def test_light_time(self):
        self.assertAlmostEqual(cd.light_time_us(cd.C_M_S / 1e6), 1.0, places=12)
        self.assertAlmostEqual(cd.light_time_us(200.0), 0.667, delta=0.001)

    def test_rejects_bad_optics(self):
        with self.assertRaises(ValueError):
            cd.rayleigh_range_m(0.0, 1.55e-6)
        with self.assertRaises(ValueError):
            cd.beam_radius_m(-1.0, 0.025, 1.55e-6)
        with self.assertRaises(ValueError):
            cd.light_time_us(-1.0)


class TestPropagation(unittest.TestCase):
    def test_sample_count_and_phase_coverage(self):
        r = cd.propagate_period(cd.build_cluster(1000.0), 650.0, steps=12)
        self.assertEqual(len(r["samples"]), 12)
        self.assertAlmostEqual(r["samples"][0]["phase_deg"], 0.0, places=12)
        self.assertAlmostEqual(r["samples"][-1]["phase_deg"], 330.0, places=12)
        self.assertAlmostEqual(r["period_min"], 97.73, delta=0.02)

    def test_latency_tracks_distance(self):
        r = cd.propagate_period(cd.build_cluster(1000.0), 650.0, steps=8)
        s = r["samples"][0]
        self.assertAlmostEqual(s["nearest_latency_us"],
                               cd.light_time_us(s["nearest_mean_m"]), places=12)
        self.assertLess(s["nearest_latency_us"], 1.0)
        self.assertLess(s["span_latency_us"], 10.0)

    def test_rejects_bad_input(self):
        with self.assertRaises(ValueError):
            cd.propagate_period([], 650.0)
        with self.assertRaises(ValueError):
            cd.propagate_period(cd.build_cluster(1000.0), 650.0, steps=1)
        with self.assertRaises(ValueError):
            cd.propagate_period(cd.build_cluster(1000.0), -5.0)


class TestReportAndCli(unittest.TestCase):
    def test_report_sections(self):
        text = cd.format_period_report(
            cd.propagate_period(cd.build_cluster(1000.0), 650.0, steps=12))
        for marker in ("时间序列", "波动幅度", "逐链路波动", "瑞利距离"):
            self.assertIn(marker, text)
        self.assertIn("刚性构型", text)

    def test_subcommands_exit_zero(self):
        cases = [
            ["period", "650", "1000"],
            ["period", "650", "1000", "--mode", "pco", "--steps", "8",
             "--aperture-cm", "1.0"],
            ["period", "650", "500", "--rings", "1", "6", "12"],
            ["link", "200"],
            ["link", "2000", "--aperture-cm", "0.5"],
            ["demo"],
            [],
        ]
        for argv in cases:
            with contextlib.redirect_stdout(io.StringIO()) as out:
                self.assertEqual(cd.main(argv), 0)
            self.assertTrue(out.getvalue().strip(), f"{argv} 无输出")


if __name__ == "__main__":
    unittest.main(verbosity=2)
