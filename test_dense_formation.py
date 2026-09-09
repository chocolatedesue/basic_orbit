"""dense_formation 的解析校验用例（纯 unittest，无第三方依赖）。"""

import contextlib
import io
import math
import unittest

import dense_formation as df
import orbit_solver as osv


class TestPacking(unittest.TestCase):
    def test_hex_packing_area_conservation(self):
        n, r = 81, 1000.0
        d = df.packing_spacing_m(n, r, "hex")
        self.assertAlmostEqual((math.sqrt(3) / 2) * d ** 2 * n,
                               math.pi * r ** 2, places=6)

    def test_suncatcher_spacing_matches_reported_range(self):
        # 公开报道的次近邻间距 100~200 m；六方密排给出同量级
        d = df.packing_spacing_m(81, 1000.0, "hex")
        self.assertAlmostEqual(d, 212.0, delta=1.0)

    def test_grid_packing(self):
        self.assertAlmostEqual(df.packing_spacing_m(81, 1000.0, "grid"), 250.0,
                               places=9)

    def test_grid_requires_perfect_square(self):
        with self.assertRaises(ValueError):
            df.packing_spacing_m(80, 1000.0, "grid")

    def test_spacing_shrinks_with_more_satellites(self):
        self.assertLess(df.packing_spacing_m(324, 1000.0),
                        df.packing_spacing_m(81, 1000.0))

    def test_rejects_bad_input(self):
        with self.assertRaises(ValueError):
            df.packing_spacing_m(1, 1000.0)
        with self.assertRaises(ValueError):
            df.packing_spacing_m(81, 0.0)
        with self.assertRaises(ValueError):
            df.packing_spacing_m(81, 1000.0, "spiral")


class TestLinkPhysics(unittest.TestCase):
    def test_light_travels_300m_per_microsecond(self):
        self.assertAlmostEqual(df.light_time_us(df.C_KM_S * 1000.0 / 1e6),
                               1.0, places=9)

    def test_light_time_is_linear(self):
        self.assertAlmostEqual(df.light_time_us(2000.0),
                               2 * df.light_time_us(1000.0), places=9)

    def test_inverse_square_power_gain(self):
        self.assertAlmostEqual(df.link_power_gain(1000.0, 100.0), 100.0, places=9)
        self.assertAlmostEqual(df.link_power_gain(500.0, 500.0), 1.0, places=9)

    def test_rejects_bad_distances(self):
        with self.assertRaises(ValueError):
            df.light_time_us(-1.0)
        with self.assertRaises(ValueError):
            df.link_power_gain(0.0, 100.0)


class TestAlongTrackDrift(unittest.TestCase):
    def test_drift_per_orbit_is_three_pi_delta_a(self):
        d = df.along_track_drift(650, 10.0)
        self.assertAlmostEqual(d["drift_per_orbit_m"], 3 * math.pi * 10.0, places=9)

    def test_drift_scales_linearly(self):
        a = df.along_track_drift(650, 1.0)["drift_per_day_m"]
        b = df.along_track_drift(650, 100.0)["drift_per_day_m"]
        self.assertAlmostEqual(b, 100 * a, places=6)

    def test_one_metre_mismatch_already_drifts_kilometres_per_year(self):
        d = df.along_track_drift(650, 1.0)
        self.assertAlmostEqual(d["drift_per_orbit_m"], 9.42, delta=0.01)
        self.assertGreater(d["drift_per_year_km"], 50.0)

    def test_tolerance_round_trips_through_drift(self):
        for days in (1.0, 30.0, 365.25):
            tol = df.semi_major_axis_tolerance_m(650, 100.0, days)
            drift = df.along_track_drift(650, tol)["drift_per_day_m"] * days
            self.assertAlmostEqual(drift, 100.0, places=6)

    def test_sub_metre_tolerance_at_650km(self):
        self.assertAlmostEqual(df.semi_major_axis_tolerance_m(650, 100.0, 1.0),
                               0.72, delta=0.01)

    def test_rejects_bad_budget(self):
        with self.assertRaises(ValueError):
            df.semi_major_axis_tolerance_m(650, 0.0)
        with self.assertRaises(ValueError):
            df.semi_major_axis_tolerance_m(650, 100.0, days=0.0)


class TestDifferentialJ2(unittest.TestCase):
    def test_raan_sensitivity_is_exactly_zero(self):
        """密集集群可行性的核心：Omega_dot 不依赖 Omega。"""
        s = df.precession_sensitivity(650, 97.99)
        self.assertEqual(s["d_rate_per_deg_raan"], 0.0)

    def test_sensitivity_matches_finite_difference(self):
        h, inc, eps = 650.0, 97.99, 1e-6
        s = df.precession_sensitivity(h, inc)
        a = osv.R_E + h

        def rate(a_km, i_deg):
            return osv.rad_s_to_deg_day(osv.nodal_precession_rate(a_km, i_deg))

        fd_a = (rate(a + eps, inc) - rate(a - eps, inc)) / (2 * eps)
        fd_i = (rate(a, inc + eps) - rate(a, inc - eps)) / (2 * eps)
        self.assertAlmostEqual(s["d_rate_per_km_altitude"], fd_a, places=6)
        self.assertAlmostEqual(s["d_rate_per_deg_inclination"], fd_i, places=6)

    def test_pure_raan_spread_has_no_secular_drift(self):
        d = df.differential_nodal_drift(650, 97.99, delta_inc_deg=0.0)
        self.assertEqual(d["differential_rate_deg_day"], 0.0)
        self.assertEqual(d["cross_track_km"], 0.0)

    def test_inclination_spread_smears_the_cluster(self):
        d_inc = df.inclination_offset_for_cross_track(650, 1000.0)
        d = df.differential_nodal_drift(650, 97.99, delta_inc_deg=d_inc)
        self.assertNotEqual(d["differential_rate_deg_day"], 0.0)
        # 1 km 横向间距若用 Δi 实现，一年散开数十公里
        self.assertAlmostEqual(d["cross_track_km"], 44.8, delta=1.0)

    def test_two_cross_track_routes_need_similar_angles(self):
        h, inc, rho = 650.0, 97.99, 1000.0
        d_inc = df.inclination_offset_for_cross_track(h, rho)
        d_raan = df.raan_offset_for_cross_track(h, inc, rho)
        self.assertAlmostEqual(d_raan * math.sin(math.radians(inc)), d_inc,
                               places=9)
        self.assertAlmostEqual(d_inc, 0.00815, delta=1e-4)

    def test_cross_track_offsets_scale_linearly(self):
        self.assertAlmostEqual(
            df.inclination_offset_for_cross_track(650, 2000.0),
            2 * df.inclination_offset_for_cross_track(650, 1000.0), places=12)

    def test_equatorial_orbit_cannot_use_raan(self):
        with self.assertRaises(ValueError):
            df.raan_offset_for_cross_track(650, 0.0, 1000.0)

    def test_altitude_mismatch_also_shifts_precession(self):
        d = df.differential_nodal_drift(650, 97.99, delta_a_m=1000.0)
        self.assertAlmostEqual(abs(d["differential_rate_deg_day"]), 4.9e-4,
                               delta=1e-5)


class TestDenseClusterAnalysis(unittest.TestCase):
    def setUp(self):
        self.inc = osv.solve_sun_synchronous_inclination(650)["inclination_deg"]
        self.r = df.analyze_dense_cluster(650, self.inc, 81, 1000.0)

    def test_reproduces_the_reference_configuration(self):
        self.assertAlmostEqual(self.r["period_min"], 97.73, delta=0.02)
        self.assertAlmostEqual(self.r["revs_per_day"], 14.735, delta=0.01)
        self.assertAlmostEqual(self.r["neighbor_spacing_m"], 212.0, delta=1.0)
        self.assertTrue(self.r["is_sun_synchronous"])

    def test_sun_synchronous_precession_tracks_the_sun(self):
        self.assertAlmostEqual(self.r["precession_deg_day"], 360 / 365.2422,
                               delta=1e-3)

    def test_relative_ellipse_closes_at_the_orbital_period(self):
        self.assertAlmostEqual(self.r["relative_ellipse_period_min"],
                               self.r["period_min"], places=12)
        ell = self.r["relative_ellipse"]
        self.assertAlmostEqual(ell["along_track_amplitude_m"], 1000.0, places=6)
        self.assertAlmostEqual(ell["radial_amplitude_m"], 500.0, places=6)

    def test_raan_route_beats_inclination_route(self):
        via_i = self.r["cross_track_via_inclination"]
        via_o = self.r["cross_track_via_raan"]
        self.assertGreater(via_i["cross_track_km"], 10.0)
        self.assertEqual(via_o["cross_track_km"], 0.0)

    def test_neighbour_latency_is_sub_microsecond(self):
        self.assertLess(self.r["neighbor_light_time_us"], 1.0)
        self.assertGreater(self.r["link_power_gain"], 1e7)

    def test_flags_a_non_sun_synchronous_choice(self):
        r = df.analyze_dense_cluster(650, 53.0, 81, 1000.0)
        self.assertFalse(r["is_sun_synchronous"])


class TestReport(unittest.TestCase):
    def test_report_covers_every_section(self):
        inc = osv.solve_sun_synchronous_inclination(650)["inclination_deg"]
        text = df.format_dense_report(
            df.analyze_dense_cluster(650, inc, 81, 1000.0))
        for marker in ("轨道基本盘", "集群几何", "Δa 沿轨漂移", "横向铺开", "链路"):
            self.assertIn(marker, text)


class TestCommandLine(unittest.TestCase):
    def test_subcommands_exit_zero(self):
        cases = [
            ["analyze", "650", "81", "1000"],
            ["analyze", "550", "49", "500", "--inclination", "53",
             "--drift-budget", "50"],
            ["tolerance", "650", "100"],
            ["tolerance", "650", "100", "--days", "365.25"],
            ["demo"],
            [],
        ]
        for argv in cases:
            with contextlib.redirect_stdout(io.StringIO()) as out:
                self.assertEqual(df.main(argv), 0)
            self.assertTrue(out.getvalue().strip(), f"{argv} 无输出")

    def test_analyze_defaults_to_sun_synchronous_inclination(self):
        with contextlib.redirect_stdout(io.StringIO()) as out:
            df.main(["analyze", "650", "81", "1000"])
        self.assertIn("97.99", out.getvalue())

    def test_rejects_altitude_without_sun_sync_solution(self):
        with self.assertRaises(ValueError):
            df.main(["analyze", "20000", "81", "1000"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
