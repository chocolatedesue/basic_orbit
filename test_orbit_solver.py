"""orbit_solver 的解析校验用例（纯 unittest，无第三方依赖）。

用公开的工程基准值做交叉验证：
  * 星链 550 km / 53° 壳：周期 ≈ 95.6 分钟，进动 ≈ -4.5 度/天
  * 800 km 太阳同步轨道：倾角 ≈ 98.6°
  * 地球同步轨道：恒星日周期反解高度 ≈ 35786 km
  * 临界倾角 63.435°：近地点幅角进动为零
"""

import math
import unittest

import orbit_solver as osv


class TestBasicElements(unittest.TestCase):
    def test_period_matches_kepler_third_law(self):
        a = osv.R_E + 550.0
        t = osv.period_from_a(a)
        self.assertAlmostEqual(t, 2 * math.pi * math.sqrt(a ** 3 / osv.MU), places=9)
        self.assertAlmostEqual(t / 60.0, 95.65, places=1)

    def test_period_and_semi_major_axis_are_inverse(self):
        for h in (300.0, 550.0, 1200.0, 20000.0, 35786.0):
            a = osv.R_E + h
            self.assertAlmostEqual(osv.a_from_period(osv.period_from_a(a)), a, places=6)

    def test_geostationary_altitude_from_sidereal_day(self):
        sidereal_day = 86164.0905
        h = osv.a_from_period(sidereal_day) - osv.R_E
        self.assertAlmostEqual(h, 35786.0, delta=1.0)

    def test_starlink_shell_precession(self):
        feat = osv.get_orbit_features(550, 53)
        self.assertAlmostEqual(feat["period_min"], 95.65, delta=0.05)
        self.assertAlmostEqual(feat["precession_deg_day"], -4.49, delta=0.02)
        self.assertAlmostEqual(feat["revs_per_day"], 15.06, delta=0.01)

    def test_prograde_regresses_retrograde_advances(self):
        self.assertLess(osv.get_orbit_features(550, 53)["precession_deg_day"], 0)
        self.assertGreater(osv.get_orbit_features(550, 120)["precession_deg_day"], 0)
        self.assertAlmostEqual(
            osv.get_orbit_features(550, 90)["precession_deg_day"], 0.0, places=12)

    def test_critical_inclination_freezes_perigee(self):
        crit = math.degrees(math.acos(math.sqrt(1.0 / 5.0)))
        self.assertAlmostEqual(crit, 63.4349, places=3)
        self.assertAlmostEqual(
            osv.get_orbit_features(800, crit)["perigee_rate_deg_day"], 0.0, places=12)

    def test_eccentricity_amplifies_precession(self):
        circ = abs(osv.nodal_precession_rate(osv.R_E + 800, 45, 0.0))
        ecc = abs(osv.nodal_precession_rate(osv.R_E + 800, 45, 0.2))
        self.assertGreater(ecc, circ)


class TestMatchingShell(unittest.TestCase):
    def test_matched_shells_share_precession_rate(self):
        m = osv.solve_matching_shell(550, 53, 1000)
        self.assertTrue(m["feasible"])
        rate1 = osv.nodal_precession_rate(osv.R_E + 550, 53)
        rate2 = osv.nodal_precession_rate(osv.R_E + 1000, m["matched_inc2_deg"])
        self.assertAlmostEqual(rate1, rate2, delta=abs(rate1) * 1e-12)

    def test_known_pairing_value(self):
        # 550 km/53° 与 1000 km 配对，闭式解为 41.40°（不是 44.13°）
        m = osv.solve_matching_shell(550, 53, 1000)
        self.assertAlmostEqual(m["matched_inc2_deg"], 41.40, delta=0.01)

    def test_same_altitude_returns_same_inclination(self):
        m = osv.solve_matching_shell(550, 53, 550)
        self.assertAlmostEqual(m["matched_inc2_deg"], 53.0, places=9)

    def test_ceiling_drives_inclination_to_zero(self):
        m = osv.solve_matching_shell(550, 53, 1000)
        at_max = osv.solve_matching_shell(550, 53, m["max_h2_km"])
        self.assertTrue(at_max["feasible"])
        self.assertAlmostEqual(at_max["matched_inc2_deg"], 0.0, places=4)

    def test_above_ceiling_is_infeasible(self):
        m = osv.solve_matching_shell(550, 53, 1000)
        too_high = osv.solve_matching_shell(550, 53, m["max_h2_km"] + 50.0)
        self.assertFalse(too_high["feasible"])
        self.assertIn("reason", too_high)

    def test_retrograde_base_has_finite_ceiling(self):
        m = osv.solve_matching_shell(780, 100, 1500)
        self.assertTrue(m["feasible"])
        self.assertGreater(m["matched_inc2_deg"], 90.0)
        self.assertTrue(math.isfinite(m["max_h2_km"]))
        rate1 = osv.nodal_precession_rate(osv.R_E + 780, 100)
        rate2 = osv.nodal_precession_rate(osv.R_E + 1500, m["matched_inc2_deg"])
        self.assertAlmostEqual(rate1, rate2, delta=abs(rate1) * 1e-12)

    def test_polar_base_has_no_ceiling(self):
        m = osv.solve_matching_shell(780, 90, 20000)
        self.assertTrue(m["feasible"])
        self.assertEqual(m["max_h2_km"], float("inf"))
        self.assertAlmostEqual(m["matched_inc2_deg"], 90.0, places=9)

    def test_near_polar_ceiling_is_high(self):
        m = osv.solve_matching_shell(780, 86.4, 2000)
        self.assertAlmostEqual(m["max_h2_km"], 9407.4, delta=1.0)


class TestSunSynchronous(unittest.TestCase):
    def test_800km_sun_synchronous_inclination(self):
        s = osv.solve_sun_synchronous_inclination(800)
        self.assertTrue(s["feasible"])
        self.assertAlmostEqual(s["inclination_deg"], 98.6, delta=0.1)

    def test_sun_sync_precession_equals_one_year(self):
        s = osv.solve_sun_synchronous_inclination(600)
        rate = osv.nodal_precession_rate(osv.R_E + 600, s["inclination_deg"])
        self.assertAlmostEqual(osv.rad_s_to_deg_day(rate), 0.9856, delta=1e-3)

    def test_too_high_has_no_solution(self):
        self.assertFalse(osv.solve_sun_synchronous_inclination(20000)["feasible"])


class TestResonantShells(unittest.TestCase):
    def test_k_revolutions_fill_the_common_period(self):
        for shell in osv.solve_resonant_shells(24, 14, 17):
            k = shell["revs_per_common_period"]
            total = osv.period_from_a(shell["semi_major_axis_km"]) * k
            self.assertAlmostEqual(total, 24 * 3600.0, places=6)

    def test_more_revolutions_means_lower_orbit(self):
        shells = osv.solve_resonant_shells(24, 14, 18)
        altitudes = [s["altitude_km"] for s in shells]
        self.assertEqual(altitudes, sorted(altitudes, reverse=True))

    def test_flags_shells_below_the_surface(self):
        shells = {s["revs_per_common_period"]: s
                  for s in osv.solve_resonant_shells(24, 1, 40)}
        self.assertTrue(shells[1]["above_surface"])
        self.assertFalse(shells[40]["above_surface"])

    def test_one_rev_per_sidereal_day_is_geo(self):
        shell = osv.solve_resonant_shells(86164.0905 / 3600.0, 1, 1)[0]
        self.assertAlmostEqual(shell["altitude_km"], 35786.0, delta=1.0)

    def test_rejects_bad_ranges(self):
        with self.assertRaises(ValueError):
            osv.solve_resonant_shells(24, 0, 5)
        with self.assertRaises(ValueError):
            osv.solve_resonant_shells(24, 10, 5)
        with self.assertRaises(ValueError):
            osv.solve_resonant_shells(-1, 1, 5)


class TestFormationGeometry(unittest.TestCase):
    def test_two_to_one_ellipse_ratio(self):
        f = osv.solve_formation_geometry(550, 1000)
        self.assertAlmostEqual(f["along_track_amplitude_m"], 1000.0, places=6)
        self.assertAlmostEqual(f["radial_amplitude_m"], 500.0, places=6)
        self.assertAlmostEqual(
            f["along_track_amplitude_m"] / f["radial_amplitude_m"], 2.0, places=12)
        self.assertAlmostEqual(f["required_eccentricity"], 7.217e-5, delta=1e-8)

    def test_offsets_stay_inside_the_ellipse(self):
        h, target = 550.0, 1000.0
        e = osv.solve_formation_geometry(h, target)["required_eccentricity"]
        period = osv.period_from_a(osv.R_E + h)
        for i in range(24):
            p = osv.formation_offset(h, e, period * i / 24.0)
            x, y = p["radial_dx_m"], p["along_track_dy_m"]
            self.assertLessEqual((x / 500.0) ** 2 + (y / 1000.0) ** 2, 1.0 + 1e-9)

    def test_relative_trajectory_closes_after_one_period(self):
        h = 550.0
        e = osv.solve_formation_geometry(h, 1000)["required_eccentricity"]
        period = osv.period_from_a(osv.R_E + h)
        start = osv.formation_offset(h, e, 0.0)
        end = osv.formation_offset(h, e, period)
        self.assertAlmostEqual(start["radial_dx_m"], end["radial_dx_m"], places=6)
        self.assertAlmostEqual(start["along_track_dy_m"], end["along_track_dy_m"],
                               places=6)

    def test_quarter_period_swaps_radial_for_along_track(self):
        h = 550.0
        e = osv.solve_formation_geometry(h, 1000)["required_eccentricity"]
        quarter = osv.period_from_a(osv.R_E + h) / 4.0
        p = osv.formation_offset(h, e, quarter)
        self.assertAlmostEqual(p["radial_dx_m"], 0.0, places=6)
        self.assertAlmostEqual(p["along_track_dy_m"], -1000.0, places=3)


class TestLineOfSight(unittest.TestCase):
    def test_matches_closed_form(self):
        los = osv.solve_los_range(550, 1000)
        r_t = osv.R_E + 80.0
        expected = (math.sqrt((osv.R_E + 550) ** 2 - r_t ** 2)
                    + math.sqrt((osv.R_E + 1000) ** 2 - r_t ** 2))
        self.assertAlmostEqual(los["max_range_km"], expected, places=9)

    def test_max_range_chord_grazes_the_margin_sphere(self):
        """最大链路的连线到地心的最近距离应恰好等于 R_E + margin。"""
        h1, h2, margin = 550.0, 1000.0, 80.0
        los = osv.solve_los_range(h1, h2, margin)
        r1, r2 = osv.R_E + h1, osv.R_E + h2
        theta = math.radians(los["max_central_angle_deg"])
        # 两星位置向量，夹角为最大地心张角
        p1 = (r1, 0.0)
        p2 = (r2 * math.cos(theta), r2 * math.sin(theta))
        dx, dy = p2[0] - p1[0], p2[1] - p1[1]
        seg_len = math.hypot(dx, dy)
        self.assertAlmostEqual(seg_len, los["max_range_km"], places=6)
        # 点到直线距离 = |cross(p1, p2)| / |p2 - p1|
        distance = abs(p1[0] * p2[1] - p1[1] * p2[0]) / seg_len
        self.assertAlmostEqual(distance, osv.R_E + margin, places=6)

    def test_higher_orbit_sees_farther(self):
        low = osv.solve_los_range(550, 550)["max_range_km"]
        high = osv.solve_los_range(1200, 1200)["max_range_km"]
        self.assertGreater(high, low)

    def test_larger_margin_shortens_the_link(self):
        self.assertLess(osv.solve_los_range(550, 550, 200)["max_range_km"],
                        osv.solve_los_range(550, 550, 0)["max_range_km"])

    def test_rejects_satellite_below_margin(self):
        with self.assertRaises(ValueError):
            osv.solve_los_range(50, 550, 80)


class TestInputValidation(unittest.TestCase):
    def test_rejects_bad_altitude(self):
        for bad in (0.0, -100.0):
            with self.assertRaises(ValueError):
                osv.get_orbit_features(bad, 53)

    def test_rejects_bad_inclination(self):
        for bad in (-1.0, 181.0):
            with self.assertRaises(ValueError):
                osv.get_orbit_features(550, bad)

    def test_rejects_bad_eccentricity(self):
        for bad in (-0.1, 1.0, 1.5):
            with self.assertRaises(ValueError):
                osv.get_orbit_features(550, 53, bad)

    def test_rejects_non_positive_formation_size(self):
        with self.assertRaises(ValueError):
            osv.solve_formation_geometry(550, 0)


class TestCommandLine(unittest.TestCase):
    def test_subcommands_exit_zero(self):
        import contextlib
        import io

        cases = [
            ["features", "550", "53"],
            ["match", "550", "53", "1000"],
            ["sunsync", "800"],
            ["resonant", "24", "15", "16"],
            ["formation", "550", "1000"],
            ["los", "550", "1000"],
            ["demo"],
            [],
        ]
        for argv in cases:
            with contextlib.redirect_stdout(io.StringIO()) as out:
                self.assertEqual(osv.main(argv), 0)
            self.assertTrue(out.getvalue().strip(), f"{argv} 无输出")

    def test_json_output_parses(self):
        import contextlib
        import io
        import json

        with contextlib.redirect_stdout(io.StringIO()) as out:
            osv.main(["--json", "match", "550", "53", "1000"])
        payload = json.loads(out.getvalue())
        self.assertTrue(payload["feasible"])
        self.assertAlmostEqual(payload["matched_inc2_deg"], 41.40, delta=0.01)


if __name__ == "__main__":
    unittest.main(verbosity=2)
