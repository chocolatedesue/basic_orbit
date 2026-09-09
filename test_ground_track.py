"""ground_track 的解析校验用例（纯 unittest，无第三方依赖）。"""

import contextlib
import io
import math
import unittest

import ground_track as gt
import orbit_solver as osv


class TestNodalRates(unittest.TestCase):
    def test_critical_inclination_zeroes_perigee_rate(self):
        crit = math.degrees(math.acos(math.sqrt(0.2)))
        self.assertAlmostEqual(gt.perigee_rate(osv.R_E + 800, crit), 0.0, places=15)

    def test_j2_lengthens_the_nodal_period_of_a_retrograde_sso(self):
        h = 650.0
        inc = osv.solve_sun_synchronous_inclination(h)["inclination_deg"]
        a = osv.R_E + h
        self.assertGreater(gt.nodal_period(a, inc), osv.period_from_a(a))
        self.assertAlmostEqual(
            gt.nodal_period(a, inc) - osv.period_from_a(a), 7.25, delta=0.1)

    def test_mean_anomaly_rate_reduces_to_n_without_j2(self):
        a = osv.R_E + 650
        n = osv.mean_motion(a)
        original_gt, original_osv = gt.J2, osv.J2
        try:
            gt.J2 = osv.J2 = 0.0
            self.assertAlmostEqual(gt.mean_anomaly_rate(a, 97.99), n, places=15)
            self.assertAlmostEqual(gt.nodal_period(a, 97.99),
                                   osv.period_from_a(a), places=9)
        finally:
            gt.J2, osv.J2 = original_gt, original_osv

    def test_sun_synchronous_nodal_day_is_one_solar_day(self):
        """SSO 的进动恰好抵消地球公转，交点日精确等于平太阳日。"""
        for h in (400.0, 650.0, 900.0):
            inc = osv.solve_sun_synchronous_inclination(h)["inclination_deg"]
            self.assertAlmostEqual(gt.nodal_day(osv.R_E + h, inc), 86400.0,
                                   delta=1.0)

    def test_non_precessing_orbit_gives_a_sidereal_day(self):
        sidereal = 2 * math.pi / gt.OMEGA_EARTH
        self.assertAlmostEqual(gt.nodal_day(osv.R_E + 650, 90.0), sidereal,
                               places=6)
        self.assertAlmostEqual(sidereal, 86164.1, delta=0.5)

    def test_revs_per_nodal_day_at_650km_sso(self):
        inc = osv.solve_sun_synchronous_inclination(650)["inclination_deg"]
        self.assertAlmostEqual(gt.revs_per_nodal_day(650, inc), 14.7166,
                               delta=0.001)

    def test_j2_correction_changes_the_rev_count_materially(self):
        """忽略 J2 会把圈数算成 14.735，足以选错重复周期。"""
        h = 650.0
        inc = osv.solve_sun_synchronous_inclination(h)["inclination_deg"]
        naive = 86400.0 / osv.period_from_a(osv.R_E + h)
        self.assertAlmostEqual(naive, 14.7347, delta=0.001)
        self.assertGreater(abs(naive - gt.revs_per_nodal_day(h, inc)), 0.015)

    def test_revs_decrease_with_altitude(self):
        inc_lo = osv.solve_sun_synchronous_inclination(500)["inclination_deg"]
        inc_hi = osv.solve_sun_synchronous_inclination(900)["inclination_deg"]
        self.assertGreater(gt.revs_per_nodal_day(500, inc_lo),
                           gt.revs_per_nodal_day(900, inc_hi))


class TestRepeatGroundTrack(unittest.TestCase):
    def test_solution_reproduces_the_target_ratio(self):
        for n, m in ((103, 7), (147, 10), (250, 17), (31, 2)):
            sol = gt.solve_repeat_altitude(n, m)
            actual = gt.revs_per_nodal_day(sol["altitude_km"],
                                           sol["inclination_deg"])
            self.assertAlmostEqual(actual, n / m, places=9)

    def test_103_over_7_sits_just_above_650km(self):
        sol = gt.solve_repeat_altitude(103, 7)
        self.assertAlmostEqual(sol["altitude_km"], 650.72, delta=0.05)
        self.assertAlmostEqual(sol["inclination_deg"], 97.989, delta=0.01)
        self.assertEqual(sol["repeat_cycle_days"], 7)

    def test_ground_track_spacing(self):
        sol = gt.solve_repeat_altitude(103, 7)
        self.assertAlmostEqual(sol["ground_track_spacing_deg"], 360 * 7 / 103,
                               places=9)

    def test_fixed_inclination_mode(self):
        sol = gt.solve_repeat_altitude(31, 2, sun_synchronous=False,
                                       inc_deg=53.0)
        self.assertAlmostEqual(sol["inclination_deg"], 53.0, places=12)
        self.assertAlmostEqual(
            gt.revs_per_nodal_day(sol["altitude_km"], 53.0), 31 / 2, places=9)

    def test_requires_inclination_when_not_sun_synchronous(self):
        with self.assertRaises(ValueError):
            gt.solve_repeat_altitude(103, 7, sun_synchronous=False)

    def test_rejects_bad_counts_and_unreachable_ratios(self):
        with self.assertRaises(ValueError):
            gt.solve_repeat_altitude(0, 7)
        with self.assertRaises(ValueError):
            gt.solve_repeat_altitude(103, 0)
        with self.assertRaises(ValueError):
            gt.solve_repeat_altitude(400, 7)   # 圈数过多，高度会在地表以下

    def test_find_options_includes_the_nearest_cycle(self):
        found = gt.find_repeat_options(650.0)
        pairs = [(o["revs_n"], o["days_m"]) for o in found["options"]]
        self.assertIn((103, 7), pairs)
        best = found["options"][0]
        self.assertEqual((best["revs_n"], best["days_m"]), (103, 7))
        self.assertLess(abs(best["altitude_shift_km"]), 1.0)

    def test_options_are_sorted_by_altitude_shift(self):
        opts = gt.find_repeat_options(650.0)["options"]
        shifts = [abs(o["altitude_shift_km"]) for o in opts]
        self.assertEqual(shifts, sorted(shifts))

    def test_options_are_in_lowest_terms(self):
        opts = gt.find_repeat_options(650.0)["options"]
        for o in opts:
            self.assertEqual(math.gcd(o["revs_n"], o["days_m"]), 1)


class TestContactWindow(unittest.TestCase):
    def test_half_angle_closed_form(self):
        h, eps = 650.0, 10.0
        w = gt.contact_window(h, eps)
        expected = math.degrees(math.acos(
            osv.R_E * math.cos(math.radians(eps)) / (osv.R_E + h))) - eps
        self.assertAlmostEqual(w["half_angle_deg"], expected, places=12)
        self.assertAlmostEqual(w["half_angle_deg"], 16.65, delta=0.02)

    def test_max_pass_duration_at_650km(self):
        self.assertAlmostEqual(gt.contact_window(650, 10)["max_pass_min"],
                               9.05, delta=0.05)

    def test_zero_mask_gives_the_horizon_geometry(self):
        w = gt.contact_window(650, 0.0)
        self.assertAlmostEqual(
            w["half_angle_deg"],
            math.degrees(math.acos(osv.R_E / (osv.R_E + 650))), places=12)

    def test_higher_mask_shortens_the_pass(self):
        self.assertLess(gt.contact_window(650, 20)["max_pass_min"],
                        gt.contact_window(650, 5)["max_pass_min"])

    def test_higher_orbit_widens_the_footprint(self):
        self.assertGreater(gt.contact_window(1200, 10)["footprint_radius_km"],
                           gt.contact_window(650, 10)["footprint_radius_km"])

    def test_rejects_bad_mask(self):
        for bad in (-1.0, 90.0, 120.0):
            with self.assertRaises(ValueError):
                gt.contact_window(650, bad)


class TestAccessEstimate(unittest.TestCase):
    def test_low_duty_cycle_for_a_single_station(self):
        acc = gt.passes_per_day_estimate(650, 60.0, 10.0)
        self.assertTrue(acc["in_coverage"])
        self.assertLess(acc["duty_cycle"], 0.05)
        self.assertAlmostEqual(acc["passes_per_day"], 5.4, delta=0.5)

    def test_high_latitude_station_sees_more_passes(self):
        low = gt.passes_per_day_estimate(650, 10.0)["passes_per_day"]
        high = gt.passes_per_day_estimate(650, 75.0)["passes_per_day"]
        self.assertGreater(high, low)

    def test_passes_are_capped_by_the_orbit_count(self):
        acc = gt.passes_per_day_estimate(650, 89.0)
        self.assertLessEqual(acc["passes_per_day"],
                             2.0 * acc["revs_per_nodal_day"] + 1e-9)

    def test_station_outside_the_reachable_band(self):
        acc = gt.passes_per_day_estimate(650, 40.0, inc_deg=5.0,
                                         sun_synchronous=False)
        self.assertFalse(acc["in_coverage"])
        self.assertEqual(acc["passes_per_day"], 0.0)

    def test_rejects_polar_station(self):
        with self.assertRaises(ValueError):
            gt.passes_per_day_estimate(650, 90.0)


class TestSpectrum(unittest.TestCase):
    def test_spectrum_lists_every_timescale(self):
        s = gt.periodicity_spectrum(650, station_latitude_deg=60)
        names = [t["name"] for t in s["timescales"]]
        self.assertTrue(any("交点周期" in n for n in names))
        self.assertTrue(any("重复星下点" in n for n in names))
        self.assertTrue(any("太阳同步年周期" in n for n in names))
        self.assertAlmostEqual(s["j2_period_correction_sec"], 7.25, delta=0.1)

    def test_spectrum_picks_the_nearest_repeat_by_default(self):
        s = gt.periodicity_spectrum(650)
        self.assertEqual((s["chosen_repeat"]["revs_n"],
                          s["chosen_repeat"]["days_m"]), (103, 7))

    def test_spectrum_honours_a_requested_repeat_cycle(self):
        s = gt.periodicity_spectrum(650, repeat_days=10)
        self.assertEqual(s["chosen_repeat"]["days_m"], 10)
        self.assertEqual(s["chosen_repeat"]["revs_n"], 147)

    def test_report_renders_every_section(self):
        text = gt.format_spectrum(gt.periodicity_spectrum(650))
        for marker in ("周期修正", "时间尺度谱", "地面接触", "可选重复星下点周期"):
            self.assertIn(marker, text)


class TestCommandLine(unittest.TestCase):
    def test_subcommands_exit_zero(self):
        cases = [
            ["spectrum", "650"],
            ["spectrum", "650", "--station-latitude", "78",
             "--elevation-mask", "5", "--repeat-days", "10"],
            ["repeat", "103", "7"],
            ["repeat", "31", "2", "--inclination", "53"],
            ["contact", "650"],
            ["contact", "650", "--elevation-mask", "20"],
            ["demo"],
            [],
        ]
        for argv in cases:
            with contextlib.redirect_stdout(io.StringIO()) as out:
                self.assertEqual(gt.main(argv), 0)
            self.assertTrue(out.getvalue().strip(), f"{argv} 无输出")

    def test_repeat_command_reports_the_altitude(self):
        with contextlib.redirect_stdout(io.StringIO()) as out:
            gt.main(["repeat", "103", "7"])
        self.assertIn("650.7", out.getvalue())


if __name__ == "__main__":
    unittest.main(verbosity=2)
