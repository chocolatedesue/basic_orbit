"""cluster_solver 的解析校验用例（纯 unittest，无第三方依赖）。"""

import contextlib
import io
import math
import unittest

import cluster_solver as cs
import orbit_solver as osv


def _precession(shell):
    return osv.rad_s_to_deg_day(
        osv.nodal_precession_rate(shell.semi_major_axis_km, shell.inclination_deg))


class TestShell(unittest.TestCase):
    def test_derived_properties(self):
        s = cs.Shell("A", 550, 53, planes=6, sats_per_plane=12, phasing_f=1)
        self.assertAlmostEqual(s.semi_major_axis_km, osv.R_E + 550)
        self.assertEqual(s.num_satellites, 72)
        self.assertEqual(s.walker_notation, "53.00°: 72/6/1")

    def test_rejects_bad_constellation_counts(self):
        with self.assertRaises(ValueError):
            cs.Shell("A", 550, 53, planes=0)
        with self.assertRaises(ValueError):
            cs.Shell("A", 550, 53, sats_per_plane=0)

    def test_rejects_out_of_range_phasing(self):
        with self.assertRaises(ValueError):
            cs.Shell("A", 550, 53, planes=6, phasing_f=6)
        with self.assertRaises(ValueError):
            cs.Shell("A", 550, 53, planes=6, phasing_f=-1)

    def test_rejects_bad_raan_span(self):
        with self.assertRaises(ValueError):
            cs.Shell("A", 550, 53, raan_span_deg=0.0)
        with self.assertRaises(ValueError):
            cs.Shell("A", 550, 53, raan_span_deg=361.0)

    def test_inherits_orbit_validation(self):
        with self.assertRaises(ValueError):
            cs.Shell("A", -10, 53)
        with self.assertRaises(ValueError):
            cs.Shell("A", 550, 200)


class TestPlaneGeometry(unittest.TestCase):
    def test_equator_separation_equals_raan_difference(self):
        for inc in (30.0, 53.0, 90.0):
            self.assertAlmostEqual(
                cs.central_angle_between_planes(inc, 25.0, 0.0), 25.0, places=9)

    def test_planes_converge_at_high_latitude(self):
        equator = cs.central_angle_between_planes(53.0, 30.0, 0.0)
        apex = cs.central_angle_between_planes(53.0, 30.0, 90.0)
        self.assertLess(apex, equator)

    def test_polar_planes_meet_at_the_pole(self):
        self.assertAlmostEqual(
            cs.central_angle_between_planes(90.0, 40.0, 90.0), 0.0, places=6)

    def test_equatorial_planes_never_converge(self):
        for u in (0.0, 30.0, 90.0):
            self.assertAlmostEqual(
                cs.central_angle_between_planes(0.0, 20.0, u), 20.0, places=9)

    def test_chord_geometry(self):
        a = 7000.0
        self.assertAlmostEqual(cs.chord_km(a, 180.0), 2 * a, places=9)
        self.assertAlmostEqual(cs.chord_km(a, 60.0), a, places=9)
        self.assertAlmostEqual(cs.chord_km(a, 0.0), 0.0, places=9)


class TestIntraShellGeometry(unittest.TestCase):
    def test_in_plane_spacing_matches_chord(self):
        s = cs.Shell("A", 550, 53, planes=24, sats_per_plane=20)
        g = cs.intra_shell_geometry(s)
        self.assertAlmostEqual(g["in_plane_spacing_deg"], 18.0, places=9)
        self.assertAlmostEqual(g["in_plane_spacing_km"],
                               cs.chord_km(s.semi_major_axis_km, 18.0), places=9)
        self.assertAlmostEqual(g["raan_step_deg"], 15.0, places=9)

    def test_single_plane_has_no_cross_plane_links(self):
        g = cs.intra_shell_geometry(cs.Shell("A", 550, 53, sats_per_plane=8))
        self.assertIsNone(g["cross_plane_max_km"])
        self.assertIsNone(g["cross_plane_min_km"])

    def test_cross_plane_min_below_max(self):
        g = cs.intra_shell_geometry(cs.Shell("A", 550, 53, planes=12,
                                             sats_per_plane=20))
        self.assertLess(g["cross_plane_min_km"], g["cross_plane_max_km"])

    def test_dense_shell_links_fit_inside_los(self):
        g = cs.intra_shell_geometry(cs.Shell("A", 550, 53, planes=24,
                                             sats_per_plane=20))
        self.assertTrue(g["all_links_within_los"])

    def test_sparse_shell_breaks_los(self):
        g = cs.intra_shell_geometry(cs.Shell("A", 550, 53, planes=3,
                                             sats_per_plane=3))
        self.assertFalse(g["all_links_within_los"])


class TestPairAnalysis(unittest.TestCase):
    def test_synodic_period_of_15_and_14_revs_is_one_day(self):
        shells = cs.design_resonant_locked_cluster(53.0, 24.0, [15, 14])["shells"]
        self.assertAlmostEqual(cs.synodic_period(*shells), 86400.0, places=3)

    def test_same_altitude_freezes_relative_phase(self):
        a = cs.Shell("A", 550, 53)
        b = cs.Shell("B", 550, 40)
        self.assertEqual(cs.synodic_period(a, b), float("inf"))

    def test_matched_pair_is_locked(self):
        base = cs.Shell("base", 550, 53)
        design = cs.design_locked_cluster(base, [1000])
        p = cs.pair_analysis(*design["shells"])
        self.assertTrue(p["planes_locked"])
        self.assertEqual(p["days_to_drift_budget"], float("inf"))
        self.assertLess(abs(p["differential_precession_deg_day"]),
                        cs.LOCK_TOLERANCE_DEG_DAY)

    def test_unmatched_pair_drifts_on_a_finite_clock(self):
        p = cs.pair_analysis(cs.Shell("A", 550, 53), cs.Shell("B", 1200, 53))
        self.assertFalse(p["planes_locked"])
        self.assertTrue(math.isfinite(p["days_to_drift_budget"]))
        # 差分进动 1.209 度/天 → 张开 1 度约 0.83 天
        self.assertAlmostEqual(p["days_to_drift_budget"], 0.83, delta=0.02)

    def test_min_range_is_the_radial_gap(self):
        p = cs.pair_analysis(cs.Shell("A", 550, 53), cs.Shell("B", 1200, 53))
        self.assertAlmostEqual(p["min_range_km"], 650.0, places=9)
        self.assertGreater(p["los_limit_km"], p["min_range_km"])


class TestCommonRepeat(unittest.TestCase):
    def test_resonant_ladder_has_zero_residual(self):
        shells = cs.design_resonant_locked_cluster(53.0, 24.0, [15, 14, 13])["shells"]
        rep = cs.check_common_repeat(shells, 24.0)
        self.assertTrue(rep["cluster_repeats"])
        for row in rep["shells"]:
            self.assertAlmostEqual(row["residual_phase_deg"], 0.0, places=6)
            self.assertEqual(row["nearest_integer_revs"],
                             int(row["name"].split("=")[1]))

    def test_arbitrary_altitudes_do_not_repeat(self):
        shells = [cs.Shell("A", 550, 53), cs.Shell("B", 800, 53)]
        rep = cs.check_common_repeat(shells, 24.0)
        self.assertFalse(rep["cluster_repeats"])
        self.assertTrue(any(abs(r["residual_phase_deg"]) > 1.0
                            for r in rep["shells"]))

    def test_residual_predicts_lap_time(self):
        rep = cs.check_common_repeat([cs.Shell("A", 550, 53)], 24.0)
        row = rep["shells"][0]
        self.assertAlmostEqual(
            row["cycles_to_full_lap"] * abs(row["residual_phase_deg"]), 360.0,
            places=6)

    def test_rejects_bad_input(self):
        with self.assertRaises(ValueError):
            cs.check_common_repeat([cs.Shell("A", 550, 53)], 0.0)
        with self.assertRaises(ValueError):
            cs.check_common_repeat([], 24.0)


class TestAnalyzeCluster(unittest.TestCase):
    def setUp(self):
        self.naive = [
            cs.Shell("低层", 550, 53, planes=24, sats_per_plane=20),
            cs.Shell("中层", 800, 53, planes=18, sats_per_plane=16),
            cs.Shell("高层", 1200, 53, planes=12, sats_per_plane=12),
        ]

    def test_counts_every_satellite(self):
        a = cs.analyze_cluster(self.naive)
        self.assertEqual(a["num_shells"], 3)
        self.assertEqual(a["num_satellites"], 480 + 288 + 144)

    def test_same_inclination_cluster_is_not_locked(self):
        a = cs.analyze_cluster(self.naive)
        self.assertFalse(a["planes_fully_locked"])
        self.assertAlmostEqual(a["precession_spread_deg_day"], 1.209, delta=0.01)
        self.assertAlmostEqual(a["worst_case_drift_days_per_deg"], 0.83, delta=0.02)

    def test_pair_count_is_combinatorial(self):
        self.assertEqual(len(cs.analyze_cluster(self.naive)["pairs"]), 3)

    def test_designed_cluster_is_locked(self):
        shells = cs.design_resonant_locked_cluster(53.0, 24.0, [15, 14, 13])["shells"]
        a = cs.analyze_cluster(shells, common_period_hours=24.0)
        self.assertTrue(a["planes_fully_locked"])
        self.assertTrue(a["common_repeat"]["cluster_repeats"])
        self.assertEqual(a["worst_case_drift_days_per_deg"], float("inf"))

    def test_single_shell_has_no_spread(self):
        a = cs.analyze_cluster([cs.Shell("A", 550, 53)])
        self.assertEqual(a["precession_spread_deg_day"], 0.0)
        self.assertTrue(a["planes_fully_locked"])
        self.assertEqual(a["pairs"], [])

    def test_common_repeat_is_optional(self):
        self.assertIsNone(cs.analyze_cluster(self.naive)["common_repeat"])

    def test_rejects_empty_and_duplicate_names(self):
        with self.assertRaises(ValueError):
            cs.analyze_cluster([])
        with self.assertRaises(ValueError):
            cs.analyze_cluster([cs.Shell("A", 550, 53), cs.Shell("A", 800, 47)])


class TestClusterDesign(unittest.TestCase):
    def test_locked_design_shares_one_precession_rate(self):
        base = cs.Shell("base", 550, 53)
        design = cs.design_locked_cluster(base, [700, 900, 1100])
        self.assertEqual(len(design["shells"]), 4)
        rates = [_precession(s) for s in design["shells"]]
        self.assertAlmostEqual(max(rates) - min(rates), 0.0, places=12)

    def test_locked_design_rejects_above_ceiling(self):
        base = cs.Shell("base", 550, 53)
        design = cs.design_locked_cluster(base, [900, 5000])
        self.assertEqual(len(design["shells"]), 2)
        self.assertEqual(len(design["rejected"]), 1)
        self.assertAlmostEqual(design["rejected"][0]["altitude_km"], 5000.0)
        self.assertAlmostEqual(design["ceiling_km"], 1631.8, delta=0.5)

    def test_resonant_design_satisfies_both_constraints(self):
        """核心断言：高度整数共振 与 倾角进动锁定 可以同时成立。"""
        design = cs.design_resonant_locked_cluster(53.0, 24.0, [15, 14, 13])
        shells = design["shells"]
        self.assertEqual(len(shells), 3)

        rates = [_precession(s) for s in shells]
        self.assertAlmostEqual(max(rates) - min(rates), 0.0, places=12)

        for shell, k in zip(shells, (15, 14, 13)):
            period = osv.period_from_a(shell.semi_major_axis_km)
            self.assertAlmostEqual(period * k, 24 * 3600.0, places=6)

    def test_resonant_design_orders_high_k_lowest(self):
        shells = cs.design_resonant_locked_cluster(53.0, 24.0, [13, 15, 14])["shells"]
        self.assertEqual([s.name for s in shells], ["k=15", "k=14", "k=13"])
        altitudes = [s.altitude_km for s in shells]
        self.assertEqual(altitudes, sorted(altitudes))

    def test_resonant_design_drops_shells_above_ceiling(self):
        design = cs.design_resonant_locked_cluster(53.0, 24.0, [15, 14, 13, 11, 9])
        kept = {s.name for s in design["shells"]}
        self.assertIn("k=15", kept)
        self.assertTrue(design["rejected"])
        for r in design["rejected"]:
            self.assertGreater(r["altitude_km"], design["ceiling_km"])

    def test_polar_base_admits_a_much_taller_cluster(self):
        design = cs.design_resonant_locked_cluster(88.0, 24.0, [15, 14, 13, 11, 9])
        self.assertEqual(design["rejected"], [])
        self.assertGreater(design["ceiling_km"], 8000.0)

    def test_carries_constellation_layout_to_every_shell(self):
        design = cs.design_resonant_locked_cluster(
            53.0, 24.0, [15, 14], planes=8, sats_per_plane=10)
        for s in design["shells"]:
            self.assertEqual(s.planes, 8)
            self.assertEqual(s.num_satellites, 80)

    def test_rejects_bad_k_values(self):
        with self.assertRaises(ValueError):
            cs.design_resonant_locked_cluster(53.0, 24.0, [])
        with self.assertRaises(ValueError):
            cs.design_resonant_locked_cluster(53.0, 24.0, [15, 15])
        with self.assertRaises(ValueError):
            cs.design_resonant_locked_cluster(53.0, 24.0, [15, 40])


class TestReport(unittest.TestCase):
    def test_report_covers_every_section(self):
        shells = cs.design_resonant_locked_cluster(53.0, 24.0, [15, 14])["shells"]
        text = cs.format_report(cs.analyze_cluster(shells, common_period_hours=24.0))
        for marker in ("壳层清单", "面共动", "壳层间关系", "链路可达", "相位复位"):
            self.assertIn(marker, text)
        self.assertIn("已锁定", text)

    def test_report_flags_a_drifting_cluster(self):
        text = cs.format_report(cs.analyze_cluster(
            [cs.Shell("A", 550, 53), cs.Shell("B", 1200, 53)]))
        self.assertIn("未锁定", text)


class TestCommandLine(unittest.TestCase):
    def test_subcommands_exit_zero(self):
        cases = [
            ["analyze", "A:550:53:24:20", "B:800:47.05:18:16"],
            ["analyze", "A:550:53", "--common-period", "24"],
            ["lock", "550", "53", "800", "1200"],
            ["design", "53", "24", "15", "14", "13"],
            ["design", "53", "24", "15", "14", "--planes", "8",
             "--sats-per-plane", "10"],
            ["demo"],
            [],
        ]
        for argv in cases:
            with contextlib.redirect_stdout(io.StringIO()) as out:
                self.assertEqual(cs.main(argv), 0)
            self.assertTrue(out.getvalue().strip(), f"{argv} 无输出")

    def test_shell_spec_parsing(self):
        s = cs._parse_shell_spec("低层:550:53:24:20")
        self.assertEqual((s.name, s.altitude_km, s.inclination_deg), ("低层", 550, 53))
        self.assertEqual(s.num_satellites, 480)

    def test_rejects_malformed_shell_specs(self):
        import argparse
        for bad in ("A:550", "A:550:53:24", "A:abc:53", "A:-5:53", "A:550:53:0:20"):
            with self.assertRaises(argparse.ArgumentTypeError, msg=bad):
                cs._parse_shell_spec(bad)


if __name__ == "__main__":
    unittest.main(verbosity=2)
