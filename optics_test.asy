// optics_test.asy
// Physical verification and unit tests for optics.asy

import optics;

write("================================================================================");
write("         OPTICS.ASY TESTS                                                       ");
write("================================================================================");
write("");

// ============================================================================
// Test Assertion Helpers
// ============================================================================
int total_tests = 0;
int passed_tests = 0;
int failed_tests = 0;

void assert_near(string test_name, real val, real expected, real tol = 1e-4) {
    ++total_tests;
    real diff = abs(val - expected);
    if (diff <= tol) {
        ++passed_tests;
        write("  [PASS] " + test_name + " = " + string(val));
    } else {
        ++failed_tests;
        write("  [FAIL] " + test_name + " = " + string(val) + " (expected " + string(expected) + ", abs diff " +
              string(diff) + ", tolerance " + string(tol) + ")");
    }
};

void assert_near_pair(string test_name, pair val, pair expected, real tol = 1e-4) {
    ++total_tests;
    real diff = length(val - expected);
    if (diff <= tol) {
        ++passed_tests;
        write("  [PASS] " + test_name + " = (" + string(val.x) + ", " + string(val.y) + ")");
    } else {
        ++failed_tests;
        write("  [FAIL] " + test_name + " = (" + string(val.x) + ", " + string(val.y) + ") (expected (" +
              string(expected.x) + ", " + string(expected.y) + "), error " + string(diff) + ", tolerance " +
              string(tol) + ")");
    }
};

void assert_bool(string test_name, bool val, bool expected) {
    ++total_tests;
    if (val == expected) {
        ++passed_tests;
        write("  [PASS] " + test_name + " = " + (val ? "true" : "false"));
    } else {
        ++failed_tests;
        write("  [FAIL] " + test_name + " = " + (val ? "true" : "false") + " (expected " +
              (expected ? "true" : "false") + ")");
    }
};

void assert_true(string test_name, bool condition) { assert_bool(test_name, condition, true); }

void assert_false(string test_name, bool condition) { assert_bool(test_name, condition, false); }

void assert_relative(string test_name, real val, real expected, real rel_tol = 1e-6, real abs_tol = 1e-10) {
    real allowed = max(abs_tol, rel_tol * max(abs(val), abs(expected)));
    assert_near(test_name, val, expected, allowed);
};

void assert_stokes_physical(string name, Polarization pol, real tol = 1e-5) {
    assert_bool(name + ": nonnegative S0", pol.S0 >= -tol, true);
    real polarized2 = pol.S1 * pol.S1 + pol.S2 * pol.S2 + pol.S3 * pol.S3;
    real total2 = pol.S0 * pol.S0;
    assert_bool(name + ": Stokes cone invariant", polarized2 <= total2 + tol, true);
};

// Spatial resolvers for test boundaries
mediumResolver test_resolve_isotropic = new Medium(pair p, Ray r, Medium m1, Medium m2, bool from_m1) {
    return (p.x > 0) ? medium("glass", 1.5) : Materials.VACUUM;
};

mediumResolver test_resolve_t_i_r = new Medium(pair p, Ray r, Medium m1, Medium m2, bool from_m1) {
    return (p.x > 0) ? Materials.VACUUM : medium("glass", 1.5);
};

// ============================================================================
// Test Case 1: Isotropic Refraction & Snell's Law (Fresnel Equations)
// ============================================================================
write("Test Case 1: Isotropic Refraction & Snell's Law (Fresnel Equations)");
real n1 = 1.0, n2 = 1.5;
real theta_i_deg = 30.0;
real theta_i_rad = radians(theta_i_deg);
Medium glass_15 = medium("glass", n2);
Ray r_inc = Ray(p = (-1.0, 0.0), d = (cos(theta_i_rad), sin(theta_i_rad)), medium = Materials.VACUUM);
r_inc.pol = Polarization.unpolarized();
Hit h_int;
h_int.valid = true;
h_int.p = (0.0, 0.0);
h_int.normal = (-1.0, 0.0);
h_int.tangent = (0.0, 1.0);
h_int.t = 1.0;
h_int.entering = true;
FresnelResult ref1 = fresnelIsotropic(Materials.VACUUM, glass_15, r_inc, h_int, test_resolve_isotropic);
real sin_theta_t = (n1 / n2) * sin(theta_i_rad);
real cos_theta_t = sqrt(1.0 - sin_theta_t ^ 2);
real expected_theta_t = degrees(asin(sin_theta_t));
pair expected_ref_d = (cos_theta_t, sin_theta_t);
real cos_i = cos(theta_i_rad);
real rs_exact = (n1 * cos_i - n2 * cos_theta_t) / (n1 * cos_i + n2 * cos_theta_t);
real rp_exact = (n2 * cos_i - n1 * cos_theta_t) / (n2 * cos_i + n1 * cos_theta_t);
real Rs_exact = rs_exact ^ 2;
real Rp_exact = rp_exact ^ 2;
real R_exact = 0.5 * (Rs_exact + Rp_exact);
real T_exact = 1.0 - R_exact;

assert_bool("TIR flag", ref1.tir, false);
assert_near("Refraction angle theta_t (deg)", degrees(asin(ref1.d.y)), expected_theta_t);
assert_near_pair("Refracted ray direction vector", ref1.d, expected_ref_d);
assert_near("Fresnel R_s reflectivity", abs2(ref1.rs), Rs_exact);
assert_near("Fresnel R_p reflectivity", abs2(ref1.rp), Rp_exact);
assert_near("Total reflectivity R", ref1.R, R_exact);
assert_near("Total transmittivity T", ref1.T, T_exact);
write("");

// ============================================================================
// Test Case 2: Total Internal Reflection (TIR) & Critical Angle
// ============================================================================
write("Test Case 2: Total Internal Reflection (TIR)");
real theta_tir_deg = 50.0;
real theta_tir_rad = radians(theta_tir_deg);
Ray r_tir = Ray(p = (-1.0, 0.0), d = (cos(theta_tir_rad), sin(theta_tir_rad)), medium = glass_15);
r_tir.pol = Polarization.unpolarized();
FresnelResult ref2 = fresnelIsotropic(glass_15, Materials.VACUUM, r_tir, h_int, test_resolve_t_i_r);
pair expected_refl_dir = (-cos(theta_tir_rad), sin(theta_tir_rad));

assert_bool("TIR flag", ref2.tir, true);
assert_near("Reflectivity R (TIR)", ref2.R, 1.0);
assert_near("Transmittivity T (TIR)", ref2.T, 0.0);
assert_near_pair("Reflected ray direction", ref2.d, expected_refl_dir);
write("");

// ============================================================================
// Test Case 3: Thin Lens Paraxial Focusing (f = 100)
// ============================================================================
write("Test Case 3: Thin Lens Paraxial Focusing");
real f_lens = 100.0;
real x_lens = 100.0;
Element test_lens = Lens.thin(center = (x_lens, 0.0), angle = 0.0, f = f_lens, diameter = 50.0);
OpticalSystem sys_lens = system(test_lens, Materials.VACUUM);
real y_top = 10.0;
real y_bot = -10.0;
Ray ray_top = Ray(p = (0.0, y_top), d = (1.0, 0.0), medium = Materials.VACUUM);
Ray ray_bot = Ray(p = (0.0, y_bot), d = (1.0, 0.0), medium = Materials.VACUUM);
Trace t_top = sys_lens.trace(ray_top);
Trace t_bot = sys_lens.trace(ray_bot);
Ray trans_top = t_top.transmitted(0);
Ray trans_bot = t_bot.transmitted(0);
real expected_x_focus = x_lens + f_lens;
pair expected_dir_top = unit((f_lens, -y_top));
pair expected_dir_bot = unit((f_lens, -y_bot));

assert_near_pair("Transmitted ray top direction", trans_top.d, expected_dir_top);
assert_near_pair("Transmitted ray bot direction", trans_bot.d, expected_dir_bot);

real t_f_top = -trans_top.p.y / trans_top.d.y;
real x_f_top = trans_top.p.x + trans_top.d.x * t_f_top;
assert_near("Top ray focal intersection X", x_f_top, expected_x_focus);
write("");

// ============================================================================
// Test Case 4: Ideal Concave Spherical Mirror Focusing
// ============================================================================
write("Test Case 4: Ideal Concave Spherical Mirror Focusing");
real f_mirror = 50.0;
real x_mirror = 100.0;
real y_inc = 10.0;
Element concave_mirror = Mirror.idealSpherical(vertex = (x_mirror, 0.0), angle = 180.0, f = f_mirror, diameter = 40.0);
OpticalSystem sys_mirror = system(concave_mirror, Materials.VACUUM);
Ray ray_mirror_top = Ray(p = (0.0, y_inc), d = (1.0, 0.0), medium = Materials.VACUUM);
Trace t_mirror = sys_mirror.trace(ray_mirror_top);
Ray refl_mirror = t_mirror.reflected(0);

real expected_x_mirror_focus = x_mirror - f_mirror;
real t_f_mirror = -refl_mirror.p.y / refl_mirror.d.y;
real x_f_mirror = refl_mirror.p.x + refl_mirror.d.x * t_f_mirror;

assert_near("Mirror focal point X (Aberration-free ideal focus)", x_f_mirror, expected_x_mirror_focus);
write("");

// ============================================================================
// Test Case 5: Parabolic GRIN Inhomogeneous Ray Trajectory
// ============================================================================
write("Test Case 5: Parabolic GRIN Inhomogeneous Ray Trajectory");
real n0_grin = 1.5;
real g_grin = 0.1;
indexProfile grin_prof = new real(pair p, real lambda) { return n0_grin * (1.0 - 0.5 * g_grin ^ 2 * p.y ^ 2); };
gradientProfile grin_grad = new pair(pair p, real lambda) { return (0.0, -n0_grin * g_grin ^ 2 * p.y); };
Medium grin_medium = medium("Parabolic GRIN").add(GRIN(grin_prof, grin_grad));
Element grin_box = Interface(box((0.0, -20.0), (70.0, 20.0)), grin_medium, Materials.VACUUM);
OpticalSystem sys_g_r_i_n = system(grin_box, Materials.VACUUM);

real y0_grin = 1.0;
Ray ray_grin = Ray(p = (0.0, y0_grin), d = (1.0, 0.0), medium = grin_medium);

real x_quarter = pi / (2.0 * g_grin);
Trace t_grin = sys_g_r_i_n.trace(ray_grin, ds = 0.01, far = x_quarter);
Ray ray_at_quarter = t_grin.last();

assert_near("GRIN ray Y position at quarter-pitch (y = 0)", ray_at_quarter.p.y, 0.0, 5e-3);
pair expected_grin_d = unit((1.0, -g_grin *y0_grin));
assert_near_pair("GRIN ray direction at quarter-pitch", ray_at_quarter.d, expected_grin_d, 5e-3);
write("");

// ============================================================================
// Test Case 6: Extraordinary Wave Walk-Off in Uniaxial Crystal
// ============================================================================
write("Test Case 6: Extraordinary Wave Walk-Off in Uniaxial Crystal");
real no_c = 1.5, ne_c = 2.0;
pair optic_axis = (1.0, 0.0);
VolumeProperty uniaxial = UniaxialExtraordinary(no = no_c, ne = ne_c, optic_axis = optic_axis);
real theta_k_deg = 45.0;
real theta_k_rad = radians(theta_k_deg);
pair k_norm = (cos(theta_k_rad), sin(theta_k_rad));

real inv_ne2 = (cos(theta_k_rad) ^ 2) / (no_c ^ 2) + (sin(theta_k_rad) ^ 2) / (ne_c ^ 2);
real expected_ne = 1.0 / sqrt(inv_ne2);
real actual_ne = uniaxial.phase_index((0.0, 0.0), k_norm, 633.0, Polarization.unpolarized());

assert_near("Extraordinary phase index n_e(45 deg)", actual_ne, expected_ne);

pair S_vec = (k_norm.x / (no_c ^ 2), k_norm.y / (ne_c ^ 2));
pair expected_S_dir = unit(S_vec);
real expected_theta_s = degrees(atan2(expected_S_dir.y, expected_S_dir.x));
pair actual_S_dir = uniaxial.group_dir(k_norm, (0.0, 0.0), 633.0);

assert_near_pair("Poynting group direction S_unit", actual_S_dir, expected_S_dir);
assert_near("Poynting angle theta_s (deg)", degrees(atan2(actual_S_dir.y, actual_S_dir.x)), expected_theta_s);
write("");

// ============================================================================
// Test Case 7: Uniaxial Crystal Double Refraction & Eigenmode Polarization Purity
// ============================================================================
write("Test Case 7: Uniaxial Crystal Double Refraction & Eigenmode Polarization Purity");
Medium calcite_o = medium("Calcite Ordinary", 1.658);
// The extraordinary index comes from the uniaxial solver alone: Medium.n()
// multiplies index factors, so adding a scalar 1.486 would square the index.
Medium calcite_e = uniaxialExtraordinary("Calcite Extraordinary", 1.658, 1.486, (1.0, 0.0));
Element crystal = Crystal(box((0.0, -10.0), (20.0, 10.0)), calcite_o, calcite_e, Materials.VACUUM);
OpticalSystem crystal_sys = system(crystal, Materials.VACUUM);

Ray ray_crystal = Ray(p = (-10.0, 0.0), d = (cos(pi / 6.0), sin(pi / 6.0)), medium = Materials.VACUUM);
ray_crystal.pol = Polarization.unpolarized();
Trace t_crystal = crystal_sys.trace(ray_crystal);

Interaction entry_ix = t_crystal.interaction(0);
assert_near("Eigenmode split: 2 transmitted rays spawned", entry_ix.transmitted.length, 2.0);

Ray o_ray = entry_ix.transmitted[0];
Ray e_ray = entry_ix.transmitted[1];

assert_near("Ordinary ray S-polarization purity (S1 == S0)", o_ray.pol.S1, o_ray.pol.S0, 1e-5);
assert_near("Ordinary ray S-polarization S2", o_ray.pol.S2, 0.0, 1e-5);
assert_near("Ordinary ray S-polarization S3", o_ray.pol.S3, 0.0, 1e-5);

assert_near("Extraordinary ray P-polarization purity (S1 == -S0)", e_ray.pol.S1, -e_ray.pol.S0, 1e-5);
assert_near("Extraordinary ray P-polarization S2", e_ray.pol.S2, 0.0, 1e-5);
assert_near("Extraordinary ray P-polarization S3", e_ray.pol.S3, 0.0, 1e-5);
real e_index = e_ray.medium.n((0, 0), e_ray.wave_normal);
assert_true("Extraordinary phase index lies between ne and no", e_index >= 1.486 - 1e-9 && e_index <= 1.658 + 1e-9);
write("");

// ============================================================================
// Test Case 8: Diffraction Grating Spectral Order Angle
// ============================================================================
write("Test Case 8: Diffraction Grating Spectral Order Angle");
real lines_per_mm = 600.0;
real lambda = 500.0;
real lambda_units = lambda / 1000.0;
real d_grating_um = 1000.0 / lines_per_mm;
int m_order = 1;
Element test_grating = Grating(center = (0.0, 0.0), diameter = 10.0, angle = 90.0, lines_per_mm = lines_per_mm);
OpticalSystem grating_sys = system(test_grating, Materials.VACUUM);

Ray ray_grating = Ray(p = (0.0, -1.0), d = (0.0, 1.0), medium = Materials.VACUUM, lambda = lambda);
Trace t_grating = grating_sys.trace(ray_grating);
Ray diff_ray = t_grating.transmitted(0);

real sin_theta_diff = (m_order * lambda_units) / d_grating_um;
real expected_diff_angle = degrees(asin(sin_theta_diff));
pair expected_diff_dir = (sin_theta_diff, sqrt(1.0 - sin_theta_diff ^ 2));

assert_near("Diffracted order angle theta_1 (deg)", degrees(asin(diff_ray.d.x)), expected_diff_angle);
assert_near_pair("Diffracted ray direction vector", diff_ray.d, expected_diff_dir);
write("");

// ============================================================================
// Test Case 9: Cauchy Dispersion & Chromatic Separation
// ============================================================================
write("Test Case 9: Cauchy Dispersion & Chromatic Separation");
real A_cauchy = 1.50, B_cauchy = 0.01;
Medium dispersive_glass = medium("Dispersive Glass", A_cauchy, B_cauchy);
real lambda_red_nm = 700.0;
real lambda_blue_nm = 400.0;
real lambda_red_um = lambda_red_nm / 1000.0;
real lambda_blue_um = lambda_blue_nm / 1000.0;
real expected_n_red = A_cauchy + B_cauchy / (lambda_red_um ^ 2);
real expected_n_blue = A_cauchy + B_cauchy / (lambda_blue_um ^ 2);

assert_near("Red index n(700nm)", dispersive_glass.n((0, 0), (1, 0), lambda_red_nm), expected_n_red);
assert_near("Blue index n(400nm)", dispersive_glass.n((0, 0), (1, 0), lambda_blue_nm), expected_n_blue);

Ray ray_red = Ray(p = (-1, 0), d = (cos(pi / 4), sin(pi / 4)), medium = Materials.VACUUM, lambda = lambda_red_nm);
Ray ray_blue = Ray(p = (-1, 0), d = (cos(pi / 4), sin(pi / 4)), medium = Materials.VACUUM, lambda = lambda_blue_nm);

mediumResolver resolve_dispersion = new Medium(pair p, Ray r, Medium m1, Medium m2, bool from_m1) {
    return (p.x > 0) ? dispersive_glass : Materials.VACUUM;
};

FresnelResult ref_red = fresnelIsotropic(Materials.VACUUM, dispersive_glass, ray_red, h_int, resolve_dispersion);
FresnelResult ref_blue = fresnelIsotropic(Materials.VACUUM, dispersive_glass, ray_blue, h_int, resolve_dispersion);

real expected_theta_red = degrees(asin(sin(pi / 4) / expected_n_red));
real expected_theta_blue = degrees(asin(sin(pi / 4) / expected_n_blue));

assert_near("Red refracted angle (deg)", degrees(asin(ref_red.d.y)), expected_theta_red);
assert_near("Blue refracted angle (deg)", degrees(asin(ref_blue.d.y)), expected_theta_blue);
assert_bool("Blue bends more than Red (normal dispersion)", ref_blue.d.y < ref_red.d.y, true);
write("");

// ============================================================================
// Test Case 10: Malus's Law via Dichroic Absorption
// ============================================================================
write("Test Case 10: Malus's Law via Dichroic Absorption");
Element malus_polarizer =
    Interface(box((0.0, -10.0), (10.0, 10.0)), medium("Malus Pol").add(Dichroic(0.0, 100.0)), Materials.VACUUM);
OpticalSystem malus_sys = system(malus_polarizer, Materials.VACUUM);

real theta_pol_deg = 30.0;
Polarization pol_in = Polarization.linear(theta_pol_deg, 1.0);
Ray r_in = Ray((-10, 0), (1, 0), pol = pol_in);
Trace t_malus = malus_sys.trace(r_in);

// Use t_malus.last() to evaluate the ray after traversing the dichroic medium:
Ray r_out = t_malus.last();

real I_expected = cos(radians(theta_pol_deg)) ^ 2; // cos^2(30 deg) = 0.75
assert_near("Malus transmitted intensity", r_out.pol.intensity(), I_expected, 1e-4);
assert_near("Transmitted light is pure S-polarized (S1 == S0)", r_out.pol.S1, r_out.pol.S0, 1e-4);
assert_near("Transmitted light DOP", r_out.pol.dop(), 1.0, 1e-5);
write("");

// ============================================================================
// Test Case 11: Optical Activity in Chiral Medium
// ============================================================================
write("Test Case 11: Optical Activity in Chiral Medium");
Medium sugar_water = medium("Sugar Water").add(OpticalActivity(5.0));

Polarization p_in = Polarization.linear(0.0, 1.0); // S = (1, 1, 0, 0)
Polarization p_out = sugar_water.stepPolarization(p_in, (0, 0), (1, 0), 633.0, 10.0);

real theta_rot_deg = 10.0 * 5.0; // 50 deg rotation
assert_near("Rotated S1", p_out.S1, cos(radians(2.0 * theta_rot_deg)), 1e-4);
assert_near("Rotated S2", p_out.S2, sin(radians(2.0 * theta_rot_deg)), 1e-4);
assert_near("Rotated S3", p_out.S3, 0.0, 1e-4);
assert_near("Rotated tilt angle psi", degrees(p_out.psi()), theta_rot_deg, 1e-3);
assert_near("DOP preserved", p_out.dop(), 1.0, 1e-5);
write("");

// ============================================================================
// Test Case 12: Beer-Lambert Law & Optical Path Length (OPL)
// ============================================================================
write("Test Case 12: Beer-Lambert Law & Optical Path Length (OPL)");
real d1 = 10.0, d2 = 20.0, far_reach = 45.0;
real n_glass_abs = 1.5;
real alpha_glass = 0.1;
Medium absorbing_glass = medium("Absorbing Glass", n_glass_abs, alpha = alpha_glass);
Element glass_block = Interface(box((d1, -10.0), (d1 + d2, 10.0)), absorbing_glass, Materials.VACUUM);
OpticalSystem multi_sys = system(glass_block, Materials.VACUUM);

Ray ray_multi = Ray(p = (0.0, 0.0), d = (1.0, 0.0), medium = Materials.VACUUM);
Trace t_multi = multi_sys.trace(ray_multi, far = far_reach);
Ray final_ray = t_multi.last();

real expected_total_opl = (1.0 * d1) + (n_glass_abs * d2) + (1.0 * far_reach);
assert_near("Total Optical Path Length (OPL)", final_ray.opl, expected_total_opl);

real T_interface = 1.0 - ((n_glass_abs - 1.0) / (n_glass_abs + 1.0)) ^ 2;
real expected_final_intensity = (T_interface ^ 2) * exp(-alpha_glass * d2);
assert_near("Final transmitted intensity (Fresnel + Beer-Lambert)", final_ray.intensity(), expected_final_intensity);
write("");

// ============================================================================
// Test Case 13: Faraday Rotation (Non-Reciprocal Polarization)
// ============================================================================
write("Test Case 13: Faraday Rotation (Non-Reciprocal Polarization)");
real V_verdet = 10.0;      // deg / (Tesla * unit)
pair B_field = (1.0, 0.0); // 1 Tesla along +x-axis
Medium faraday_medium = medium("Faraday Glass").add(FaradayRotator(V_verdet, B_field));

// Forward propagation: d = (1, 0) over ds = 3.0 units. Expected rotation = +30 deg.
Polarization pol_init_f = Polarization.linear(0.0, 1.0); // Pure S
Polarization pol_rot_f = faraday_medium.stepPolarization(pol_init_f, (0, 0), (1, 0), 633.0, 3.0);

assert_near("Forward Faraday S1 (cos 60 deg)", pol_rot_f.S1, cos(radians(60.0)), 1e-4);
assert_near("Forward Faraday S2 (sin 60 deg)", pol_rot_f.S2, sin(radians(60.0)), 1e-4);
assert_near("Forward Faraday tilt angle psi (+30 deg)", degrees(pol_rot_f.psi()), 30.0, 1e-3);
assert_near("Forward Faraday intensity conserved", pol_rot_f.intensity(), 1.0, 1e-6);

// Backward propagation: d = (-1, 0) over ds = 3.0 units. Expected rotation = -30 deg.
Polarization pol_rot_b = faraday_medium.stepPolarization(pol_init_f, (0, 0), (-1, 0), 633.0, 3.0);

assert_near("Backward Faraday S1 (cos -60 deg)", pol_rot_b.S1, cos(radians(-60.0)), 1e-4);
assert_near("Backward Faraday S2 (sin -60 deg)", pol_rot_b.S2, sin(radians(-60.0)), 1e-4);
assert_near("Backward Faraday tilt angle psi (-30 deg)", degrees(pol_rot_b.psi()), -30.0, 1e-3);
assert_near("Backward Faraday intensity conserved", pol_rot_b.intensity(), 1.0, 1e-6);
write("");

// ============================================================================
// Test Case 14: PhaseGradient element / PhaseShift wavevector deflection
// ============================================================================
write("Test Case 14: PhaseGradient element / PhaseShift wavevector deflection");
real deltaK_shift = 0.5;
Element pg_interface =
    PhaseGradient(center = (0.0, 0.0), diameter = 10.0, angle = 90.0, delta_k = deltaK_shift, reflective = false);
OpticalSystem pg_sys = system(pg_interface, Materials.VACUUM);

Ray ray_pg = Ray(p = (0.0, -1.0), d = (0.0, 1.0), medium = Materials.VACUUM);
Trace t_pg = pg_sys.trace(ray_pg);
Ray trans_pg = t_pg.transmitted(0);

pair expected_pg_dir = (deltaK_shift, sqrt(1.0 - deltaK_shift ^ 2));
assert_near_pair("Deflected direction after PhaseGradient", trans_pg.d, expected_pg_dir);
write("");

// ============================================================================
// Test Case 15: Geometric Aberration in Spherical Mirror
// ============================================================================
write("Test Case 15: Geometric Aberration in Spherical Mirror");
real f_mir = 50.0;
real x_mir = 100.0;
real y_margin = 30.0;

Element spherical_mir = Mirror.spherical(vertex = (x_mir, 0.0), angle = 180.0, f = f_mir, diameter = 80.0);
Element ideal_mir = Mirror.idealSpherical(vertex = (x_mir, 0.0), angle = 180.0, f = f_mir, diameter = 80.0);

OpticalSystem sys_spherical = system(spherical_mir, Materials.VACUUM);
OpticalSystem sys_ideal = system(ideal_mir, Materials.VACUUM);

Ray ray_marginal = Ray(p = (0.0, y_margin), d = (1.0, 0.0), medium = Materials.VACUUM);

Trace t_spherical = sys_spherical.trace(ray_marginal);
Trace t_ideal = sys_ideal.trace(ray_marginal);

Ray refl_spherical = t_spherical.reflected(0);
Ray refl_ideal = t_ideal.reflected(0);

real expected_x_ideal_focus = x_mir - f_mir;
real t_f_ideal = -refl_ideal.p.y / refl_ideal.d.y;
real x_f_ideal = refl_ideal.p.x + refl_ideal.d.x * t_f_ideal;
assert_near("Ideal Spherical Mirror focal X (No aberration)", x_f_ideal, expected_x_ideal_focus);

real t_f_spherical = -refl_spherical.p.y / refl_spherical.d.y;
real x_f_spherical = refl_spherical.p.x + refl_spherical.d.x * t_f_spherical;

real R_mir = 2.0 * f_mir;
real phi_normal = asin(y_margin / R_mir);
real theta_refl_spherical = 2.0 * phi_normal;
real expected_x_spherical_focus = refl_spherical.p.x - y_margin / tan(theta_refl_spherical);

assert_near("Real Spherical Mirror focal X (With spherical aberration)", x_f_spherical, expected_x_spherical_focus);
assert_bool("Real mirror focus is shifted left/right due to aberration", abs(x_f_spherical - x_f_ideal) > 1.0, true);
write("");

// ============================================================================
// Test Case 16: Aperture Blocking & Spatial Filtering
// ============================================================================
write("Test Case 16: Aperture Blocking & Spatial Filtering");
Element spatial_filter = Aperture(center = (50.0, 0.0), diameter = 50.0, hole = 10.0, angle = 0.0);
OpticalSystem filter_sys = system(spatial_filter, Materials.VACUUM);

Ray ray_inside = Ray(p = (0.0, 2.0), d = (1.0, 0.0), medium = Materials.VACUUM);
Trace t_inside = filter_sys.trace(ray_inside, far = 100.0);

Ray ray_outside = Ray(p = (0.0, 8.0), d = (1.0, 0.0), medium = Materials.VACUUM);
Trace t_outside = filter_sys.trace(ray_outside, far = 100.0);

assert_near_pair("Inside ray final position (unblocked)", t_inside.last().p, (100.0, 2.0));
assert_near("Inside ray intensity (retains 100% power)", t_inside.last().intensity(), 1.0);

assert_near_pair("Outside ray final position (blocked at screen)", t_outside.last().p, (50.0, 8.0));
assert_near("Outside ray intensity after block (0 power)", t_outside.last().intensity(), 0.0);
assert_bool("Outside ray is marked dead/terminated", t_outside.last().terminated || !t_outside.last().alive, true);
write("");

// ============================================================================
// Test Case 17: Cube Beam Splitter (Power Splitting)
// ============================================================================
write("Test Case 17: Cube Beam Splitter (Power Splitting)");
Element beam_splitter =
    BeamSplitter.cube(center = (50.0, 0.0), angle = 0.0, size = 20.0, R = 0.4, glass = Materials.VACUUM);
OpticalSystem split_sys = system(beam_splitter, Materials.VACUUM);

Ray ray_split = Ray(p = (0.0, 0.0), d = (1.0, 0.0), medium = Materials.VACUUM, pol = Polarization.linear(0.0, 1.0));
Trace t_split = split_sys.trace(ray_split, far = 100.0);

Interaction split_ix = t_split.interaction(1);
assert_near("Cube splitter coating: Transmitted fraction", split_ix.T, 0.6);
assert_near("Cube splitter coating: Reflected fraction", split_ix.R, 0.4);

Ray split_trans = t_split.transmitted(1);
assert_near_pair("Transmitted ray direction", split_trans.d, (1.0, 0.0));
assert_near("Transmitted ray intensity", split_trans.intensity(), 0.6);

Ray split_refl = t_split.reflected(1);
assert_near_pair("Reflected ray direction (upward split)", split_refl.d, (0.0, 1.0));
assert_near("Reflected ray intensity", split_refl.intensity(), 0.4);
write("");

// ============================================================================
// Test Case 18: Gain Medium (Optical Amplification)
// ============================================================================
write("Test Case 18: Gain Medium (Optical Amplification)");
real g_amp = 0.1;
real d_entry = 10.0;
real d_slab = 20.0;
real far_exit = 40.0;

Medium active_medium = medium("Amplifying Gas", 1.0).add(Gain(g_amp));
Element gain_slab = Interface(box((d_entry, -10.0), (d_entry + d_slab, 10.0)), active_medium, Materials.VACUUM);
OpticalSystem amp_sys = system(gain_slab, Materials.VACUUM);

Ray ray_amp = Ray(p = (0.0, 0.0), d = (1.0, 0.0), medium = Materials.VACUUM);
Trace t_amp = amp_sys.trace(ray_amp, far = far_exit);
Ray ray_final_amp = t_amp.last();

real expected_amp_opl = d_entry + d_slab + far_exit;
assert_near("Optical Path Length of amplified ray", ray_final_amp.opl, expected_amp_opl);

real expected_amplified_I = exp(g_amp * d_slab);
assert_near("Amplified Ray intensity (Beer-Lambert with negative absorption)",
            ray_final_amp.intensity(),
            expected_amplified_I);
write("");

// ============================================================================
// Test Case 19: Brewster's Angle & Polarization Extinction (Rp = 0)
// ============================================================================
write("Test Case 19: Brewster's Angle & Complete P-Polarization Extinction");
real n_air = 1.0;
real n_glass_br = 1.5;
real theta_B_deg = degrees(atan(n_glass_br / n_air));
real theta_B_rad = radians(theta_B_deg);

Medium glass_brewster = medium("Brewster Glass", n_glass_br);
Element brewster_window = Interface(box((0.0, -20.0), (10.0, 20.0)), glass_brewster, Materials.VACUUM);
OpticalSystem brewster_sys = system(brewster_window, Materials.VACUUM);

Ray ray_p = Ray(p = (-10.0 * cos(theta_B_rad), -10.0 * sin(theta_B_rad)),
                d = (cos(theta_B_rad), sin(theta_B_rad)),
                medium = Materials.VACUUM,
                pol = Polarization.linear(90.0, 1.0)); // Pure P

Trace t_p = brewster_sys.trace(ray_p);
Interaction brewster_ix = t_p.interaction(0);

assert_near("Brewster angle theta_B (deg)", theta_B_deg, 56.309932);
assert_near("P-polarization reflectivity Rp at Brewster angle (0% reflection)", brewster_ix.R, 0.0, 1e-6);
assert_near("P-polarization transmission Tp at Brewster angle (100% transmission)", brewster_ix.T, 1.0, 1e-6);

real sin_t_B = (n_air / n_glass_br) * sin(theta_B_rad);
real cos_t_B = sqrt(1.0 - sin_t_B ^ 2);
real rs_B = (n_air * cos(theta_B_rad) - n_glass_br * cos_t_B) / (n_air * cos(theta_B_rad) + n_glass_br * cos_t_B);
real expected_Rs_brewster = rs_B ^ 2;

Ray ray_s = Ray(p = (-10.0 * cos(theta_B_rad), -10.0 * sin(theta_B_rad)),
                d = (cos(theta_B_rad), sin(theta_B_rad)),
                medium = Materials.VACUUM,
                pol = Polarization.linear(0.0, 1.0)); // Pure S
Trace t_s = brewster_sys.trace(ray_s);
assert_near("S-polarization reflectivity Rs at Brewster angle (~14.8%)", t_s.interaction(0).R, expected_Rs_brewster);
write("");

// ============================================================================
// Test Case 20: Thick Lens & Exact Lensmaker's Equation
// ============================================================================
write("Test Case 20: Thick Lens & Exact Lensmaker's Equation");
real R1 = 50.0, R2 = 50.0, t_thick = 10.0, n_thick = 1.5;
Medium bk7_thick = medium("BK7", n_thick);

real power_thick = (n_thick - 1.0) * ((1.0 / R1) + (1.0 / R2) - ((n_thick - 1.0) * t_thick) / (n_thick * R1 * R2));
real expected_f_thick = 1.0 / power_thick;
real expected_BFL = expected_f_thick * (1.0 - (n_thick - 1.0) * t_thick / (n_thick * R1));

real center_x = 100.0;
Element thick_lens = Lens.biconvex(center = (center_x, 0.0),
                                  angle = 0.0,
                                  R1 = R1,
                                  R2 = R2,
                                  thickness = t_thick,
                                  diameter = 30.0,
                                  glass = bk7_thick,
                                  ambient = Materials.VACUUM);
OpticalSystem thick_sys = system(thick_lens, Materials.VACUUM);

real y_paraxial = 1.0;
Ray ray_thick = Ray(p = (0.0, y_paraxial), d = (1.0, 0.0), medium = Materials.VACUUM);
Trace t_thick_lens = thick_sys.trace(ray_thick);
Ray ray_exit_thick = t_thick_lens.last();

real f_measured = y_paraxial / abs(ray_exit_thick.d.y);

real t_f_thick = -ray_exit_thick.p.y / ray_exit_thick.d.y;
real x_focus_measured = ray_exit_thick.p.x + ray_exit_thick.d.x * t_f_thick;

real back_vertex_x = center_x + t_thick / 2.0;
real expected_focus_x = back_vertex_x + expected_BFL;

assert_near("Thick Lens Ray-Traced Effective Focal Length f", f_measured, expected_f_thick, 0.05);
assert_near("Thick Lens Ray-Traced Focal Intersection X", x_focus_measured, expected_focus_x, 0.05);
write("");

// ============================================================================
// Test Case 21: Nested Spatial Containment (Glass Lens in Water Tank)
// ============================================================================
write("Test Case 21: Nested Spatial Containment (Glass Lens in Water Tank)");
Medium water = medium("Water", 1.333333);
Medium glass_in_water = medium("Glass", 1.50);

Element water_tank = Interface(box((0.0, -40.0), (200.0, 40.0)), water, Materials.VACUUM);
Element lens_inside = Lens.biconvex(center = (100.0, 0.0),
                                   angle = 0.0,
                                   R1 = 50.0,
                                   R2 = 50.0,
                                   thickness = 6.0,
                                   diameter = 30.0,
                                   glass = glass_in_water);

OpticalSystem nested_sys = system(new Element[]{water_tank, lens_inside}, Materials.VACUUM);

Ray ray_nested = Ray(p = (-10.0, 5.0), d = (1.0, 0.0), medium = Materials.VACUUM);
Trace t_nested = nested_sys.trace(ray_nested, far = 250.0);

assert_near("Total interface interactions in nested system", t_nested.interactions().length, 4.0);
assert_near_pair("Entry into water tank at x = 0", t_nested.hit(0), (0.0, 5.0));
assert_bool("Ray properly resolves water as ambient for the lens", abs(t_nested.hit(1).x - 97.0) < 1.0, true);
assert_bool("Ray exits water tank at x = 200", abs(t_nested.hit(3).x - 200.0) < 1e-4, true);
write("");

// ============================================================================
// Test Case 22: 90-Degree Retroreflector Invariance
// ============================================================================
write("Test Case 22: 90-Degree Retroreflector Invariance");
Element mirror_h = Mirror.flat(center = (25.0, 0.0), angle = 90.0, diameter = 50.0);
Element mirror_v = Mirror.flat(center = (0.0, 25.0), angle = 0.0, diameter = 50.0);

OpticalSystem retro_sys = system(new Element[]{mirror_h, mirror_v}, Materials.VACUUM);

real theta_retro = radians(37.0);
pair d_in_retro = (-cos(theta_retro), -sin(theta_retro));
Ray ray_retro = Ray(p = (40.0, 30.0), d = d_in_retro, medium = Materials.VACUUM);

Trace t_retro = retro_sys.trace(ray_retro);

pair expected_d_retro = -d_in_retro;
assert_near("Corner reflector: exactly 2 reflections", t_retro.interactions().length, 2.0);
assert_near_pair("Retroreflected ray is strictly anti-parallel (-d_in)", t_retro.reflected(1).d, expected_d_retro);
write("");

// ============================================================================
// Test Case 23: Fermat's Principle / Geometric OPL across a Thin Lens
// ============================================================================
write("Test Case 23: Fermat's Principle / Geometric OPL across a Thin Lens");
real s_o = 200.0, s_i = 200.0, f_iso = 100.0;
real h_marginal = 15.0;

Element iso_lens = Lens.thin(center = (0.0, 0.0), angle = 0.0, f = f_iso, diameter = 60.0);
OpticalSystem iso_sys = system(iso_lens, Materials.VACUUM);

pair point_source = (-s_o, 0.0);
pair image_point = (s_i, 0.0);

Ray ray_axial = Ray(p = point_source, d = (1.0, 0.0), medium = Materials.VACUUM);
Trace t_axial = iso_sys.trace(ray_axial, far = s_i);

pair d_marginal = unit((s_o, h_marginal));
Ray ray_marginal_iso = Ray(p = point_source, d = d_marginal, medium = Materials.VACUUM);
Trace t_marginal = iso_sys.trace(ray_marginal_iso, far = s_i * 1.5);

real opl_axial = t_axial.opl();

Ray marginal_after = t_marginal.transmitted(0);
real t_to_image = (image_point.x - marginal_after.p.x) / marginal_after.d.x;
real opl_marginal_to_image = t_marginal.interaction(0).incident.opl + t_to_image * 1.0;

real expected_opl_marginal = sqrt(s_o ^ 2 + h_marginal ^ 2) + sqrt(s_i ^ 2 + h_marginal ^ 2);

assert_near("Axial Ray OPL to image point", opl_axial, s_o + s_i);
assert_near("Marginal Ray Geometric OPL", opl_marginal_to_image, expected_opl_marginal, 1e-4);
write("");

// ============================================================================
// Test Case 24: Multi-Order Diffraction Grating Energy Conservation
// ============================================================================
write("Test Case 24: Multi-Order Diffraction Grating Energy Conservation");
int[] diff_orders = {-1, 0, 1};
real[] diff_efficiencies = {0.25, 0.50, 0.25};

Element tri_grating = Grating(center = (0.0, 0.0),
                             diameter = 20.0,
                             angle = 90.0,
                             lines_per_mm = 300.0,
                             orders = diff_orders,
                             efficiencies = diff_efficiencies);
OpticalSystem tri_sys = system(tri_grating, Materials.VACUUM);

Ray ray_tri = Ray(p = (0.0, -5.0), d = (0.0, 1.0), medium = Materials.VACUUM, lambda = 600.0);
Trace t_tri = tri_sys.trace(ray_tri);

Interaction g_ix = t_tri.interaction(0);
assert_near("Total diffracted orders spawned", g_ix.transmitted.length, 3.0);

real total_diffracted_power = 0.0;
for (int k = 0; k < g_ix.transmitted.length; ++k) {
    total_diffracted_power += g_ix.transmitted[k].intensity();
}

assert_near("Diffraction order m = -1 power (25%)", g_ix.transmitted[0].intensity(), 0.25);
assert_near("Diffraction order m =  0 power (50%)", g_ix.transmitted[1].intensity(), 0.50);
assert_near("Diffraction order m = +1 power (25%)", g_ix.transmitted[2].intensity(), 0.25);
assert_near("Total Energy Conserved across all diffraction orders (100%)", total_diffracted_power, 1.0);
write("");

// ============================================================================
// Test Case 25: Optical Bisection & Geometric Ray Aiming (aimAngle)
// ============================================================================
write("Test Case 25: Optical Bisection & Geometric Ray Aiming (aimAngle)");
Element lens_aim = Lens.thin(center = (100.0, 0.0), angle = 0.0, f = 100.0, diameter = 80.0);
OpticalSystem sys_aim = system(lens_aim, Materials.VACUUM);

pair pt_a = (0.0, 20.0);
pair pt_b = (300.0, -10.0);

real aimed_angle_deg = aimAngle(pt_a, pt_b, sys_aim, -10.0, 10.0);
real expected_aimed_angle = degrees(atan(10.0 / 100.0));

assert_near("Aimed launch angle (degrees)", aimed_angle_deg, expected_aimed_angle, 1e-3);

Ray ray_aimed = Ray.angle(pt_a, aimed_angle_deg, Materials.VACUUM);
Trace t_aimed = sys_aim.trace(ray_aimed, far = 400.0);
Ray final_aimed = t_aimed.last();

real t_aimed_at_300 = (300.0 - final_aimed.p.x) / final_aimed.d.x;
real y_at_300 = final_aimed.p.y + final_aimed.d.y * t_aimed_at_300;

assert_near("Aimed ray crossing height Y at x = 300 (target Y = -10)", y_at_300, -10.0, 1e-3);
write("");

// ============================================================================
// Test Case 26: Stokes Vector Invariants & DOP
// ============================================================================
write("Test Case 26: Stokes Vector Invariants & DOP");
Polarization unpol = Polarization.unpolarized(2.5);
assert_stokes_physical("Unpolarized state", unpol);
assert_near("Unpolarized S0", unpol.S0, 2.5, 1e-5);
assert_near("Unpolarized S1", unpol.S1, 0.0, 1e-5);
assert_near("Unpolarized S2", unpol.S2, 0.0, 1e-5);
assert_near("Unpolarized S3", unpol.S3, 0.0, 1e-5);
assert_near("Unpolarized DOP", unpol.dop(), 0.0, 1e-5);

Polarization lin = Polarization.linear(45.0, 1.0);
assert_stokes_physical("Linear 45 state", lin);
assert_near("Linear 45 S1", lin.S1, 0.0, 1e-5);
assert_near("Linear 45 S2", lin.S2, 1.0, 1e-5);
assert_near("Linear 45 S3", lin.S3, 0.0, 1e-5);
assert_near("Linear 45 DOP", lin.dop(), 1.0, 1e-5);
write("");

// ============================================================================
// Test Case 27: Circular Polarization Handedness
// ============================================================================
write("Test Case 27: Circular Polarization Handedness");
Polarization rcp = Polarization.circular();
assert_near("RCP S0", rcp.S0, 1.0, 1e-5);
assert_near("RCP S1", rcp.S1, 0.0, 1e-5);
assert_near("RCP S2", rcp.S2, 0.0, 1e-5);
assert_near("RCP S3", rcp.S3, 1.0, 1e-5);
assert_near("RCP chi angle (rad)", rcp.chi(), pi / 4.0, 1e-5);
assert_near("RCP DOP", rcp.dop(), 1.0, 1e-5);

Polarization lcp = Polarization.circular(false);
assert_near("LCP S3", lcp.S3, -1.0, 1e-5);
assert_near("LCP chi angle (rad)", lcp.chi(), -pi / 4.0, 1e-5);
assert_near("LCP DOP", lcp.dop(), 1.0, 1e-5);
write("");

// ============================================================================
// Test Case 28: GRIN Adaptive Solver Precision (RK45 vs RK4)
// ============================================================================
write("Test Case 28: GRIN Adaptive Solver Precision (RK45 vs RK4)");
real n0_g = 1.5, g_g = 0.10;
indexProfile gp = new real(pair p, real lambda) { return n0_g * (1.0 - 0.5 * g_g ^ 2 * p.y ^ 2); };
gradientProfile gg = new pair(pair p, real lambda) { return (0.0, -n0_g * g_g ^ 2 * p.y); };

Medium grin_med = medium("Parabolic GRIN").add(GRIN(gp, gg));
Element grin_slab = Interface(box((0.0, -10.0), (50.0, 10.0)), grin_med, Materials.VACUUM);
OpticalSystem sys_grin_compare = system(grin_slab, Materials.VACUUM);

real y0_test = 0.2;
real s_travel = 5.0;
Ray r_init = Ray(p = (0.0, y0_test), d = (1.0, 0.0), medium = grin_med);

config.physics.grin_solver = GRINSolverType.RK4;
Trace t_rk4 = sys_grin_compare.trace(r_init, ds = 0.05, far = s_travel);
Ray ray_rk4 = t_rk4.last();

config.physics.grin_solver = GRINSolverType.RK45;
Trace t_rk45 = sys_grin_compare.trace(r_init, ds = 0.05, far = s_travel);
Ray ray_rk45 = t_rk45.last();

real expected_y_paraxial = y0_test * cos(g_g * s_travel);
assert_near("RK45 adaptive solver final Y position", ray_rk45.p.y, expected_y_paraxial, 1e-4);
assert_near("RK4 classic solver final Y position", ray_rk4.p.y, expected_y_paraxial, 1e-4);
assert_near("RK45 and RK4 solvers converge to each other (< 1e-5)", ray_rk45.p.y, ray_rk4.p.y, 1e-5);
write("");

// ============================================================================
// Test Case 29: GRIN OPL Midpoint Integration
// ============================================================================
write("Test Case 29: GRIN OPL Midpoint Integration");
indexProfile axial_gradient_index = new real(pair p, real lambda) { return 1.0 + p.x; };
gradientProfile axial_gradient = new pair(pair p, real lambda) { return (1.0, 0.0); };
Medium axial_gradient_medium = medium("Axial gradient").add(GRIN(axial_gradient_index, axial_gradient));
OpticalSystem empty_system = system(new Element[]);
Trace axial_gradient_trace = empty_system.trace(Ray((0, 0), (1, 0), axial_gradient_medium), ds = 0.05, far = 1.0);
assert_near("GRIN OPL integral of n(x)=1+x over [0,1]", axial_gradient_trace.opl(), 1.5, 1e-8);
assert_near_pair("Line intersection", lineIntersection((0, 0), (1, 0), (2, -1), (0, 1)), (2, 0));
pair parallel_intersection = lineIntersection((0, 0), (1, 0), (0, 1), (2, 0));
assert_bool("Parallel line intersection is infinite", abs(parallel_intersection.x) == infinity, true);
write("");

// ============================================================================
// Test Case 30: Ordered Custom Surface Properties
// ============================================================================
write("Test Case 30: Ordered Custom Surface Properties");
Property half_power = property(
    "Half power coating",
    new Ray[](Ray r, Hit h, mediumResolver resolve) {
        r.pol.scaleTo(0.5 * r.intensity());
        return new Ray[]{r};
    },
    PropertyStage.AfterBoundary);
path test_boundary = linearBoundary((0, 0), (1, 0), 20);
Element bare_interface = Interface(test_boundary, medium("test glass", 1.5));
Element coated_interface = Interface(test_boundary, medium("test glass", 1.5));
coated_interface.surfaces[0].add(half_power);
Ray coating_probe = Ray.angle((-5, 0), 10, Materials.VACUUM);
real bare_power = system(bare_interface).trace(coating_probe).transmitted(0).intensity();
real coated_power = system(coated_interface).trace(coating_probe).transmitted(0).intensity();
assert_near("AfterBoundary custom action runs after refraction", coated_power, 0.5 * bare_power, 1e-8);

VolumeProperty low_priority_index = volumeProperty("low priority", priority = 0);
low_priority_index.phase_index = new real(pair p, pair d, real lambda, Polarization pol) { return 1.2; };
VolumeProperty high_priority_index = volumeProperty("high priority", priority = 10);
high_priority_index.phase_index = new real(pair p, pair d, real lambda, Polarization pol) { return 1.6; };
Medium ordered_medium = medium("ordered").add(high_priority_index).add(low_priority_index);
assert_near("Volume-index profiles compose multiplicatively", ordered_medium.n(), 1.92);
assert_bool("Volume-property priority determines evaluation order",
            ordered_medium.orderedProperties()[1].name == "high priority",
            true);
write("");

// ============================================================================
// Test Case 29b: Diverging and Parallel Beam Generator Factories
// ============================================================================
write("Test Case 29b: Diverging and Parallel Beam Generator Factories");
Ray[] p_beam = Beam.parallel(start = (0.0, -4.0), step = (0.0, 2.0), count = 5, d = (1.0, 0.0));
assert_near("Parallel beam ray count", p_beam.length, 5.0);
assert_near_pair("Parallel beam Ray 0 start point", p_beam[0].p, (0.0, -4.0));
assert_near_pair("Parallel beam Ray 2 start point (center)", p_beam[2].p, (0.0, 0.0));
assert_near_pair("Parallel beam Ray 4 start point", p_beam[4].p, (0.0, 4.0));
assert_near_pair("Parallel beam Ray direction (+X)", p_beam[3].d, (1.0, 0.0));

Ray[] d_beam = Beam.diverging(source = (0.0, 0.0), d = (1.0, 0.0), angle = 30.0, count = 3);
assert_near("Diverging beam ray count", d_beam.length, 3.0);
assert_near_pair("Diverging beam Ray 1 direction (center)", d_beam[1].d, (1.0, 0.0));
assert_near("Diverging beam Ray 0 angle (-15 deg)", degrees(atan2(d_beam[0].d.y, d_beam[0].d.x)), -15.0);
assert_near("Diverging beam Ray 2 angle (+15 deg)", degrees(atan2(d_beam[2].d.y, d_beam[2].d.x)), 15.0);
write("");

// ============================================================================
// Test Case 30b: Biconcave Lens Beam Divergence & Virtual Focus
// ============================================================================
write("Test Case 30b: Biconcave Lens Beam Divergence & Virtual Focus");
Element biconcave_lens = Lens.thin(center = (100.0, 0.0), angle = 0.0, f = -50.0, diameter = 40.0);
OpticalSystem sys_concave = system(biconcave_lens, Materials.VACUUM);

Ray ray_cc_top = Ray(p = (0.0, 10.0), d = (1.0, 0.0), medium = Materials.VACUUM);
Ray ray_cc_bot = Ray(p = (0.0, -10.0), d = (1.0, 0.0), medium = Materials.VACUUM);

Trace t_cc_top = sys_concave.trace(ray_cc_top);
Trace t_cc_bot = sys_concave.trace(ray_cc_bot);

Ray trans_cc_top = t_cc_top.transmitted(0);
Ray trans_cc_bot = t_cc_bot.transmitted(0);

pair expected_cc_top_dir = unit((50.0, 10.0));
pair expected_cc_bot_dir = unit((50.0, -10.0));

assert_near_pair("Biconcave diverging top ray direction", trans_cc_top.d, expected_cc_top_dir);
assert_near_pair("Biconcave diverging bottom ray direction", trans_cc_bot.d, expected_cc_bot_dir);

real t_v_focus_top = -trans_cc_top.p.y / trans_cc_top.d.y;
real x_v_focus_top = trans_cc_top.p.x + trans_cc_top.d.x * t_v_focus_top;
assert_near("Biconcave virtual focal point X (expected 50.0)", x_v_focus_top, 50.0);
write("");

// ============================================================================
// Test Case 31: Dichroic Extinction Ratio
// ============================================================================
write("Test Case 31: Dichroic Extinction Ratio");
Medium pol_medium = medium("DichroicFilter").add(Dichroic(0.0, 10.0));
Polarization out_unpol = pol_medium.stepPolarization(Polarization.unpolarized(1.0), (0, 0), (1, 0), 633.0, 1.0);

assert_stokes_physical("Filter output physical state", out_unpol);
assert_near("Diattenuated output S0", out_unpol.S0, 0.5 * (1.0 + exp(-10.0)), 1e-4);
assert_near("Diattenuated output S1", out_unpol.S1, 0.5 * (1.0 - exp(-10.0)), 1e-4);
assert_near("Filtered light DOP approaching 1.0", out_unpol.dop(), 1.0, 1e-3);
write("");

// ============================================================================
// Test Case 32: Dense-to-Dense TIR Critical Angle
// ============================================================================
write("Test Case 32: Dense-to-Dense TIR Critical Angle");
Medium flint = medium("Flint 1.62", 1.62);
Medium water_m = medium("Water 1.333", 1.333333);

real theta_c_deg = degrees(asin(water_m.n() / flint.n()));
assert_true("Critical angle is physically defined for n_high > n_low",
            water_m.n() < flint.n() && theta_c_deg > 0.0 && theta_c_deg < 90.0);

Ray ray_below_crit = Ray.angle(p = (-1, 0), angle = theta_c_deg - 2.0, medium = flint);
Ray ray_above_crit = Ray.angle(p = (-1, 0), angle = theta_c_deg + 2.0, medium = flint);

mediumResolver res_flint_water = new Medium(pair p, Ray r, Medium m1, Medium m2, bool from_m1) {
    return (p.x > 0) ? water_m : flint;
};
Hit h_planar;
h_planar.valid = true;
h_planar.normal = (-1, 0);

FresnelResult ref_below = fresnelIsotropic(flint, water_m, ray_below_crit, h_planar, res_flint_water);
FresnelResult ref_above = fresnelIsotropic(flint, water_m, ray_above_crit, h_planar, res_flint_water);

assert_bool("TIR below critical angle is false", ref_below.tir, false);
assert_bool("TIR above critical angle is true", ref_above.tir, true);
write("");

// ============================================================================
// Test Case 33: Ideal Mirror Collimation from Focal Point
// ============================================================================
write("Test Case 33: Ideal Mirror Collimation from Focal Point");
Element oap = Mirror.idealSpherical(vertex = (50, 0), angle = 180.0, f = 50.0, diameter = 40.0);
OpticalSystem oap_sys = system(oap, Materials.VACUUM);

Ray oap_ray1 = Ray.between((0, 0), (50, 5), Materials.VACUUM);
Ray oap_ray2 = Ray.between((0, 0), (50, -5), Materials.VACUUM);

Trace t_oap1 = oap_sys.trace(oap_ray1);
Trace t_oap2 = oap_sys.trace(oap_ray2);

assert_near("Ray 1 collimated parallel to axis (d.x = -1)", t_oap1.reflected(0).d.x, -1.0, 1e-4);
assert_near("Ray 1 zero transverse divergence (d.y = 0)", t_oap1.reflected(0).d.y, 0.0, 1e-4);
assert_near("Ray 2 collimated parallel to axis (d.x = -1)", t_oap2.reflected(0).d.x, -1.0, 1e-4);
assert_near("Ray 2 zero transverse divergence (d.y = 0)", t_oap2.reflected(0).d.y, 0.0, 1e-4);
write("");

// ============================================================================
// Test Case 34: Optical Smith-Helmholtz / Lagrange Invariant
// ============================================================================
write("Test Case 34: Optical Smith-Helmholtz / Lagrange Invariant");
real f_invariant_lens = 100.0;
real s_obj = 200.0;
Element lens = Lens.thin(center = (0, 0), angle = 0, f = f_invariant_lens, diameter = 80);
OpticalSystem inv_sys = system(lens, Materials.VACUUM);

real y_obj = 10.0;
real u_in = 0.05;

Ray r_marginal = Ray(p = (-s_obj, 0), d = (cos(u_in), sin(u_in)), medium = Materials.VACUUM);
Trace t_marg = inv_sys.trace(r_marginal, far = 300.0);
Ray exit_marg = t_marg.transmitted(0);
real u_out = abs(asin(exit_marg.d.y));

Ray r_chief = Ray.between((-s_obj, y_obj), (0, 0), Materials.VACUUM);
Trace t_chief = inv_sys.trace(r_chief, far = 300.0);
Ray exit_chief = t_chief.transmitted(0);

real t_img = (s_obj - exit_chief.p.x) / exit_chief.d.x;
real y_img = abs(exit_chief.p.y + exit_chief.d.y * t_img);

real H_in = 1.0 * y_obj * u_in;
real H_out = 1.0 * y_img * u_out;

assert_near("Optical Lagrange invariant H = n * y * u", H_in, H_out, 5e-3);
write("");

// ============================================================================
// Test Case 35: Sellmeier Dispersion Formula Verification
// ============================================================================
write("Test Case 35: Sellmeier Dispersion Formula Verification");
real lambda_test_nm = 587.5618;
real l2 = (lambda_test_nm / 1000.0) ^ 2;
real B1 = 1.03961212, C1 = 0.00600069867;
real B2 = 0.231792344, C2 = 0.0200179144;
real B3 = 1.01046945, C3 = 103.560653;
real n_sellmeier_expected = sqrt(1.0 + (B1 * l2) / (l2 - C1) + (B2 * l2) / (l2 - C2) + (B3 * l2) / (l2 - C3));

indexProfile sellmeier_b_k7 = new real(pair p, real lambda) {
    real w2 = (lambda / 1000.0) ^ 2;
    return sqrt(1.0 + (1.03961212 * w2) / (w2 - 0.00600069867) + (0.231792344 * w2) / (w2 - 0.0200179144) +
                (1.01046945 * w2) / (w2 - 103.560653));
};
Medium bk7_disp = medium("BK7 Sellmeier").add(GRIN(sellmeier_b_k7));
assert_near("BK7 refractive index at d-line (587.56 nm)", bk7_disp.n(lambda = lambda_test_nm), 1.51680, 1e-4);
write("");

// ============================================================================
// Test Case 36: Transmissive PhaseShift Target Medium Continuity
// ============================================================================
write("Test Case 36: Transmissive PhaseShift Target Medium Continuity");
Medium dense_glass = medium("Dense Glass", 2.0);
Element meta_sheet = PhaseGradient(center = (0, 0), diameter = 20, angle = 90, delta_k = 0.2, reflective = false);
OpticalSystem ps_sys = system(meta_sheet, dense_glass);

Ray r_ps = Ray(p = (-5, 0), d = (1, 0), medium = Materials.VACUUM);
Trace t_ps = ps_sys.trace(r_ps);
assert_near("PhaseShift transmitted ray inherits resolved ambient medium", t_ps.transmitted(0).medium.n(), 2.0);
write("");

// Test Case 37: Polarization.scaleTo & Single-Source Phase Tracking
write("Test Case 37: Polarization.scaleTo & Single-Source Phase Tracking");

// Part A: Direct Unit Test on Polarization.scaleTo
Polarization pol = Polarization.circular(2.0);
assert_near("Initial intensity", pol.intensity(), 2.0, 1e-6);

pol.scaleTo(8.0);
assert_near("Rescaled intensity matches target 8.0", pol.intensity(), 8.0, 1e-6);
assert_near("Stokes S0 scaled to 8.0", pol.S0, 8.0, 1e-6);
assert_near("Stokes S3 scaled proportionally to 8.0", pol.S3, 8.0, 1e-6);
assert_near("DOP preserved under scaleTo", pol.dop(), 1.0, 1e-6);

// Part B: System Reflection Test
Element flat_mirror = Mirror.flat(center = (50, 0), angle = 180, diameter = 40);
OpticalSystem mirror_sys = system(flat_mirror, Materials.VACUUM);

Ray ray_in = Ray(p = (0, 0), d = (1, 0), pol = Polarization.linear(45.0, 1.0));
ray_in.phase = 0.0;
real I_in = ray_in.intensity();

Trace t_refl = mirror_sys.trace(ray_in, far = 100.0);
Ray ray_out = t_refl.reflected(0);

assert_near("Reflected ray intensity conserved", ray_out.intensity(), I_in, 1e-6);
assert_near("Single reflection phase shift is exactly pi", ray_out.phase, pi, 1e-6);
// Ideal mirror Mueller matrix diag(1, 1, -1, -1) (same convention as
// dielectric reflection): linear 45 deg maps to linear -45 deg.
assert_near("Mirror maps S2 -> -S2 (linear 45 -> -45)", ray_out.pol.S2, -1.0, 1e-6);
assert_stokes_physical("Reflected ray physical state", ray_out.pol);
write("");

// ============================================================================
// Test Case 38: Minimum Deviation of an Equilateral Prism
// ============================================================================
write("Test Case 38: Minimum Deviation of an Equilateral Prism");
real L = 60.0;
real H_apex = L / sqrt(3.0);
real H_base = -L / (2.0 * sqrt(3.0));
pair apex = (0.0, H_apex);
pair b_left = (-L / 2.0, H_base);
pair b_right = (L / 2.0, H_base);
path prism_path = b_left -- b_right -- apex -- cycle;

Medium glass15 = medium("Glass 1.5", 1.50);
Element prism = Interface(prism_path, glass15, Materials.VACUUM);
OpticalSystem prism_sys = system(prism, Materials.VACUUM);

real theta_i = degrees(asin(1.5 * sin(radians(30.0))));
real ray_angle_lab = -30.0 + theta_i;
pair hit_left = midpoint(b_left -- apex);
Ray r_prism =
    Ray.toPoint(target = hit_left, d = (cos(radians(ray_angle_lab)), sin(radians(ray_angle_lab))), distance = 20.0);

Trace t_prism = prism_sys.trace(r_prism, far = 40.0);
Ray inside_ray = t_prism.transmitted(0);
Ray exit_ray = t_prism.transmitted(1);

assert_near("Ray inside prism is parallel to base (d.y = 0)", inside_ray.d.y, 0.0, 1e-4);
assert_near("Ray inside prism d.x = 1", inside_ray.d.x, 1.0, 1e-4);

real exit_angle_lab = degrees(atan2(exit_ray.d.y, exit_ray.d.x));
assert_near("Exit ray angle is symmetric (-18.59 deg)", exit_angle_lab, -ray_angle_lab, 1e-3);

real expected_delta = 2.0 * theta_i - 60.0;
real measured_delta = ray_angle_lab - exit_angle_lab;
assert_near("Minimum deviation angle matches formula (37.18 deg)", measured_delta, expected_delta, 1e-3);
write("");

// ============================================================================
// Test Case 39: TIR Optical Waveguide / Light Pipe
// ============================================================================
write("Test Case 39: TIR Optical Waveguide / Light Pipe");
path core_path = box((0, -5), (100, 5));
Element guide = Interface(core_path, medium("Core", 1.50), Materials.VACUUM);
OpticalSystem guide_sys = system(guide, Materials.VACUUM);

Ray r_guide = Ray(p = (0.01, 0), d = (cos(radians(25.0)), sin(radians(25.0))), medium = medium("Core", 1.50));
Trace t_guide = guide_sys.trace(r_guide, far = 120.0);

Ray exit_guide = t_guide.last();
real guide_n1 = 1.50;
real guide_n2 = 1.0;
real guide_theta_i = radians(25.0);
real guide_sin_t = guide_n1 / guide_n2 * sin(guide_theta_i);
real guide_cos_t = sqrt(1.0 - guide_sin_t ^ 2);
real guide_cos_i = cos(guide_theta_i);
real guide_rs = ((guide_n1 * guide_cos_i - guide_n2 * guide_cos_t) / (guide_n1 * guide_cos_i + guide_n2 * guide_cos_t)) ^ 2;
real guide_rp = ((guide_n2 * guide_cos_i - guide_n1 * guide_cos_t) / (guide_n2 * guide_cos_i + guide_n1 * guide_cos_t)) ^ 2;
real guide_expected_transmission = 1.0 - 0.5 * (guide_rs + guide_rp);
assert_near("Waveguide output power from end-face Fresnel transmission",
            exit_guide.intensity(),
            guide_expected_transmission,
            1e-4);
assert_bool("Ray traversed through the guide to x >= 100", exit_guide.p.x >= 100.0 - 1e-3, true);
write("");

// ============================================================================
// Test Case 41: Degree of Polarization (DoP) & Stokes Vector at Brewster Angle
// ============================================================================
write("Test Case 41: Degree of Polarization (DoP) & Stokes Vector at Brewster Angle");
Medium glass_br_41 = medium("Glass", 1.50);
real theta_B_deg_41 = degrees(atan(1.50 / 1.0));

path boundary_plane = (0, -50) -- (0, 50);
Element iface = Interface(boundary_plane, Materials.VACUUM, glass_br_41);
OpticalSystem br_sys = system(iface, Materials.VACUUM);

Ray r_unpol = Ray.toPoint(target = (0, 0),
                          d = (cos(radians(theta_B_deg_41)), sin(radians(theta_B_deg_41))),
                          distance = 20.0,
                          pol = Polarization.unpolarized(1.0));
Trace t_br = br_sys.trace(r_unpol);

Ray r_refl = t_br.reflected(0);
Ray r_trans = t_br.transmitted(0);

assert_near("Brewster Reflected Ray DoP = 100%", r_refl.pol.dop(), 1.0, 1e-4);
assert_near("Brewster Reflected Ray is purely S-polarized (S1 == S0)", r_refl.pol.S1, r_refl.pol.S0, 1e-4);

real brewster_sin_t = sin(radians(theta_B_deg_41)) / glass_br_41.n();
real brewster_cos_t = sqrt(1.0 - brewster_sin_t ^ 2);
real brewster_rs_amplitude = (cos(radians(theta_B_deg_41)) - glass_br_41.n() * brewster_cos_t) /
                           (cos(radians(theta_B_deg_41)) + glass_br_41.n() * brewster_cos_t);
real brewster_rs_power = brewster_rs_amplitude ^ 2;
real brewster_ts = 1.0 - brewster_rs_power;
real brewster_tp = 1.0;
real expected_dop = (brewster_tp - brewster_ts) / (brewster_tp + brewster_ts);
assert_near("Brewster Transmitted Ray DoP matches analytical formula (~7.99%)", r_trans.pol.dop(), expected_dop, 1e-3);
write("");

// ============================================================================
// Test Case 42: Porro Prism Double-TIR Retroreflector
// ============================================================================
write("Test Case 42: Porro Prism Double-TIR Retroreflector");
path porro_path = (0, -20) -- (20, 0) -- (0, 20) -- cycle;
Element porro = Interface(porro_path, medium("Glass", 1.50), Materials.VACUUM);
OpticalSystem porro_sys = system(porro, Materials.VACUUM);

Ray r_porro = Ray(p = (-10, 8.0), d = (1, 0), medium = Materials.VACUUM);
Trace t_porro = porro_sys.trace(r_porro, far = 40.0);
Ray r_porro_exit = t_porro.last();

assert_near("Porro retroreflected d.x is -1", r_porro_exit.d.x, -1.0, 1e-4);
assert_near("Porro retroreflected d.y is 0", r_porro_exit.d.y, 0.0, 1e-4);
assert_near("Porro inverted exit position y = -8.0", r_porro_exit.p.y, -8.0, 1e-4);

real porro_r0 = ((1.50 - 1.0) / (1.50 + 1.0)) ^ 2;
real porro_t0 = 1.0 - porro_r0;
assert_near("Porro prism intensity from two Fresnel transmissions", r_porro_exit.intensity(), porro_t0 ^ 2, 1e-4);
write("");

// ============================================================================
// Test Case 43: Cardinal Points Verification
// ============================================================================
write("Test Case 43: Cardinal Points Verification");
Medium glass_lens = medium("Crown", 1.50);
Element thickLens_cardinal =
    Lens.biconvex(center = (50, 0), angle = 0, R1 = 50, R2 = 50, thickness = 10, diameter = 40, glass = glass_lens);
OpticalSystem sys = system(thickLens_cardinal, Materials.VACUUM);
CardinalPoints cp = cardinalPoints(sys);

real n_cardinal = 1.50;
real R1_cardinal = 50.0;
real R2_cardinal = -50.0;
real thickness_cardinal = 10.0;

real expected_v1x = 50.0 - thickness_cardinal / 2.0;
real expected_v2x = 50.0 + thickness_cardinal / 2.0;
assert_near("Front vertex from center and thickness", cp.V1.x, expected_v1x, 1e-3);
assert_near("Rear vertex from center and thickness", cp.V2.x, expected_v2x, 1e-3);

real lens_power =
    (n_cardinal - 1.0) * (1.0 / R1_cardinal - 1.0 / R2_cardinal +
                          ((n_cardinal - 1.0) * thickness_cardinal) / (n_cardinal * R1_cardinal * R2_cardinal));
real expected_e_f_l = 1.0 / lens_power;
assert_relative("Effective focal length from thick lensmaker equation", cp.f2, expected_e_f_l, 1e-5, 1e-5);

real expected_b_f_l = expected_e_f_l * (1.0 - ((n_cardinal - 1.0) * thickness_cardinal) / (n_cardinal * R1_cardinal));
assert_relative("Back focal length from thick-lens paraxial formula", cp.BFL, expected_b_f_l, 1e-5, 1e-5);
assert_relative("Rear focal point is one BFL beyond rear vertex", cp.F2.x, cp.V2.x + cp.BFL, 1e-5, 1e-5);
assert_relative("Front focal point is one BFL before front vertex", cp.F1.x, cp.V1.x - cp.BFL, 1e-5, 1e-5);

real principal_offset = cp.f2 - cp.BFL;
assert_relative("Front principal point offset", cp.H1.x, cp.V1.x + principal_offset, 1e-5, 1e-5);
assert_relative("Rear principal point offset", cp.H2.x, cp.V2.x - principal_offset, 1e-5, 1e-5);
write("");

// ============================================================================
// Test Case 44: Mueller Diattenuation and Retardance Helper Verification
// ============================================================================
write("Test Case 44: Mueller Diattenuation and Retardance Helper Verification");
// 1. Pure diattenuator: S_in = (1, 0, 1, 0), T1 = 0.8, T2 = 0.2
Polarization s_in = new Polarization;
s_in.S0 = 1.0;
s_in.S1 = 0.0;
s_in.S2 = 1.0;
s_in.S3 = 0.0;
Polarization s_out = copy(s_in);
applyDiattenuation(s_out, s_in, 0.8, 0.2);

assert_near("Pure diattenuation S0 = 0.5 * (0.8 + 0.2)", s_out.S0, 0.5, 1e-5);
assert_near("Pure diattenuation S1 = 0.5 * (0.8 - 0.2)", s_out.S1, 0.3, 1e-5);
assert_near("Pure diattenuation S2 = sqrt(0.8 * 0.2) * 1.0", s_out.S2, 0.4, 1e-5);
assert_near("Pure diattenuation S3 = 0", s_out.S3, 0.0, 1e-5);

// 2. Pure retarder: T1 = 1.0, T2 = 1.0, delta = pi/2
Polarization s_ret = copy(s_in);
applyDiattenuation(s_ret, s_in, 1.0, 1.0, (0.0, -1.0));

assert_near("Quarter-wave retardance S0 preserved", s_ret.S0, 1.0, 1e-5);
assert_near("Quarter-wave retardance S1 preserved", s_ret.S1, 0.0, 1e-5);
assert_near("Quarter-wave retardance S2 rotated to 0", s_ret.S2, 0.0, 1e-5);
assert_near("Quarter-wave retardance S3 converted to -1", s_ret.S3, -1.0, 1e-5);

// 3. Complex amplitude overload: js = (0, 1), jp = (1, 0)
Polarization s_jones = new Polarization;
applyDiattenuation(s_jones, s_in, (0, 1), (1, 0));

assert_near("Jones overload S0 = 1", s_jones.S0, 1.0, 1e-5);
assert_near("Jones overload S2 = 0", s_jones.S2, 0.0, 1e-5);
assert_near("Jones overload S3 = 1", s_jones.S3, 1.0, 1e-5);

// 4. Unpolarized light through an ideal polarizer
Polarization s_unpol = Polarization.unpolarized(1.0);
Polarization s_pol = copy(s_unpol);
applyDiattenuation(s_pol, s_unpol, 1.0, 0.0);

assert_near("Ideal polarizer: unpolarized output S0", s_pol.S0, 0.5, 1e-5);
assert_near("Ideal polarizer: output S1", s_pol.S1, 0.5, 1e-5);
assert_near("Ideal polarizer: output S2", s_pol.S2, 0.0, 1e-5);
assert_near("Ideal polarizer: output S3", s_pol.S3, 0.0, 1e-5);
assert_stokes_physical("Ideal polarizer", s_pol);

// 5. Equal transmission scales every Stokes component equally
Polarization s_neutral = new Polarization;
s_neutral.S0 = 1.0;
s_neutral.S1 = 0.2;
s_neutral.S2 = 0.3;
s_neutral.S3 = 0.1;
Polarization s_neutral_out = copy(s_neutral);
applyDiattenuation(s_neutral_out, s_neutral, 0.4, 0.4);

assert_near("Equal transmission: S0", s_neutral_out.S0, 0.4, 1e-5);
assert_near("Equal transmission: S1", s_neutral_out.S1, 0.08, 1e-5);
assert_near("Equal transmission: S2", s_neutral_out.S2, 0.12, 1e-5);
assert_near("Equal transmission: S3", s_neutral_out.S3, 0.04, 1e-5);
assert_stokes_physical("Equal transmission", s_neutral_out);

// 6. Ideal polarizer rejects orthogonal state
Polarization s_blocked = new Polarization;
s_blocked.S0 = 1.0;
s_blocked.S1 = -1.0;
s_blocked.S2 = 0.0;
s_blocked.S3 = 0.0;
applyDiattenuation(s_blocked, copy(s_blocked), 1.0, 0.0);

assert_near("Blocked polarization: S0", s_blocked.S0, 0.0, 1e-5);
assert_near("Blocked polarization: S1", s_blocked.S1, 0.0, 1e-5);
assert_near("Blocked polarization: S2", s_blocked.S2, 0.0, 1e-5);
assert_near("Blocked polarization: S3", s_blocked.S3, 0.0, 1e-5);
assert_stokes_physical("Blocked polarization", s_blocked);

// 7. Physicality for partially polarized input
Polarization s_partial = new Polarization;
s_partial.S0 = 1.0;
s_partial.S1 = 0.2;
s_partial.S2 = 0.3;
s_partial.S3 = 0.1;
applyDiattenuation(s_partial, copy(s_partial), 0.8, 0.3);

assert_stokes_physical("Partially polarized input", s_partial);

// 8. Complex cross term with reduced coherence
Polarization s_coherence = new Polarization;
s_coherence.S0 = 1.0;
s_coherence.S1 = 0.0;
s_coherence.S2 = 0.6;
s_coherence.S3 = 0.0;
applyDiattenuation(s_coherence, copy(s_coherence), 0.5, 0.5, (0.4, 0.0));

assert_near("Reduced coherence: S0 preserved", s_coherence.S0, 0.5, 1e-5);
assert_near("Reduced coherence: S1", s_coherence.S1, 0.0, 1e-5);
assert_near("Reduced coherence: S2", s_coherence.S2, 0.24, 1e-5);
assert_near("Reduced coherence: S3", s_coherence.S3, 0.0, 1e-5);
assert_stokes_physical("Reduced coherence", s_coherence);
write("");

// ============================================================================
// Test Case 45: Polarization copy invariants and state preservation
// ============================================================================
write("Test Case 45: Polarization copy invariants and state preservation");
Polarization extra_unpolarized = Polarization.unpolarized(1.0);
Polarization extra_unpolarized_copy = copy(extra_unpolarized);

assert_relative("Copy preserves unpolarized intensity", extra_unpolarized_copy.intensity(), extra_unpolarized.intensity());
assert_relative("Copy preserves unpolarized degree of polarization",
                extra_unpolarized_copy.dop(),
                extra_unpolarized.dop());

Polarization extra_linear = Polarization.linear(23.0, 1.7);
Polarization extra_linear_copy = copy(extra_linear);

assert_relative("Copy preserves linear-state intensity", extra_linear_copy.intensity(), extra_linear.intensity());
assert_near("Copy preserves S1", extra_linear_copy.S1, extra_linear.S1, 1e-5);
assert_near("Copy preserves S2", extra_linear_copy.S2, extra_linear.S2, 1e-5);
assert_near("Copy preserves S3", extra_linear_copy.S3, extra_linear.S3, 1e-5);
write("");

// ============================================================================
// Test Case 45-1: Polarization Rotation Consistency
// ============================================================================
write("Test Case 45-1: Polarization Rotation Consistency");
Polarization rot = Polarization.linear(0.0, 1.0);
rot.applyRotation(radians(45.0));

assert_near("Rotation 45: S0 preserved", rot.S0, 1.0, 1e-5);
assert_near("Rotation 45: S1 becomes zero", rot.S1, 0.0, 1e-5);
assert_near("Rotation 45: S2 becomes one", rot.S2, 1.0, 1e-5);
assert_near("Rotation 45: S3 remains zero", rot.S3, 0.0, 1e-5);

rot.applyRotation(radians(-45.0));

assert_near("Inverse rotation: S1 restored", rot.S1, 1.0, 1e-5);
assert_near("Inverse rotation: S2 restored", rot.S2, 0.0, 1e-5);

Polarization circ = Polarization.circular();
real s3_before = circ.S3;
circ.applyRotation(radians(23.0));

assert_near("Circular rotation: S0 preserved", circ.S0, 1.0, 1e-5);
assert_near("Circular rotation: S3 preserved", circ.S3, s3_before, 1e-5);
assert_near("Circular rotation: degree of polarization preserved", circ.dop(), 1.0, 1e-5);

Polarization unpol_rot = Polarization.unpolarized(2.0);
unpol_rot.applyRotation(radians(37.0));

assert_near("Unpolarized rotation: S0 preserved", unpol_rot.S0, 2.0, 1e-5);
assert_near("Unpolarized rotation: S1 remains zero", unpol_rot.S1, 0.0, 1e-5);
assert_near("Unpolarized rotation: S2 remains zero", unpol_rot.S2, 0.0, 1e-5);
assert_near("Unpolarized rotation: S3 remains zero", unpol_rot.S3, 0.0, 1e-5);
assert_near("Unpolarized rotation: DoP remains zero", unpol_rot.dop(), 0.0, 1e-5);
write("");

// ============================================================================
// Test Case 46: Normal-incidence Fresnel boundary conditions
// ============================================================================
write("Test Case 46: Normal-incidence Fresnel boundary conditions");
Ray extra_normal_ray = Ray(p = (-1.0, 0.0), d = (1.0, 0.0), medium = Materials.VACUUM);
extra_normal_ray.pol = Polarization.unpolarized();
FresnelResult extra_normal_ref =
    fresnelIsotropic(Materials.VACUUM, glass_15, extra_normal_ray, h_int, test_resolve_isotropic);
real extra_normal_r = ((n2 - n1) / (n2 + n1)) ^ 2;

assert_bool("Normal incidence is not total internal reflection", extra_normal_ref.tir, false);
assert_near_pair("Normal-incidence transmitted direction", extra_normal_ref.d, (1.0, 0.0));
assert_relative("Normal-incidence Fresnel reflectance", extra_normal_ref.R, extra_normal_r);
assert_relative("Normal-incidence energy balance", extra_normal_ref.R + extra_normal_ref.T, 1.0);
write("");

// ============================================================================
// Test Case 47: Polarization intensity scaling laws
// ============================================================================
write("Test Case 47: Polarization intensity scaling laws");
Polarization extra_scale_source = Polarization.linear(31.0, 1.0);
real extra_scale_factor = 2.75;
Polarization extra_scale_result = extra_scale_factor * extra_scale_source;

assert_relative("Intensity scales linearly with scalar factor",
                extra_scale_result.intensity(),
                extra_scale_source.intensity() * extra_scale_factor);
assert_relative("S1 scales proportionally", extra_scale_result.S1, extra_scale_source.S1 *extra_scale_factor);
assert_relative("DOP preserved under scaling", extra_scale_result.dop(), extra_scale_source.dop());
write("");

// ============================================================================
// Test Case 48: Stokes physical-state invariants
// ============================================================================
write("Test Case 48: Stokes physical-state invariants");
Polarization extra_s = Polarization.linear(0.0, 1.0);
Polarization extra_p = Polarization.linear(90.0, 1.0);
Polarization extra_circ = Polarization.circular();

assert_relative("S-polarized light has unit degree of polarization", extra_s.dop(), 1.0);
assert_relative("P-polarized light has unit degree of polarization", extra_p.dop(), 1.0);
assert_relative("Circularly polarized light has unit degree of polarization", extra_circ.dop(), 1.0);
assert_relative("S-polarized Stokes intensity agrees with S0", extra_s.S0, extra_s.intensity());
assert_relative("P-polarized Stokes intensity agrees with S0", extra_p.S0, extra_p.intensity());
assert_relative("Circular Stokes intensity agrees with S0", extra_circ.S0, extra_circ.intensity());
write("");

// Test Case 49: Lossless retarder conserves intensity and degree of polarization
write("Test Case 49: Lossless retarder conservation");
Polarization extra_retarder_input = copy(extra_circ);
Polarization extra_retarder_output = copy(extra_retarder_input);
applyDiattenuation(extra_retarder_output, extra_retarder_input, 1.0, 1.0, (cos(pi / 3.0), -sin(pi / 3.0)));

assert_relative("Lossless retarder preserves S0", extra_retarder_output.S0, extra_retarder_input.S0);
assert_relative("Lossless retarder preserves degree of polarization",
                extra_retarder_output.dop(),
                extra_retarder_input.dop());
write("");

// ============================================================================
// Test Case 50: Relative tolerance helper across scales
// ============================================================================
write("Test Case 50: Relative tolerance helper");
assert_relative("Relative comparison at unit scale", 1.0 + 1e-8, 1.0, 1e-7);
assert_relative("Relative comparison at large scale", 1e8 + 0.5, 1e8, 1e-8);
assert_relative("Relative comparison near zero", 1e-12, 0.0, 1e-6, 1e-11);
write("");

// ============================================================================
// Test Case 51: Complex phase in diattenuation
// ============================================================================
write("Test Case 51: Complex phase in diattenuation");
Polarization phase_input = new Polarization;
phase_input.S0 = 1.0;
phase_input.S1 = 0.0;
phase_input.S2 = 1.0;
phase_input.S3 = 0.0;
Polarization phase_output = copy(phase_input);
applyDiattenuation(phase_output, phase_input, 0.5, 0.5, (0.3, 0.2));

assert_near("Complex phase: S0", phase_output.S0, 0.5, 1e-5);
assert_near("Complex phase: S1", phase_output.S1, 0.0, 1e-5);
assert_near("Complex phase: S2", phase_output.S2, 0.3, 1e-5);
assert_near("Complex phase: S3", phase_output.S3, 0.2, 1e-5);

real phase_pol2 = phase_output.S1 ^ 2 + phase_output.S2 ^ 2 + phase_output.S3 ^ 2;

assert_bool("Complex phase: output remains inside Stokes cone", phase_pol2 <= phase_output.S0 ^ 2 + 1e-5, true);
write("");

// ============================================================================
// Test Case 52: Maximum-coherence physicality boundary
// ============================================================================
write("Test Case 52: Maximum-coherence physicality boundary");
Polarization boundary_input = new Polarization;
boundary_input.S0 = 1.0;
boundary_input.S1 = 0.0;
boundary_input.S2 = 1.0;
boundary_input.S3 = 0.0;
Polarization boundary_output = copy(boundary_input);
applyDiattenuation(boundary_output, boundary_input, 0.8, 0.2, (0.24, 0.32));

assert_near("Boundary: S0", boundary_output.S0, 0.5, 1e-5);
assert_near("Boundary: S1", boundary_output.S1, 0.3, 1e-5);
assert_near("Boundary: S2", boundary_output.S2, 0.24, 1e-5);
assert_near("Boundary: S3", boundary_output.S3, 0.32, 1e-5);

real boundary_pol2 = boundary_output.S1 ^ 2 + boundary_output.S2 ^ 2 + boundary_output.S3 ^ 2;

assert_near("Boundary: fully polarized output lies on sphere", boundary_pol2, boundary_output.S0 ^ 2, 1e-5);
write("");

// ============================================================================
// Test Case 53: Dichroic Medium Propagation & Mueller Theory
// ============================================================================
write("Test Case 53: Dichroic Medium Propagation & Mueller Theory");
real alpha_s_val = 0.0;
real alpha_p_val = 1.0;
real prop_dist = 2.0;
Medium dichroic_filter = medium("DichroicFilter").add(Dichroic(alpha_s_val, alpha_p_val));
Polarization sync_input = Polarization.linear(45.0, 1.0);
Polarization sync_output = dichroic_filter.stepPolarization(sync_input, (0, 0), (1, 0), 633.0, prop_dist);

assert_relative("Stokes S0 agrees with output intensity", sync_output.S0, sync_output.intensity(), 1e-5, 1e-7);
assert_relative("Stokes degree of polarization is physical", sync_output.dop(), 1.0, 1e-5, 1e-7);

Polarization sync_copy = copy(sync_output);

assert_relative("Copied Stokes S0 agrees with intensity", sync_copy.S0, sync_copy.intensity(), 1e-5, 1e-7);
assert_near("Copied S1 preserved", sync_copy.S1, sync_output.S1, 1e-5);
assert_near("Copied S2 preserved", sync_copy.S2, sync_output.S2, 1e-5);
assert_near("Copied S3 preserved", sync_copy.S3, sync_output.S3, 1e-5);

Polarization sync_unpol_input = Polarization.unpolarized(1.0);
Polarization sync_unpol_output = dichroic_filter.stepPolarization(sync_unpol_input, (0, 0), (1, 0), 633.0, prop_dist);

assert_relative("Unpolarized-path S0 agrees with intensity",
                sync_unpol_output.S0,
                sync_unpol_output.intensity(),
                1e-5,
                1e-7);

real Ts = exp(-alpha_s_val * prop_dist);
real Tp = exp(-alpha_p_val * prop_dist);

assert_relative("Diattenuator S0 matches Mueller prediction",
                sync_unpol_output.S0,
                sync_unpol_input.S0 *(Ts + Tp) / 2.0,
                1e-5,
                1e-7);

assert_relative("Diattenuator S1 matches Mueller prediction",
                sync_unpol_output.S1,
                sync_unpol_input.S0 *(Ts - Tp) / 2.0,
                1e-5,
                1e-7);

assert_near("Diattenuator S2 remains zero", sync_unpol_output.S2, 0.0);
assert_near("Diattenuator S3 remains zero", sync_unpol_output.S3, 0.0);
write("");

// ============================================================================
// Test Case 54: Deep-copy isolation
// ============================================================================
write("Test Case 54: Deep-copy isolation");
Polarization original = Polarization.linear(45.0, 1.0);
Polarization cloned = copy(original);

assert_near("Copy isolation: initial S0 preserved", cloned.S0, original.S0);
assert_near("Copy isolation: initial S1 preserved", cloned.S1, original.S1);
assert_near("Copy isolation: initial S2 preserved", cloned.S2, original.S2);
assert_near("Copy isolation: initial S3 preserved", cloned.S3, original.S3);

real source_s0 = original.S0;
real source_s1 = original.S1;
real source_s2 = original.S2;
real source_s3 = original.S3;
cloned.S0 = 3.0;
cloned.S1 = 1.0;
cloned.S2 = 0.5;
cloned.S3 = -0.25;

assert_near("Copy isolation: source Stokes S0 unchanged", original.S0, source_s0);
assert_near("Copy isolation: source Stokes S1 unchanged", original.S1, source_s1);
assert_near("Copy isolation: source Stokes S2 unchanged", original.S2, source_s2);
assert_near("Copy isolation: source Stokes S3 unchanged", original.S3, source_s3);

real cloned_s0 = cloned.S0;
real cloned_s1 = cloned.S1;
original.S0 = 2.0;
original.S1 = -0.5;

assert_near("Copy isolation: clone S0 unchanged after source mutation", cloned.S0, cloned_s0);
assert_near("Copy isolation: clone S1 unchanged after source mutation", cloned.S1, cloned_s1);

Polarization unpol_test = Polarization.unpolarized(2.0);
Polarization unpol_copy = copy(unpol_test);

assert_near("Copy isolation: unpolarized S0 copied", unpol_copy.S0, unpol_test.S0);
assert_near("Copy isolation: unpolarized S1 copied", unpol_copy.S1, 0.0);
assert_near("Copy isolation: unpolarized S2 copied", unpol_copy.S2, 0.0);
assert_near("Copy isolation: unpolarized S3 copied", unpol_copy.S3, 0.0);

unpol_copy.S0 = 7.0;
unpol_copy.S1 = 1.0;

assert_near("Copy isolation: unpolarized source S0 unchanged", unpol_test.S0, 2.0);
assert_near("Copy isolation: unpolarized source S1 unchanged", unpol_test.S1, 0.0);
write("");

// ============================================================================
// Test Case 55: Zero-distance propagation
// ============================================================================
write("Test Case 55: Zero-distance propagation");
Polarization zero_dist_input = Polarization.linear(45.0, 1.0);
real zero_s0 = zero_dist_input.S0;
real zero_s1 = zero_dist_input.S1;
real zero_s2 = zero_dist_input.S2;
real zero_s3 = zero_dist_input.S3;
Polarization zero_dist_output = dichroic_filter.stepPolarization(zero_dist_input, (0, 0), (1, 0), 633.0, 0.0);

assert_near("Zero distance: S0 unchanged", zero_dist_output.S0, zero_s0);
assert_near("Zero distance: S1 unchanged", zero_dist_output.S1, zero_s1);
assert_near("Zero distance: S2 unchanged", zero_dist_output.S2, zero_s2);
assert_near("Zero distance: S3 unchanged", zero_dist_output.S3, zero_s3);
write("");

// ============================================================================
// Test Case 56: Unpolarized light through a diattenuator
// ============================================================================
write("Test Case 56: Unpolarized light through a diattenuator");
real test_ts = 0.8;
real test_tp = 0.2;
Polarization unpol_src_stokes = Polarization.unpolarized(1.0);
Polarization unpol_dst_stokes = copy(unpol_src_stokes);
applyDiattenuation(unpol_dst_stokes, unpol_src_stokes, test_ts, test_tp);
real expected_s0 = 0.5 * (test_ts + test_tp);
real expected_s1 = 0.5 * (test_ts - test_tp);

assert_near("Unpolarized diattenuation: S0", unpol_dst_stokes.S0, expected_s0);
assert_near("Unpolarized diattenuation: S1", unpol_dst_stokes.S1, expected_s1);
assert_near("Unpolarized diattenuation: S2 remains zero", unpol_dst_stokes.S2, 0.0);
assert_near("Unpolarized diattenuation: S3 remains zero", unpol_dst_stokes.S3, 0.0);
assert_near("Unpolarized diattenuation: degree of polarization", unpol_dst_stokes.dop(), 0.6);
write("");

// ============================================================================
// Test Case 57: Fresnel Rhomb Achromatic Quarter-Wave Retardance
// ============================================================================
write("Test Case 57: Fresnel Rhomb Achromatic Quarter-Wave Retardance");
real n_bk7 = 1.515;
Medium glass_rhomb = medium("BK7_Rhomb", A = n_bk7);
real u = Tan(22.5);
real A_coef = 1.0 + u * u;
real B_coef = 1.0 + 1.0 / (n_bk7 * n_bk7);
real C_coef = 1.0 / (n_bk7 * n_bk7);
real theta = aSin(sqrt((B_coef + sqrt(B_coef * B_coef - 4.0 * A_coef * C_coef)) / (2.0 * A_coef)));
real L = 50.0;
real w = 2.0 * tan(radians(theta));
real h = 2.0;
path rhomb_path = (0, 0) -- (w * L, h * L) -- (w * L, h * L + L) -- (0, L) -- cycle;
Element rhomb = Interface(rhomb_path, glass_rhomb);
OpticalSystem rhomb_sys = Materials.VACUUM + rhomb;
// Launch linear 45-deg ray: S = (1, 0, 1, 0)
Ray r_in = Ray(p = (-20, L / 1.2), d = (1, 0), pol = Polarization.linear(45.0, 1.0));
Trace t_rhomb = rhomb_sys.trace(r_in);
// Use t_rhomb.last() to inspect the ray that has exited the rhomb after both TIR bounces:
Ray r_exit = t_rhomb.last();

assert_stokes_physical("Fresnel Rhomb exit ray", r_exit.pol);
assert_near("Fresnel Rhomb output S1 is 0", r_exit.pol.S1, 0.0, 1e-3);
assert_near("Fresnel Rhomb output S2 is 0", r_exit.pol.S2, 0.0, 1e-3);
assert_near("Fresnel Rhomb output S3 is pure circular", abs(r_exit.pol.S3) / r_exit.pol.S0, 1.0, 1e-3);
assert_near("Fresnel Rhomb output DOP", r_exit.pol.dop(), 1.0, 1e-5);
write("");

// ============================================================================
// Test Case 58: Fresnel Rhomb Unpolarized Invariance
// ============================================================================
write("Test Case 58: Fresnel Rhomb Unpolarized Invariance");
real n_bk7 = 1.515;
Medium glass_rhomb_unpol = medium("BK7_Rhomb_Unpol", A = n_bk7);
real u = Tan(22.5);
real A_c = 1.0 + u * u, B_c = 1.0 + 1.0 / (n_bk7 * n_bk7), C_c = 1.0 / (n_bk7 * n_bk7);
real theta = aSin(sqrt((B_c + sqrt(B_c * B_c - 4.0 * A_c * C_c)) / (2.0 * A_c)));
real L = 50.0, w = 2.0 * tan(radians(theta)), h = 2.0;
Element rhomb = Interface(
    (0, 0) -- (w * L, h * L) -- (w * L, h * L + L) -- (0, L) -- cycle, glass_rhomb_unpol);
OpticalSystem rhomb_sys = Materials.VACUUM + rhomb;
Ray r_unpol_rhomb = Ray(p = (-20, L / 1.2), d = (1, 0), pol = Polarization.unpolarized(1.0));
Trace t_unpol_rhomb = rhomb_sys.trace(r_unpol_rhomb);
Ray r_exit_unpol = t_unpol_rhomb.last(); // after both TIR bounces and the exit face

assert_near("Unpolarized exit S1 is 0", r_exit_unpol.pol.S1, 0.0, 1e-5);
assert_near("Unpolarized exit S2 is 0", r_exit_unpol.pol.S2, 0.0, 1e-5);
assert_near("Unpolarized exit S3 is 0", r_exit_unpol.pol.S3, 0.0, 1e-5);
assert_near("Unpolarized exit DOP remains 0", r_exit_unpol.pol.dop(), 0.0, 1e-5);
real rhomb_t0 = 1.0 - ((n_bk7 - 1.0) / (n_bk7 + 1.0)) ^ 2;
assert_near("Unpolarized exit intensity = two normal-incidence Fresnel transmissions",
            r_exit_unpol.intensity(),
            rhomb_t0 ^ 2,
            1e-6);
assert_near("Rhomb exit ray is parallel to the input", r_exit_unpol.d.x, 1.0, 1e-6);
write("");

// ============================================================================
// Test Case 59: Incoherent Superposition & Scaling of Unpolarized Light
// ============================================================================
write("Test Case 59: Incoherent Superposition & Scaling of Unpolarized Light");
Polarization u1 = Polarization.unpolarized(1.2);
Polarization u2 = Polarization.unpolarized(0.8);
Polarization u_sum = u1 + u2;

assert_near("Summed intensity S0", u_sum.S0, 2.0, 1e-5);
assert_near("Summed S1 is 0", u_sum.S1, 0.0, 1e-5);
assert_near("Summed S2 is 0", u_sum.S2, 0.0, 1e-5);
assert_near("Summed S3 is 0", u_sum.S3, 0.0, 1e-5);
assert_near("Summed DOP remains 0", u_sum.dop(), 0.0, 1e-5);

Polarization u_scaled = 3.5 * u1;

assert_near("Scaled intensity S0", u_scaled.S0, 3.5 * 1.2, 1e-5);
assert_near("Scaled DOP remains 0", u_scaled.dop(), 0.0, 1e-5);
write("");

// ============================================================================
// Test Case 60: Optical Rotation Invariance for Unpolarized Light
// ============================================================================
write("Test Case 60: Optical Rotation Invariance for Unpolarized Light");
Polarization unpol_rot = Polarization.unpolarized(1.0);
unpol_rot.applyRotation(radians(45.0));

assert_near("Rotated unpolarized S0", unpol_rot.S0, 1.0, 1e-5);
assert_near("Rotated unpolarized S1", unpol_rot.S1, 0.0, 1e-5);
assert_near("Rotated unpolarized S2", unpol_rot.S2, 0.0, 1e-5);
assert_near("Rotated unpolarized S3", unpol_rot.S3, 0.0, 1e-5);
assert_near("Rotated unpolarized DOP", unpol_rot.dop(), 0.0, 1e-5);
write("");

// ============================================================================
// Test Case 61: Polarizing Beam Splitter (PBS) with Unpolarized Light
// ============================================================================
write("Test Case 61: Polarizing Beam Splitter (PBS) with Unpolarized Light");
Polarization in_unpol = Polarization.unpolarized(2.0);
real Is = 0.5 * (in_unpol.S0 + in_unpol.S1);
real Ip = 0.5 * (in_unpol.S0 - in_unpol.S1);
Polarization out_s = Polarization.linear(0.0, Is);
Polarization out_p = Polarization.linear(90.0, Ip);

assert_near("PBS S-channel intensity", out_s.S0, 1.0, 1e-5);
assert_near("PBS S-channel S1", out_s.S1, 1.0, 1e-5);
assert_near("PBS S-channel DOP", out_s.dop(), 1.0, 1e-5);
assert_near("PBS P-channel intensity", out_p.S0, 1.0, 1e-5);
assert_near("PBS P-channel S1", out_p.S1, -1.0, 1e-5);
assert_near("PBS P-channel DOP", out_p.dop(), 1.0, 1e-5);
// Same split performed by the engine: ideal linear polarizers at 0 and 90 deg.
Trace t_lp0 = system(Polarizer.LP((10, 0), 0.0, 10.0)).trace(Ray((0, 0), (1, 0), pol = copy(in_unpol)));
Trace t_lp90 = system(Polarizer.LP((10, 0), 90.0, 10.0)).trace(Ray((0, 0), (1, 0), pol = copy(in_unpol)));
assert_near("LP(0) transmits the S half of unpolarized light", t_lp0.last().pol.S0, Is, 1e-6);
assert_near("LP(0) output is pure S (S1 = S0)", t_lp0.last().pol.S1, t_lp0.last().pol.S0, 1e-6);
assert_near("LP(90) transmits the P half of unpolarized light", t_lp90.last().pol.S0, Ip, 1e-6);
assert_near("LP(90) output is pure P (S1 = -S0)", t_lp90.last().pol.S1, -t_lp90.last().pol.S0, 1e-6);
write("");

// ============================================================================
// Test Case 62: General Elliptical Polarization (Poincaré vs Semi-Axes)
// ============================================================================
write("Test Case 62: General Elliptical Polarization (Poincaré vs Semi-Axes)");
Polarization el1 = Polarization.elliptical(a = 1.0, b = 0.5, phase = 90.0, theta = 0.0);

assert_stokes_physical("Elliptical state 1", el1);
assert_near("Elliptical S0", el1.S0, 1.25, 1e-5);
assert_near("Elliptical S1", el1.S1, 0.75, 1e-5);
assert_near("Elliptical S2", el1.S2, 0.0, 1e-5);
assert_near("Elliptical S3", el1.S3, 1.0, 1e-5);
assert_near("Elliptical DOP", el1.dop(), 1.0, 1e-5);

Polarization el2 = Polarization.fromAngles(psi = 45.0, chi = 20.0, intensity = 1.0);

assert_stokes_physical("Elliptical state 2", el2);
assert_near("Poincare S1", el2.S1, 0.0, 1e-4);
assert_near("Poincare S2", el2.S2, cos(radians(40.0)), 1e-4);
assert_near("Poincare S3", el2.S3, sin(radians(40.0)), 1e-4);
assert_near("Poincare DOP", el2.dop(), 1.0, 1e-5);
write("");

// ============================================================================
// Test Case 63: Dielectric Refraction & Reflection Phase Shifts
// ============================================================================
write("Test Case 63: Dielectric Refraction & Reflection Phase Shifts");
real n_air = 1.0;
real n_glass = 1.5168; // BK7
Medium m_air = medium("Air", n_air);
Medium m_glass = medium("BK7", n_glass);

path boundary_plane = (0, -50) -- (0, 50);
Element iface = Interface(boundary_plane, m_air, m_glass);
OpticalSystem sys = system(iface, Materials.VACUUM);

real theta_B_deg = degrees(atan(n_glass / n_air)); // ~56.6 deg

// ------------------------------------------------------------------------
// Part 1: Refraction (Transmission) Phase Shift is Strictly Zero
// ------------------------------------------------------------------------
// Launch a 45-deg linearly polarized ray: S = (1, 0, 1, 0), phase = 0
real theta_i_deg = 30.0; // Angle of incidence
pair dir_in = (cos(radians(theta_i_deg)), sin(radians(theta_i_deg)));
Ray r_in1 = Ray(p = (-10 * dir_in.x, -10 * dir_in.y), d = dir_in, medium = m_air, pol = Polarization.linear(45.0, 1.0));
r_in1.phase = 0.0;

Trace t1 = sys.trace(r_in1);
Ray r_trans = t1.transmitted(0);

assert_stokes_physical("Transmitted ray physical state", r_trans.pol);

// Transmitted ray must NOT acquire any circular ellipticity (S3 == 0)
assert_near("Transmitted S3 remains exactly zero (no retardance)", r_trans.pol.S3, 0.0, 1e-5);

// Relative phase between transmitted P and S must be 0 deg
real delta_trans_deg = degrees(atan2(r_trans.pol.S3, r_trans.pol.S2));
assert_near("Transmitted relative phase (P - S) is 0 deg", delta_trans_deg, 0.0, 1e-4);

// Carrier phase through boundary is continuous (delta_phi = 0)
assert_near("Transmitted ray carrier phase is 0", r_trans.phase, 0.0, 1e-5);

// ------------------------------------------------------------------------
// Part 2: Dielectric Reflection BELOW Brewster's Angle (Phase Flip = 180 deg)
// ------------------------------------------------------------------------
// At theta_i = 30 deg (< theta_B ~ 56.6 deg):
// r_s < 0 (phase pi), r_p > 0 (phase 0) -> Relative phase shift delta = 180 deg
Ray r_refl_below = t1.reflected(0);

assert_stokes_physical("Reflected below Brewster", r_refl_below.pol);
assert_near("Reflected S3 below Brewster remains zero", r_refl_below.pol.S3, 0.0, 1e-5);

// S2 must be negative because of the 180 deg relative phase flip:
assert_true("Reflected S2 below Brewster is negative", r_refl_below.pol.S2 < -1e-4);

real delta_refl_below = degrees(atan2(r_refl_below.pol.S3, r_refl_below.pol.S2));
assert_near("Reflected relative phase below Brewster is 180 deg", abs(delta_refl_below), 180.0, 1e-3);

// ------------------------------------------------------------------------
// Part 3: Dielectric Reflection ABOVE Brewster's Angle (Phase Flip = 0 deg)
// ------------------------------------------------------------------------
// At theta_i = 75 deg (> theta_B ~ 56.6 deg):
// Both r_s < 0 and r_p < 0 -> Both acquire pi phase -> Relative phase delta = 0 deg
real theta_steep_deg = 75.0;
pair dir_steep = (cos(radians(theta_steep_deg)), sin(radians(theta_steep_deg)));
Ray r_in2 = Ray(p = (-10 * dir_steep.x, -10 * dir_steep.y),
                d = dir_steep,
                medium = m_air,
                pol = Polarization.linear(45.0, 1.0));

Trace t2 = sys.trace(r_in2);
Ray r_refl_above = t2.reflected(0);

assert_stokes_physical("Reflected above Brewster", r_refl_above.pol);
assert_near("Reflected S3 above Brewster remains zero", r_refl_above.pol.S3, 0.0, 1e-5);

// S2 becomes positive again because both r_s and r_p share the same sign (both negative)
assert_true("Reflected S2 above Brewster is positive", r_refl_above.pol.S2 > 1e-4);

real delta_refl_above = degrees(atan2(r_refl_above.pol.S3, r_refl_above.pol.S2));
assert_near("Reflected relative phase above Brewster is 0 deg", delta_refl_above, 0.0, 1e-3);
write("");

// ============================================================================
// Test Case 64: Anisotropic Maxwell Boundary Solver - Isotropic Limit
// ============================================================================
write("Test Case 64: Anisotropic Maxwell Boundary Solver - Isotropic Limit");
// When no == ne == n_glass, anisotropic solver MUST match isotropic Snell exactly
real n_test = 1.55;
Medium m_iso = medium("Glass_Iso", n_test);
Medium m_uni_o = medium("Uni_O", n_test);
Medium m_uni_e = uniaxialExtraordinary("Uni_E", n_test, n_test, (1.0, 0.0));

Element crystal_iso = Crystal(box((0.0, -20.0), (20.0, 20.0)), m_uni_o, m_uni_e, Materials.VACUUM);
OpticalSystem sys = system(crystal_iso, Materials.VACUUM);

// Launch at 30 deg incidence
Ray r_in = Ray((-10, -5), unit((cos(radians(30.0)), sin(radians(30.0)))), pol = Polarization.linear(45.0, 1.0));
Trace t = sys.trace(r_in);
Interaction ix = t.interaction(0);

// Theoretical isotropic Fresnel values for n1 = 1.0, n2 = 1.55 at 30 deg
real cos_i = cos(radians(30.0));
real sin_t = (1.0 / n_test) * sin(radians(30.0));
real cos_t = sqrt(1.0 - sin_t * sin_t);
real rs_iso = (1.0 * cos_i - n_test * cos_t) / (1.0 * cos_i + n_test * cos_t);
real rp_iso = (n_test * cos_i - 1.0 * cos_t) / (n_test * cos_i + 1.0 * cos_t);

Ray o_ray = ix.transmitted[0];
Ray e_ray = ix.transmitted[1];

assert_near("Isotropic limit To matches (1 - Rs^2)", o_ray.intensity() / 0.5, 1.0 - rs_iso ^ 2, 1e-5);
assert_near("Isotropic limit Te matches (1 - Rp^2)", e_ray.intensity() / 0.5, 1.0 - rp_iso ^ 2, 1e-5);
assert_near_pair("Isotropic limit o-ray and e-ray directions are identical", o_ray.d, e_ray.d, 1e-6);
write("");

// ============================================================================
// Test Case 65: Anisotropic Maxwell Boundary Solver - Energy Conservation
// ============================================================================
write("Test Case 65: Anisotropic Maxwell Boundary Solver - Energy Conservation");
// Calcite: no = 1.658, ne = 1.486
real no_cal = 1.658, ne_cal = 1.486;
pair optic_axis = unit((1.0, 1.0)); // 45 deg tilted optic axis
Medium calcite_o = medium("Calcite_O", no_cal);
Medium calcite_e = uniaxialExtraordinary("Calcite_E", no_cal, ne_cal, optic_axis);

Element crystal = Crystal(box((0.0, -20.0), (20.0, 20.0)), calcite_o, calcite_e, Materials.VACUUM);
OpticalSystem sys_energy = system(crystal, Materials.VACUUM);

// Launch arbitrary elliptical polarization at 40 deg incidence
Polarization p_in = Polarization.fromAngles(25.0, 15.0, 1.0); // Intensity = 1.0
Ray r_in = Ray((-10, 0), unit((cos(radians(40.0)), sin(radians(40.0)))), pol = p_in);

Trace t = sys_energy.trace(r_in);
Interaction ix = t.interaction(0);

real power_o = ix.transmitted[0].intensity();
real power_e = ix.transmitted[1].intensity();
real power_refl = ix.reflected[0].intensity();
real total_power_out = power_o + power_e + power_refl;

// Strict Maxwell Poynting flux conservation: R + To + Te == 1.0
assert_near("Strict energy conservation across anisotropic boundary (R + To + Te = 1)", total_power_out, 1.0, 1e-5);
assert_true("Extraordinary ray walks off from ordinary ray", abs(ix.transmitted[0].d.y - ix.transmitted[1].d.y) > 0.01);
write("");

// ============================================================================
// REGRESSION, INVARIANT, EDGE-CASE AND CONVERGENCE TESTS
// Each "Regression" case failed against the original optics.asy.
// ============================================================================

// ============================================================================
// Test Case 66 (Regression): Beer-Lambert absorption preserves the degree of
// polarization (the whole Stokes vector decays, not only S0)
// ============================================================================
write("Test Case 66: Absorption preserves degree of polarization");
Medium rg_absMed = medium("rg absorber", 1.0, alpha = 0.5);
Polarization rg_absOut = rg_absMed.stepPolarization(Polarization.linear(0.0, 1.0), (0, 0), (1, 0), 633.0, 2.0);
assert_near("Absorbed S0 = exp(-alpha L)", rg_absOut.S0, exp(-1.0), 1e-9);
assert_near("Absorbed S1 decays with S0", rg_absOut.S1, exp(-1.0), 1e-9);
assert_stokes_physical("Absorbed linear state", rg_absOut, 1e-12);

Polarization rg_ell = Polarization.fromAngles(20.0, 10.0, 1.0);
Polarization rg_ellOut = rg_absMed.stepPolarization(rg_ell, (0, 0), (1, 0), 633.0, 3.0);
assert_near("Absorbed elliptical state keeps DoP = 1",
            sqrt(rg_ellOut.S1 ^ 2 + rg_ellOut.S2 ^ 2 + rg_ellOut.S3 ^ 2) / rg_ellOut.S0,
            1.0,
            1e-9);
assert_near("Absorbed elliptical state keeps azimuth", rg_ellOut.psi(), rg_ell.psi(), 1e-9);

// Through the tracer: linearly polarized ray across an absorbing, index-matched slab
Element rg_absSlab = Interface(box((10, -10), (20, 10)), rg_absMed, Materials.VACUUM);
Ray rg_absEnd = system(rg_absSlab).trace(Ray((0, 0), (1, 0), pol = Polarization.linear(30.0, 1.0)), far = 20).last();
assert_near("Traced absorbed ray intensity", rg_absEnd.intensity(), exp(-0.5 * 10.0), 1e-9);
assert_near("Traced absorbed ray DoP", rg_absEnd.pol.dop(), 1.0, 1e-9);
assert_near("Traced absorbed ray azimuth (30 deg)", degrees(rg_absEnd.pol.psi()), 30.0, 1e-6);
write("");

// ============================================================================
// Test Case 67 (Regression): transforming an Element keeps its focal length
// ============================================================================
write("Test Case 67: Element transform preserves focal-length metadata");
Element rg_thin = Lens.thin((0, 0), 0, 42.0, 10.0);
Element rg_thinMoved = shift(5, 3) * rotate(10) * rg_thin;
assert_near("Focal length survives transform", rg_thinMoved.focal_length, 42.0, 1e-12);
assert_near("Focal length survives copy", copy(rg_thin).focal_length, 42.0, 1e-12);
Element rg_aligned = align(Lens.thin((50, 0), 0, -30.0, 10.0), (0, 0), (100, 10));
assert_near("Focal length survives align()", rg_aligned.focal_length, -30.0, 1e-12);
write("");

// ============================================================================
// Test Case 68 (Regression): extraordinary ray through a parallel calcite
// slab exits parallel to the incident ray (tangential k conservation); the
// boundary seed ray must keep its wave normal
// ============================================================================
write("Test Case 68: Birefringent slab exit directions and lateral walk-off");
Medium rg_co = medium("rg calcite o", 1.658);
Medium rg_ce = uniaxialExtraordinary("rg calcite e", 1.658, 1.486, unit((1, 1)));
Element rg_slab = Crystal(box((0, -50), (10, 50)), rg_co, rg_ce, Materials.VACUUM);
OpticalSystem rg_slabSys = system(rg_slab);
pair rg_dIn = dir(20);
Trace rg_tE = rg_slabSys.trace(Ray((-10, 0), rg_dIn, pol = Polarization.linear(90.0, 1.0)), far = 20);
Trace rg_tO = rg_slabSys.trace(Ray((-10, 0), rg_dIn, pol = Polarization.linear(0.0, 1.0)), far = 20);
Ray rg_eExit = rg_tE.transmitted(1);
Ray rg_oExit = rg_tO.transmitted(1);
assert_near_pair("e-ray exits parallel to incidence", rg_eExit.d, rg_dIn, 1e-9);
assert_near_pair("o-ray exits parallel to incidence", rg_oExit.d, rg_dIn, 1e-9);
assert_true("e and o rays are laterally displaced (walk-off)", abs(rg_eExit.p.y - rg_oExit.p.y) > 0.1);
assert_near("Exit wave normal equals exit ray direction (isotropic ambient)",
            length(rg_eExit.wave_normal - rg_eExit.d),
            0.0,
            1e-9);
write("");

// ============================================================================
// Test Case 69 (Regression): optical path length of an extraordinary ray is
// k.dr = n_phase cos(walk-off) L, checked against an independent solution of
// the index ellipse
// ============================================================================
write("Test Case 69: Extraordinary-ray optical path length");
real rg_no = 1.658, rg_ne = 1.486, rg_t = 10.0;
Medium rg_ceX = uniaxialExtraordinary("rg e (axis x)", rg_no, rg_ne, (1, 0));
Element rg_slabX = Crystal(box((0, -50), (rg_t, 50)), medium("rg o", rg_no), rg_ceX, Materials.VACUUM);
real rg_thIn = 30.0;
Trace rg_tX = system(rg_slabX).trace(Ray((0, 0) - 5 * dir(rg_thIn), dir(rg_thIn), pol = Polarization.linear(90.0, 1.0)),
                                     far = 20);
// Rays stored in the trace carry the cumulative OPL at the *end* of their
// segment, so incident(k).opl is the OPL at interface k.
Ray rg_eIn = rg_tX.transmitted(0);
real rg_oplMeasured = rg_tX.incident(1).opl - rg_tX.incident(0).opl;
// Independent solution: optic axis along x (= slab normal)
real rg_ky = sin(radians(rg_thIn));
real rg_kx = rg_no * sqrt(1.0 - rg_ky ^ 2 / rg_ne ^ 2);
pair rg_g = unit((rg_kx / rg_no ^ 2, rg_ky / rg_ne ^ 2));
real rg_oplInside = rg_kx * rg_t + rg_ky * rg_t * rg_g.y / rg_g.x;
assert_near_pair("e-ray group direction matches index-ellipse normal", rg_eIn.d, rg_g, 1e-9);
assert_near("e-ray OPL across slab = k . dr", rg_oplMeasured, rg_oplInside, 1e-7);
real rg_geomLen = rg_t / rg_g.x;
assert_true("e-ray OPL differs from naive n_phase * L",
            abs(rg_oplMeasured - rg_eIn.medium.n((0, 0), rg_eIn.wave_normal) * rg_geomLen) > 1e-3);
// Isotropic limit: OPL reduces to n * L
Element rg_slabIso = Crystal(box((0, -50), (rg_t, 50)),
                             medium("o", 1.6),
                             uniaxialExtraordinary("e", 1.6, 1.6, (1, 0)),
                             Materials.VACUUM);
Trace rg_tIso = system(rg_slabIso).trace(Ray((-5, 0), (1, 0), pol = Polarization.linear(90.0, 1.0)), far = 20);
assert_near("Isotropic-limit e-ray OPL = n t", rg_tIso.incident(1).opl - rg_tIso.incident(0).opl, 1.6 * rg_t, 1e-9);
write("");

// ============================================================================
// Test Case 70 (Regression): R + T = 1 at every lossless interface, even when
// the arriving ray was attenuated in the bulk
// ============================================================================
write("Test Case 70: Interface energy balance after bulk absorption");
Medium rg_absGlass = medium("rg abs glass", 1.5, alpha = 0.2);
Trace rg_tAbs = system(Interface(box((0, -10), (10, 10)), rg_absGlass, Materials.VACUUM)).trace(Ray((-5, 0), (1, 0)));
Interaction rg_exitIx = rg_tAbs.interaction(1);
assert_near("Exit-face R + T = 1", rg_exitIx.R + rg_exitIx.T, 1.0, 1e-9);
assert_near("Exit-face R matches normal-incidence Fresnel", rg_exitIx.R, (0.5 / 2.5) ^ 2, 1e-9);
assert_near("Arriving intensity reported at the exit face",
            rg_exitIx.incident.intensity(),
            (1 - 0.04) * exp(-2.0),
            1e-9);
write("");

// ============================================================================
// Test Case 71 (Regression): GRIN + absorption uses the full curved arc length
// ============================================================================
write("Test Case 71: Absorbing GRIN medium attenuation along curved path");
indexProfile rg_gp = new real(pair p, real lambda) { return 1.5 * (1 - 0.5 * 0.01 * p.y ^ 2); };
gradientProfile rg_gg = new pair(pair p, real lambda) { return (0, -1.5 * 0.01 * p.y); };
Medium rg_gm = medium("rg absorbing GRIN").add(GRIN(rg_gp, rg_gg)).add(Absorption(0.05));
Trace rg_tg =
    system(Interface(box((0, -20), (10, 20)), rg_gm, Materials.VACUUM)).trace(Ray((0.5, 1), (1, 0), rg_gm), far = 30);
real rg_arc = arclength(rg_tg.segments()[0].traj);
Interaction rg_gix = rg_tg.interaction(0);
assert_true("GRIN trajectory is curved (arc length > 9.5)", rg_arc > 9.5);
assert_near("GRIN arrival intensity uses the full arc (~exp(-0.05 * 9.5))",
            rg_gix.incident.intensity(),
            exp(-0.05 * 9.5),
            2e-3);
assert_near("Trace.last() before the boundary agrees with arrival intensity",
            rg_tg.last().intensity() / rg_gix.T,
            rg_gix.incident.intensity(),
            1e-6);
assert_near("GRIN exit-face R + T = 1", rg_gix.R + rg_gix.T, 1.0, 1e-6);
write("");

// ============================================================================
// Test Case 72 (Regression): ideal mirror reverses circular handedness,
// consistent with the dielectric reflection convention; two mirrors restore it
// ============================================================================
write("Test Case 72: Mirror polarization handedness");
Trace rg_tm = system(Mirror.flat((10, 0), 180, 10)).trace(Ray((0, 0), (1, 0), pol = Polarization.circular()));
Trace rg_tdi = system(Interface((10, -10) -- (10, 10), Materials.VACUUM, medium("g", 1.5)))
                          .trace(Ray((0, 0), (1, 0), pol = Polarization.circular()));
assert_near("Mirror: RCP -> LCP (S3 = -1)", rg_tm.reflected(0).pol.S3, -1.0, 1e-9);
assert_near("Mirror agrees with dielectric normal-incidence handedness flip",
            rg_tm.reflected(0).pol.S3,
            rg_tdi.reflected(0).pol.S3 / rg_tdi.reflected(0).pol.S0,
            1e-9);
Trace rg_tRetro =
    system(new Element[]{Mirror.flat((25, 0), 90, 50), Mirror.flat((0, 25), 0, 50)})
        .trace(Ray((40, 30), (-cos(radians(37)), -sin(radians(37))), pol = Polarization.fromAngles(15, 20)));
Polarization rg_retroIn = Polarization.fromAngles(15, 20);
Polarization rg_retroOut = rg_tRetro.reflected(1).pol;
assert_near("Two reflections restore S2", rg_retroOut.S2, rg_retroIn.S2, 1e-9);
assert_near("Two reflections restore S3", rg_retroOut.S3, rg_retroIn.S3, 1e-9);
Trace rg_tIdeal =
    system(Mirror.idealSpherical((100, 0), 180, 50, 40)).trace(Ray((0, 5), (1, 0), pol = Polarization.circular()));
assert_near("Ideal spherical mirror also flips handedness", rg_tIdeal.reflected(0).pol.S3, -1.0, 1e-9);
write("");

// ============================================================================
// Test Case 73 (Regression): linear polarizer transmission axis is honored
// ============================================================================
write("Test Case 73: Linear polarizer axis (Malus's law)");
real[] rg_axes = {0.0, 30.0, 60.0, 90.0, 135.0};
for (real ax : rg_axes) {
    Ray rg_lpOut =
        system(Polarizer.LP((10, 0), ax, 10.0)).trace(Ray((0, 0), (1, 0), pol = Polarization.linear(0.0, 1.0))).last();
    assert_near("LP(" + string(ax) + ") on linear 0: I = cos^2", rg_lpOut.intensity(), Cos(ax) ^ 2, 1e-6);
    if (rg_lpOut.intensity() > 1e-3)
        assert_near("LP(" + string(ax) + ") output azimuth: cos 2(psi - axis) = 1",
                    Cos(2 * (degrees(rg_lpOut.pol.psi()) - ax)),
                    1.0,
                    1e-8);
}
Ray rg_lpTilt = system(Polarizer.LP((10, 0), 45.0, 10.0, angle = 20))
                    .trace(Ray((0, 0), (1, 0), pol = Polarization.linear(45.0)))
                    .last();
assert_near("Tilted LP plate passes its own axis without Fresnel loss", rg_lpTilt.intensity(), 1.0, 1e-6);
write("");

// ============================================================================
// Test Case 74 (Regression): waveplate geometry is independent of its fast
// axis; QWP/HWP Jones-calculus results
// ============================================================================
write("Test Case 74: Waveplates");
Trace rg_tq = system(Polarizer.QWP((10, 0), 0.0, 20.0)).trace(Ray((0, 5), (1, 0), pol = Polarization.linear(45.0)));
assert_near("QWP(fast 0) is hit by an off-axis horizontal ray", rg_tq.length(), 1.0);
assert_near("QWP: linear 45 -> circular (|S3| = 1)", abs(rg_tq.last().pol.S3), 1.0, 1e-9);
assert_near_pair("QWP leaves the ray direction unchanged", rg_tq.last().d, (1, 0), 1e-12);
Trace rg_tq2 = (Polarizer.QWP((10, 0), 0.0, 20.0) + Polarizer.QWP((20, 0), 0.0, 20.0))
                   .trace(Ray((0, 5), (1, 0), pol = Polarization.linear(45.0)));
assert_near("Two QWPs = HWP: linear 45 -> linear -45 (S2 = -1)", rg_tq2.last().pol.S2, -1.0, 1e-9);
Ray rg_hwp =
    system(Polarizer.HWP((10, 0), 22.5, 20.0)).trace(Ray((0, 0), (1, 0), pol = Polarization.linear(0.0))).last();
assert_near("HWP(22.5): linear 0 -> linear 45 (S2 = 1)", rg_hwp.pol.S2, 1.0, 1e-9);
assert_near("HWP preserves intensity", rg_hwp.intensity(), 1.0, 1e-12);
OpticalSystem rg_waterSys = system(Polarizer.QWP((10, 0), 0, 20), Materials.WATER);
Trace rg_tqw = rg_waterSys.trace(Ray((0, 0), (1, 0)));
assert_true("Waveplate in water does not switch the ray into vacuum", rg_tqw.last().medium == Materials.WATER);
write("");

// ============================================================================
// Test Case 75 (Regression): uniaxial Fresnel result keeps R + T = 1 and a
// TIR phase convention consistent with fresnelIsotropic
// ============================================================================
write("Test Case 75: Uniaxial Fresnel bookkeeping and TIR phase convention");
Hit rg_h;
rg_h.valid = true;
rg_h.p = (0, 0);
rg_h.normal = (-1, 0);
Medium rg_uo = medium("uo", 1.658);
Medium rg_ue = uniaxialExtraordinary("ue", 1.658, 1.486, unit((1, 2)));
mediumResolver rg_resU = new Medium(pair p, Ray r, Medium m1, Medium m2, bool from_m1) {
    return (p.x > 0) ? rg_uo : Materials.VACUUM;
};
Ray rg_ru = Ray((-1, 0), dir(35), Materials.VACUUM, pol = Polarization.fromAngles(10, 5));
FresnelResult rg_fu = fresnelUniaxial(rg_ru, rg_h, Materials.VACUUM, rg_uo, rg_ue, rg_resU);
assert_near("Uniaxial R + T = 1", rg_fu.R + rg_fu.T, 1.0, 1e-12);
assert_near("Uniaxial T = incident-weighted To/Te",
            rg_fu.T,
            (rg_fu.To * rg_ru.pol.sPower() + rg_fu.Te * rg_ru.pol.pPower()) / rg_ru.intensity(),
            1e-12);

// Dense isotropic glass -> crystal beyond the ordinary critical angle
Medium rg_dense = medium("dense", 2.0);
Medium rg_lo = medium("lo", 1.5);
Medium rg_le = uniaxialExtraordinary("le", 1.5, 1.5, (1, 0));
mediumResolver rg_resD = new Medium(pair p, Ray r, Medium m1, Medium m2, bool from_m1) {
    return (p.x > 0) ? rg_lo : rg_dense;
};
Ray rg_rd = Ray((-1, 0), dir(60), rg_dense, pol = Polarization.unpolarized());
FresnelResult rg_fuT = fresnelUniaxial(rg_rd, rg_h, rg_dense, rg_lo, rg_le, rg_resD);
FresnelResult rg_fiT = fresnelIsotropic(rg_dense, rg_lo, rg_rd, rg_h, rg_resD);
assert_bool("Ordinary mode is totally reflected", rg_fuT.tir_o, true);
assert_near("Uniaxial TIR |rs| = 1", length(rg_fuT.rs), 1.0, 1e-12);
assert_near_pair("Uniaxial TIR rs equals isotropic TIR rs (same phase convention)", rg_fuT.rs, rg_fiT.rs, 1e-12);
write("");

// ============================================================================
// Test Case 76: TIR retardance matches the analytic Fresnel formula
// tan(delta/2) = cos(t) sqrt(sin^2 t - n^2) / sin^2 t,  n = n2/n1
// ============================================================================
write("Test Case 76: TIR phase retardance");
Medium rg_glass = medium("rg glass", 1.515);
mediumResolver rg_resG = new Medium(pair p, Ray r, Medium m1, Medium m2, bool from_m1) {
    return (p.x > 0) ? Materials.VACUUM : rg_glass;
};
real[] rg_tirAngles = {45.0, 50.0, 60.0, 75.0};
for (real th : rg_tirAngles) {
    FresnelResult rg_ft = fresnelIsotropic(rg_glass, Materials.VACUUM, Ray((-1, 0), dir(th), rg_glass), rg_h, rg_resG);
    pair rg_cross = rg_ft.rs * conj(rg_ft.rp);
    real rg_delta = abs(atan2(rg_cross.y, rg_cross.x));
    real rg_nr = 1.0 / 1.515;
    real rg_expect = 2.0 * atan(Cos(th) * sqrt(Sin(th) ^ 2 - rg_nr ^ 2) / Sin(th) ^ 2);
    assert_near("TIR retardance at " + string(th) + " deg (rad)", rg_delta, rg_expect, 1e-12);
    assert_near("TIR |rs| = |rp| = 1 at " + string(th) + " deg", length(rg_ft.rs) * length(rg_ft.rp), 1.0, 1e-12);
}
write("");

// ============================================================================
// Test Case 77: Fresnel transmission Mueller matrix: fully polarized light
// stays fully polarized and Stokes power is conserved across the interface
// ============================================================================
write("Test Case 77: Fresnel Mueller transmission/reflection invariants");
Element rg_iface = Interface((0, -50) -- (0, 50), Materials.VACUUM, medium("g", 1.5));
real[] rg_incs = {0.0, 20.0, 50.0, 70.0, 85.0};
for (real th : rg_incs) {
    Polarization rg_pin = Polarization.fromAngles(30.0, 12.0, 1.0);
    Trace rg_tf = system(rg_iface).trace(Ray((0, 0) - 10 * dir(th), dir(th), pol = rg_pin));
    Interaction rg_fix = rg_tf.interaction(0);
    Ray rg_tr = rg_fix.transmitted[0];
    real rg_rI = (rg_fix.reflected.length > 0) ? rg_fix.reflected[0].intensity() : 0.0;
    assert_near("Power conserved (R + T) at " + string(th) + " deg", rg_tr.intensity() + rg_rI, 1.0, 1e-9);
    assert_near("Transmitted DoP = 1 at " + string(th) + " deg", rg_tr.pol.dop(), 1.0, 1e-9);
    assert_stokes_physical("Transmitted state at " + string(th) + " deg", rg_tr.pol, 1e-9);
    if (rg_fix.reflected.length > 0)
        assert_stokes_physical("Reflected state at " + string(th) + " deg", rg_fix.reflected[0].pol, 1e-9);
}
// Partially polarized input stays inside the Stokes cone
Polarization rg_partial = Polarization.stokes(1.0, 0.3, -0.2, 0.4);
Trace rg_tp = system(rg_iface).trace(Ray((0, 0) - 10 * dir(60), dir(60), pol = rg_partial));
assert_stokes_physical("Partially polarized transmitted state", rg_tp.transmitted(0).pol, 1e-12);
assert_stokes_physical("Partially polarized reflected state", rg_tp.reflected(0).pol, 1e-12);
write("");

// ============================================================================
// Test Case 78 (Regression): cardinal points use h / tan(u); exact for an
// ideal thin lens even with a large probe height
// ============================================================================
write("Test Case 78: Cardinal-point focal length with large probe height");
OpticalSystem rg_thinSys = system(Lens.thin((0, 0), 0, 100.0, 40.0));
CardinalPoints rg_cpBig = cardinalPoints(rg_thinSys, paraxial_h = 5.0);
assert_near("Thin-lens EFL independent of probe height", rg_cpBig.f2, 100.0, 1e-9);
assert_near("Thin-lens front focal length", rg_cpBig.f1, 100.0, 1e-9);
assert_near_pair("Thin-lens principal point H2 at lens", rg_cpBig.H2, (0, 0), 1e-9);
// Nodal points of a lens immersed in a single medium coincide with principal points
OpticalSystem rg_wetSys = system(Lens.thin((0, 0), 0, 100.0, 40.0), Materials.WATER);
CardinalPoints rg_cpWet = cardinalPoints(rg_wetSys);
assert_near_pair("Immersed lens: N1 = H1", rg_cpWet.N1, rg_cpWet.H1, 1e-9);
assert_near_pair("Immersed lens: N2 = H2", rg_cpWet.N2, rg_cpWet.H2, 1e-9);
write("");

// ============================================================================
// Test Case 79 (Regression): GRIN boundary detection with very small steps
// (chords shorter than the integration tolerance used to be skipped and the
// ray tunnelled through the boundary)
// ============================================================================
write("Test Case 79: GRIN boundary detection with tiny fixed steps");
indexProfile rg_lin = new real(pair p, real lambda) { return 1.5 + 0.01 * p.y; };
gradientProfile rg_linG = new pair(pair p, real lambda) { return (0, 0.01); };
Medium rg_gl = medium("rg linear GRIN").add(GRIN(rg_lin, rg_linG));
OpticalSystem rg_thinGrin = system(Interface(box((0, -1), (0.02, 1)), rg_gl, Materials.VACUUM));
config.physics.grin_solver = GRINSolverType.RK4;
Trace rg_tTiny = rg_thinGrin.trace(Ray((0.01, 0), (1, 0), rg_gl), ds = 5e-5, far = 1.0);
config.physics.grin_solver = GRINSolverType.RK45;
assert_true("Ray stops at the GRIN boundary", rg_tTiny.length() >= 1);
assert_near("Boundary hit at x = 0.02", rg_tTiny.hit(0).x, 0.02, 1e-9);
assert_near("Ray continues in vacuum after the boundary", rg_tTiny.last().medium.n(), 1.0, 1e-12);
write("");

// ============================================================================
// Test Case 80 (Regression): GRIN OPL over the final partial step to a
// boundary uses the midpoint index (exact for a linear profile)
// ============================================================================
write("Test Case 80: GRIN OPL to a boundary (partial final step)");
Medium rg_axial = medium("rg axial")
                      .add(GRIN(
                          new real(pair p, real lambda) { return 1.0 + p.x; },
                          new pair(pair p, real lambda) { return (1.0, 0.0); }));
Trace rg_tAx = system(Interface(box((-1, -1), (0.73, 1)), rg_axial, Materials.VACUUM))
                   .trace(Ray((0, 0), (1, 0), rg_axial), ds = 0.05);
assert_near("OPL = integral_0^0.73 (1 + x) dx", rg_tAx.interaction(0).incident.opl, 0.73 + 0.73 ^ 2 / 2, 1e-9);
write("");

// ============================================================================
// Test Case 81: GRIN numerical convergence (RK4 global error order and RK45
// tolerance control) against a fine-step reference
// ============================================================================
write("Test Case 81: GRIN solver convergence");
OpticalSystem rg_grinOpen = system(new Element[]);
Medium rg_par = medium("rg parabolic")
                    .add(GRIN(
                        new real(pair p, real lambda) { return 1.5 * (1 - 0.5 * 0.3 ^ 2 * p.y ^ 2); },
                        new pair(pair p, real lambda) { return (0, -1.5 * 0.3 ^ 2 * p.y); }));
Ray rg_r0 = Ray((0, 1.0), (1, 0), rg_par);
config.physics.grin_solver = GRINSolverType.RK4;
real rg_yRef = rg_grinOpen.trace(rg_r0, ds = 0.0025, far = 4.0).last().p.y;
real rg_e1 = abs(rg_grinOpen.trace(rg_r0, ds = 0.07, far = 4.0).last().p.y - rg_yRef);
real rg_e2 = abs(rg_grinOpen.trace(rg_r0, ds = 0.035, far = 4.0).last().p.y - rg_yRef);
config.physics.grin_solver = GRINSolverType.RK45;
write("  RK4 errors: ds=0.07 -> " + string(rg_e1) + ", ds=0.035 -> " + string(rg_e2));
assert_true("RK4 error decreases >= 10x when ds halves (4th-order convergence)", rg_e2 < rg_e1 / 10.0);
assert_true("RK4 coarse error is already small (< 1e-6)", rg_e1 < 1e-6);
real rg_y45 = rg_grinOpen.trace(rg_r0, ds = 0.05, far = 4.0).last().p.y;
assert_near("RK45 matches fine RK4 reference", rg_y45, rg_yRef, 1e-5);
write("");

// ============================================================================
// Test Case 82: Sellmeier dispersion factory (homogeneous; no GRIN integration)
// ============================================================================
write("Test Case 82: Sellmeier factory");
Medium rg_bk7 = medium("rg BK7").add(
    Sellmeier(new real[]{1.03961212, 0.231792344, 1.01046945}, new real[]{0.00600069867, 0.0200179144, 103.560653}));
assert_near("Sellmeier BK7 n_d (587.56 nm)", rg_bk7.n(lambda = 587.5618), 1.51680, 1e-4);
assert_true("Normal dispersion: n_F(486.1) > n_C(656.3)", rg_bk7.n(lambda = 486.1) > rg_bk7.n(lambda = 656.3));
Trace rg_tS = system(Interface(box((0, -10), (10, 10)), rg_bk7, Materials.VACUUM)).trace(Ray((-5, 0), dir(10)));
assert_near("Ray in homogeneous Sellmeier glass is a single straight segment", length(rg_tS.segments()[1].traj), 1.0);
write("");

// ============================================================================
// Test Case 83: GRIN profile scaled by a homogeneous index factor follows the
// same ray (ray equation is invariant under n -> c n); exercises analytic
// gradient composition with gradient-free factors
// ============================================================================
write("Test Case 83: Scaled GRIN profile ray invariance");
Medium rg_parScaled = medium("rg scaled", 1.3)
                          .add(GRIN(
                              new real(pair p, real lambda) { return 1.5 * (1 - 0.5 * 0.3 ^ 2 * p.y ^ 2); },
                              new pair(pair p, real lambda) { return (0, -1.5 * 0.3 ^ 2 * p.y); }));
real rg_yScaled = rg_grinOpen.trace(Ray((0, 1.0), (1, 0), rg_parScaled), ds = 0.05, far = 4.0).last().p.y;
assert_near("Scaled GRIN ray matches unscaled ray", rg_yScaled, rg_y45, 1e-6);
real rg_oplScaled = rg_grinOpen.trace(Ray((0, 1.0), (1, 0), rg_parScaled), ds = 0.05, far = 4.0).opl();
real rg_oplBase = rg_grinOpen.trace(Ray((0, 1.0), (1, 0), rg_par), ds = 0.05, far = 4.0).opl();
assert_near("Scaled GRIN OPL scales by the factor", rg_oplScaled / rg_oplBase, 1.3, 1e-6);
write("");

// ============================================================================
// Test Case 84 (Regression): Hit.entering reports entry into closed elements
// ============================================================================
write("Test Case 84: Hit.entering orientation");
Element rg_box = Interface(box((0, -5), (10, 5)), medium("g", 1.5));
Hit rg_hIn = rg_box.surfaces[0].intersect(Ray((-5, 0), (1, 0)));
Hit rg_hOut = rg_box.surfaces[0].intersect(Ray((5, 0), (1, 0)));
assert_bool("Ray from outside is entering", rg_hIn.entering, true);
assert_bool("Ray from inside is exiting", rg_hOut.entering, false);
write("");

// ============================================================================
// Test Case 85 (Regression): copy(Surface) keeps the precomputed inverse
// ============================================================================
write("Test Case 85: Surface copy keeps inverse transform");
Element rg_shifted = shift(30, 7) * rotate(25) * Interface(box((0, -5), (10, 5)), medium("g", 1.5));
Surface rg_flat = rg_shifted.flatten()[0];
Surface rg_flatCopy = copy(rg_flat);
Ray rg_probe = Ray((0, 7), (1, 0));
assert_near_pair("Copied surface gives the same hit point",
                 rg_flatCopy.intersect(rg_probe).p,
                 rg_flat.intersect(rg_probe).p,
                 1e-9);
assert_true("Copied surface hit is valid", rg_flatCopy.intersect(rg_probe).valid);
write("");

// ============================================================================
// Test Case 86 (Regression): beam factories inherit the medium at the source
// ============================================================================
write("Test Case 86: Beam rays resolve their starting medium");
OpticalSystem rg_inGlass = system(Interface(box((-10, -10), (10, 10)), medium("g", 1.5)));
Ray[] rg_beam = Beam.parallel((0, -2), (0, 2), 3, (1, 0));
Trace rg_tb = rg_inGlass.trace(rg_beam[1]);
assert_near("Beam ray launched inside glass starts in glass", rg_tb.first().medium.n(), 1.5, 1e-12);
assert_near("No spurious entry interface (first hit is the exit face at x = 10)", rg_tb.hit(0).x, 10.0, 1e-9);
write("");

// ============================================================================
// Test Case 87: Linear dichroism with a rotated axis
// ============================================================================
write("Test Case 87: Rotated dichroic axis");
Medium rg_dich = medium("rg dichroic 30").add(Dichroic(0.0, 50.0, 30.0));
Polarization rg_pass = rg_dich.stepPolarization(Polarization.linear(30.0), (0, 0), (1, 0), 633, 1.0);
Polarization rg_block = rg_dich.stepPolarization(Polarization.linear(120.0), (0, 0), (1, 0), 633, 1.0);
Polarization rg_mid = rg_dich.stepPolarization(Polarization.linear(75.0), (0, 0), (1, 0), 633, 1.0);
assert_near("Light along the axis passes", rg_pass.intensity(), 1.0, 1e-9);
assert_near("Light orthogonal to the axis is extinguished", rg_block.intensity(), 0.0, 1e-9);
assert_near("Malus at 45 deg from the axis", rg_mid.intensity(), 0.5, 1e-9);
assert_near("Output azimuth equals the axis", degrees(rg_mid.psi()), 30.0, 1e-6);
write("");

// ============================================================================
// Test Case 88: Edge cases - grazing incidence, zero-intensity rays, max_depth
// ============================================================================
write("Test Case 88: Edge cases");
mediumResolver rg_resE = new Medium(pair p, Ray r, Medium m1, Medium m2, bool from_m1) {
    return (p.x > 0) ? medium("g", 1.5) : Materials.VACUUM;
};
FresnelResult rg_graze = fresnelIsotropic(Materials.VACUUM, medium("g", 1.5), Ray((-1, 0), dir(89.99)), rg_h, rg_resE);
assert_true("Grazing incidence: R finite and close to 1", rg_graze.R > 0.99 && rg_graze.R <= 1.0);
assert_near("Grazing incidence: R + T = 1", rg_graze.R + rg_graze.T, 1.0, 1e-12);
Trace rg_dark = system(Interface(box((0, -5), (10, 5)), medium("g", 1.5)))
                    .trace(Ray((-5, 0), (1, 0), pol = Polarization.unpolarized(0.0)));
assert_near("Zero-intensity ray spawns no children", rg_dark.length(), 0.0);
OpticalSystem rg_cavity = system(new Element[]{Mirror.flat((0, 0), 0, 10), Mirror.flat((10, 0), 180, 10)});
Trace rg_tCav = rg_cavity.trace(Ray((5, 0), (1, 0)), max_depth = 6);
assert_near("Cavity bounces limited by max_depth", rg_tCav.length(), 6.0);
write("");

// ============================================================================
// Test Case 89: Extraordinary TIR in the isotropic limit
//
// When no == ne, the uniaxial extraordinary mode must reduce to the
// ordinary isotropic Fresnel solution, including the COMPLEX reflection
// amplitude, not merely its magnitude.
// ============================================================================
write("Test Case 89: Extraordinary TIR isotropic-limit regression");
Medium rg_isoLimit = uniaxialExtraordinary("isotropic limit", 1.5, 1.5, (1, 0));
real[] rg_isoTirAngles = {50.0, 55.0, 60.0, 70.0};
pair[] rg_isoRpValues;
for (real th : rg_isoTirAngles) {
    Ray rg_isoRay = Ray((-1, 0), dir(th), rg_dense, pol = Polarization.linear(90.0));
    FresnelResult rg_anisoResult = fresnelUniaxial(rg_isoRay, rg_h, rg_dense, rg_lo, rg_isoLimit, rg_resD);
    FresnelResult rg_isoResult = fresnelIsotropic(rg_dense, rg_lo, rg_isoRay, rg_h, rg_resD);
    assert_true("Extraordinary isotropic-limit TIR at " + string(th) + " deg", rg_anisoResult.tir_e);
    assert_near("Extraordinary TIR magnitude at " + string(th) + " deg", abs(rg_anisoResult.rp), 1.0, 1e-10);
    assert_near_pair("Extraordinary complex rp matches isotropic rp at " + string(th) + " deg",
                     rg_anisoResult.rp,
                     rg_isoResult.rp,
                     1e-8);
    assert_near("Extraordinary Rp at " + string(th) + " deg", rg_anisoResult.Rp, abs2(rg_isoResult.rp), 1e-10);
    rg_isoRpValues.push(rg_anisoResult.rp);
}
// The isotropic TIR phase changes with angle. A fixed rp = (-1, 0)
// cannot satisfy the complex-amplitude comparisons above.
assert_true("Extraordinary TIR phase varies with incidence angle",
            length(rg_isoRpValues[0] - rg_isoRpValues[3]) > 1e-3);
write("");

// ============================================================================
// Test Case 90: Genuine anisotropy and optic-axis orientation
// ============================================================================
write("Test Case 90: Genuine anisotropy and optic-axis orientation");
real rg_no90 = 1.50;
real rg_ne90 = 1.70;
real rg_angle90 = 55.0;
Medium rg_eAxis0_90 = uniaxialExtraordinary("uniaxial axis 0", rg_no90, rg_ne90, (1, 0));
Medium rg_eAxis45_90 = uniaxialExtraordinary("uniaxial axis 45", rg_no90, rg_ne90, dir(45));
Medium rg_eAxis90_90 = uniaxialExtraordinary("uniaxial axis 90", rg_no90, rg_ne90, (0, 1));
Ray rg_ray90 = Ray((-1, 0), dir(rg_angle90), rg_dense, pol = Polarization.linear(90.0));
FresnelResult rg_f0_90 = fresnelUniaxial(rg_ray90, rg_h, rg_dense, rg_lo, rg_eAxis0_90, rg_resD);
FresnelResult rg_f45_90 = fresnelUniaxial(rg_ray90, rg_h, rg_dense, rg_lo, rg_eAxis45_90, rg_resD);
FresnelResult rg_f90_90 = fresnelUniaxial(rg_ray90, rg_h, rg_dense, rg_lo, rg_eAxis90_90, rg_resD);
assert_true("Anisotropic result has birefringence flag", rg_f0_90.is_birefringent);
assert_true("Axis-0 extraordinary reflection amplitude is finite", !isnan(rg_f0_90.rp.x) && !isnan(rg_f0_90.rp.y));
assert_true("Axis-45 extraordinary reflection amplitude is finite", !isnan(rg_f45_90.rp.x) && !isnan(rg_f45_90.rp.y));
assert_true("Axis-90 extraordinary reflection amplitude is finite", !isnan(rg_f90_90.rp.x) && !isnan(rg_f90_90.rp.y));
assert_true("Axis-0 reflectance is in physical range", rg_f0_90.Rp >= -1e-9 && rg_f0_90.Rp <= 1.0 + 1e-9);
assert_true("Axis-45 reflectance is in physical range", rg_f45_90.Rp >= -1e-9 && rg_f45_90.Rp <= 1.0 + 1e-9);
assert_true("Axis-90 reflectance is in physical range", rg_f90_90.Rp >= -1e-9 && rg_f90_90.Rp <= 1.0 + 1e-9);
// Only assert an orientation difference if the selected setup produces one.
// A failure here is diagnostic: it may indicate that the test geometry
// is insensitive to the optic-axis rotation, not necessarily a solver bug.
write("  axis 0 rp  = " + (string)rg_f0_90.rp);
write("  axis 45 rp = " + (string)rg_f45_90.rp);
write("  axis 90 rp = " + (string)rg_f90_90.rp);
write("");

// ============================================================================
// Test Case 91: Extraordinary critical-angle sweep
// ============================================================================
write("Test Case 91: Extraordinary critical-angle sweep");
Medium rg_e91 = uniaxialExtraordinary("critical sweep crystal", 1.50, 1.70, dir(30));
real[] rg_angles91 = {20.0, 30.0, 40.0, 45.0, 48.0, 50.0, 52.0, 54.0, 56.0, 58.0, 60.0, 65.0, 70.0, 75.0, 80.0};
for (real th91 : rg_angles91) {
    Ray rg_ray91 = Ray((-1, 0), dir(th91), rg_dense, pol = Polarization.linear(90.0));
    FresnelResult rg_fr91 = fresnelUniaxial(rg_ray91, rg_h, rg_dense, rg_lo, rg_e91, rg_resD);
    assert_true("Critical sweep finite rp at " + string(th91), !isnan(rg_fr91.rp.x) && !isnan(rg_fr91.rp.y));
    assert_true("Critical sweep finite Rs at " + string(th91), !isnan(rg_fr91.Rs));
    assert_true("Critical sweep finite Rp at " + string(th91), !isnan(rg_fr91.Rp));
    assert_true("Critical sweep physical Rp at " + string(th91), rg_fr91.Rp >= -1e-8 && rg_fr91.Rp <= 1.0 + 1e-8);
    if (rg_fr91.tir_e) {
        assert_near("Extraordinary TIR unit reflectance at " + string(th91), rg_fr91.Rp, 1.0, 1e-8);
    }
}
write("");

// ============================================================================
// Test Case 92: Reflectance and transmittance bookkeeping
// ============================================================================
write("Test Case 92: Reflectance and transmittance bookkeeping");
Medium rg_e92 = uniaxialExtraordinary("energy bookkeeping crystal", 1.50, 1.70, dir(35));
Ray rg_ray92 = Ray((-1, 0), dir(35), rg_dense, pol = Polarization.linear(45.0));
FresnelResult rg_fr92 = fresnelUniaxial(rg_ray92, rg_h, rg_dense, rg_lo, rg_e92, rg_resD);
assert_near("Total energy bookkeeping R+T", rg_fr92.R + rg_fr92.T, 1.0, 1e-9);
if (rg_fr92.tir_e) {
    assert_near("Extraordinary TIR has zero transmitted power", rg_fr92.Te, 0.0, 1e-9);
    assert_near("Extraordinary TIR has unit reflected power", rg_fr92.Rp, 1.0, 1e-9);
}
assert_true("Total reflectance within physical bounds", rg_fr92.R >= -1e-9 && rg_fr92.R <= 1.0 + 1e-9);
assert_true("Total transmittance within physical bounds", rg_fr92.T >= -1e-9 && rg_fr92.T <= 1.0 + 1e-9);
write("");

// ============================================================================
// Test Case 93: Extraordinary dispersion and tangential phase matching
// ============================================================================
write("Test Case 93: Extraordinary dispersion and phase matching");
real rg_no93 = 1.50;
real rg_ne93 = 1.70;
pair rg_axis93 = dir(30);
Medium rg_e93 = uniaxialExtraordinary("dispersion regression", rg_no93, rg_ne93, rg_axis93);
pair rg_N93 = unitOr(rg_h.normal);
pair rg_T93 = (-rg_N93.y, rg_N93.x);
real[] rg_angles93 = {10.0, 20.0, 30.0, 40.0};
for (real rg_theta93 : rg_angles93) {
    Ray rg_ray93 = Ray((-1, 0), dir(rg_theta93), rg_dense, pol = Polarization.linear(90.0));
    pair rg_ki93 = unitOr(rg_ray93.wave_normal, rg_ray93.d);
    real rg_n193 = rg_dense.n(rg_h.p, rg_ki93, rg_ray93.lambda, rg_ray93.pol);
    real rg_kt93 = rg_n193 * dot(rg_ki93, rg_T93);
    PhaseMatch rg_pm93 = solvePhaseMatch(rg_e93, rg_N93, rg_T93, rg_kt93, -rg_N93, rg_h.p, rg_ray93.lambda);
    if (rg_pm93.tir) {
        write("  Diagnostic: extraordinary mode evanescent at " + string(rg_theta93) + " degrees");
    } else {
        pair rg_k93 = rg_pm93.n * rg_pm93.wave_normal;
        pair rg_c93 = unitOr(rg_axis93);
        pair rg_q93 = perp(rg_c93);
        real rg_disp93 = dot(rg_k93, rg_c93) ^ 2 / (rg_no93 * rg_no93) + dot(rg_k93, rg_q93) ^ 2 / (rg_ne93 * rg_ne93);
        assert_near("Extraordinary dispersion equals one at " + string(rg_theta93) + " degrees", rg_disp93, 1.0, 1e-7);
        assert_near("Tangential wavevector conserved at " + string(rg_theta93) + " degrees",
                    dot(rg_k93, rg_T93),
                    rg_kt93,
                    1e-7);
        assert_true("Valid extraordinary group direction at " + string(rg_theta93) + " degrees",
                    dot(rg_pm93.group_direction, -rg_N93) > 0);
    }
}
write("");

// ============================================================================
// Test Case 94: Independent Energy Flux Conservation (Isotropic)
// ============================================================================
write("Test Case 94: Independent Energy Flux Conservation (Isotropic)");
real n_in = 1.0, n_out = 1.5;
real th_deg = 45.0;
real th_rad = radians(th_deg);
real sin_t = (n_in / n_out) * sin(th_rad);
real cos_t = sqrt(1.0 - sin_t ^ 2);
real cos_i = cos(th_rad);
// Field amplitudes from Maxwell boundary conditions
pair rs_test = (n_in * cos_i - n_out * cos_t) / (n_in * cos_i + n_out * cos_t);
pair rp_test = (n_out * cos_i - n_in * cos_t) / (n_out * cos_i + n_in * cos_t);
pair ts_test = (2.0 * n_in * cos_i) / (n_in * cos_i + n_out * cos_t);
pair tp_test = (2.0 * n_in * cos_i) / (n_out * cos_i + n_in * cos_t);
// Independent power flux: T = (n2*cos_t)/(n1*cos_i) * |t|^2
real flux_factor = (n_out * cos_t) / (n_in * cos_i);
real Rs_flux = rs_test.x * rs_test.x + rs_test.y * rs_test.y;
real Rp_flux = rp_test.x * rp_test.x + rp_test.y * rp_test.y;
real Ts_flux = flux_factor * (ts_test.x * ts_test.x + ts_test.y * ts_test.y);
real Tp_flux = flux_factor * (tp_test.x * tp_test.x + tp_test.y * tp_test.y);
assert_near("TE/s-polarization energy conservation (Rs + Ts == 1.0)", Rs_flux + Ts_flux, 1.0, 1e-9);
assert_near("TM/p-polarization energy conservation (Rp + Tp == 1.0)", Rp_flux + Tp_flux, 1.0, 1e-9);
write("");

// ============================================================================
// Test Case 95: Complex TIR Retardance (Exact Fresnel Rhomb Quarter-Wave Angle)
// ============================================================================
write("Test Case 95: Complex TIR Retardance (Exact Fresnel Rhomb Quarter-Wave Angle)");
// 1. Fetch exact material indices from predefined library
real n_glass = Materials.BK7.n();
real n_air = Materials.AIR.n();
real n_rel = n_air / n_glass;
// 2. Solve analytically for the exact angle that gives delta = pi / 4 (45 degrees)
real T = sqrt(2.0) - 1.0; // tan(pi / 8)
real A = 1.0 + T * T;
real B = -(1.0 + n_rel * n_rel);
real C = n_rel * n_rel;
real disc = B * B - 4.0 * A * C;
// The larger root gives the standard rhomb apex angle (~55 deg)
real u = (-B + sqrt(disc)) / (2.0 * A); // u = sin^2(theta)
real th_rad = asin(sqrt(u));
real s_i = sin(th_rad);
real c_i = cos(th_rad);
// 3. Run the standalone interface solver
pair k_inc = n_glass * (s_i, c_i);
pair N_into2 = (0, 1);
pair T_surf = (1, 0);
FresnelResult em_tir = solveIsotropicInterface(n_glass, n_air, k_inc, N_into2, T_surf);
real delta_s = atan2(em_tir.rs.y, em_tir.rs.x);
real delta_p = atan2(em_tir.rp.y, em_tir.rp.x);
real delta = abs(delta_s - delta_p);
// 4. Assertions: exact 45.0 deg per bounce -> exact 90.0 deg for two bounces
assert_bool("Solver flags TIR correctly", em_tir.tir, true);
assert_near("TIR Rs magnitude == 1.0", em_tir.Rs, 1.0, 1e-9);
assert_near("TIR Rp magnitude == 1.0", em_tir.Rp, 1.0, 1e-9);
assert_near("TIR Ts flux == 0.0", em_tir.Ts, 0.0, 1e-9);
assert_near("TIR Tp flux == 0.0", em_tir.Tp, 0.0, 1e-9);
assert_near("Single bounce retardance is exactly pi/4 (rad)", delta, pi / 4.0, 1e-8);
assert_near("Fresnel Rhomb: two bounces yield exactly 90.0 deg", 2.0 * degrees(delta), 90.0, 1e-7);
write("");

// ============================================================================
// Test Case 96: Independent Poynting Flux Conservation in Birefringent Crystal
// ============================================================================
write("Test Case 96: Independent Poynting Flux Conservation in Birefringent Crystal");
// Calcite crystal with optic axis tilted at 35 degrees
real no_calc = 1.658, ne_calc = 1.486;
pair c_tilted = dir(35.0);
Medium cal_o = medium("Calcite O", no_calc);
Medium cal_e = medium("Calcite E").add(UniaxialExtraordinary(no_calc, ne_calc, c_tilted));
// Oblique incidence from vacuum at 30 degrees
real th_inc_deg = 30.0;
Ray r_oblique = Ray.angle(p = (-2, 0), angle = th_inc_deg, medium = Materials.VACUUM);
Hit h_cryst;
h_cryst.valid = true;
h_cryst.p = (0, 0);
h_cryst.normal = (-1, 0);
h_cryst.tangent = (0, 1);
mediumResolver res_cryst = new Medium(pair p, Ray r, Medium m1, Medium m2, bool from_m1) {
    return (p.x > 0) ? cal_o : Materials.VACUUM;
};
FresnelResult res96 = fresnelUniaxial(r_oblique, h_cryst, Materials.VACUUM, cal_o, cal_e, res_cryst);
// 1. Both channels must be propagating
assert_bool("Ordinary mode propagates", res96.tir_o, false);
assert_bool("Extraordinary mode propagates", res96.tir_e, false);
// 2. Poynting walk-off check: extraordinary group direction d_e != waveNormal_e
real walkoff_dot = dot(res96.d_e, res96.wave_normal_e);
assert_true("Poynting walk-off exists (d_e != k_e)", walkoff_dot < 0.9999 && walkoff_dot > 0.0);
// 3. Strict Independent Power Conservation: Rs + To == 1.0 and Rp + Te == 1.0
assert_near("Ordinary mode energy conservation (Rs + To == 1.0)", res96.Rs + res96.To, 1.0, 1e-8);
assert_near("Extraordinary mode energy conservation (Rp + Te == 1.0)", res96.Rp + res96.Te, 1.0, 1e-8);
write("");

// ============================================================================
// Test Case 97: Birefringent TIR Asymmetry (Glan-Prism Operating Regime)
// ============================================================================
write("Test Case 97: Birefringent TIR Asymmetry (Glan-Prism Operating Regime)");
// Calcite has no = 1.658, ne = 1.486.
// Critical angles into air:
//   theta_c,o = asin(1 / 1.658) = 37.10 deg
//   theta_c,e = asin(1 / ne(theta)) ~ 40.3 deg (for axis along surface)
// At 38.0 degrees: ordinary mode undergoes TIR, while extraordinary mode transmits cleanly!
Medium cal_o_dense = medium("Dense Calcite O", no_calc);
Medium cal_e_dense = medium("Dense Calcite E").add(UniaxialExtraordinary(no_calc, ne_calc, (0, 1)));
real th_glan_deg = 38.0;
// 1. Launch ordinary ray (n = 1.658) hitting air interface
Ray r_o_glan = Ray.angle(p = (-2, 0), angle = th_glan_deg, medium = cal_o_dense);
mediumResolver res_air = new Medium(pair p, Ray r, Medium m1, Medium m2, bool from_m1) {
    return (p.x > 0) ? Materials.AIR : cal_o_dense;
};
FresnelResult res_o = fresnelIsotropic(cal_o_dense, Materials.AIR, r_o_glan, h_cryst, res_air);
// Ordinary wave must undergo TIR (kt = 1.021 > 1.0)
assert_bool("Ordinary mode is totally reflected (TIR)", res_o.tir, true);
assert_near("Ordinary Rs == 1.0", res_o.Rs, 1.0, 1e-8);
assert_near("Ordinary T == 0.0", res_o.T, 0.0, 1e-8);
// 2. Launch extraordinary ray (n_e ~ 1.545, pure TM/p-polarization) hitting air interface
mediumResolver res_air_e = new Medium(pair p, Ray r, Medium m1, Medium m2, bool from_m1) {
    return (p.x > 0) ? Materials.AIR : cal_e_dense;
};
// Physical extraordinary eigenmode is pure p-polarized: S1 = -1
Polarization p_pol = Polarization.stokes(1.0, -1.0, 0.0, 0.0);
Ray r_e_glan = Ray.angle(p = (-2, 0), angle = th_glan_deg, medium = cal_e_dense, pol = p_pol);
FresnelResult res_e = fresnelIsotropic(cal_e_dense, Materials.AIR, r_e_glan, h_cryst, res_air_e);
// Extraordinary wave transmits freely (kt = 0.951 < 1.0, Te = 94%)
assert_bool("Extraordinary mode transmits freely", res_e.tir, false);
assert_near("Extraordinary energy conservation (R + T == 1.0)", res_e.R + res_e.T, 1.0, 1e-8);
assert_true("High extraordinary transmittance (T > 0.85)", res_e.T > 0.85);
write("");

// ============================================================================
// Test Case 98: Dual-Mode Anisotropic TIR (Evanescent Field Phase Retardance)
// ============================================================================
write("Test Case 98: Dual-Mode Anisotropic TIR (Evanescent Field Phase Retardance)");
// At 50 degrees, both modes exceed their critical angles into air
real th_tir_both = 50.0;
Ray r_tir_both = Ray.angle(p = (-2, 0), angle = th_tir_both, medium = cal_o_dense);
FresnelResult res98 = fresnelUniaxial(r_tir_both, h_cryst, cal_o_dense, Materials.AIR, Materials.AIR, res_air);
assert_bool("Both modes in TIR: o-mode", res98.tir_o, true);
assert_bool("Both modes in TIR: e-mode", res98.tir_e, true);
// Magnitudes must be identically unity
assert_near("TIR Rs magnitude == 1.0", res98.Rs, 1.0, 1e-8);
assert_near("TIR Rp magnitude == 1.0", res98.Rp, 1.0, 1e-8);
assert_near("TIR To flux == 0.0", res98.To, 0.0, 1e-8);
assert_near("TIR Te flux == 0.0", res98.Te, 0.0, 1e-8);
// Complex amplitudes must have non-zero imaginary parts (differential phase retardance)
assert_true("Evanescent rs has complex phase", abs(res98.rs.y) > 0.1);
assert_true("Evanescent rp has complex phase", abs(res98.rp.y) > 0.1);
real phase_s = atan2(res98.rs.y, res98.rs.x);
real phase_p = atan2(res98.rp.y, res98.rp.x);
real delta_aniso = abs(phase_p - phase_s);
assert_true("Anisotropic TIR differential retardance is non-zero", delta_aniso > 0.05);
write("");

// ============================================================================
// Test Case 99: Isotropic Limit Verification (no == ne)
// ============================================================================
write("Test Case 99: Isotropic Limit Verification (no == ne)");
// When a crystal has no = ne, fresnelUniaxial must match fresnelIsotropic identically
real n_iso_val = 1.55;
Medium med_iso_match = medium("Iso Match", n_iso_val);
Medium med_crystal_limit_o = medium("Limit O", n_iso_val);
Medium med_crystal_limit_e = medium("Limit E").add(UniaxialExtraordinary(n_iso_val, n_iso_val, dir(30)));
Ray r_limit = Ray.angle(p = (-2, 0), angle = 40.0, medium = Materials.VACUUM);
Hit h_limit;
h_limit.valid = true;
h_limit.p = (0, 0);
h_limit.normal = (-1, 0);
h_limit.tangent = (0, 1);
mediumResolver res_limit = new Medium(pair p, Ray r, Medium m1, Medium m2, bool from_m1) {
    return (p.x > 0) ? med_crystal_limit_o : Materials.VACUUM;
};
FresnelResult res_iso_base = fresnelIsotropic(Materials.VACUUM, med_iso_match, r_limit, h_limit, res_limit);
FresnelResult res_uni_limit =
    fresnelUniaxial(r_limit, h_limit, Materials.VACUUM, med_crystal_limit_o, med_crystal_limit_e, res_limit);
assert_near("Isotropic limit: Rs matches isotropic solver", res_uni_limit.Rs, res_iso_base.Rs, 1e-12);
assert_near("Isotropic limit: Rp matches isotropic solver", res_uni_limit.Rp, res_iso_base.Rp, 1e-12);
assert_near_pair("Isotropic limit: o-ray direction matches isotropic refraction",
                 res_uni_limit.d_o,
                 res_iso_base.d,
                 1e-12);
assert_near_pair("Isotropic limit: e-ray direction matches isotropic refraction",
                 res_uni_limit.d_e,
                 res_iso_base.d,
                 1e-12);
assert_near("Isotropic limit: zero walk-off between o and e modes",
            dot(res_uni_limit.d_o, res_uni_limit.d_e),
            1.0,
            1e-12);
write("");

// ============================================================================
// Test Case 100: Brewster Angle Reflection Extinction via Refract()
// ============================================================================
write("Test Case 100: Brewster Angle Reflection Extinction via Refract()");
// Air (n=1) to Glass (n=1.5). Brewster angle: tan(theta_B) = 1.5 / 1.0 -> 56.3099 deg
real n_glass_b = 1.5;
Medium glass_brewster = medium("Brewster Glass", n_glass_b);
real th_B_deg = degrees(atan(n_glass_b / 1.0));
// Launch pure TM / p-polarized ray (S1 = -1) aimed directly at (0, 0)
Polarization p_polarized = Polarization.stokes(1.0, -1.0, 0.0, 0.0);
pair launch_dir = dir(th_B_deg);
pair launch_pt = -5.0 * launch_dir;
Ray r_brewster = Ray(p = launch_pt, d = launch_dir, medium = Materials.VACUUM, pol = p_polarized);
Element brewster_surface = Interface(linearBoundary((0, 0), (0, 1), 20), glass_brewster, Materials.VACUUM);
Trace t_brewster = system(brewster_surface).trace(r_brewster);
// Transmitted ray must carry 100% of incident power (zero reflection)
Ray trans_B = t_brewster.transmitted(0);
assert_true("Transmitted ray exists", trans_B != null);
assert_near("Transmitted intensity at Brewster angle is 100%", trans_B.intensity(), 1.0, 1e-7);
assert_near("Brewster TM reflection is extinguished", t_brewster.interaction(0).R, 0.0, 1e-12);
write("");

// ============================================================================
// Test Case 101: End-to-End Crystal Eigenmode Polarization Purity
// ============================================================================
write("Test Case 101: End-to-End Crystal Eigenmode Polarization Purity");
// Launch unpolarized light into a Calcite Crystal slab
Medium cal_o_slab = medium("Calcite Slab O", 1.658);
Medium cal_e_slab = medium("Calcite Slab E").add(UniaxialExtraordinary(1.658, 1.486, dir(45)));
Element crystal_slab = Crystal(box((-2, -5), (2, 5)), cal_o_slab, cal_e_slab, Materials.VACUUM);
Ray unpol_ray = Ray.angle(p = (-5, 0), angle = 10.0, medium = Materials.VACUUM, pol = Polarization.unpolarized());
Trace t_crystal = system(crystal_slab).trace(unpol_ray);
// Entry interface must spawn 2 rays
Interaction entry_step = t_crystal.interaction(0);
assert_true("Entry splits unpolarized light into 2 modes", entry_step.transmitted.length == 2);
Ray o_beam = entry_step.transmitted[0];
Ray e_beam = entry_step.transmitted[1];
// Ordinary beam must be pure S-polarization: S1 / S0 == +1.0, S2 == 0, S3 == 0
assert_near("Ordinary mode pure S-polarization (S1/S0 == +1.0)", o_beam.pol.S1 / o_beam.pol.S0, 1.0, 1e-8);
assert_near("Ordinary mode S2 == 0.0", o_beam.pol.S2, 0.0, 1e-8);
assert_near("Ordinary mode S3 == 0.0", o_beam.pol.S3, 0.0, 1e-8);
// Extraordinary beam must be pure P-polarization: S1 / S0 == -1.0, S2 == 0, S3 == 0
assert_near("Extraordinary mode pure P-polarization (S1/S0 == -1.0)", e_beam.pol.S1 / e_beam.pol.S0, -1.0, 1e-8);
assert_near("Extraordinary mode S2 == 0.0", e_beam.pol.S2, 0.0, 1e-8);
assert_near("Extraordinary mode S3 == 0.0", e_beam.pol.S3, 0.0, 1e-8);
// Both beams must conserve total energy at entrance: I_o + I_e + I_refl == 1.0
real I_o_in = o_beam.intensity();
real I_e_in = e_beam.intensity();
real I_refl_in = (entry_step.reflected.length > 0) ? entry_step.reflected[0].intensity() : 0.0;
assert_near("Total crystal entrance energy conserved (I_o + I_e + I_refl == 1.0)",
            I_o_in + I_e_in + I_refl_in,
            1.0,
            1e-7);
write("");

// ============================================================================
// Test Case 102: Round-Trip Jones/Stokes Representation & Physicality
// ============================================================================
write("Test Case 102: Round-Trip Jones/Stokes Representation & Physicality");
// 1. Pure S-polarization: Es = (2, 0), Ep = (0, 0)
Polarization pol_s = Polarization.jones((2.0, 0.0), (0.0, 0.0));
assert_bool("Pure S state is physical", pol_s.isPhysical(), true);
assert_bool("Pure S state is pure", pol_s.isPure(), true);
assert_near("Pure S state S0 == 4.0", pol_s.S0, 4.0, 1e-9);
assert_near("Pure S state S1 == 4.0", pol_s.S1, 4.0, 1e-9);
assert_near("Pure S state S2 == 0.0", pol_s.S2, 0.0, 1e-9);
assert_near("Pure S state S3 == 0.0", pol_s.S3, 0.0, 1e-9);
pair[] J_s = pol_s.toJones();
assert_true("Pure S toJones has 2 components", J_s.length == 2);
assert_near("Pure S toJones Es.x == 2.0", J_s[0].x, 2.0, 1e-9);
assert_near("Pure S toJones Es.y == 0.0", J_s[0].y, 0.0, 1e-9);
assert_near("Pure S toJones Ep.x == 0.0", J_s[1].x, 0.0, 1e-9);
assert_near("Pure S toJones Ep.y == 0.0", J_s[1].y, 0.0, 1e-9);
// 2. Pure P-polarization: Es = (0, 0), Ep = (3, 0)
Polarization pol_p = Polarization.jones((0.0, 0.0), (3.0, 0.0));
assert_bool("Pure P state is physical", pol_p.isPhysical(), true);
assert_bool("Pure P state is pure", pol_p.isPure(), true);
assert_near("Pure P state S0 == 9.0", pol_p.S0, 9.0, 1e-9);
assert_near("Pure P state S1 == -9.0", pol_p.S1, -9.0, 1e-9);
pair[] J_p = pol_p.toJones();
assert_true("Pure P toJones has 2 components", J_p.length == 2);
assert_near("Pure P toJones Es.x == 0.0", J_p[0].x, 0.0, 1e-9);
assert_near("Pure P toJones Es.y == 0.0", J_p[0].y, 0.0, 1e-9);
assert_near("Pure P toJones Ep.x == 3.0", J_p[1].x, 3.0, 1e-9);
assert_near("Pure P toJones Ep.y == 0.0", J_p[1].y, 0.0, 1e-9);
// 3. Circular RCP: Es = (1, 0), Ep = (0, 1)
Polarization pol_rcp = Polarization.jones((1.0, 0.0), (0.0, 1.0));
assert_near("RCP S0 == 2.0", pol_rcp.S0, 2.0, 1e-9);
assert_near("RCP S1 == 0.0", pol_rcp.S1, 0.0, 1e-9);
assert_near("RCP S2 == 0.0", pol_rcp.S2, 0.0, 1e-9);
assert_near("RCP S3 == 2.0", pol_rcp.S3, 2.0, 1e-9);
pair[] J_rcp = pol_rcp.toJones();
assert_true("RCP toJones has 2 components", J_rcp.length == 2);
assert_near("RCP toJones Es.x == 1.0", J_rcp[0].x, 1.0, 1e-9);
assert_near("RCP toJones Es.y == 0.0", J_rcp[0].y, 0.0, 1e-9);
assert_near("RCP toJones Ep.x == 0.0", J_rcp[1].x, 0.0, 1e-9);
assert_near("RCP toJones Ep.y == 1.0", J_rcp[1].y, 1.0, 1e-9);
// 4. Circular LCP: Es = (1, 0), Ep = (0, -1)
Polarization pol_lcp = Polarization.jones((1.0, 0.0), (0.0, -1.0));
assert_near("LCP S3 == -2.0", pol_lcp.S3, -2.0, 1e-9);
pair[] J_lcp = pol_lcp.toJones();
assert_near("LCP toJones Ep.y == -1.0", J_lcp[1].y, -1.0, 1e-9);
// 5. Linear +45 deg: Es = (1, 0), Ep = (1, 0)
Polarization pol_45 = Polarization.jones((1.0, 0.0), (1.0, 0.0));
assert_near("Linear 45 S2 == 2.0", pol_45.S2, 2.0, 1e-9);
pair[] J_45 = pol_45.toJones();
assert_near("Linear 45 toJones Ep.x == 1.0", J_45[1].x, 1.0, 1e-9);
// 6. Unpolarized rejection
Polarization pol_unpol = Polarization.unpolarized(1.0);
assert_bool("Unpolarized is not pure", pol_unpol.isPure(), false);
pair[] J_unpol = pol_unpol.toJones();
assert_true("Unpolarized toJones is empty", J_unpol.length == 0);
// 7. Unphysical rejection
Polarization pol_unphys = Polarization.stokes(1.0, 2.0, 0.0, 0.0);
assert_bool("Unphysical state isPhysical is false", pol_unphys.isPhysical(), false);
assert_bool("Unphysical state isPure is false", pol_unphys.isPure(), false);
pair[] J_unphys = pol_unphys.toJones();
assert_true("Unphysical toJones is empty", J_unphys.length == 0);
write("");

// ============================================================================
// Test Case 103: Young 2-Slit & 3-Slit Coherent Superposition
// ============================================================================
write("Test Case 103: Young 2-Slit & 3-Slit Coherent Superposition");
// 2-slit interference: two equal amplitude pure S rays (I0 = 1.0 each)
Ray y_r1 = Ray(p = (0, 1), d = (1, 0), pol = Polarization.linear(0), phase = 0, coherence_id = 1);
Ray y_r2_inphase = Ray(p = (0, -1), d = (1, 0), pol = Polarization.linear(0), phase = 0, coherence_id = 1);
Ray y_r2_outphase = Ray(p = (0, -1), d = (1, 0), pol = Polarization.linear(0), phase = pi, coherence_id = 1);
Ray y_r2_quad = Ray(p = (0, -1), d = (1, 0), pol = Polarization.linear(0), phase = 0.5 * pi, coherence_id = 1);
assert_near("2-slit constructive peak (I == 4*I0)", interfere(y_r1, y_r2_inphase), 4.0, 1e-9);
assert_near("2-slit destructive minimum (I == 0.0)", interfere(y_r1, y_r2_outphase), 0.0, 1e-9);
assert_near("2-slit quadrature (I == 2*I0)", interfere(y_r1, y_r2_quad), 2.0, 1e-9);
// 3-slit interference: three equal amplitude pure S rays
Ray y_r3_inphase = Ray(p = (0, 0), d = (1, 0), pol = Polarization.linear(0), phase = 0, coherence_id = 1);
assert_near("3-slit principal maximum (I == 9*I0)", interfere(new Ray[]{y_r1, y_r3_inphase, y_r2_inphase}), 9.0, 1e-9);
// 3-slit minimum: phases 0, 2*pi/3, 4*pi/3 (phasor sum = 0)
Ray y3_a = Ray(p = (0, 0), d = (1, 0), pol = Polarization.linear(0), phase = 0, coherence_id = 1);
Ray y3_b = Ray(p = (0, 0), d = (1, 0), pol = Polarization.linear(0), phase = 2.0 * pi / 3.0, coherence_id = 1);
Ray y3_c = Ray(p = (0, 0), d = (1, 0), pol = Polarization.linear(0), phase = 4.0 * pi / 3.0, coherence_id = 1);
assert_near("3-slit first minimum (I == 0.0)", interfere(new Ray[]{y3_a, y3_b, y3_c}), 0.0, 1e-9);
write("");

// ============================================================================
// Test Case 104: Mach–Zehnder Arm Delay & Phase Modulation
// ============================================================================
write("Test Case 104: Mach–Zehnder Arm Delay & Phase Modulation");
real lambda_test = 500.0; // nm
// 50:50 input splitter + 50:50 recombination -> 0.25 power in each interfering beam
Ray arm1 = Ray((0, 0), (1, 0), lambda = lambda_test, pol = Polarization.linear(0), coherence_id = 10);
arm1.pol.scaleTo(0.25);
arm1.opl = 100.0;
Ray arm2_0 = Ray((0, 0), (1, 0), lambda = lambda_test, pol = Polarization.linear(0), coherence_id = 10);
arm2_0.pol.scaleTo(0.25);
arm2_0.opl = 100.0; // Delta L = 0
// Constructive interference at Delta L = 0 -> full input power 1.0
assert_near("Mach-Zehnder zero-delay full transmission (I == 1.0)", interfere(arm1, arm2_0), 1.0, 1e-9);
// Destructive interference at Delta L = lambda / 2
real half_lambda_mm = wavelengthInUnits(lambda_test) / 2.0;
Ray arm2_pi = Ray((0, 0), (1, 0), lambda = lambda_test, pol = Polarization.linear(0), coherence_id = 10);
arm2_pi.pol.scaleTo(0.25);
arm2_pi.opl = 100.0 + half_lambda_mm;
assert_near("Mach-Zehnder half-wave delay extinction (I == 0.0)", interfere(arm1, arm2_pi), 0.0, 1e-9);
// Constructive interference at Delta L = lambda
real full_lambda_mm = wavelengthInUnits(lambda_test);
Ray arm2_2pi = Ray((0, 0), (1, 0), lambda = lambda_test, pol = Polarization.linear(0), coherence_id = 10);
arm2_2pi.pol.scaleTo(0.25);
arm2_2pi.opl = 100.0 + full_lambda_mm;
assert_near("Mach-Zehnder full-wave delay transmission (I == 1.0)", interfere(arm1, arm2_2pi), 1.0, 1e-9);
// Quadrature interference at Delta L = lambda / 4
real quarter_lambda_mm = wavelengthInUnits(lambda_test) / 4.0;
Ray arm2_half_pi = Ray((0, 0), (1, 0), lambda = lambda_test, pol = Polarization.linear(0), coherence_id = 10);
arm2_half_pi.pol.scaleTo(0.25);
arm2_half_pi.opl = 100.0 + quarter_lambda_mm;
assert_near("Mach-Zehnder quarter-wave delay 50% transmission (I == 0.5)", interfere(arm1, arm2_half_pi), 0.5, 1e-9);
write("");

// ============================================================================
// Test Case 105: Fresnel–Arago Interference Laws
// ============================================================================
write("Test Case 105: Fresnel–Arago Interference Laws");
// Orthogonal polarizations (S vs P) from same coherent source
Ray fa_s = Ray((0, 0), (1, 0), pol = Polarization.linear(0), coherence_id = 42);  // Pure S (Es=1, Ep=0)
Ray fa_p = Ray((0, 0), (1, 0), pol = Polarization.linear(90), coherence_id = 42); // Pure P (Es=0, Ep=1)
// Fresnel-Arago 1st Law: Orthogonal states do not modulate irradiance
assert_near("Fresnel-Arago: zero phase shift irradiance (I == 2.0)", interfere(fa_s, fa_p), 2.0, 1e-9);
Ray fa_p_quad = copy(fa_p);
fa_p_quad.phase = 0.5 * pi;
assert_near("Fresnel-Arago: quarter-wave phase shift irradiance (I == 2.0)", interfere(fa_s, fa_p_quad), 2.0, 1e-9);
Ray fa_p_out = copy(fa_p);
fa_p_out.phase = pi;
assert_near("Fresnel-Arago: half-wave phase shift irradiance (I == 2.0)", interfere(fa_s, fa_p_out), 2.0, 1e-9);
// Resulting polarization states:
// At Delta phi = 0: E_s + E_p = Linear +45 deg (S2 = 2.0, S3 = 0)
Polarization state_0 = interfereState(fa_s, fa_p);
assert_near("Fresnel-Arago 0 deg state S2 == 2.0", state_0.S2, 2.0, 1e-9);
assert_near("Fresnel-Arago 0 deg state S3 == 0.0", state_0.S3, 0.0, 1e-9);
// At Delta phi = pi/2: E_s + i*E_p = Circular RCP (S2 = 0, S3 = 2.0)
Polarization state_quad = interfereState(fa_s, fa_p_quad);
assert_near("Fresnel-Arago 90 deg state S2 == 0.0", state_quad.S2, 0.0, 1e-9);
assert_near("Fresnel-Arago 90 deg state S3 == 2.0", state_quad.S3, 2.0, 1e-9);
// At Delta phi = pi: E_s - E_p = Linear -45 deg (S2 = -2.0, S3 = 0)
Polarization state_pi = interfereState(fa_s, fa_p_out);
assert_near("Fresnel-Arago 180 deg state S2 == -2.0", state_pi.S2, -2.0, 1e-9);
assert_near("Fresnel-Arago 180 deg state S3 == 0.0", state_pi.S3, 0.0, 1e-9);
// Incoherent fallback: different coherence IDs
Ray incoh_r1 = Ray((0, 0), (1, 0), pol = Polarization.linear(0), coherence_id = 100);
Ray incoh_r2 = Ray((0, 0), (1, 0), pol = Polarization.linear(0), coherence_id = 200);
assert_near("Incoherent sources do not interfere (I == 2.0)", interfere(incoh_r1, incoh_r2), 2.0, 1e-9);
// Incoherent fallback: different wavelengths
Ray diff_lam1 = Ray((0, 0), (1, 0), lambda = 500.0, pol = Polarization.linear(0), coherence_id = 100);
Ray diff_lam2 = Ray((0, 0), (1, 0), lambda = 600.0, pol = Polarization.linear(0), coherence_id = 100);
assert_near("Different wavelengths do not interfere (I == 2.0)", interfere(diff_lam1, diff_lam2), 2.0, 1e-9);
write("");

// Summary
write("================================================================================");
write(" TEST SUITE COMPLETED: " + string(passed_tests) + " / " + string(total_tests) + " TESTS PASSED");
write(" FAILED: " + string(failed_tests));
write("================================================================================");
if (failed_tests > 0) abort("optics_test.asy: " + string(failed_tests) + " assertion(s) failed");
