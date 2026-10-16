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

void assert_near(string testName, real val, real expected, real tol = 1e-4) {
    ++total_tests;
    real diff = abs(val - expected);
    if (diff <= tol) {
        ++passed_tests;
        write("  [PASS] " + testName + " = " + string(val));
    } else {
        ++failed_tests;
        write("  [FAIL] " + testName + " = " + string(val) + " (expected " + string(expected) + ", abs diff " +
              string(diff) + ", tolerance " + string(tol) + ")");
    }
};

void assert_near_pair(string testName, pair val, pair expected, real tol = 1e-4) {
    ++total_tests;
    real diff = length(val - expected);
    if (diff <= tol) {
        ++passed_tests;
        write("  [PASS] " + testName + " = (" + string(val.x) + ", " + string(val.y) + ")");
    } else {
        ++failed_tests;
        write("  [FAIL] " + testName + " = (" + string(val.x) + ", " + string(val.y) + ") (expected (" +
              string(expected.x) + ", " + string(expected.y) + "), error " + string(diff) + ", tolerance " +
              string(tol) + ")");
    }
};

void assert_bool(string testName, bool val, bool expected) {
    ++total_tests;
    if (val == expected) {
        ++passed_tests;
        write("  [PASS] " + testName + " = " + (val ? "true" : "false"));
    } else {
        ++failed_tests;
        write("  [FAIL] " + testName + " = " + (val ? "true" : "false") + " (expected " +
              (expected ? "true" : "false") + ")");
    }
};

void assert_true(string testName, bool condition) { assert_bool(testName, condition, true); }

void assert_false(string testName, bool condition) { assert_bool(testName, condition, false); }

void assert_relative(string testName, real val, real expected, real relTol = 1e-6, real absTol = 1e-10) {
    real allowed = max(absTol, relTol * max(abs(val), abs(expected)));
    assert_near(testName, val, expected, allowed);
};

void assert_stokes_physical(string name, Polarization pol, real tol = 1e-5) {
    assert_bool(name + ": nonnegative S0", pol.S0 >= -tol, true);
    real polarized2 = pol.S1 * pol.S1 + pol.S2 * pol.S2 + pol.S3 * pol.S3;
    real total2 = pol.S0 * pol.S0;
    assert_bool(name + ": Stokes cone invariant", polarized2 <= total2 + tol, true);
};

// Spatial resolvers for test boundaries
mediumResolver testResolveIsotropic = new Medium(pair p, Ray r, Medium m1, Medium m2, bool fromM1) {
    return (p.x > 0) ? medium("glass", 1.5) : Materials.VACUUM;
};

mediumResolver testResolveTIR = new Medium(pair p, Ray r, Medium m1, Medium m2, bool fromM1) {
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
FresnelResult ref1 = fresnelIsotropic(Materials.VACUUM, glass_15, r_inc, h_int, testResolveIsotropic);
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

assert_bool("TIR flag", ref1.totalInternalReflection, false);
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
FresnelResult ref2 = fresnelIsotropic(glass_15, Materials.VACUUM, r_tir, h_int, testResolveTIR);
pair expected_refl_dir = (-cos(theta_tir_rad), sin(theta_tir_rad));

assert_bool("TIR flag", ref2.totalInternalReflection, true);
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
Element testLens = Lens.thin(center = (x_lens, 0.0), angle = 0.0, f = f_lens, diameter = 50.0);
OpticalSystem sysLens = system(testLens, Materials.VACUUM);
real y_top = 10.0;
real y_bot = -10.0;
Ray ray_top = Ray(p = (0.0, y_top), d = (1.0, 0.0), medium = Materials.VACUUM);
Ray ray_bot = Ray(p = (0.0, y_bot), d = (1.0, 0.0), medium = Materials.VACUUM);
Trace t_top = sysLens.trace(ray_top);
Trace t_bot = sysLens.trace(ray_bot);
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
Element concaveMirror = Mirror.idealSpherical(vertex = (x_mirror, 0.0), angle = 180.0, f = f_mirror, diameter = 40.0);
OpticalSystem sysMirror = system(concaveMirror, Materials.VACUUM);
Ray ray_mirror_top = Ray(p = (0.0, y_inc), d = (1.0, 0.0), medium = Materials.VACUUM);
Trace t_mirror = sysMirror.trace(ray_mirror_top);
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
indexProfile grinProf = new real(pair p, real lambda) { return n0_grin * (1.0 - 0.5 * g_grin ^ 2 * p.y ^ 2); };
gradientProfile grinGrad = new pair(pair p, real lambda) { return (0.0, -n0_grin * g_grin ^ 2 * p.y); };
Medium grin_medium = medium("Parabolic GRIN").add(GRIN(grinProf, grinGrad));
Element grinBox = Interface(box((0.0, -20.0), (70.0, 20.0)), grin_medium, Materials.VACUUM);
OpticalSystem sysGRIN = system(grinBox, Materials.VACUUM);

real y0_grin = 1.0;
Ray ray_grin = Ray(p = (0.0, y0_grin), d = (1.0, 0.0), medium = grin_medium);

real x_quarter = pi / (2.0 * g_grin);
Trace t_grin = sysGRIN.trace(ray_grin, ds = 0.01, far = x_quarter);
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
VolumeProperty uniaxial = UniaxialExtraordinary(no = no_c, ne = ne_c, opticAxis = optic_axis);
real theta_k_deg = 45.0;
real theta_k_rad = radians(theta_k_deg);
pair k_norm = (cos(theta_k_rad), sin(theta_k_rad));

real inv_ne2 = (cos(theta_k_rad) ^ 2) / (no_c ^ 2) + (sin(theta_k_rad) ^ 2) / (ne_c ^ 2);
real expected_ne = 1.0 / sqrt(inv_ne2);
real actual_ne = uniaxial.phaseIndex((0.0, 0.0), k_norm, 633.0, Polarization.unpolarized());

assert_near("Extraordinary phase index n_e(45 deg)", actual_ne, expected_ne);

pair S_vec = (k_norm.x / (no_c ^ 2), k_norm.y / (ne_c ^ 2));
pair expected_S_dir = unit(S_vec);
real expected_theta_s = degrees(atan2(expected_S_dir.y, expected_S_dir.x));
pair actual_S_dir = uniaxial.groupDir(k_norm, (0.0, 0.0), 633.0);

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
OpticalSystem crystalSys = system(crystal, Materials.VACUUM);

Ray ray_crystal = Ray(p = (-10.0, 0.0), d = (cos(pi / 6.0), sin(pi / 6.0)), medium = Materials.VACUUM);
ray_crystal.pol = Polarization.unpolarized();
Trace t_crystal = crystalSys.trace(ray_crystal);

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
real e_index = e_ray.medium.n((0, 0), e_ray.waveNormal);
assert_true("Extraordinary phase index lies between ne and no", e_index >= 1.486 - 1e-9 && e_index <= 1.658 + 1e-9);
write("");

// ============================================================================
// Test Case 8: Diffraction Grating Spectral Order Angle
// ============================================================================
write("Test Case 8: Diffraction Grating Spectral Order Angle");
real lines_per_mm = 600.0;
real lambda_nm = 500.0;
real lambda_um = lambda_nm / 1000.0;
real d_grating_um = 1000.0 / lines_per_mm;
int m_order = 1;
Element testGrating = Grating(center = (0.0, 0.0), diameter = 10.0, angle = 90.0, lines_per_mm = lines_per_mm);
OpticalSystem gratingSys = system(testGrating, Materials.VACUUM);

Ray ray_grating = Ray(p = (0.0, -1.0), d = (0.0, 1.0), medium = Materials.VACUUM, lambda = lambda_nm);
Trace t_grating = gratingSys.trace(ray_grating);
Ray diff_ray = t_grating.transmitted(0);

real sin_theta_diff = (m_order * lambda_um) / d_grating_um;
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
Medium dispersiveGlass = medium("Dispersive Glass", A_cauchy, B_cauchy);
real lambda_red_nm = 700.0;
real lambda_blue_nm = 400.0;
real lambda_red_um = lambda_red_nm / 1000.0;
real lambda_blue_um = lambda_blue_nm / 1000.0;
real expected_n_red = A_cauchy + B_cauchy / (lambda_red_um ^ 2);
real expected_n_blue = A_cauchy + B_cauchy / (lambda_blue_um ^ 2);

assert_near("Red index n(700nm)", dispersiveGlass.n((0, 0), (1, 0), lambda_red_nm), expected_n_red);
assert_near("Blue index n(400nm)", dispersiveGlass.n((0, 0), (1, 0), lambda_blue_nm), expected_n_blue);

Ray ray_red = Ray(p = (-1, 0), d = (cos(pi / 4), sin(pi / 4)), medium = Materials.VACUUM, lambda = lambda_red_nm);
Ray ray_blue = Ray(p = (-1, 0), d = (cos(pi / 4), sin(pi / 4)), medium = Materials.VACUUM, lambda = lambda_blue_nm);

mediumResolver resolveDispersion = new Medium(pair p, Ray r, Medium m1, Medium m2, bool fromM1) {
    return (p.x > 0) ? dispersiveGlass : Materials.VACUUM;
};

FresnelResult ref_red = fresnelIsotropic(Materials.VACUUM, dispersiveGlass, ray_red, h_int, resolveDispersion);
FresnelResult ref_blue = fresnelIsotropic(Materials.VACUUM, dispersiveGlass, ray_blue, h_int, resolveDispersion);

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
Element malusPolarizer =
    Interface(box((0.0, -10.0), (10.0, 10.0)), medium("Malus Pol").add(Dichroic(0.0, 100.0)), Materials.VACUUM);
OpticalSystem malusSys = system(malusPolarizer, Materials.VACUUM);

real theta_pol_deg = 30.0;
Polarization pol_in = Polarization.linear(theta_pol_deg, 1.0);
Ray r_in = Ray((-10, 0), (1, 0), pol = pol_in);
Trace t_malus = malusSys.trace(r_in);

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
Medium sugarWater = medium("Sugar Water").add(OpticalActivity(5.0));

Polarization p_in = Polarization.linear(0.0, 1.0); // S = (1, 1, 0, 0)
Polarization p_out = sugarWater.stepPolarization(p_in, (0, 0), (1, 0), 633.0, 10.0);

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
Medium absorbingGlass = medium("Absorbing Glass", n_glass_abs, alpha = alpha_glass);
Element glassBlock = Interface(box((d1, -10.0), (d1 + d2, 10.0)), absorbingGlass, Materials.VACUUM);
OpticalSystem multiSys = system(glassBlock, Materials.VACUUM);

Ray ray_multi = Ray(p = (0.0, 0.0), d = (1.0, 0.0), medium = Materials.VACUUM);
Trace t_multi = multiSys.trace(ray_multi, far = far_reach);
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
Medium faradayMedium = medium("Faraday Glass").add(FaradayRotator(V_verdet, B_field));

// Forward propagation: d = (1, 0) over ds = 3.0 units. Expected rotation = +30 deg.
Polarization pol_init_f = Polarization.linear(0.0, 1.0); // Pure S
Polarization pol_rot_f = faradayMedium.stepPolarization(pol_init_f, (0, 0), (1, 0), 633.0, 3.0);

assert_near("Forward Faraday S1 (cos 60 deg)", pol_rot_f.S1, cos(radians(60.0)), 1e-4);
assert_near("Forward Faraday S2 (sin 60 deg)", pol_rot_f.S2, sin(radians(60.0)), 1e-4);
assert_near("Forward Faraday tilt angle psi (+30 deg)", degrees(pol_rot_f.psi()), 30.0, 1e-3);
assert_near("Forward Faraday intensity conserved", pol_rot_f.intensity(), 1.0, 1e-6);

// Backward propagation: d = (-1, 0) over ds = 3.0 units. Expected rotation = -30 deg.
Polarization pol_rot_b = faradayMedium.stepPolarization(pol_init_f, (0, 0), (-1, 0), 633.0, 3.0);

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
Element pgInterface =
    PhaseGradient(center = (0.0, 0.0), diameter = 10.0, angle = 90.0, deltaK = deltaK_shift, reflective = false);
OpticalSystem pgSys = system(pgInterface, Materials.VACUUM);

Ray ray_pg = Ray(p = (0.0, -1.0), d = (0.0, 1.0), medium = Materials.VACUUM);
Trace t_pg = pgSys.trace(ray_pg);
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

Element sphericalMir = Mirror.spherical(vertex = (x_mir, 0.0), angle = 180.0, f = f_mir, diameter = 80.0);
Element idealMir = Mirror.idealSpherical(vertex = (x_mir, 0.0), angle = 180.0, f = f_mir, diameter = 80.0);

OpticalSystem sysSpherical = system(sphericalMir, Materials.VACUUM);
OpticalSystem sysIdeal = system(idealMir, Materials.VACUUM);

Ray ray_marginal = Ray(p = (0.0, y_margin), d = (1.0, 0.0), medium = Materials.VACUUM);

Trace t_spherical = sysSpherical.trace(ray_marginal);
Trace t_ideal = sysIdeal.trace(ray_marginal);

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
Element spatialFilter = Aperture(center = (50.0, 0.0), diameter = 50.0, hole = 10.0, angle = 0.0);
OpticalSystem filterSys = system(spatialFilter, Materials.VACUUM);

Ray ray_inside = Ray(p = (0.0, 2.0), d = (1.0, 0.0), medium = Materials.VACUUM);
Trace t_inside = filterSys.trace(ray_inside, far = 100.0);

Ray ray_outside = Ray(p = (0.0, 8.0), d = (1.0, 0.0), medium = Materials.VACUUM);
Trace t_outside = filterSys.trace(ray_outside, far = 100.0);

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
Element beamSplitter =
    BeamSplitter.cube(center = (50.0, 0.0), angle = 0.0, size = 20.0, R = 0.4, glass = Materials.VACUUM);
OpticalSystem splitSys = system(beamSplitter, Materials.VACUUM);

Ray ray_split = Ray(p = (0.0, 0.0), d = (1.0, 0.0), medium = Materials.VACUUM, pol = Polarization.linear(0.0, 1.0));
Trace t_split = splitSys.trace(ray_split, far = 100.0);

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

Medium activeMedium = medium("Amplifying Gas", 1.0).add(Gain(g_amp));
Element gainSlab = Interface(box((d_entry, -10.0), (d_entry + d_slab, 10.0)), activeMedium, Materials.VACUUM);
OpticalSystem ampSys = system(gainSlab, Materials.VACUUM);

Ray ray_amp = Ray(p = (0.0, 0.0), d = (1.0, 0.0), medium = Materials.VACUUM);
Trace t_amp = ampSys.trace(ray_amp, far = far_exit);
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
Element brewsterWindow = Interface(box((0.0, -20.0), (10.0, 20.0)), glass_brewster, Materials.VACUUM);
OpticalSystem brewsterSys = system(brewsterWindow, Materials.VACUUM);

Ray ray_p = Ray(p = (-10.0 * cos(theta_B_rad), -10.0 * sin(theta_B_rad)),
                d = (cos(theta_B_rad), sin(theta_B_rad)),
                medium = Materials.VACUUM,
                pol = Polarization.linear(90.0, 1.0)); // Pure P

Trace t_p = brewsterSys.trace(ray_p);
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
Trace t_s = brewsterSys.trace(ray_s);
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
Element thickLens = Lens.biconvex(center = (center_x, 0.0),
                                  angle = 0.0,
                                  R1 = R1,
                                  R2 = R2,
                                  thickness = t_thick,
                                  diameter = 30.0,
                                  glass = bk7_thick,
                                  ambient = Materials.VACUUM);
OpticalSystem thickSys = system(thickLens, Materials.VACUUM);

real y_paraxial = 1.0;
Ray ray_thick = Ray(p = (0.0, y_paraxial), d = (1.0, 0.0), medium = Materials.VACUUM);
Trace t_thick_lens = thickSys.trace(ray_thick);
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

Element waterTank = Interface(box((0.0, -40.0), (200.0, 40.0)), water, Materials.VACUUM);
Element lensInside = Lens.biconvex(center = (100.0, 0.0),
                                   angle = 0.0,
                                   R1 = 50.0,
                                   R2 = 50.0,
                                   thickness = 6.0,
                                   diameter = 30.0,
                                   glass = glass_in_water);

OpticalSystem nestedSys = system(new Element[]{waterTank, lensInside}, Materials.VACUUM);

Ray ray_nested = Ray(p = (-10.0, 5.0), d = (1.0, 0.0), medium = Materials.VACUUM);
Trace t_nested = nestedSys.trace(ray_nested, far = 250.0);

assert_near("Total interface interactions in nested system", t_nested.interactions().length, 4.0);
assert_near_pair("Entry into water tank at x = 0", t_nested.hit(0), (0.0, 5.0));
assert_bool("Ray properly resolves water as ambient for the lens", abs(t_nested.hit(1).x - 97.0) < 1.0, true);
assert_bool("Ray exits water tank at x = 200", abs(t_nested.hit(3).x - 200.0) < 1e-4, true);
write("");

// ============================================================================
// Test Case 22: 90-Degree Retroreflector Invariance
// ============================================================================
write("Test Case 22: 90-Degree Retroreflector Invariance");
Element mirrorH = Mirror.flat(center = (25.0, 0.0), angle = 90.0, diameter = 50.0);
Element mirrorV = Mirror.flat(center = (0.0, 25.0), angle = 0.0, diameter = 50.0);

OpticalSystem retroSys = system(new Element[]{mirrorH, mirrorV}, Materials.VACUUM);

real theta_retro = radians(37.0);
pair d_in_retro = (-cos(theta_retro), -sin(theta_retro));
Ray ray_retro = Ray(p = (40.0, 30.0), d = d_in_retro, medium = Materials.VACUUM);

Trace t_retro = retroSys.trace(ray_retro);

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

Element isoLens = Lens.thin(center = (0.0, 0.0), angle = 0.0, f = f_iso, diameter = 60.0);
OpticalSystem isoSys = system(isoLens, Materials.VACUUM);

pair point_source = (-s_o, 0.0);
pair image_point = (s_i, 0.0);

Ray ray_axial = Ray(p = point_source, d = (1.0, 0.0), medium = Materials.VACUUM);
Trace t_axial = isoSys.trace(ray_axial, far = s_i);

pair d_marginal = unit((s_o, h_marginal));
Ray ray_marginal_iso = Ray(p = point_source, d = d_marginal, medium = Materials.VACUUM);
Trace t_marginal = isoSys.trace(ray_marginal_iso, far = s_i * 1.5);

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

Element triGrating = Grating(center = (0.0, 0.0),
                             diameter = 20.0,
                             angle = 90.0,
                             lines_per_mm = 300.0,
                             orders = diff_orders,
                             efficiencies = diff_efficiencies);
OpticalSystem triSys = system(triGrating, Materials.VACUUM);

Ray ray_tri = Ray(p = (0.0, -5.0), d = (0.0, 1.0), medium = Materials.VACUUM, lambda = 600.0);
Trace t_tri = triSys.trace(ray_tri);

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
Element lensAim = Lens.thin(center = (100.0, 0.0), angle = 0.0, f = 100.0, diameter = 80.0);
OpticalSystem sysAim = system(lensAim, Materials.VACUUM);

pair ptA = (0.0, 20.0);
pair ptB = (300.0, -10.0);

real aimed_angle_deg = aimAngle(ptA, ptB, sysAim, -10.0, 10.0);
real expected_aimed_angle = degrees(atan(10.0 / 100.0));

assert_near("Aimed launch angle (degrees)", aimed_angle_deg, expected_aimed_angle, 1e-3);

Ray ray_aimed = Ray.angle(ptA, aimed_angle_deg, Materials.VACUUM);
Trace t_aimed = sysAim.trace(ray_aimed, far = 400.0);
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

Medium grinMed = medium("Parabolic GRIN").add(GRIN(gp, gg));
Element grinSlab = Interface(box((0.0, -10.0), (50.0, 10.0)), grinMed, Materials.VACUUM);
OpticalSystem sysGrinCompare = system(grinSlab, Materials.VACUUM);

real y0_test = 0.2;
real s_travel = 5.0;
Ray r_init = Ray(p = (0.0, y0_test), d = (1.0, 0.0), medium = grinMed);

config.physics.grinSolver = GRINSolverType.RK4;
Trace t_rk4 = sysGrinCompare.trace(r_init, ds = 0.05, far = s_travel);
Ray ray_rk4 = t_rk4.last();

config.physics.grinSolver = GRINSolverType.RK45;
Trace t_rk45 = sysGrinCompare.trace(r_init, ds = 0.05, far = s_travel);
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
indexProfile axialGradientIndex = new real(pair p, real lambda) { return 1.0 + p.x; };
gradientProfile axialGradient = new pair(pair p, real lambda) { return (1.0, 0.0); };
Medium axialGradientMedium = medium("Axial gradient").add(GRIN(axialGradientIndex, axialGradient));
OpticalSystem emptySystem = system(new Element[]);
Trace axialGradientTrace = emptySystem.trace(Ray((0, 0), (1, 0), axialGradientMedium), ds = 0.05, far = 1.0);
assert_near("GRIN OPL integral of n(x)=1+x over [0,1]", axialGradientTrace.opl(), 1.5, 1e-8);
assert_near_pair("Line intersection", lineIntersection((0, 0), (1, 0), (2, -1), (0, 1)), (2, 0));
pair parallelIntersection = lineIntersection((0, 0), (1, 0), (0, 1), (2, 0));
assert_bool("Parallel line intersection is infinite", abs(parallelIntersection.x) == infinity, true);
write("");

// ============================================================================
// Test Case 30: Ordered Custom Surface Properties
// ============================================================================
write("Test Case 30: Ordered Custom Surface Properties");
Property halfPower = property(
    "Half power coating",
    new Ray[](Ray r, Hit h, mediumResolver resolve) {
        r.pol.scaleTo(0.5 * r.intensity());
        return new Ray[]{r};
    },
    PropertyStage.AfterBoundary);
path testBoundary = linearBoundary((0, 0), (1, 0), 20);
Element bareInterface = Interface(testBoundary, medium("test glass", 1.5));
Element coatedInterface = Interface(testBoundary, medium("test glass", 1.5));
coatedInterface.surfaces[0].add(halfPower);
Ray coatingProbe = Ray.angle((-5, 0), 10, Materials.VACUUM);
real barePower = system(bareInterface).trace(coatingProbe).transmitted(0).intensity();
real coatedPower = system(coatedInterface).trace(coatingProbe).transmitted(0).intensity();
assert_near("AfterBoundary custom action runs after refraction", coatedPower, 0.5 * barePower, 1e-8);

VolumeProperty lowPriorityIndex = volumeProperty("low priority", priority = 0);
lowPriorityIndex.phaseIndex = new real(pair p, pair d, real lambda, Polarization pol) { return 1.2; };
VolumeProperty highPriorityIndex = volumeProperty("high priority", priority = 10);
highPriorityIndex.phaseIndex = new real(pair p, pair d, real lambda, Polarization pol) { return 1.6; };
Medium orderedMedium = medium("ordered").add(highPriorityIndex).add(lowPriorityIndex);
assert_near("Volume-index profiles compose multiplicatively", orderedMedium.n(), 1.92);
assert_bool("Volume-property priority determines evaluation order",
            orderedMedium.orderedProperties()[1].name == "high priority",
            true);
write("");

// ============================================================================
// Test Case 29b: Diverging and Parallel Beam Generator Factories
// ============================================================================
write("Test Case 29b: Diverging and Parallel Beam Generator Factories");
Ray[] pBeam = Beam.parallel(start = (0.0, -4.0), step = (0.0, 2.0), count = 5, d = (1.0, 0.0));
assert_near("Parallel beam ray count", pBeam.length, 5.0);
assert_near_pair("Parallel beam Ray 0 start point", pBeam[0].p, (0.0, -4.0));
assert_near_pair("Parallel beam Ray 2 start point (center)", pBeam[2].p, (0.0, 0.0));
assert_near_pair("Parallel beam Ray 4 start point", pBeam[4].p, (0.0, 4.0));
assert_near_pair("Parallel beam Ray direction (+X)", pBeam[3].d, (1.0, 0.0));

Ray[] dBeam = Beam.diverging(source = (0.0, 0.0), d = (1.0, 0.0), angle = 30.0, count = 3);
assert_near("Diverging beam ray count", dBeam.length, 3.0);
assert_near_pair("Diverging beam Ray 1 direction (center)", dBeam[1].d, (1.0, 0.0));
assert_near("Diverging beam Ray 0 angle (-15 deg)", degrees(atan2(dBeam[0].d.y, dBeam[0].d.x)), -15.0);
assert_near("Diverging beam Ray 2 angle (+15 deg)", degrees(atan2(dBeam[2].d.y, dBeam[2].d.x)), 15.0);
write("");

// ============================================================================
// Test Case 30b: Biconcave Lens Beam Divergence & Virtual Focus
// ============================================================================
write("Test Case 30b: Biconcave Lens Beam Divergence & Virtual Focus");
Element biconcaveLens = Lens.thin(center = (100.0, 0.0), angle = 0.0, f = -50.0, diameter = 40.0);
OpticalSystem sysConcave = system(biconcaveLens, Materials.VACUUM);

Ray ray_cc_top = Ray(p = (0.0, 10.0), d = (1.0, 0.0), medium = Materials.VACUUM);
Ray ray_cc_bot = Ray(p = (0.0, -10.0), d = (1.0, 0.0), medium = Materials.VACUUM);

Trace t_cc_top = sysConcave.trace(ray_cc_top);
Trace t_cc_bot = sysConcave.trace(ray_cc_bot);

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
Medium polMedium = medium("DichroicFilter").add(Dichroic(0.0, 10.0));
Polarization out_unpol = polMedium.stepPolarization(Polarization.unpolarized(1.0), (0, 0), (1, 0), 633.0, 1.0);

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

mediumResolver resFlintWater = new Medium(pair p, Ray r, Medium m1, Medium m2, bool fromM1) {
    return (p.x > 0) ? water_m : flint;
};
Hit h_planar;
h_planar.valid = true;
h_planar.normal = (-1, 0);

FresnelResult ref_below = fresnelIsotropic(flint, water_m, ray_below_crit, h_planar, resFlintWater);
FresnelResult ref_above = fresnelIsotropic(flint, water_m, ray_above_crit, h_planar, resFlintWater);

assert_bool("TIR below critical angle is false", ref_below.totalInternalReflection, false);
assert_bool("TIR above critical angle is true", ref_above.totalInternalReflection, true);
write("");

// ============================================================================
// Test Case 33: Ideal Mirror Collimation from Focal Point
// ============================================================================
write("Test Case 33: Ideal Mirror Collimation from Focal Point");
Element oap = Mirror.idealSpherical(vertex = (50, 0), angle = 180.0, f = 50.0, diameter = 40.0);
OpticalSystem oapSys = system(oap, Materials.VACUUM);

Ray oap_ray1 = Ray.between((0, 0), (50, 5), Materials.VACUUM);
Ray oap_ray2 = Ray.between((0, 0), (50, -5), Materials.VACUUM);

Trace t_oap1 = oapSys.trace(oap_ray1);
Trace t_oap2 = oapSys.trace(oap_ray2);

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
OpticalSystem invSys = system(lens, Materials.VACUUM);

real y_obj = 10.0;
real u_in = 0.05;

Ray r_marginal = Ray(p = (-s_obj, 0), d = (cos(u_in), sin(u_in)), medium = Materials.VACUUM);
Trace t_marg = invSys.trace(r_marginal, far = 300.0);
Ray exit_marg = t_marg.transmitted(0);
real u_out = abs(asin(exit_marg.d.y));

Ray r_chief = Ray.between((-s_obj, y_obj), (0, 0), Materials.VACUUM);
Trace t_chief = invSys.trace(r_chief, far = 300.0);
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
real lambda_nm_test = 587.5618;
real l2 = (lambda_nm_test / 1000.0) ^ 2;
real B1 = 1.03961212, C1 = 0.00600069867;
real B2 = 0.231792344, C2 = 0.0200179144;
real B3 = 1.01046945, C3 = 103.560653;
real n_sellmeier_expected = sqrt(1.0 + (B1 * l2) / (l2 - C1) + (B2 * l2) / (l2 - C2) + (B3 * l2) / (l2 - C3));

indexProfile sellmeierBK7 = new real(pair p, real lambda) {
    real w2 = (lambda / 1000.0) ^ 2;
    return sqrt(1.0 + (1.03961212 * w2) / (w2 - 0.00600069867) + (0.231792344 * w2) / (w2 - 0.0200179144) +
                (1.01046945 * w2) / (w2 - 103.560653));
};
Medium bk7_disp = medium("BK7 Sellmeier").add(GRIN(sellmeierBK7));
assert_near("BK7 refractive index at d-line (587.56 nm)", bk7_disp.n(lambda = lambda_nm_test), 1.51680, 1e-4);
write("");

// ============================================================================
// Test Case 36: Transmissive PhaseShift Target Medium Continuity
// ============================================================================
write("Test Case 36: Transmissive PhaseShift Target Medium Continuity");
Medium denseGlass = medium("Dense Glass", 2.0);
Element metaSheet = PhaseGradient(center = (0, 0), diameter = 20, angle = 90, deltaK = 0.2, reflective = false);
OpticalSystem psSys = system(metaSheet, denseGlass);

Ray r_ps = Ray(p = (-5, 0), d = (1, 0), medium = Materials.VACUUM);
Trace t_ps = psSys.trace(r_ps);
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
Element flatMirror = Mirror.flat(center = (50, 0), angle = 180, diameter = 40);
OpticalSystem mirrorSys = system(flatMirror, Materials.VACUUM);

Ray ray_in = Ray(p = (0, 0), d = (1, 0), pol = Polarization.linear(45.0, 1.0));
ray_in.phase = 0.0;
real I_in = ray_in.intensity();

Trace t_refl = mirrorSys.trace(ray_in, far = 100.0);
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
path prismPath = b_left -- b_right -- apex -- cycle;

Medium glass15 = medium("Glass 1.5", 1.50);
Element prism = Interface(prismPath, glass15, Materials.VACUUM);
OpticalSystem prismSys = system(prism, Materials.VACUUM);

real theta_i = degrees(asin(1.5 * sin(radians(30.0))));
real ray_angle_lab = -30.0 + theta_i;
pair hit_left = midpoint(b_left -- apex);
Ray r_prism =
    Ray.toPoint(target = hit_left, d = (cos(radians(ray_angle_lab)), sin(radians(ray_angle_lab))), distance = 20.0);

Trace t_prism = prismSys.trace(r_prism, far = 40.0);
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
path corePath = box((0, -5), (100, 5));
Element guide = Interface(corePath, medium("Core", 1.50), Materials.VACUUM);
OpticalSystem guideSys = system(guide, Materials.VACUUM);

Ray r_guide = Ray(p = (0.01, 0), d = (cos(radians(25.0)), sin(radians(25.0))), medium = medium("Core", 1.50));
Trace t_guide = guideSys.trace(r_guide, far = 120.0);

Ray exit_guide = t_guide.last();
real guideN1 = 1.50;
real guideN2 = 1.0;
real guideThetaI = radians(25.0);
real guideSinT = guideN1 / guideN2 * sin(guideThetaI);
real guideCosT = sqrt(1.0 - guideSinT ^ 2);
real guideCosI = cos(guideThetaI);
real guideRs = ((guideN1 * guideCosI - guideN2 * guideCosT) / (guideN1 * guideCosI + guideN2 * guideCosT)) ^ 2;
real guideRp = ((guideN2 * guideCosI - guideN1 * guideCosT) / (guideN2 * guideCosI + guideN1 * guideCosT)) ^ 2;
real guideExpectedTransmission = 1.0 - 0.5 * (guideRs + guideRp);
assert_near("Waveguide output power from end-face Fresnel transmission",
            exit_guide.intensity(),
            guideExpectedTransmission,
            1e-4);
assert_bool("Ray traversed through the guide to x >= 100", exit_guide.p.x >= 100.0 - 1e-3, true);
write("");

// ============================================================================
// Test Case 41: Degree of Polarization (DoP) & Stokes Vector at Brewster Angle
// ============================================================================
write("Test Case 41: Degree of Polarization (DoP) & Stokes Vector at Brewster Angle");
Medium glass_br_41 = medium("Glass", 1.50);
real theta_B_deg_41 = degrees(atan(1.50 / 1.0));

path boundaryPlane = (0, -50) -- (0, 50);
Element iface = Interface(boundaryPlane, Materials.VACUUM, glass_br_41);
OpticalSystem brSys = system(iface, Materials.VACUUM);

Ray r_unpol = Ray.toPoint(target = (0, 0),
                          d = (cos(radians(theta_B_deg_41)), sin(radians(theta_B_deg_41))),
                          distance = 20.0,
                          pol = Polarization.unpolarized(1.0));
Trace t_br = brSys.trace(r_unpol);

Ray r_refl = t_br.reflected(0);
Ray r_trans = t_br.transmitted(0);

assert_near("Brewster Reflected Ray DoP = 100%", r_refl.pol.dop(), 1.0, 1e-4);
assert_near("Brewster Reflected Ray is purely S-polarized (S1 == S0)", r_refl.pol.S1, r_refl.pol.S0, 1e-4);

real brewsterSinT = sin(radians(theta_B_deg_41)) / glass_br_41.n();
real brewsterCosT = sqrt(1.0 - brewsterSinT ^ 2);
real brewsterRsAmplitude = (cos(radians(theta_B_deg_41)) - glass_br_41.n() * brewsterCosT) /
                           (cos(radians(theta_B_deg_41)) + glass_br_41.n() * brewsterCosT);
real brewsterRsPower = brewsterRsAmplitude ^ 2;
real brewsterTs = 1.0 - brewsterRsPower;
real brewsterTp = 1.0;
real expected_dop = (brewsterTp - brewsterTs) / (brewsterTp + brewsterTs);
assert_near("Brewster Transmitted Ray DoP matches analytical formula (~7.99%)", r_trans.pol.dop(), expected_dop, 1e-3);
write("");

// ============================================================================
// Test Case 42: Porro Prism Double-TIR Retroreflector
// ============================================================================
write("Test Case 42: Porro Prism Double-TIR Retroreflector");
path porroPath = (0, -20) -- (20, 0) -- (0, 20) -- cycle;
Element porro = Interface(porroPath, medium("Glass", 1.50), Materials.VACUUM);
OpticalSystem porroSys = system(porro, Materials.VACUUM);

Ray r_porro = Ray(p = (-10, 8.0), d = (1, 0), medium = Materials.VACUUM);
Trace t_porro = porroSys.trace(r_porro, far = 40.0);
Ray r_porro_exit = t_porro.last();

assert_near("Porro retroreflected d.x is -1", r_porro_exit.d.x, -1.0, 1e-4);
assert_near("Porro retroreflected d.y is 0", r_porro_exit.d.y, 0.0, 1e-4);
assert_near("Porro inverted exit position y = -8.0", r_porro_exit.p.y, -8.0, 1e-4);

real porroR0 = ((1.50 - 1.0) / (1.50 + 1.0)) ^ 2;
real porroT0 = 1.0 - porroR0;
assert_near("Porro prism intensity from two Fresnel transmissions", r_porro_exit.intensity(), porroT0 ^ 2, 1e-4);
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

real expectedV1x = 50.0 - thickness_cardinal / 2.0;
real expectedV2x = 50.0 + thickness_cardinal / 2.0;
assert_near("Front vertex from center and thickness", cp.V1.x, expectedV1x, 1e-3);
assert_near("Rear vertex from center and thickness", cp.V2.x, expectedV2x, 1e-3);

real lensPower =
    (n_cardinal - 1.0) * (1.0 / R1_cardinal - 1.0 / R2_cardinal +
                          ((n_cardinal - 1.0) * thickness_cardinal) / (n_cardinal * R1_cardinal * R2_cardinal));
real expectedEFL = 1.0 / lensPower;
assert_relative("Effective focal length from thick lensmaker equation", cp.f2, expectedEFL, 1e-5, 1e-5);

real expectedBFL = expectedEFL * (1.0 - ((n_cardinal - 1.0) * thickness_cardinal) / (n_cardinal * R1_cardinal));
assert_relative("Back focal length from thick-lens paraxial formula", cp.BFL, expectedBFL, 1e-5, 1e-5);
assert_relative("Rear focal point is one BFL beyond rear vertex", cp.F2.x, cp.V2.x + cp.BFL, 1e-5, 1e-5);
assert_relative("Front focal point is one BFL before front vertex", cp.F1.x, cp.V1.x - cp.BFL, 1e-5, 1e-5);

real principalOffset = cp.f2 - cp.BFL;
assert_relative("Front principal point offset", cp.H1.x, cp.V1.x + principalOffset, 1e-5, 1e-5);
assert_relative("Rear principal point offset", cp.H2.x, cp.V2.x - principalOffset, 1e-5, 1e-5);
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
Polarization extraUnpolarized = Polarization.unpolarized(1.0);
Polarization extraUnpolarizedCopy = copy(extraUnpolarized);

assert_relative("Copy preserves unpolarized intensity", extraUnpolarizedCopy.intensity(), extraUnpolarized.intensity());
assert_relative("Copy preserves unpolarized degree of polarization",
                extraUnpolarizedCopy.dop(),
                extraUnpolarized.dop());

Polarization extraLinear = Polarization.linear(23.0, 1.7);
Polarization extraLinearCopy = copy(extraLinear);

assert_relative("Copy preserves linear-state intensity", extraLinearCopy.intensity(), extraLinear.intensity());
assert_near("Copy preserves S1", extraLinearCopy.S1, extraLinear.S1, 1e-5);
assert_near("Copy preserves S2", extraLinearCopy.S2, extraLinear.S2, 1e-5);
assert_near("Copy preserves S3", extraLinearCopy.S3, extraLinear.S3, 1e-5);
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

Polarization unpolRot = Polarization.unpolarized(2.0);
unpolRot.applyRotation(radians(37.0));

assert_near("Unpolarized rotation: S0 preserved", unpolRot.S0, 2.0, 1e-5);
assert_near("Unpolarized rotation: S1 remains zero", unpolRot.S1, 0.0, 1e-5);
assert_near("Unpolarized rotation: S2 remains zero", unpolRot.S2, 0.0, 1e-5);
assert_near("Unpolarized rotation: S3 remains zero", unpolRot.S3, 0.0, 1e-5);
assert_near("Unpolarized rotation: DoP remains zero", unpolRot.dop(), 0.0, 1e-5);
write("");

// ============================================================================
// Test Case 46: Normal-incidence Fresnel boundary conditions
// ============================================================================
write("Test Case 46: Normal-incidence Fresnel boundary conditions");
Ray extraNormalRay = Ray(p = (-1.0, 0.0), d = (1.0, 0.0), medium = Materials.VACUUM);
extraNormalRay.pol = Polarization.unpolarized();
FresnelResult extraNormalRef =
    fresnelIsotropic(Materials.VACUUM, glass_15, extraNormalRay, h_int, testResolveIsotropic);
real extraNormalR = ((n2 - n1) / (n2 + n1)) ^ 2;

assert_bool("Normal incidence is not total internal reflection", extraNormalRef.totalInternalReflection, false);
assert_near_pair("Normal-incidence transmitted direction", extraNormalRef.d, (1.0, 0.0));
assert_relative("Normal-incidence Fresnel reflectance", extraNormalRef.R, extraNormalR);
assert_relative("Normal-incidence energy balance", extraNormalRef.R + extraNormalRef.T, 1.0);
write("");

// ============================================================================
// Test Case 47: Polarization intensity scaling laws
// ============================================================================
write("Test Case 47: Polarization intensity scaling laws");
Polarization extraScaleSource = Polarization.linear(31.0, 1.0);
real extraScaleFactor = 2.75;
Polarization extraScaleResult = extraScaleFactor * extraScaleSource;

assert_relative("Intensity scales linearly with scalar factor",
                extraScaleResult.intensity(),
                extraScaleSource.intensity() * extraScaleFactor);
assert_relative("S1 scales proportionally", extraScaleResult.S1, extraScaleSource.S1 *extraScaleFactor);
assert_relative("DOP preserved under scaling", extraScaleResult.dop(), extraScaleSource.dop());
write("");

// ============================================================================
// Test Case 48: Stokes physical-state invariants
// ============================================================================
write("Test Case 48: Stokes physical-state invariants");
Polarization extraS = Polarization.linear(0.0, 1.0);
Polarization extraP = Polarization.linear(90.0, 1.0);
Polarization extraCirc = Polarization.circular();

assert_relative("S-polarized light has unit degree of polarization", extraS.dop(), 1.0);
assert_relative("P-polarized light has unit degree of polarization", extraP.dop(), 1.0);
assert_relative("Circularly polarized light has unit degree of polarization", extraCirc.dop(), 1.0);
assert_relative("S-polarized Stokes intensity agrees with S0", extraS.S0, extraS.intensity());
assert_relative("P-polarized Stokes intensity agrees with S0", extraP.S0, extraP.intensity());
assert_relative("Circular Stokes intensity agrees with S0", extraCirc.S0, extraCirc.intensity());
write("");

// Test Case 49: Lossless retarder conserves intensity and degree of polarization
write("Test Case 49: Lossless retarder conservation");
Polarization extraRetarderInput = copy(extraCirc);
Polarization extraRetarderOutput = copy(extraRetarderInput);
applyDiattenuation(extraRetarderOutput, extraRetarderInput, 1.0, 1.0, (cos(pi / 3.0), -sin(pi / 3.0)));

assert_relative("Lossless retarder preserves S0", extraRetarderOutput.S0, extraRetarderInput.S0);
assert_relative("Lossless retarder preserves degree of polarization",
                extraRetarderOutput.dop(),
                extraRetarderInput.dop());
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
Polarization phaseInput = new Polarization;
phaseInput.S0 = 1.0;
phaseInput.S1 = 0.0;
phaseInput.S2 = 1.0;
phaseInput.S3 = 0.0;
Polarization phaseOutput = copy(phaseInput);
applyDiattenuation(phaseOutput, phaseInput, 0.5, 0.5, (0.3, 0.2));

assert_near("Complex phase: S0", phaseOutput.S0, 0.5, 1e-5);
assert_near("Complex phase: S1", phaseOutput.S1, 0.0, 1e-5);
assert_near("Complex phase: S2", phaseOutput.S2, 0.3, 1e-5);
assert_near("Complex phase: S3", phaseOutput.S3, 0.2, 1e-5);

real phasePol2 = phaseOutput.S1 ^ 2 + phaseOutput.S2 ^ 2 + phaseOutput.S3 ^ 2;

assert_bool("Complex phase: output remains inside Stokes cone", phasePol2 <= phaseOutput.S0 ^ 2 + 1e-5, true);
write("");

// ============================================================================
// Test Case 52: Maximum-coherence physicality boundary
// ============================================================================
write("Test Case 52: Maximum-coherence physicality boundary");
Polarization boundaryInput = new Polarization;
boundaryInput.S0 = 1.0;
boundaryInput.S1 = 0.0;
boundaryInput.S2 = 1.0;
boundaryInput.S3 = 0.0;
Polarization boundaryOutput = copy(boundaryInput);
applyDiattenuation(boundaryOutput, boundaryInput, 0.8, 0.2, (0.24, 0.32));

assert_near("Boundary: S0", boundaryOutput.S0, 0.5, 1e-5);
assert_near("Boundary: S1", boundaryOutput.S1, 0.3, 1e-5);
assert_near("Boundary: S2", boundaryOutput.S2, 0.24, 1e-5);
assert_near("Boundary: S3", boundaryOutput.S3, 0.32, 1e-5);

real boundaryPol2 = boundaryOutput.S1 ^ 2 + boundaryOutput.S2 ^ 2 + boundaryOutput.S3 ^ 2;

assert_near("Boundary: fully polarized output lies on sphere", boundaryPol2, boundaryOutput.S0 ^ 2, 1e-5);
write("");

// ============================================================================
// Test Case 53: Dichroic Medium Propagation & Mueller Theory
// ============================================================================
write("Test Case 53: Dichroic Medium Propagation & Mueller Theory");
real alpha_s_val = 0.0;
real alpha_p_val = 1.0;
real prop_dist = 2.0;
Medium dichroicFilter = medium("DichroicFilter").add(Dichroic(alpha_s_val, alpha_p_val));
Polarization syncInput = Polarization.linear(45.0, 1.0);
Polarization syncOutput = dichroicFilter.stepPolarization(syncInput, (0, 0), (1, 0), 633.0, prop_dist);

assert_relative("Stokes S0 agrees with output intensity", syncOutput.S0, syncOutput.intensity(), 1e-5, 1e-7);
assert_relative("Stokes degree of polarization is physical", syncOutput.dop(), 1.0, 1e-5, 1e-7);

Polarization syncCopy = copy(syncOutput);

assert_relative("Copied Stokes S0 agrees with intensity", syncCopy.S0, syncCopy.intensity(), 1e-5, 1e-7);
assert_near("Copied S1 preserved", syncCopy.S1, syncOutput.S1, 1e-5);
assert_near("Copied S2 preserved", syncCopy.S2, syncOutput.S2, 1e-5);
assert_near("Copied S3 preserved", syncCopy.S3, syncOutput.S3, 1e-5);

Polarization syncUnpolInput = Polarization.unpolarized(1.0);
Polarization syncUnpolOutput = dichroicFilter.stepPolarization(syncUnpolInput, (0, 0), (1, 0), 633.0, prop_dist);

assert_relative("Unpolarized-path S0 agrees with intensity",
                syncUnpolOutput.S0,
                syncUnpolOutput.intensity(),
                1e-5,
                1e-7);

real Ts = exp(-alpha_s_val * prop_dist);
real Tp = exp(-alpha_p_val * prop_dist);

assert_relative("Diattenuator S0 matches Mueller prediction",
                syncUnpolOutput.S0,
                syncUnpolInput.S0 *(Ts + Tp) / 2.0,
                1e-5,
                1e-7);

assert_relative("Diattenuator S1 matches Mueller prediction",
                syncUnpolOutput.S1,
                syncUnpolInput.S0 *(Ts - Tp) / 2.0,
                1e-5,
                1e-7);

assert_near("Diattenuator S2 remains zero", syncUnpolOutput.S2, 0.0);
assert_near("Diattenuator S3 remains zero", syncUnpolOutput.S3, 0.0);
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

real sourceS0 = original.S0;
real sourceS1 = original.S1;
real sourceS2 = original.S2;
real sourceS3 = original.S3;
cloned.S0 = 3.0;
cloned.S1 = 1.0;
cloned.S2 = 0.5;
cloned.S3 = -0.25;

assert_near("Copy isolation: source Stokes S0 unchanged", original.S0, sourceS0);
assert_near("Copy isolation: source Stokes S1 unchanged", original.S1, sourceS1);
assert_near("Copy isolation: source Stokes S2 unchanged", original.S2, sourceS2);
assert_near("Copy isolation: source Stokes S3 unchanged", original.S3, sourceS3);

real clonedS0 = cloned.S0;
real clonedS1 = cloned.S1;
original.S0 = 2.0;
original.S1 = -0.5;

assert_near("Copy isolation: clone S0 unchanged after source mutation", cloned.S0, clonedS0);
assert_near("Copy isolation: clone S1 unchanged after source mutation", cloned.S1, clonedS1);

Polarization unpol_test = Polarization.unpolarized(2.0);
Polarization unpolCopy = copy(unpol_test);

assert_near("Copy isolation: unpolarized S0 copied", unpolCopy.S0, unpol_test.S0);
assert_near("Copy isolation: unpolarized S1 copied", unpolCopy.S1, 0.0);
assert_near("Copy isolation: unpolarized S2 copied", unpolCopy.S2, 0.0);
assert_near("Copy isolation: unpolarized S3 copied", unpolCopy.S3, 0.0);

unpolCopy.S0 = 7.0;
unpolCopy.S1 = 1.0;

assert_near("Copy isolation: unpolarized source S0 unchanged", unpol_test.S0, 2.0);
assert_near("Copy isolation: unpolarized source S1 unchanged", unpol_test.S1, 0.0);
write("");

// ============================================================================
// Test Case 55: Zero-distance propagation
// ============================================================================
write("Test Case 55: Zero-distance propagation");
Polarization zeroDistInput = Polarization.linear(45.0, 1.0);
real zeroS0 = zeroDistInput.S0;
real zeroS1 = zeroDistInput.S1;
real zeroS2 = zeroDistInput.S2;
real zeroS3 = zeroDistInput.S3;
Polarization zeroDistOutput = dichroicFilter.stepPolarization(zeroDistInput, (0, 0), (1, 0), 633.0, 0.0);

assert_near("Zero distance: S0 unchanged", zeroDistOutput.S0, zeroS0);
assert_near("Zero distance: S1 unchanged", zeroDistOutput.S1, zeroS1);
assert_near("Zero distance: S2 unchanged", zeroDistOutput.S2, zeroS2);
assert_near("Zero distance: S3 unchanged", zeroDistOutput.S3, zeroS3);
write("");

// ============================================================================
// Test Case 56: Unpolarized light through a diattenuator
// ============================================================================
write("Test Case 56: Unpolarized light through a diattenuator");
real testTs = 0.8;
real testTp = 0.2;
Polarization unpolSrcStokes = Polarization.unpolarized(1.0);
Polarization unpolDstStokes = copy(unpolSrcStokes);
applyDiattenuation(unpolDstStokes, unpolSrcStokes, testTs, testTp);
real expectedS0 = 0.5 * (testTs + testTp);
real expectedS1 = 0.5 * (testTs - testTp);

assert_near("Unpolarized diattenuation: S0", unpolDstStokes.S0, expectedS0);
assert_near("Unpolarized diattenuation: S1", unpolDstStokes.S1, expectedS1);
assert_near("Unpolarized diattenuation: S2 remains zero", unpolDstStokes.S2, 0.0);
assert_near("Unpolarized diattenuation: S3 remains zero", unpolDstStokes.S3, 0.0);
assert_near("Unpolarized diattenuation: degree of polarization", unpolDstStokes.dop(), 0.6);
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
path rhombPath = (0, 0) -- (w * L, h * L) -- (w * L, h * L + L) -- (0, L) -- cycle;
Element rhomb = Interface(rhombPath, glass_rhomb);
OpticalSystem rhombSys = Materials.VACUUM + rhomb;
// Launch linear 45-deg ray: S = (1, 0, 1, 0)
Ray r_in = Ray(p = (-20, L / 1.2), d = (1, 0), pol = Polarization.linear(45.0, 1.0));
Trace t_rhomb = rhombSys.trace(r_in);
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
OpticalSystem rhombSys = Materials.VACUUM + rhomb;
Ray r_unpol_rhomb = Ray(p = (-20, L / 1.2), d = (1, 0), pol = Polarization.unpolarized(1.0));
Trace t_unpol_rhomb = rhombSys.trace(r_unpol_rhomb);
Ray r_exit_unpol = t_unpol_rhomb.last(); // after both TIR bounces and the exit face

assert_near("Unpolarized exit S1 is 0", r_exit_unpol.pol.S1, 0.0, 1e-5);
assert_near("Unpolarized exit S2 is 0", r_exit_unpol.pol.S2, 0.0, 1e-5);
assert_near("Unpolarized exit S3 is 0", r_exit_unpol.pol.S3, 0.0, 1e-5);
assert_near("Unpolarized exit DOP remains 0", r_exit_unpol.pol.dop(), 0.0, 1e-5);
real rhombT0 = 1.0 - ((n_bk7 - 1.0) / (n_bk7 + 1.0)) ^ 2;
assert_near("Unpolarized exit intensity = two normal-incidence Fresnel transmissions",
            r_exit_unpol.intensity(),
            rhombT0 ^ 2,
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
Polarization inUnpol = Polarization.unpolarized(2.0);
real Is = 0.5 * (inUnpol.S0 + inUnpol.S1);
real Ip = 0.5 * (inUnpol.S0 - inUnpol.S1);
Polarization out_s = Polarization.linear(0.0, Is);
Polarization out_p = Polarization.linear(90.0, Ip);

assert_near("PBS S-channel intensity", out_s.S0, 1.0, 1e-5);
assert_near("PBS S-channel S1", out_s.S1, 1.0, 1e-5);
assert_near("PBS S-channel DOP", out_s.dop(), 1.0, 1e-5);
assert_near("PBS P-channel intensity", out_p.S0, 1.0, 1e-5);
assert_near("PBS P-channel S1", out_p.S1, -1.0, 1e-5);
assert_near("PBS P-channel DOP", out_p.dop(), 1.0, 1e-5);
// Same split performed by the engine: ideal linear polarizers at 0 and 90 deg.
Trace t_lp0 = system(Polarizer.LP((10, 0), 0.0, 10.0)).trace(Ray((0, 0), (1, 0), pol = copy(inUnpol)));
Trace t_lp90 = system(Polarizer.LP((10, 0), 90.0, 10.0)).trace(Ray((0, 0), (1, 0), pol = copy(inUnpol)));
assert_near("LP(0) transmits the S half of unpolarized light", t_lp0.last().pol.S0, Is, 1e-6);
assert_near("LP(0) output is pure S (S1 = S0)", t_lp0.last().pol.S1, t_lp0.last().pol.S0, 1e-6);
assert_near("LP(90) transmits the P half of unpolarized light", t_lp90.last().pol.S0, Ip, 1e-6);
assert_near("LP(90) output is pure P (S1 = -S0)", t_lp90.last().pol.S1, -t_lp90.last().pol.S0, 1e-6);
write("");

// ============================================================================
// Test Case 62: General Elliptical Polarization (Poincaré vs Semi-Axes)
// ============================================================================
write("Test Case 62: General Elliptical Polarization (Poincaré vs Semi-Axes)");
Polarization el1 = Polarization.elliptical(a = 1.0, b = 0.5, phase_deg = 90.0, theta_deg = 0.0);

assert_stokes_physical("Elliptical state 1", el1);
assert_near("Elliptical S0", el1.S0, 1.25, 1e-5);
assert_near("Elliptical S1", el1.S1, 0.75, 1e-5);
assert_near("Elliptical S2", el1.S2, 0.0, 1e-5);
assert_near("Elliptical S3", el1.S3, 1.0, 1e-5);
assert_near("Elliptical DOP", el1.dop(), 1.0, 1e-5);

Polarization el2 = Polarization.fromAngles(psi_deg = 45.0, chi_deg = 20.0, intensity = 1.0);

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

path boundaryPlane = (0, -50) -- (0, 50);
Element iface = Interface(boundaryPlane, m_air, m_glass);
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
pair c_axis = unit((1.0, 1.0)); // 45 deg tilted optic axis
Medium calcite_o = medium("Calcite_O", no_cal);
Medium calcite_e = uniaxialExtraordinary("Calcite_E", no_cal, ne_cal, c_axis);

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
assert_near("Focal length survives transform", rg_thinMoved.focalLength, 42.0, 1e-12);
assert_near("Focal length survives copy", copy(rg_thin).focalLength, 42.0, 1e-12);
Element rg_aligned = align(Lens.thin((50, 0), 0, -30.0, 10.0), (0, 0), (100, 10));
assert_near("Focal length survives align()", rg_aligned.focalLength, -30.0, 1e-12);
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
            length(rg_eExit.waveNormal - rg_eExit.d),
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
            abs(rg_oplMeasured - rg_eIn.medium.n((0, 0), rg_eIn.waveNormal) * rg_geomLen) > 1e-3);
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
mediumResolver rg_resU = new Medium(pair p, Ray r, Medium m1, Medium m2, bool fromM1) {
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
mediumResolver rg_resD = new Medium(pair p, Ray r, Medium m1, Medium m2, bool fromM1) {
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
mediumResolver rg_resG = new Medium(pair p, Ray r, Medium m1, Medium m2, bool fromM1) {
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
CardinalPoints rg_cpBig = cardinalPoints(rg_thinSys, paraxialH = 5.0);
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
config.physics.grinSolver = GRINSolverType.RK4;
Trace rg_tTiny = rg_thinGrin.trace(Ray((0.01, 0), (1, 0), rg_gl), ds = 5e-5, far = 1.0);
config.physics.grinSolver = GRINSolverType.RK45;
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
config.physics.grinSolver = GRINSolverType.RK4;
real rg_yRef = rg_grinOpen.trace(rg_r0, ds = 0.0025, far = 4.0).last().p.y;
real rg_e1 = abs(rg_grinOpen.trace(rg_r0, ds = 0.07, far = 4.0).last().p.y - rg_yRef);
real rg_e2 = abs(rg_grinOpen.trace(rg_r0, ds = 0.035, far = 4.0).last().p.y - rg_yRef);
config.physics.grinSolver = GRINSolverType.RK45;
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
// Test Case 88: Edge cases - grazing incidence, zero-intensity rays, maxDepth
// ============================================================================
write("Test Case 88: Edge cases");
mediumResolver rg_resE = new Medium(pair p, Ray r, Medium m1, Medium m2, bool fromM1) {
    return (p.x > 0) ? medium("g", 1.5) : Materials.VACUUM;
};
FresnelResult rg_graze = fresnelIsotropic(Materials.VACUUM, medium("g", 1.5), Ray((-1, 0), dir(89.99)), rg_h, rg_resE);
assert_true("Grazing incidence: R finite and close to 1", rg_graze.R > 0.99 && rg_graze.R <= 1.0);
assert_near("Grazing incidence: R + T = 1", rg_graze.R + rg_graze.T, 1.0, 1e-12);
Trace rg_dark = system(Interface(box((0, -5), (10, 5)), medium("g", 1.5)))
                    .trace(Ray((-5, 0), (1, 0), pol = Polarization.unpolarized(0.0)));
assert_near("Zero-intensity ray spawns no children", rg_dark.length(), 0.0);
OpticalSystem rg_cavity = system(new Element[]{Mirror.flat((0, 0), 0, 10), Mirror.flat((10, 0), 180, 10)});
Trace rg_tCav = rg_cavity.trace(Ray((5, 0), (1, 0)), maxDepth = 6);
assert_near("Cavity bounces limited by maxDepth", rg_tCav.length(), 6.0);
write("");

// Summary
write("================================================================================");
write(" TEST SUITE COMPLETED: " + string(passed_tests) + " / " + string(total_tests) + " TESTS PASSED");
write(" FAILED: " + string(failed_tests));
write("================================================================================");
if (failed_tests > 0) abort("optics_test.asy: " + string(failed_tests) + " assertion(s) failed");