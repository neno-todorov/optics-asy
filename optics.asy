// optics.asy
//
// A 2D optical ray tracing engine for Asymptote (requires Asymptote >= 3.11
// for the Trace subscript operator; older versions use Trace.get(i)).
//
// Naming convention:
//   * PascalCase       - structs, public optical constructors and factories
//   * camelCase        - operations, calculations, queries, and helpers
//   * snake_case       - variables, struct fields, parameters, and local values
//   * UPPER_SNAKE_CASE - constants and enum-like values
// Naming details:
//   * semantic qualifiers precede units: lambda_eff_nm, theta_inc_deg;
//   * unit suffixes are last: *_nm, *_um, *_mm, *_deg, *_rad;
//   * pen variables use the pen_* form: pen_ray, pen_laser, pen_line;
//   * standard optical notation (S0, rs, rp, n1, n2) is preserved.
//
// Units and physical scale:
//   * Public wavelength arguments are in nm; material dispersion models use µm
//     internally.  Wavelengths are never converted to geometry units for
//     physical calculations.
//   * Points, paths, ray positions, intersections, and element dimensions use
//     native Asymptote geometry units and are independent of rendering size.
//   * config.physics.meters_per_unit defines the physical length represented by
//     one geometry unit.  Physical distances are converted to µm before they
//     are used for OPL, phase, absorption, or polarization evolution.
//   * Ray.opl is stored in µm.  vacuumWaveNumber() returns rad/µm, so
//     vacuumWaveNumber(lambda) * Ray.opl is a dimensionless phase.
//   * Length-dependent material coefficients (for example absorption and
//     rotation coefficients) are interpreted per µm.
//
// Polarization and Fresnel conventions:
//   * Stokes vectors use one fixed global basis: s (= z, TE, out of the
//     drawing plane) and p (= in-plane, TM).
//   * S1 = Is - Ip; S2 and S3 describe the relative phase between the s and p
//     field components.
//   * Fresnel amplitudes use rs = (n1 ci - n2 ct)/(...) and
//     rp = (n2 ci - n1 ct)/(...); mirrors use the n2 -> infinity limit.
//   * In anisotropic media, d is the energy/group direction and wave_normal
//     is the phase-propagation direction.
// Run the test suite with:  asy -noV optics_test.asy

import geometry;
import graph;
import math;
import version;

usepackage("fontenc", "T2A");
usepackage("amsmath");
usepackage("siunitx");
usepackage("icomma");

// Soft Color Palette
pen base_red = rgb(0.96, 0.50, 0.50);    // Soft red
pen base_blue = rgb(0.40, 0.65, 0.90);   // Soft blue
pen base_green = rgb(0.45, 0.80, 0.45);  // Soft green
pen base_yellow = rgb(0.98, 0.85, 0.35); // Soft yellow
pen base_brown = rgb(0.75, 0.55, 0.40);  // Soft brown
pen dark_red = rgb(0.62, 0.32, 0.32);    // Dark red
pen dark_blue = rgb(0.26, 0.42, 0.58);   // Dark blue
pen dark_green = rgb(0.29, 0.52, 0.29);  // Dark green
pen dark_yellow = rgb(0.64, 0.55, 0.23); // Dark yellow
pen dark_brown = rgb(0.49, 0.36, 0.26);  // Dark brown

pen Red(real x = 1.0) { return interp(white, base_red, x); }
pen Blue(real x = 1.0) { return interp(white, base_blue, x); }
pen Green(real x = 1.0) { return interp(white, base_green, x); }
pen Yellow(real x = 1.0) { return interp(white, base_yellow, x); }
pen Brown(real x = 1.0) { return interp(white, base_brown, x); }
pen DarkRed(real x = 1.0) { return interp(white, dark_red, x); }
pen DarkBlue(real x = 1.0) { return interp(white, dark_blue, x); }
pen DarkGreen(real x = 1.0) { return interp(white, dark_green, x); }
pen DarkYellow(real x = 1.0) { return interp(white, dark_yellow, x); }
pen DarkBrown(real x = 1.0) { return interp(white, dark_brown, x); }

// ============================================================================
// PART I: CORE ENGINE
// ============================================================================
struct GRINSolverType {
    static int RK4 = 0;  // Classic Fixed-Step
    static int RK45 = 1; // Adaptive Cash-Karp
}

// Layered Settings & Config
struct MathConfig {
    real eps = 1e-8;
    real stokes_intensity_eps = 1e-12;
    real solver_denominator_eps = 1e-12;
}

struct PhysicsConfig {
    real lambda = 633.0;
    // Physical scale: metres represented by one native Asymptote geometry unit.
    // The historical convention is one geometry unit = one micrometre.
    real meters_per_unit = 1e-6;
    real min_intensity = 1e-3;
    real min_refractive_index = 1e-1;
    real eps_integration = 1e-4;
    real wavevector_alignment_tolerance = 1e-6;
    real aiming_direction_eps = 1e-7;
    real aiming_height_tolerance = 1e-5;
    real far = 1e3;
    real ds = 1e-2;
    real paraxial_height = 5e-2;
    real differential_angle = 5e-1;
    int grin_solver = GRINSolverType.RK45;
    int max_depth = 20;
    int max_branches = 128;
    bool trace_reflections = true;
    bool trace_diffraction = true;
}

struct GeometryConfig {
    real eps = 1e-8;
    real eps_intersection = 1e-8;
    real self_intersection_factor = 1e2;
    real near_intersection_factor = 1e1;
    real ray_to_path_extension = 1e3;
    real hemisphere_margin = 0.999;
    int max_bisection_steps = 25;
    int arc_subdivisions = 64;
}

struct RenderingConfig {
    real vis_lambda = nan;
    bool use_intensity_opacity = true;
    bool draw_polarization = false;
    bool draw_wavefronts = false;
    bool draw_reflections = false;
    bool draw_absorption_fading = true;
    pen pen_ray = black;
    pen pen_beam = black;
    pen pen_mirror = black;
    pen pen_absorber = black;
    pen pen_lens = black;
    real clip_velocity_eps = 1e-9;
    real min_lens_head = 10pt;
    real min_drawable_length = 1e-10;
    int fade_subpath_steps = 32;
    int grating_tick_count = 12;
    real grating_tick_angle = 20;
    real grating_tick_length = 0.12;
    real polarization_spacing = 0.5;
    real min_polarization_marker_with_arrow = 0.5;
    real min_polarization_segment_length = 0.40;    // Segments shorter than this show 0 markers
    real min_polarization_marker_separation = 0.25; // Minimum distance between consecutive markers
    real arrow_size = 5pt;
    real arrow_pol_size = arrow_size / 3;
    arrowhead arrow_head = HookHead;
    arrowhead arrow_pol_head = arrow_head;
    real polarization_inset_min_arrow_fraction = 0.18;
    real polarization_inset_unpolarized_threshold = 1e-4;
    int polarization_inset_unpolarized_spokes = 4;
    int polarization_inset_ellipse_steps = 36;
    real wavefront_max_dist = 120.0;
    int wavefront_max_orders = 3;
    int wavefront_max_markers = 2000;
}

struct Config {
    MathConfig math;
    PhysicsConfig physics;
    GeometryConfig geometry;
    RenderingConfig rendering;
}

Config
defaultConfig() {
    Config c = new Config;
    c.math = new MathConfig;
    c.physics = new PhysicsConfig;
    c.geometry = new GeometryConfig;
    c.rendering = new RenderingConfig;
    return c;
}

Config config = defaultConfig();

// Wavelengths are expressed in nm at the public API; geometric lengths use um.
// Unit-conversion ratios: multiply to convert between nm and um.
real nm_to_um = 1e-3; // multiply a value in nm to get um
real um_to_nm = 1e3;  // multiply a value in um to get nm

// Canonical physical length unit for optical calculations is micrometres.
real metersToUm(real metres) { return metres * 1e6; }
real geometryToUm(real distance) { return distance * metersToUm(config.physics.meters_per_unit); }
real geometryWavelengthUm(real lambda_nm) {
    return (lambda_nm * nm_to_um) / metersToUm(config.physics.meters_per_unit);
}
// resolveEffectiveLambda() is a rendering quantity: it returns geometry units.

// local_vis_lambda is an optional diagram wavelength in geometric units;
// src_lambda is a physical wavelength in nanometers.
real resolveEffectiveLambda(real local_vis_lambda, real src_lambda) {
    if (!isnan(local_vis_lambda) && local_vis_lambda > 0) return local_vis_lambda;
    if (!isnan(config.rendering.vis_lambda) && config.rendering.vis_lambda > 0) return config.rendering.vis_lambda;
    return geometryWavelengthUm(src_lambda);
}

real vacuumWaveNumber(real lambda_nm) {
    real lambda = lambda_nm * nm_to_um;
    return (lambda > config.math.eps) ? 2.0 * pi / lambda : 0.0;
}

// String Case Helpers
string LowerCase(string s) {
    string result = "";
    for (int i = 0; i < length(s); ++i) {
        int c = ascii(s[i]);
        if (c >= ascii("A") && c <= ascii("Z")) c += ascii("a") - ascii("A");
        result += format("%c", c);
    }
    return result;
}

string UpperCase(string s) {
    string result = "";
    for (int i = 0; i < length(s); ++i) {
        int c = ascii(s[i]);
        if (c >= ascii("a") && c <= ascii("z")) c -= ascii("a") - ascii("A");
        result += format("%c", c);
    }
    return result;
}

// Geometric Math Utilities
real clamp(real x, real a, real b) { return max(a, min(b, x)); }

// private
pair unitOr(pair v, pair fallback = (1, 0)) {
    if (length(v) > config.math.eps) return unit(v);
    if (length(fallback) > config.math.eps) return unit(fallback);
    return (1, 0);
}

// private
pair perp(pair v) {
    pair u = unitOr(v);
    return (-ypart(u), xpart(u));
}

// https://en.wikipedia.org/wiki/Shoelace_formula
private
real signedArea(path p) {
    real area = 0;
    int n = length(p);
    for (int i = 0; i < n; ++i) {
        pair p1 = point(p, i);
        pair p2 = point(p, i + 1);
        area += (p1.x * p2.y - p2.x * p1.y);
    }
    return area / 2.0;
}

private
path ensureCCW(path p) {
    if (!cyclic(p)) return p;
    return (signedArea(p) < 0) ? reverse(p) : p;
}

// --- General Path Geometry & Intersection Helpers ---
pair lineIntersection(pair p1, pair d1, pair p2, pair d2) {
    real scale = length(d1) * length(d2);
    if (scale <= config.math.eps) return (infinity, infinity);
    real denom = cross(d1, d2);
    if (abs(denom) <= config.math.eps * scale) return (infinity, infinity);
    return p1 + (cross(p2 - p1, d2) / denom) * d1;
}

pair intersection(path a, path b) {
    if (length(a) == 0 || length(b) == 0) {
        return (infinity, infinity);
    }
    if (length(a) == 1 && length(b) == 1 && straight(a, 0) && straight(b, 0)) {
        return lineIntersection(point(a, 0), dir(a, 0), point(b, 0), dir(b, 0));
    }
    pair[] pts = intersectionpoints(a, b);
    if (pts.length > 0) {
        return pts[0];
    }
    return (infinity, infinity);
}

// ============================================================================
// Polarization (Pure Stokes Representation)
// ============================================================================
struct Polarization {
    // Stokes parameters (Single Source of Truth)
    real S0 = 1.0;
    real S1 = 0.0;
    real S2 = 0.0;
    real S3 = 0.0;

    // --- State Queries ---
    real intensity() { return S0; }

    // Power carried by the s (TE, out-of-plane) and p (TM, in-plane)
    // components.  In this 2D engine the s axis is always the global z axis,
    // so S1 is referenced to the same s/p basis at every interface.
    real sPower() { return max(0.0, 0.5 * (S0 + S1)); }
    real pPower() { return max(0.0, 0.5 * (S0 - S1)); }

    real dop() {
        if (S0 <= config.math.eps) return 0.0;
        return clamp(sqrt(S1 * S1 + S2 * S2 + S3 * S3) / S0, 0.0, 1.0);
    }

    // Physicality check: S0 >= 0 and S0^2 >= S1^2 + S2^2 + S3^2
    bool isPhysical(real tol = config.math.eps) {
        if (S0 < -tol) return false;
        real p2 = S1 * S1 + S2 * S2 + S3 * S3;
        real max_s0 = S0 * (1.0 + tol);
        return p2 <= max_s0 * max_s0 + tol;
    }

    // Purity check based on Degree of Polarization (DOP)
    bool isPure(real tol = config.math.eps) {
        if (S0 <= config.math.eps) return false;
        if (!isPhysical(tol)) return false;
        return dop() >= 1.0 - tol;
    }

    // On-demand Jones vector extraction: returns [Es, Ep] (or empty array if not pure)
    pair[] toJones(real tol = config.math.eps) {
        if (!isPure(tol)) return new pair[]{};

        real s_pow = 0.5 * (S0 + S1);
        real p_pow = 0.5 * (S0 - S1);
        real Es_mag = (s_pow > config.math.eps) ? sqrt(s_pow) : 0.0;
        real Ep_mag = (p_pow > config.math.eps) ? sqrt(p_pow) : 0.0;

        pair Es;
        pair Ep;

        if (Es_mag <= config.math.eps) {
            // Pure P-polarized mode: Es is 0, Ep carries all amplitude
            Es = (0.0, 0.0);
            Ep = (Ep_mag, 0.0);
        } else if (Ep_mag <= config.math.eps) {
            // Pure S-polarized mode: Ep is 0, Es carries all amplitude
            Es = (Es_mag, 0.0);
            Ep = (0.0, 0.0);
        } else {
            // Both components nonzero: Es is real and positive, Ep carries relative phase
            real delta = atan2(S3, S2);
            Es = (Es_mag, 0.0);
            Ep = (Ep_mag * cos(delta), Ep_mag * sin(delta));
        }

        return new pair[]{Es, Ep};
    }

    // Azimuth / tilt of the major axis: [-pi/2, pi/2]
    real psi() { return 0.5 * atan2(S2, S1); }

    // Ellipticity angle: [-pi/4 (LCP) to +pi/4 (RCP)]
    real chi() {
        real p = dop();
        if (p <= config.math.eps) return 0.0;
        return 0.5 * asin(clamp(S3 / (S0 * p), -1.0, 1.0));
    }

    // --- In-place Transformations ---
    void scaleTo(real new_intensity) {
        if (this.S0 > config.math.eps) {
            real factor = new_intensity / this.S0;
            this.S0 = new_intensity;
            this.S1 *= factor;
            this.S2 *= factor;
            this.S3 *= factor;
        } else {
            this.S0 = new_intensity;
        }
    }

    // Polarization-independent attenuation/gain: scales the whole Stokes
    // vector, so the degree of polarization is preserved.
    void scale(real factor) {
        this.S0 *= factor;
        this.S1 *= factor;
        this.S2 *= factor;
        this.S3 *= factor;
    }

    // Phase retardance between p and s components (Fresnel Rhomb / Waveplates)
    void applyRetardance(real delta) {
        real c = cos(delta);
        real s = sin(delta);
        real new_s2 = S2 * c - S3 * s;
        real new_s3 = S2 * s + S3 * c;
        S2 = new_s2;
        S3 = new_s3;
    }

    // Rotation of the polarization plane by theta (Optical Activity / Faraday)
    void applyRotation(real theta) {
        real c2 = cos(2.0 * theta);
        real s2 = sin(2.0 * theta);
        real new_s1 = S1 * c2 - S2 * s2;
        real new_s2 = S1 * s2 + S2 * c2;
        S1 = new_s1;
        S2 = new_s2;
    }

    // ========================================================================
    // Static Factories
    // ========================================================================
    // 1. Unpolarized light
    static Polarization unpolarized(real intensity = 1.0) {
        Polarization p = new Polarization;
        p.S0 = intensity;
        p.S1 = 0.0;
        p.S2 = 0.0;
        p.S3 = 0.0;
        return p;
    }

    // 2. Linear polarization at angle theta
    static Polarization linear(real theta, real intensity = 1.0) {
        Polarization p = new Polarization;
        real th = radians(theta);
        p.S0 = intensity;
        p.S1 = intensity * cos(2.0 * th);
        p.S2 = intensity * sin(2.0 * th);
        p.S3 = 0.0;
        return p;
    }

    // 3. Circular polarization (+1: RCP, -1: LCP)
    static Polarization circular(bool ccw = true, real intensity = 1.0) {
        Polarization p = new Polarization;
        p.S0 = intensity;
        p.S1 = 0.0;
        p.S2 = 0.0;
        p.S3 = (ccw ? 1.0 : -1.0) * intensity;
        return p;
    }

    // 4. Elliptical polarization via Poincaré angles (psi, chi)
    static Polarization fromAngles(real psi, real chi, real intensity = 1.0) {
        Polarization p = new Polarization;
        real psi_rad = radians(psi);
        real chi_rad = radians(chi);
        p.S0 = intensity;
        p.S1 = intensity * cos(2.0 * chi_rad) * cos(2.0 * psi_rad);
        p.S2 = intensity * cos(2.0 * chi_rad) * sin(2.0 * psi_rad);
        p.S3 = intensity * sin(2.0 * chi_rad);
        return p;
    }

    // 5. Elliptical polarization via semi-axes (a, b), phase shift, and tilt
    static Polarization elliptical(real a = nan,
                                   real b = nan,
                                   real phase = 90.0,
                                   real theta = 0.0,
                                   real azimuth = nan,
                                   real ellipticity = nan,
                                   real intensity = 1.0,
                                   real amplitude = 1.0) {
        // 1. If explicit azimuth / ellipticity are given:
        if (!isnan(azimuth) && !isnan(ellipticity)) {
            real I = (!isnan(intensity) && intensity != 1.0) ? intensity : (amplitude * amplitude);
            return fromAngles(azimuth, ellipticity, I);
        }

        // 2. Safeguard: If first two positional args are large (likely azimuth & ellipticity passed positionally):
        if (!isnan(a) && !isnan(b) && (abs(a) > 2.0 || abs(b) > 2.0)) {
            real I = (!isnan(intensity) && intensity != 1.0) ? intensity : 1.0;
            return fromAngles(a, b, I);
        }

        // 3. Otherwise: semi-axes (a, b) with phase and tilt
        real a_val = isnan(a) ? 1.0 : a;
        real b_val = isnan(b) ? 0.5 : b;
        real delta_rad = radians(phase);
        real theta_rad = radians(theta);

        real s0 = a_val * a_val + b_val * b_val;
        real s1 = a_val * a_val - b_val * b_val;
        real s2 = 2.0 * a_val * b_val * cos(delta_rad);
        real s3 = 2.0 * a_val * b_val * sin(delta_rad);

        real c2 = cos(2.0 * theta);
        real s2_th = sin(2.0 * theta_rad);

        Polarization p = new Polarization;
        p.S0 = s0;
        p.S1 = s1 * c2 - s2 * s2_th;
        p.S2 = s1 * s2_th + s2 * c2;
        p.S3 = s3;
        return p;
    }

    // Convenient direct Stokes constructor:
    static Polarization stokes(real s0, real s1, real s2 = 0.0, real s3 = 0.0) {
        Polarization p = new Polarization;
        p.S0 = s0;
        p.S1 = s1;
        p.S2 = s2;
        p.S3 = s3;
        return p;
    }

    // Factory constructor from complex Jones amplitudes (parallel to .stokes())
    static Polarization jones(pair Es, pair Ep) {
        Polarization pol = new Polarization;
        real Is = abs2(Es);
        real Ip = abs2(Ep);
        pair cross = 2.0 * conj(Es) * Ep;

        pol.S0 = Is + Ip;
        pol.S1 = Is - Ip;
        pol.S2 = cross.x;
        pol.S3 = cross.y;
        return pol;
    }
};

pair[] toJones(Polarization pol, real tol = config.math.eps) {
    if (pol == null) return new pair[]{};
    return pol.toJones(tol);
}

// ============================================================================
// Standalone Functions & Operators
// ============================================================================

Polarization copy(Polarization p) {
    if (p == null) return null;
    Polarization c = new Polarization;
    c.S0 = p.S0;
    c.S1 = p.S1;
    c.S2 = p.S2;
    c.S3 = p.S3;
    return c;
}

// Incoherent ray superposition (Stokes addition)
Polarization operator+(Polarization a, Polarization b) {
    Polarization out = new Polarization;
    out.S0 = a.S0 + b.S0;
    out.S1 = a.S1 + b.S1;
    out.S2 = a.S2 + b.S2;
    out.S3 = a.S3 + b.S3;
    return out;
}

// Intensity scaling
Polarization operator*(real scalar, Polarization a) {
    Polarization out = new Polarization;
    out.S0 = a.S0 * scalar;
    out.S1 = a.S1 * scalar;
    out.S2 = a.S2 * scalar;
    out.S3 = a.S3 * scalar;
    return out;
}

Polarization operator*(Polarization a, real scalar) { return scalar * a; }

// Applies the 4x4 Mueller Diattenuation Matrix (with optional complex cross-coupling / retardance C)
// from a source Stokes vector to a target Stokes vector (or in-place when target == src),
// given orthogonal intensity transmittances/reflectances T1 and T2 (e.g. Ts, Tp or Rs, Rp).
// If C is omitted, C defaults to (sqrt(max(0.0, T1 * T2)), 0.0) (pure zero-retardance diattenuation).
// 1. General Mueller Diattenuation with optional phase cross-term C
void applyDiattenuation(Polarization target, Polarization src, real T1, real T2, pair C = (nan, nan)) {
    if (target == null || src == null) return;

    real half_sum = 0.5 * (T1 + T2);
    real half_diff = 0.5 * (T1 - T2);
    pair cross = isnan(C.x) ? (sqrt(max(0.0, T1 * T2)), 0.0) : C;

    real S0_new = half_sum * src.S0 + half_diff * src.S1;
    real S1_new = half_diff * src.S0 + half_sum * src.S1;
    real S2_new = cross.x * src.S2 - cross.y * src.S3;
    real S3_new = cross.y * src.S2 + cross.x * src.S3;

    target.S0 = S0_new;
    target.S1 = S1_new;
    target.S2 = S2_new;
    target.S3 = S3_new;
}

// 2. Overload for complex Fresnel amplitude coefficients (js, jp)
void applyDiattenuation(Polarization target, Polarization src, pair js, pair jp) {
    if (target == null || src == null) return;

    real T1 = js.x * js.x + js.y * js.y; // |js|^2
    real T2 = jp.x * jp.x + jp.y * jp.y; // |jp|^2

    // Cross-term C = js * conj(jp) = (x1 + i*y1) * (x2 - i*y2)
    pair cross = (js.x * jp.x + js.y * jp.y, js.y * jp.x - js.x * jp.y);

    applyDiattenuation(target, src, T1, T2, cross);
}

// Volume Properties (Medium Components)
// Anisotropic Boundary Packet & Function Types
struct AnisotropicTransmission {
    bool valid = false;
    pair wave_normal;
    pair group_direction;
    real phase_index;
}

struct PhaseMatch {
    bool tir = false;
    pair wave_normal;
    pair group_direction;
    real n;
}

// Function signatures
typedef real indexProfile(pair p, real lambda);
// Defines a base index generator: n = f(p, d, lambda, pol)
typedef real volumeIndexFunc(pair p, pair d, real lambda, Polarization pol);
// NEW: Defines an index modifier/transform: n_out = T(n_base, p, d, lambda, pol)
typedef real volumeModifierFunc(real n_base, pair p, pair d, real lambda, Polarization pol);
typedef real volumeAbsorptionFunc(pair p, pair d, real lambda, Polarization pol);
typedef pair gradientProfile(pair p, real lambda);
typedef pair volumeGradientFunc(pair p, pair d, real lambda);
typedef pair groupDirectionFunc(pair wave_normal, pair p, real lambda);
typedef pair waveNormalFunc(pair group_dir, pair p, real lambda);
typedef AnisotropicTransmission anisotropicSolverFunc(pair N, pair T, real kt, pair forward, pair p, real lambda);
typedef Polarization volumeStepFunc(Polarization pol, pair p, pair d, real lambda, real ds);

struct VolumeProperty {
    string name = "";
    // Lower values run first; equal priorities retain add() order.
    int priority = 0;
    bool spatial_index = false;
    volumeIndexFunc phase_index = null;       // Set for base models
    volumeModifierFunc phase_modifier = null; // Set for modifiers/perturbations
    volumeGradientFunc phase_gradient = null;
    volumeAbsorptionFunc attenuation = null;
    volumeStepFunc polarization_step = null;
    groupDirectionFunc group_dir = null;
    waveNormalFunc wave_normal = null;
    anisotropicSolverFunc anisotropic_solver = null;
}

VolumeProperty
volumeProperty(string name, int priority = 0) {
    VolumeProperty vp = new VolumeProperty;
    vp.name = name;
    vp.priority = priority;
    return vp;
}

// Composable Factories
VolumeProperty
RefractiveIndex(real A = 1.0, real B = 0.0, indexProfile n_profile = null, gradientProfile grad_profile = null) {
    VolumeProperty vp = new VolumeProperty;
    vp.name = "RefractiveIndex";
    vp.spatial_index = (n_profile != null);
    if (n_profile != null) {
        vp.phase_index = new real(pair p, pair d, real lambda, Polarization pol) {
            real l_um = lambda * nm_to_um;
            real base_n = A + (l_um > 0 ? B / (l_um * l_um) : 0);
            return base_n * n_profile(p, lambda);
        };
    } else if (A != 1.0 || B != 0.0) {
        vp.phase_index = new real(pair p, pair d, real lambda, Polarization pol) {
            real l_um = lambda * nm_to_um;
            return A + (l_um > 0 ? B / (l_um * l_um) : 0);
        };
    }
    if (grad_profile != null) {
        vp.phase_gradient = new pair(pair p, pair d, real lambda) {
            real l_um = lambda * nm_to_um;
            real base_n = A + (l_um > 0 ? B / (l_um * l_um) : 0);
            return base_n * grad_profile(p, lambda);
        };
    }
    return vp;
}

VolumeProperty GRIN(indexProfile n_profile, gradientProfile grad_profile = null) {
    VolumeProperty vp = volumeProperty("GRIN", 10);
    vp.spatial_index = true;
    vp.phase_modifier = new real(real n0, pair p, pair d, real lambda, Polarization pol) {
        return (n_profile != null) ? (n0 * n_profile(p, lambda)) : n0;
    };

    // Adapt the 2-argument gradientProfile(p, lambda) to the 3-argument volumeGradFunc(p, d, lambda)
    if (grad_profile != null) {
        vp.phase_gradient = new pair(pair p, pair d, real lambda) { return grad_profile(p, lambda); };
    } else {
        vp.phase_gradient = null;
    }

    return vp;
}

VolumeProperty Absorption(real alpha) {
    VolumeProperty vp = new VolumeProperty;
    vp.name = "Absorption";
    vp.attenuation = new real(pair p, pair d, real lambda, Polarization pol) { return alpha; };
    return vp;
}

VolumeProperty Gain(real g) {
    VolumeProperty vp = new VolumeProperty;
    vp.name = "Gain";
    vp.attenuation = new real(pair p, pair d, real lambda, Polarization pol) { return -g; };
    return vp;
}

// Linear dichroism.  alpha_s acts on light polarized along `axis`
// (measured in the Stokes frame, 0 = s/TE), alpha_p on the orthogonal state.
VolumeProperty Dichroic(real alpha_s, real alpha_p, real axis = 0.0) {
    VolumeProperty vp = new VolumeProperty;
    vp.name = "Dichroic";
    real axis = radians(axis);

    vp.polarization_step = new Polarization(Polarization in_pol, pair p, pair d, real lambda, real ds) {
        if (in_pol == null || ds <= 0) return in_pol;

        Polarization out = copy(in_pol);
        real ds_um = geometryToUm(ds);
        real Ts = exp(-alpha_s * ds_um);
        real Tp = exp(-alpha_p * ds_um);
        // Rotate into the dichroic eigenframe, apply the Mueller
        // diattenuator, and rotate back.
        if (axis != 0) out.applyRotation(-axis);
        applyDiattenuation(out, out, Ts, Tp);
        if (axis != 0) out.applyRotation(axis);
        return out;
    };
    return vp;
}

// Sellmeier dispersion n^2(lambda) = 1 + sum B_i lambda^2 / (lambda^2 - C_i),
// with lambda in um and C_i in um^2.  Homogeneous: no GRIN integration.
VolumeProperty Sellmeier(real[] B, real[] C) {
    if (B.length != C.length) abort("Sellmeier: B and C must have equal length.");
    VolumeProperty vp = new VolumeProperty;
    vp.name = "Sellmeier";
    real[] Bc = copy(B), Cc = copy(C);
    vp.phase_index = new real(pair p, pair d, real lambda, Polarization pol) {
        real l2 = (lambda * nm_to_um) ^ 2;
        real n2 = 1.0;
        for (int i = 0; i < Bc.length; ++i) n2 += Bc[i] * l2 / (l2 - Cc[i]);
        return sqrt(max(0.0, n2));
    };
    return vp;
}

VolumeProperty UniaxialExtraordinary(real no, real ne, pair optic_axis = (1, 0)) {
    if (no <= 0 || ne <= 0) abort("Uniaxial refractive indices must be positive.");
    VolumeProperty vp = new VolumeProperty;
    vp.name = "Uniaxial Extraordinary";
    pair optic_axis = unitOr(optic_axis);
    pair q_axis = perp(optic_axis);
    vp.phase_index = new real(pair p, pair wave_normal, real lambda, Polarization pol) {
        pair u = unitOr(wave_normal);
        real inv_n2 = dot(u, optic_axis) ^ 2 / (no * no) + dot(u, q_axis) ^ 2 / (ne * ne);
        return 1.0 / sqrt(max(config.math.eps, inv_n2));
    };
    vp.group_dir = new pair(pair wave_normal, pair p, real lambda) {
        pair k = unitOr(wave_normal);
        return unitOr(dot(k, optic_axis) * optic_axis / (no * no) + dot(k, q_axis) * q_axis / (ne * ne), k);
    };
    vp.wave_normal = new pair(pair group_dir, pair p, real lambda) {
        pair d = unitOr(group_dir);
        return unitOr((dot(d, optic_axis) * no * no) * optic_axis + (dot(d, q_axis) * ne * ne) * q_axis, d);
    };
    vp.anisotropic_solver = new AnisotropicTransmission(pair N, pair T, real kt, pair forward, pair p, real lambda) {
        AnisotropicTransmission out = new AnisotropicTransmission;

        real inv_no2 = 1.0 / (no * no);
        real inv_ne2 = 1.0 / (ne * ne);

        real eta_tt = dot(T, optic_axis) * dot(T, optic_axis) * inv_no2 + dot(T, q_axis) * dot(T, q_axis) * inv_ne2;
        real eta_nn = dot(N, optic_axis) * dot(N, optic_axis) * inv_no2 + dot(N, q_axis) * dot(N, q_axis) * inv_ne2;
        real eta_tn = (dot(N, optic_axis) * dot(T, optic_axis)) * inv_no2 + (dot(N, q_axis) * dot(T, q_axis)) * inv_ne2;

        real a = eta_nn;
        real b = 2.0 * kt * eta_tn;
        real cc = kt * kt * eta_tt - 1.0;
        real discriminant = b * b - 4.0 * a * cc;

        if (a <= config.math.eps) return out;

        // Phase 4: Exact Anisotropic TIR (no isotropic fallback!)
        if (discriminant < -config.math.eps) {
            out.valid = false;
            return out;
        }

        // Phase 3: Propagating Extraordinary Mode
        real root = sqrt(max(0.0, discriminant));
        real candidates[] = {(-b + root) / (2.0 * a), (-b - root) / (2.0 * a)};
        real best_forward = -infinity;
        real best_kn = candidates[0];

        for (real kn : candidates) {
            pair k = kn * N + kt * T;
            pair wave_normal = unitOr(k);
            pair group = unitOr(dot(wave_normal, optic_axis) * optic_axis * inv_no2 +
                                    dot(wave_normal, q_axis) * q_axis * inv_ne2,
                                wave_normal);
            real score = dot(group, forward);
            if (score > best_forward) {
                best_forward = score;
                best_kn = kn;
                out.valid = (score > config.math.eps);
                out.wave_normal = wave_normal;
                out.group_direction = group;
                out.phase_index = length(k);
            }
        }

        return out;
    };
    return vp;
}
// Optical Activity / Natural Polarization Rotation
VolumeProperty OpticalActivity(real rot_per_unit) {
    VolumeProperty vp = new VolumeProperty;
    vp.name = "OpticalActivity";

    vp.polarization_step = new Polarization(Polarization in_pol, pair p, pair d, real lambda, real ds) {
        if (in_pol == null || ds <= 0) return in_pol;
        Polarization out = copy(in_pol);

        real theta = radians(rot_per_unit * geometryToUm(ds));
        out.applyRotation(theta);

        return out;
    };
    return vp;
}

// Faraday Rotation: theta = V * (B . d) * ds (Non-reciprocal magnetic rotation)
VolumeProperty FaradayRotator(real V, pair B = (1, 0)) {
    VolumeProperty vp = new VolumeProperty;
    vp.name = "FaradayRotator";

    vp.polarization_step = new Polarization(Polarization in_pol, pair p, pair d, real lambda, real ds) {
        if (in_pol == null || ds <= 0) return in_pol;
        Polarization out = copy(in_pol);

        real theta = radians(V * dot(B, d) * geometryToUm(ds));
        out.applyRotation(theta);

        return out;
    };
    return vp;
}

// Modular Medium Struct
// Bulk Volume Physics:
// Phase Velocity / Eikonal Equation   -->  n(p, d, lambda)
// Attenuation Coefficient             -->  absorption(p, d, lambda)
// Continuous Field Evolution (ds)     -->  stepPolarization(pol, ds)
// Constitutive Relation (k <-> S)     -->  group_direction(), waveNormalFromGroup()
struct Medium {
    string name = "";
    VolumeProperty[] props;

    Medium add(VolumeProperty prop) {
        if (prop != null) this.props.push(prop);
        return this;
    }

    VolumeProperty[] orderedProperties() {
        VolumeProperty[] ordered = copy(this.props);
        for (int i = 1; i < ordered.length; ++i) {
            VolumeProperty key = ordered[i];
            int j = i - 1;
            while (j >= 0 && ordered[j].priority > key.priority) {
                ordered[j + 1] = ordered[j]; -- j;
            }
            ordered[j + 1] = key;
        }
        return ordered;
    }

    real n(pair p = (0, 0), pair d = (1, 0), real lambda = config.physics.lambda, Polarization pol = null) {
        real val = 1.0;
        bool has_index = false;

        for (VolumeProperty vp : orderedProperties()) {
            // Absolute index generators compose multiplicatively in priority order
            if (vp.phase_index != null) {
                real cur_n = vp.phase_index(p, d, lambda, pol);
                val = has_index ? (val * cur_n) : cur_n;
                has_index = true;
            }
            // Modifiers scale/perturb the current running index
            if (vp.phase_modifier != null) {
                val = vp.phase_modifier(val, p, d, lambda, pol);
                has_index = true;
            }
        }
        return val;
    }

    // Cumulative Absorption/Gain
    // The Beer–Lambert absorption coefficient alpha.
    real absorption(pair p = (0, 0), pair d = (1, 0), real lambda = config.physics.lambda, Polarization pol = null) {
        real total_alpha = 0.0;
        for (VolumeProperty vp : orderedProperties()) {
            if (vp.attenuation != null) total_alpha += vp.attenuation(p, d, lambda, pol);
        }
        return total_alpha;
    }
    // Polarization & Loss Stepping (Chains all steps + Beer-Lambert decay)
    // Evolves the Stokes vector over a distance ds (applying dichroism, optical activity/rotation, and
    // continuous absorption decay).
    Polarization stepPolarization(Polarization in_pol, pair p, pair d, real lambda, real ds) {
        if (in_pol == null || ds <= 0) return in_pol;
        Polarization curr_pol = copy(in_pol);
        for (VolumeProperty vp : orderedProperties()) {
            if (vp.polarization_step != null) {
                curr_pol = vp.polarization_step(curr_pol, p, d, lambda, ds);
            }
        }
        real a = absorption(p, d, lambda, in_pol);
        // Beer-Lambert: isotropic loss/gain scales the whole Stokes vector
        // by exp(-alpha*ds).  Scaling only S0 would leave |S_123| > S0.
        if (a != 0) curr_pol.scale(exp(-a * geometryToUm(ds)));

        return curr_pol;
    }
    // Anisotropic Direction & Normal Conversion (Pure delegation)
    // Maps wavevector k to Poynting vector S (and vice versa).
    pair groupDirection(pair k, pair p = (0, 0), real lambda = config.physics.lambda) {
        pair out = unitOr(k);
        for (VolumeProperty vp : orderedProperties()) {
            if (vp.group_dir != null) out = vp.group_dir(k, p, lambda);
        }
        return out;
    }
    pair waveNormalFromGroup(pair d, pair p = (0, 0), real lambda = config.physics.lambda) {
        pair out = unitOr(d);
        for (VolumeProperty vp : orderedProperties()) {
            if (vp.wave_normal != null) out = vp.wave_normal(d, p, lambda);
        }
        return out;
    }
}

Medium
medium(string name = "",
       real A = 1,
       real B = 0,
       real alpha = 0.0,
       indexProfile n_profile = null,
       gradientProfile grad_profile = null) {
    Medium m = new Medium;
    m.name = name;
    if (A != 1.0 || B != 0.0 || n_profile != null) m.add(RefractiveIndex(A, B, n_profile, grad_profile));
    if (alpha != 0) m.add(Absorption(alpha));
    return m;
}

Medium uniaxialExtraordinary(string name, real no, real ne, pair optic_axis = (1, 0)) {
    Medium m = new Medium;
    m.name = name;
    m.add(UniaxialExtraordinary(no, ne, optic_axis));
    return m;
}

struct Materials {
    static Medium VACUUM = medium("vacuum", 1.0);
    static Medium AIR = medium("air", 1.000293);
    static Medium WATER = medium("water", 1.324, 0.00312);
    static Medium BK7 = medium("BK7", 1.5046, 0.00420);
    static Medium FLINT = medium("flint", 1.62, 0.0125);
    static Medium SILICA = medium("silica", 1.4580, 0.00354);
}

// private
struct FresnelResult {
    // -------------------------------------------------------------------------
    // Backward-Compatible Primary Fields (Used by standard lenses, mirrors, tests)
    // -------------------------------------------------------------------------
    pair d;
    pair wave_normal;
    Medium medium;
    real T = 1.0;
    real R = 0.0;
    bool tir = false;

    // Boundary Context
    pair normal;
    Medium source_medium;
    real n1 = 1.0;
    real n2 = 1.0;
    real cosi = 1.0;
    real cost = 1.0;

    // Reflection
    pair reflected_direction;
    pair reflected_wave_normal;
    pair rs = (0, 0);
    pair rp = (0, 0);
    real Rs = 0.0;
    real Rp = 0.0;
    real Ts = 0.0;
    real Tp = 0.0;

    // -------------------------------------------------------------------------
    // Mode-Specific Fields (For Anisotropic Crystals / Birefringence)
    // -------------------------------------------------------------------------
    bool is_birefringent = false;

    // Channel 1: Ordinary / TE / S
    pair d_o;
    pair wave_normal_o;
    pair ts = (0, 0); // or 'to'
    real To = 0.0;
    bool tir_o = false;
    Medium medium_o;

    // Channel 2: Extraordinary / TM / P
    pair d_e;
    pair wave_normal_e;
    pair tp = (0, 0); // or 'te'
    real Te = 0.0;
    bool tir_e = false;
    Medium medium_e;

    // -------------------------------------------------------------------------
    // Synchronizer: Call at the end of solver to guarantee all aliases match
    // -------------------------------------------------------------------------
    void sync() {
        if (!this.is_birefringent) {
            // In isotropic media: both modes share the primary path
            this.d_o = this.d;
            this.d_e = this.d;
            this.wave_normal_o = this.wave_normal;
            this.wave_normal_e = this.wave_normal;
            this.medium_o = this.medium;
            this.medium_e = this.medium;
            this.To = this.Ts;
            this.Te = this.Tp;
            this.tir_o = this.tir;
            this.tir_e = this.tir;
        } else {
            // In crystals: primary aliases default to the ordinary mode
            this.d = this.d_o;
            this.wave_normal = this.wave_normal_o;
            this.medium = this.medium_o;
            this.tir = (this.tir_o && this.tir_e);
            this.Ts = this.To;
            this.Tp = this.Te;
            // To and Te are per-mode transmittances; the total transmittance
            // is their incident-power-weighted mean, i.e. T = 1 - R (already
            // set by the solver).  Summing To + Te would exceed unity.
        }
    }
}

// Standalone Maxwell boundary condition solver for isotropic interfaces
FresnelResult solveIsotropicInterface(real n1, real n2, pair inc_wave_normal, pair N_into2, pair T_surf) {
    FresnelResult res = new FresnelResult;
    res.is_birefringent = false;

    pair u = unit(inc_wave_normal);
    real cosi = clamp(dot(u, N_into2), 0.0, 1.0);
    real sini = dot(u, T_surf);

    real kt = n1 * sini;
    real kn1 = n1 * cosi;
    res.cosi = cosi;
    res.n1 = n1;
    res.n2 = n2;

    real sint_val = kt / n2;
    real sin2_t = sint_val * sint_val;
    res.tir = (sin2_t > 1.0);

    // Reflected and transmitted ray directions:
    res.reflected_direction = unit(sini * T_surf - cosi * N_into2);
    res.reflected_wave_normal = res.reflected_direction;

    if (res.tir) {
        real k_im = sqrt(sin2_t - 1.0);
        res.cost = 0.0;
        pair denom_s = (kn1, -n2 * k_im);
        pair denom_p = (n2 * cosi, -n1 * k_im);
        res.rs = (kn1, n2 * k_im) / denom_s;
        res.rp = (n2 * cosi, n1 * k_im) / denom_p;
        res.ts = (0.0, 0.0);
        res.tp = (0.0, 0.0);
        res.Rs = 1.0;
        res.Rp = 1.0;
        res.Ts = 0.0;
        res.Tp = 0.0;
        res.d = res.reflected_direction;
        res.wave_normal = res.reflected_wave_normal;
        res.R = 1.0;
        res.T = 0.0;
    } else {
        real cost_val = sqrt(max(0.0, 1.0 - sin2_t));
        res.cost = cost_val;
        real denom_s = kn1 + n2 * cost_val;
        real denom_p = n2 * cosi + n1 * cost_val;
        res.rs = ((kn1 - n2 * cost_val) / denom_s, 0.0);
        res.rp = ((n2 * cosi - n1 * cost_val) / denom_p, 0.0);
        res.ts = ((2.0 * kn1) / denom_s, 0.0);
        res.tp = ((2.0 * kn1) / denom_p, 0.0);
        res.Rs = res.rs.x * res.rs.x;
        res.Rp = res.rp.x * res.rp.x;
        real flux_ratio = (n2 * cost_val) / max(config.math.eps, kn1);
        res.Ts = flux_ratio * (res.ts.x * res.ts.x);
        res.Tp = flux_ratio * (res.tp.x * res.tp.x);
        res.d = unit(sint_val * T_surf + cost_val * N_into2);
        res.wave_normal = res.d;
        real R_unpol = 0.5 * (res.Rs + res.Rp);
        res.R = R_unpol;
        res.T = 1.0 - R_unpol;
    }

    res.sync();
    return res;
}

// Collision Packet
// private
struct Hit {
    bool valid = false;
    real t = infinity;
    pair p;
    pair normal;
    pair tangent;
    pair center;
    pair axis;
    bool entering;
    int flat_index = -1;
}

private Hit
copy(Hit src) {
    if (src == null) return null;
    Hit out = new Hit;
    out.valid = src.valid;
    out.t = src.t;
    out.p = src.p;
    out.normal = src.normal;
    out.tangent = src.tangent;
    out.center = src.center;
    out.axis = src.axis;
    out.entering = src.entering;
    out.flat_index = src.flat_index;
    return out;
}

struct Ray {
    int id = -1;
    int parent_id = -1;
    int generation = 0;
    int[] children_ids = new int[];
    int default_child_index = 0;
    pair p;
    // Wave-normal direction; `d` is the geometric/group direction.  They are
    // identical in isotropic media and differ for extraordinary rays.
    pair wave_normal;
    pair d;
    real lambda;
    real opl;
    real phase;
    Medium medium;
    int coherence_id = 0;
    bool alive = true;
    bool terminated = false;
    int[] ignored_flat_indices = new int[];
    Hit hit = null;
    path trajectory = nullpath;
    // REPLACED: Es, Ep, and unpolarized are now encapsulated!
    Polarization pol;
    static int Incident = 1;
    static int Reflected = 2;
    static int Transmitted = 3;
    pair point(real t) { return p + t * d; }
    real intensity() { return pol.intensity(); } // Delegates to pol!
    void initTree(int new_id, int new_parent_id, int new_generation) {
        this.id = new_id;
        this.parent_id = new_parent_id;
        this.generation = new_generation;
        this.children_ids = new int[];
        this.default_child_index = 0;
    }
    void operator init(pair p = (0, 0),
                       pair d = (1, 0),
                       Medium medium = null,
                       real lambda = config.physics.lambda,
                       Polarization pol = Polarization.unpolarized(),
                       real phase = 0,
                       int coherence_id = 0) {
        this.id = -1;
        this.parent_id = -1;
        this.generation = 0;
        this.children_ids = new int[];
        this.default_child_index = 0;
        this.p = p;
        this.d = unitOr(d);
        this.wave_normal = this.d;
        this.medium = medium;
        this.lambda = lambda;
        this.phase = phase;
        this.coherence_id = coherence_id;
        // Directly assign the incoming polarization state
        this.pol = copy(pol);
        this.opl = 0;
        this.alive = true;
        this.terminated = false;
        this.ignored_flat_indices = new int[];
        this.hit = null;
        this.trajectory = nullpath;
    }
    // STATIC FACTORY NAMESPACES
    static Ray angle(pair p,
                     real angle,
                     Medium medium = null,
                     real lambda = config.physics.lambda,
                     Polarization pol = Polarization.unpolarized(),
                     real phase = 0,
                     int coherence_id = 0) {
        return Ray(p, dir(angle), medium, lambda, pol, phase, coherence_id);
    }
    static Ray between(pair a,
                       pair b,
                       Medium medium = null,
                       real lambda = config.physics.lambda,
                       Polarization pol = Polarization.unpolarized(),
                       real phase = 0,
                       int coherence_id = 0) {
        return Ray(a, b - a, medium, lambda, pol, phase, coherence_id);
    }
    static Ray toPoint(pair target,
                       pair d,
                       real distance,
                       Medium medium = null,
                       real lambda = config.physics.lambda,
                       Polarization pol = Polarization.unpolarized(),
                       real phase = 0,
                       int coherence_id = 0) {
        pair d_hat = unitOr(d);
        return Ray(target - distance * d_hat, d_hat, medium, lambda, pol, phase, coherence_id);
    }
    static Ray toPoint(pair target,
                       real angle,
                       real distance,
                       Medium medium = null,
                       real lambda = config.physics.lambda,
                       Polarization pol = Polarization.unpolarized(),
                       real phase = 0,
                       int coherence_id = 0) {
        return toPoint(target, dir(angle), distance, medium, lambda, pol, phase, coherence_id);
    }
    // Added Phase Methods
    real totalPhase() {
        real k0 = vacuumWaveNumber(this.lambda);
        return k0 * this.opl + this.phase;
    }
} 

pair point(Ray r, real t) {
    return r.point(t);
}

// Optical path per unit geometric length along the ray (group) direction:
// k.d = n_phase * cos(walk-off).  Equals the phase index in isotropic media.
real rayIndex(Medium m, pair p, pair wave_normal, pair d, real lambda, Polarization pol = null) {
    if (m == null) return 1.0;
    pair u = unitOr(wave_normal, d);
    pair g = unitOr(d, u);
    return m.n(p, u, lambda, pol) * max(0.0, dot(u, g));
}

real rayIndex(Ray r, pair p) { return rayIndex(r.medium, p, r.wave_normal, r.d, r.lambda, r.pol); }
real rayIndex(Ray r) { return rayIndex(r, r.p); }

pair projection(pair p, Ray r) { return r.p + dot(p - r.p, r.d) * unit(r.d); }
pair projection(Ray r, pair p) { return projection(p, r); }

// Dynamic Polarization Evaluator
private
Polarization resolvePolarizationAt(Ray r, real s) {
    if (r == null) return new Polarization;
    if (s <= 0 || r.medium == null) return copy(r.pol);
    return r.medium.stepPolarization(r.pol, r.p, r.d, r.lambda, s);
}

Ray copy(Ray src) {
    if (src == null) return null;
    Ray r = new Ray;
    r.id = src.id;
    r.parent_id = src.parent_id;
    r.generation = src.generation;
    r.children_ids = copy(src.children_ids);
    r.default_child_index = src.default_child_index;
    r.p = src.p;
    r.wave_normal = src.wave_normal;
    r.d = src.d;
    r.medium = src.medium;
    r.lambda = src.lambda;
    r.phase = src.phase;
    r.opl = src.opl;
    r.pol = copy(src.pol);
    r.coherence_id = src.coherence_id;
    r.alive = src.alive;
    r.terminated = src.terminated;
    r.ignored_flat_indices = copy(src.ignored_flat_indices);
    if (src.hit != null) r.hit = copy(src.hit);
    r.trajectory = src.trajectory;
    return r;
}

Ray parallel(Ray r, pair p) {
    Ray out = copy(r);
    out.p = p;
    return out;
}

Ray reverse(Ray r, real distance = 0) {
    Ray out = copy(r);
    out.p = r.p + distance * r.d;
    out.d = -r.d;
    out.wave_normal = -r.wave_normal;
    return out;
}

// Advances ray geometrically along its direction by dist, accumulating OPL
Ray propagate(Ray r, real dist) {
    if (r == null) return null;
    Ray out = copy(r);
    real n_eff = rayIndex(r.medium, r.p, r.wave_normal, r.d, r.lambda, r.pol);
    out.p = r.p + dist * r.d;
    out.opl = r.opl + geometryToUm(dist) * n_eff;
    return out;
}

// Advances ray geometrically to a target point, accumulating OPL
Ray propagateTo(Ray r, pair target) {
    if (r == null) return null;
    real dist = length(target - r.p);
    return propagate(r, dist);
}

Ray spawnRay(Ray src,
             Hit h,
             pair direction,
             real intensity,
             Medium medium,
             real phase_shift = 0.0,
             pair wave_normal = (0, 0)) {
    Ray out = copy(src);
    out.initTree(-1, src.id, src.generation + 1);
    out.d = unitOr(direction, src.d);
    out.wave_normal = unitOr(wave_normal, out.d);
    out.p = h.p;
    out.medium = medium;
    // Scale polarization state
    out.pol.scaleTo(intensity);
    // Boundary phase shift
    out.phase += phase_shift;
    out.terminated = false;
    // KEEP EXISTING flat_index LOGIC EXACTLY AS IT IS:
    if (h.t > config.geometry.eps_intersection) {
        if (h.flat_index >= 0)
            out.ignored_flat_indices = new int[]{h.flat_index};
        else
            out.ignored_flat_indices = new int[];
    } else {
        out.ignored_flat_indices = copy(src.ignored_flat_indices);
        if (h.flat_index >= 0) out.ignored_flat_indices.push(h.flat_index);
    }
    out.hit = null;
    out.trajectory = nullpath;
    return out;
}

// ============================================================================
// Multi-Ray Coherent / Incoherent Interference
// ============================================================================

// Returns combined polarization state from interference at a detector point
Polarization interfereState(Ray[] rays) {
    if (rays.length == 0) return Polarization.unpolarized(0.0);
    if (rays.length == 1) return copy(rays[0].pol);

    int ref_coh = rays[0].coherence_id;
    real ref_lambda = rays[0].lambda;
    bool all_coherent = true;

    for (Ray r : rays) {
        if (r == null || r.pol == null || !r.pol.isPure() || r.coherence_id != ref_coh ||
            abs(r.lambda - ref_lambda) > config.math.eps) {
            all_coherent = false;
            break;
        }
    }

    if (!all_coherent) {
        real s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
        for (Ray r : rays) {
            if (r != null && r.pol != null) {
                s0 += r.pol.S0;
                s1 += r.pol.S1;
                s2 += r.pol.S2;
                s3 += r.pol.S3;
            }
        }
        return Polarization.stokes(s0, s1, s2, s3);
    }

    pair Es_tot = (0.0, 0.0);
    pair Ep_tot = (0.0, 0.0);

    for (Ray r : rays) {
        pair[] J = r.pol.toJones();
        real phi = r.totalPhase();
        pair phasor = (cos(phi), sin(phi));
        Es_tot += J[0] * phasor;
        Ep_tot += J[1] * phasor;
    }

    return Polarization.jones(Es_tot, Ep_tot);
}

// Computes detected irradiance from any number of rays meeting at a screen/detector point
real interfere(Ray[] rays) {
    if (rays.length == 0) return 0.0;
    if (rays.length == 1) return rays[0].intensity();
    return interfereState(rays).intensity();
}

// Varargs overloads allowing arbitrary positional rays: interfere(r1, r2, ...), interfereState(r1, r2, ...)
real interfere( ... Ray[] rays) { return interfere(rays); }

Polarization interfereState( ... Ray[] rays) { return interfereState(rays); }

// Pluggable ODE Integrator Engine
private
struct GRINStep {
    Ray ray;
    real error;
}

typedef GRINStep
GRINStepFunction(Ray r, real step);

struct IndexAndGradient {
    real n = 1.0;
    pair grad = (0, 0);
}

private IndexAndGradient
mediumIndexAndGradient(Medium m, pair p, pair d, real lambda, Polarization pol = null) {
    IndexAndGradient res = new IndexAndGradient;
    if (m == null) return res;

    real current_n = 1.0;
    pair gradient = (0, 0);
    bool has_index = false;
    bool has_complete_analytic_gradient = true;

    for (VolumeProperty vp : m.orderedProperties()) {
        if (vp.phase_index != null) {
            real local_index = vp.phase_index(p, d, lambda, pol);
            if (vp.phase_gradient == null) {
                if (vp.spatial_index)
                    has_complete_analytic_gradient = false;
                else
                    gradient = local_index * gradient;
            } else {
                gradient = local_index * gradient + current_n * vp.phase_gradient(p, d, lambda);
            }
            current_n = has_index ? (current_n * local_index) : local_index;
            has_index = true;
        }

        if (vp.phase_modifier != null) {
            real mod_n = vp.phase_modifier(current_n, p, d, lambda, pol);
            real factor = (current_n > config.math.eps) ? (mod_n / current_n) : 1.0;

            if (vp.phase_gradient != null) {
                gradient = factor * gradient + current_n * vp.phase_gradient(p, d, lambda);
            } else {
                if (vp.spatial_index)
                    has_complete_analytic_gradient = false;
                else
                    gradient = factor * gradient;
            }
            current_n = mod_n;
            has_index = true;
        }
    }

    res.n = has_index ? current_n : m.n(p, d, lambda, pol);

    if (!has_complete_analytic_gradient) {
        real eps = config.physics.eps_integration;
        real dx = (m.n(p + (eps, 0), d, lambda, pol) - m.n(p - (eps, 0), d, lambda, pol)) / (2.0 * eps);
        real dy = (m.n(p + (0, eps), d, lambda, pol) - m.n(p - (0, eps), d, lambda, pol)) / (2.0 * eps);
        res.grad = (dx, dy);
    } else {
        res.grad = gradient;
    }

    return res;
}

private
pair grinDirectionDerivative(Ray r, pair p, pair d) {
    pair u = (length(d) > config.math.eps) ? unit(d) : r.d;
    IndexAndGradient ig = mediumIndexAndGradient(r.medium, p, u, r.lambda, r.pol);
    if (ig.n <= config.physics.min_refractive_index) return (0, 0);
    return (ig.grad - dot(ig.grad, u) * u) / ig.n;
}

private
GRINStep grinAdvanceRK4(Ray r, real step) {
    GRINStep res = new GRINStep;
    res.ray = copy(r);
    if (step <= 0) {
        res.error = 0;
        return res;
    }
    pair p1 = r.p;
    pair d1 = unitOr(r.d);
    pair k1p = d1;
    pair k1d = grinDirectionDerivative(r, p1, d1);
    pair p2 = r.p + 0.5 * step * k1p;
    pair d2 = unitOr(r.d + 0.5 * step * k1d, d1);
    pair k2p = d2;
    pair k2d = grinDirectionDerivative(r, p2, d2);
    pair p3 = r.p + 0.5 * step * k2p;
    pair d3 = unitOr(r.d + 0.5 * step * k2d, d2);
    pair k3p = d3;
    pair k3d = grinDirectionDerivative(r, p3, d3);
    pair p4 = r.p + step * k3p;
    pair d4 = unitOr(r.d + step * k3d, d3);
    pair k4p = d4;
    pair k4d = grinDirectionDerivative(r, p4, d4);
    res.ray.p = r.p + (step / 6.0) * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);
    res.ray.d = unitOr(r.d + (step / 6.0) * (k1d + 2.0 * k2d + 2.0 * k3d + k4d), d4);
    res.ray.wave_normal = r.medium.waveNormalFromGroup(res.ray.d, res.ray.p, res.ray.lambda);
    res.error = 0.0;
    return res;
}

private
GRINStep grinAdvanceRK45(Ray r, real step) {
    GRINStep res = new GRINStep;
    res.ray = copy(r);
    if (step <= 0) {
        res.error = 0;
        return res;
    }
    pair p0 = r.p;
    pair d0 = unitOr(r.d);
    real c21 = 1.0 / 5.0;
    real c31 = 3.0 / 40.0, c32 = 9.0 / 40.0;
    real c41 = 3.0 / 10.0, c42 = -9.0 / 10.0, c43 = 1.2;
    real c51 = -11.0 / 54.0, c52 = 2.5, c53 = -70.0 / 27.0, c54 = 35.0 / 27.0;
    real c61 = 1631.0 / 55296.0, c62 = 175.0 / 512.0, c63 = 575.0 / 13824.0, c64 = 44275.0 / 110592.0,
         c65 = 253.0 / 4096.0;
    pair k1p = d0;
    pair k1d = grinDirectionDerivative(r, p0, d0);
    pair p2 = p0 + step * (c21 * k1p);
    pair d2 = unitOr(d0 + step * (c21 * k1d));
    pair k2d = grinDirectionDerivative(r, p2, d2);
    pair p3 = p0 + step * (c31 * k1p + c32 * d2);
    pair d3 = unitOr(d0 + step * (c31 * k1d + c32 * k2d));
    pair k3d = grinDirectionDerivative(r, p3, d3);
    pair p4 = p0 + step * (c41 * k1p + c42 * d2 + c43 * d3);
    pair d4 = unitOr(d0 + step * (c41 * k1d + c42 * k2d + c43 * k3d));
    pair k4d = grinDirectionDerivative(r, p4, d4);
    pair p5 = p0 + step * (c51 * k1p + c52 * d2 + c53 * d3 + c54 * d4);
    pair d5 = unitOr(d0 + step * (c51 * k1d + c52 * k2d + c53 * k3d + c54 * k4d));
    pair k5d = grinDirectionDerivative(r, p5, d5);
    pair p6 = p0 + step * (c61 * k1p + c62 * d2 + c63 * d3 + c64 * d4 + c65 * d5);
    pair d6 = unitOr(d0 + step * (c61 * k1d + c62 * k2d + c63 * k3d + c64 * k4d + c65 * k5d));
    pair k6d = grinDirectionDerivative(r, p6, d6);
    pair p_5th =
        p0 + step * ((37.0 / 378.0) * k1p + (250.0 / 621.0) * d3 + (125.0 / 594.0) * d4 + (512.0 / 1771.0) * d6);
    pair d_5th = unitOr(
        d0 + step * ((37.0 / 378.0) * k1d + (250.0 / 621.0) * k3d + (125.0 / 594.0) * k4d + (512.0 / 1771.0) * k6d));
    pair p_4th = p0 + step * ((2825.0 / 27648.0) * k1p + (18575.0 / 48384.0) * d3 + (13525.0 / 55296.0) * d4 +
                              (277.0 / 14336.0) * d5 + 0.25 * d6);
    pair d_4th = unitOr(d0 + step * ((2825.0 / 27648.0) * k1d + (18575.0 / 48384.0) * k3d + (13525.0 / 55296.0) * k4d +
                                     (277.0 / 14336.0) * k5d + 0.25 * k6d));
    res.ray.p = p_5th;
    res.ray.d = d_5th;
    res.ray.wave_normal = r.medium.waveNormalFromGroup(res.ray.d, res.ray.p, res.ray.lambda);
    res.error = length(p_5th - p_4th) + length(d_5th - d_4th);
    return res;
}

private
struct Integrator {
    string name;
    bool is_adaptive;
    GRINStepFunction advance;
    static Integrator RK4() {
        Integrator ig = new Integrator;
        ig.name = "RK4";
        ig.is_adaptive = false;
        ig.advance = grinAdvanceRK4;
        return ig;
    }
    static Integrator RK45() {
        Integrator ig = new Integrator;
        ig.name = "RK45";
        ig.is_adaptive = true;
        ig.advance = grinAdvanceRK45;
        return ig;
    }
}

// Collision Packet & Plugin Results
private struct InteractionResult {
    Ray[] rays;
    int default_index;
    bool absorbed;
}

private InteractionResult
interactionResult(Ray[] rays = new Ray[], int default_index = 0, bool absorbed = false) {
    InteractionResult res = new InteractionResult;
    res.rays = rays;
    res.default_index = default_index;
    res.absorbed = absorbed;
    return res;
}

private
struct PropagationResult {
    Ray ray;
    bool hit_surface;
    Hit hit;
    path trajectory;
}

private PropagationResult
propagationResult(Ray ray, bool hit_surface = false, Hit hit = null, path trajectory = nullpath) {
    PropagationResult res = new PropagationResult;
    res.ray = ray;
    res.hit_surface = hit_surface;
    res.hit = hit;
    res.trajectory = trajectory;
    return res;
}

struct Interaction {
    Hit hit;
    Ray incident;
    Ray[] reflected;
    Ray[] transmitted;
    real R = 0.0;
    real T = 0.0;
    pair normal() {
        if (incident == null) return hit.normal;
        return dot(incident.d, hit.normal) > 0 ? -hit.normal : hit.normal;
    }
}

private Interaction
copy(Interaction src) {
    if (src == null) return null;
    Interaction out = new Interaction;
    out.hit = copy(src.hit);
    out.R = src.R;
    out.T = src.T;
    out.incident = copy(src.incident);
    for (Ray r : src.reflected) out.reflected.push(copy(r));
    for (Ray r : src.transmitted) out.transmitted.push(copy(r));
    return out;
}

// Pure Physics Helper Functions (Decoupled!)
pair phaseGradientNormal(Ray r, Hit h) {
    pair N = unitOr(h.normal);
    return dot(r.d, N) > 0 ? -N : N;
}

private
struct InterfaceFrame {
    pair normal;
    real cosi;
    Medium source_medium;
    Medium target_medium;
}

typedef Medium
mediumResolver(pair, Ray, Medium, Medium, bool);

private
InterfaceFrame interfaceFrame(Medium m1, Medium m2, Ray r, Hit h, mediumResolver resolve);

// Incident-power-weighted reflectance for the global s/p basis.
private
real weightedReflectance(Polarization pol, real Rs, real Rp) {
    real I = (pol != null) ? pol.intensity() : 0.0;
    if (I > config.math.eps) return (Rs * pol.sPower() + Rp * pol.pPower()) / I;
    return 0.5 * (Rs + Rp);
}

private
pair thinLensRefractDir(Ray r, Hit h, real f) {
    if (abs(f) <= config.math.eps) return r.d;
    pair A = unitOr(h.axis);
    pair T = (-A.y, A.x);
    real height = dot(h.p - h.center, T);
    // Exact paraxial deflection: abs(dot(r.d, A)) handles oblique rays & bidirectional passes
    return unitOr(r.d - (height / f) * abs(dot(r.d, A)) * T, r.d);
}

private
pair diffractionOutputDir(pair N, pair T, real k_out, bool reflective) {
    real k_normal = sqrt(max(0.0, 1.0 - k_out * k_out));
    return k_out * T + (reflective ? k_normal * N : -k_normal * N);
}

typedef Ray[] rayGenerator(Ray r, Hit h, mediumResolver resolve);
typedef pair propertyAlignFunc(pair in_dir, pair target_dir, real lambda, int order_index);

// Surface actions run in stage/priority/declaration order.  Built-ins retain
// their legacy callbacks; `apply` is the preferred extension API.
struct PropertyStage {
    static int BeforeBoundary = 0;
    static int Deflect = 10;
    static int Interact = 20;
    static int AfterBoundary = 30;
}

struct Property {
    string name = "";
    int stage = PropertyStage.Interact;
    int priority = 0;
    rayGenerator generate = null; // Used by Refract, Reflect, Focus, etc.
    propertyAlignFunc align = null;
}

// Generic factory for custom surface properties
Property property(string name, rayGenerator action, int stage = PropertyStage.Interact, int priority = 0) {
    Property p = new Property;
    p.name = name;
    p.stage = stage;
    p.priority = priority;
    p.generate = action;
    return p;
}

Property copy(Property src) {
    if (src == null) return null;
    Property p = new Property;
    p.name = src.name;
    p.stage = src.stage;
    p.priority = src.priority;
    p.generate = src.generate;
    p.align = src.align;
    return p;
}

// Override a built-in's default placement without mutating the original value.
Property ordered(Property src, int stage = PropertyStage.Interact, int priority = 0) {
    Property p = copy(src);
    p.stage = stage;
    p.priority = priority;
    return p;
}

Property[] properties( ... Property[] props) { return props; }
// Surface Boundary & Physics Struct
struct Surface {
    string name;
    Property[] props; // Composable Properties
    transform T = identity();
    transform inv_T = identity();
    path p;
    pen stroke;
    pen fill = nullpen;
    pair center;
    pair axis;
    Surface add(Property p) {
        props.push(p);
        return this;
    }
    Surface properties( ... Property[] pList) {
        for (Property p : pList) props.push(p);
        return this;
    }
    Hit intersect(Ray r, real max_distance = infinity, bool ignore_origin = false) {
        Hit h = new Hit;
        h.valid = false;
        h.t = infinity;
        if (length(this.p) == 0) return h;
        // Uses the precomputed inverse! No inverse() computation during tracing!
        real ray_reach = (max_distance < infinity) ? max_distance : config.physics.far;
        pair local_p = this.inv_T * r.p;
        pair local_end = this.inv_T * (r.p + ray_reach * r.d);
        real[] hits = intersections(this.p, local_p, local_end);
        real best_time = -1;
        real best_distance = infinity;
        for (real u : hits) {
            pair pt = point(this.p, u);
            pair global_pt = this.T * pt;
            pair global_q = global_pt - r.p;
            real d = length(global_q);
            if (ignore_origin && d <= config.geometry.self_intersection_factor * config.geometry.eps_intersection) {
                continue;
            }
            if (d <= max_distance + config.geometry.eps_intersection && d < best_distance) {
                if (d <= config.geometry.near_intersection_factor * config.geometry.eps_intersection ||
                    dot(global_q, r.d) > config.geometry.eps) {
                    best_distance = d;
                    best_time = u;
                }
            }
        }
        if (best_time < 0) return h;
        h.valid = true;
        pair local_hit_pt = point(this.p, best_time);
        pair local_tangent = dir(this.p, best_time);
        pair local_normal = (local_tangent.y, -local_tangent.x);
        h.p = this.T * local_hit_pt;
        h.tangent = unitOr(this.T * (local_hit_pt + local_tangent) - h.p, local_tangent);
        h.normal = unitOr(this.T * (local_hit_pt + local_normal) - h.p, local_normal);
        h.center = this.center;
        h.axis = this.axis;
        h.t = best_distance;
        // Boundaries are CCW, so `normal` is the outward normal: the ray
        // enters the enclosed region when it travels against it.
        h.entering = dot(r.d, h.normal) < 0;
        return h;
    }
    // Ordered, composable interaction pipeline. Actions are stable-sorted by
    // stage and priority; equal keys retain declaration order.
    InteractionResult interact(Ray r, Hit h, mediumResolver resolve) {
        // If surface has no properties, ray is absorbed by default:
        if (this.props.length == 0) return interactionResult(new Ray[], 0, true);

        // Sort properties by stage and priority:
        Property[] ordered = copy(this.props);
        for (int i = 1; i < ordered.length; ++i) {
            Property key = ordered[i];
            int j = i - 1;
            while (j >= 0 && (ordered[j].stage > key.stage ||
                              (ordered[j].stage == key.stage && ordered[j].priority > key.priority))) {
                ordered[j + 1] = ordered[j]; -- j;
            }
            ordered[j + 1] = key;
        }

        // Seed the ray at the boundary.  The wave normal must be carried
        // over explicitly: for extraordinary rays it differs from r.d and
        // determines the conserved tangential wavevector.
        Ray[] active = {spawnRay(r, h, r.d, r.intensity(), r.medium, 0.0, r.wave_normal)};
        for (Property p : ordered) {
            if (active.length == 0) break;
            if (p.generate == null) continue;

            Ray[] next;
            for (Ray ray : active) {
                Ray[] produced = p.generate(ray, h, resolve);
                next.append(produced);
            }
            active = next;
        }

        // Return clean InteractionResult packet:
        // - active: the spawned child rays
        // - default_index: 0 (default ray to follow in tree)
        // - absorbed: true if all rays were consumed/extinguished
        return interactionResult(active, 0, active.length == 0);
    }
}

Surface
copy(Surface src) {
    Surface s = new Surface;
    s.name = src.name;
    s.props = copy(src.props);
    s.T = src.T;
    s.inv_T = src.inv_T;
    s.p = src.p;
    s.stroke = src.stroke;
    s.fill = src.fill;
    s.center = src.center;
    s.axis = src.axis;
    return s;
}

private
pair transformedDirection(transform T, pair origin, pair direction) {
    if (abs(direction) <= config.geometry.eps) return direction;
    return unitOr(T * (origin + direction) - T * origin, direction);
}

private
Surface transformedSurface(Surface src, transform parent_T) {
    Surface out = copy(src);
    out.T = parent_T * src.T;
    out.inv_T = inverse(out.T); // Precalculated ONCE here!
    out.center = parent_T * src.center;
    out.axis = transformedDirection(parent_T, src.center, src.axis);
    return out;
}

// Structural Laying Components (Element)
struct Element {
    string name;
    transform T = identity();
    Surface[] surfaces;
    Element[] children;
    pair axis;
    pair center;
    path boundary;
    Medium medium = null;
    Medium med_e = null;
    real focal_length = nan; // Optional visual metadata for lenses/mirrors
    pen stroke;
    pen fill = nullpen;
    bool show_info;
    void custom_draw(picture pic, Element e) = null;
    void add(Surface s) { surfaces.push(s); }
    void add(Element e) { children.push(e); }
    Element add(Property p) {
        for (Surface s : surfaces) s.add(p);
        for (Element child : children) child.add(p);
        return this;
    }
    Element properties( ... Property[] pList) {
        for (Property p : pList) this.add(p);
        return this;
    }
    // Flattens the element hierarchy using accumulated transform and native append()
    Surface[] flatten(transform parent_T = identity()) {
        Surface[] result;
        transform total_T = parent_T * this.T;
        // Transform and add own surfaces
        for (Surface s : this.surfaces) {
            result.push(transformedSurface(s, total_T));
        }
        // Recursively flatten children and append in 1 line
        for (Element child : this.children) {
            result.append(child.flatten(total_T)); // <-- Native append()
        }
        return result;
    }
}

Element
copy(Element src) {
    Element e = new Element;
    e.name = src.name;
    e.T = src.T;
    e.axis = src.axis;
    e.center = src.center;
    e.stroke = src.stroke;
    e.fill = src.fill;
    e.show_info = src.show_info;
    e.boundary = src.boundary;
    e.medium = src.medium;
    e.med_e = src.med_e;
    e.focal_length = src.focal_length;
    e.custom_draw = src.custom_draw;
    for (Surface s : src.surfaces) {
        e.surfaces.push(copy(s));
    }
    for (Element child : src.children) {
        e.children.push(copy(child));
    }
    return e;
}

Element element(path boundary,
                Property[] props,
                pen stroke = black,
                pen fill = nullpen,
                string name = "",
                pair center = (0, 0),
                pair axis = (1, 0),
                Medium medium = null,
                Medium med_e = null) {
    Surface s = new Surface;
    s.name = name;
    s.p = ensureCCW(boundary);
    s.stroke = stroke;
    s.fill = fill;
    s.props = copy(props);
    s.center = center;
    s.axis = unitOr(axis);
    Element e = new Element;
    e.name = name;
    e.stroke = stroke;
    e.fill = fill;
    e.center = center;
    e.axis = s.axis;
    e.boundary = s.p;
    e.medium = medium; // Explicitly assigned!
    e.med_e = med_e;   // Optional second mode for crystals!
    e.add(s);
    return e;
}

private
struct ElementMatch {
    Element elem = null;
    real area = infinity;
}

private void
searchEnclosing(Element[] list, pair p, transform current_T, ElementMatch best) {
    for (int i = 0; i < list.length; ++i) {
        Element e = list[i];
        transform total_T = current_T * e.T;
        // Only cyclic (closed) boundaries define a volumetric medium
        if (length(e.boundary) > 0 && cyclic(e.boundary)) {
            path b = total_T * e.boundary;
            if (inside(b, p)) {
                real a = abs(signedArea(b));
                if (a < best.area) {
                    best.area = a;
                    best.elem = e;
                }
            }
        }
        // Recurse into children with accumulated transform
        if (e.children.length > 0) {
            searchEnclosing(e.children, p, total_T, best);
        }
    }
}

private
Element findInnermostElement(Element[] elements, pair p) {
    if (elements.length == 0) return null;
    ElementMatch best = new ElementMatch;
    searchEnclosing(elements, p, identity(), best);
    return best.elem;
}

private
Medium resolveMediumAt(Element[] elements, Medium ambient, pair p, Ray r, Medium m1, Medium m2, bool from_m1) {
    // Explicit target override:
    // If the surface explicitly hardcodes an outside medium (m2 != null), use it.
    // Otherwise, null means "dynamically resolve the medium at point p".
    Medium explicit_target = from_m1 ? m2 : m1;
    if (explicit_target != null) return explicit_target;
    // System Ambient Fallback
    Medium sys_ambient = (ambient != null) ? ambient : Materials.VACUUM;
    // Dynamic spatial query
    Element inner = findInnermostElement(elements, p);
    if (inner != null && inner.medium != null) {
        // If the ray is currently in an extraordinary mode (med_e), preserve extraordinary state
        if (inner.med_e != null && r != null && r.medium == inner.med_e) {
            return inner.med_e;
        }
        return inner.medium;
    }
    // Ray is in open space -> return system ambient
    return sys_ambient;
}

InterfaceFrame interfaceFrame(Medium m1, Medium m2, Ray r, Hit h, mediumResolver resolve) {
    InterfaceFrame f;
    pair N = unitOr(h.normal);
    real cos_forward = dot(r.d, N);
    bool from_m1 = cos_forward >= 0;
    f.normal = from_m1 ? -N : N;
    f.cosi = abs(dot(unitOr(r.wave_normal, r.d), N));
    // Source medium is ALWAYS the ray's current physical medium
    f.source_medium = (r.medium != null) ? r.medium : Materials.VACUUM;
    // Target medium is resolved spatially in transmission direction
    pair target_point = h.p + config.geometry.eps_intersection * r.d;
    f.target_medium = resolve(target_point, r, m1, m2, from_m1);
    // Failsafe
    if (f.source_medium == null) f.source_medium = Materials.VACUUM;
    if (f.target_medium == null) f.target_medium = Materials.VACUUM;
    return f;
};

// private
PhaseMatch solvePhaseMatch(Medium m, pair N, pair T, real kt, pair forward, pair p, real lambda) {
    PhaseMatch pm;
    real local_n = m.n(p, forward, lambda);
    pm.n = local_n;
    // A constitutive anisotropic solver defines the dispersion relation; two
    // independent solvers cannot be composed without an explicit model.
    anisotropicSolverFunc solver = null;
    for (VolumeProperty vp : m.orderedProperties()) {
        if (vp.anisotropic_solver != null) {
            if (solver != null) abort("A Medium may have one anisotropic constitutive solver.");
            solver = vp.anisotropic_solver;
        }
    }
    if (solver != null) {
        AnisotropicTransmission at = solver(N, T, kt, forward, p, lambda);
        pm.tir = !at.valid;
        pm.wave_normal = at.wave_normal;
        pm.group_direction = at.group_direction;
        pm.n = (at.phase_index > 0) ? at.phase_index : local_n;
        return pm;
    }
    // Isotropic dispersion sphere
    real kn2 = local_n * local_n - kt * kt;
    if (kn2 < 0) {
        pm.tir = true;
        return pm;
    }
    real kn = sqrt(max(0.0, kn2));
    real sign = (dot(N, forward) >= 0) ? 1.0 : -1.0;
    pair k_vec = kt * T + (sign * kn) * N;
    pm.tir = false;
    pm.wave_normal = unitOr(k_vec, forward);
    pm.group_direction = pm.wave_normal;
    return pm;
}

// ============================================================================
// Maxwell Boundary Solver for Isotropic <-> Uniaxial Planar Interfaces
// ============================================================================
FresnelResult fresnelIsotropic(Medium m1, Medium m2, Ray r, Hit h, mediumResolver resolve, real reflectance = nan) {
    InterfaceFrame f = interfaceFrame(m1, m2, r, h, resolve);
    pair inc_wave_normal = unitOr(r.wave_normal, r.d);

    pair N_into2 = (dot(inc_wave_normal, f.normal) < 0) ? -f.normal : f.normal;
    pair T_surf = perp(N_into2);

    real n1 = f.source_medium.n(h.p, inc_wave_normal, r.lambda, r.pol);
    real n2 = max(config.physics.min_refractive_index,
                  (f.target_medium != null) ? f.target_medium.n(h.p, inc_wave_normal, r.lambda, r.pol) : 1.0);

    FresnelResult res = solveIsotropicInterface(n1, n2, inc_wave_normal, N_into2, T_surf);

    res.normal = f.normal;
    res.source_medium = f.source_medium;
    res.medium = res.tir ? f.source_medium : f.target_medium;

    res.R = isnan(reflectance) ? weightedReflectance(r.pol, res.Rs, res.Rp) : clamp(reflectance, 0.0, 1.0);
    res.T = 1.0 - res.R;

    res.sync();
    return res;
}

FresnelResult fresnelUniaxial(Ray r, Hit h, Medium med_iso, Medium med_o, Medium med_e, mediumResolver resolve) {
    FresnelResult res;
    res.is_birefringent = true;

    InterfaceFrame f = interfaceFrame(med_iso, med_o, r, h, resolve);
    pair N = f.normal;
    pair T = (-N.y, N.x);

    pair inc_wave_normal = unitOr(r.wave_normal, r.d);
    real n1 = f.source_medium.n(h.p, inc_wave_normal, r.lambda, r.pol);
    real kt = n1 * dot(inc_wave_normal, T);
    real kn_in = abs(n1 * dot(inc_wave_normal, N));

    res.normal = f.normal;
    res.source_medium = f.source_medium;
    res.n1 = n1;
    res.cosi = f.cosi;

    // Reflected kinematics
    PhaseMatch refl_match = solvePhaseMatch(f.source_medium, N, T, kt, N, h.p, r.lambda);
    res.reflected_direction = refl_match.group_direction;
    res.reflected_wave_normal = refl_match.wave_normal;

    // Transmitted kinematics
    PhaseMatch trans_o = solvePhaseMatch(med_o, N, T, kt, -N, h.p, r.lambda);
    PhaseMatch trans_e = solvePhaseMatch(med_e, N, T, kt, -N, h.p, r.lambda);

    // 1. Ordinary channel (TE)
    res.d_o = trans_o.group_direction;
    res.wave_normal_o = trans_o.wave_normal;
    res.tir_o = trans_o.tir;
    res.medium_o = med_o;

    if (!res.tir_o) {
        real qo = abs(trans_o.n * dot(trans_o.wave_normal, N));
        real denom_s = kn_in + qo;
        res.rs = ((kn_in - qo) / denom_s, 0.0);
        res.ts = ((2.0 * kn_in) / denom_s, 0.0);
        res.Rs = res.rs.x * res.rs.x;
        res.To = 1.0 - res.Rs;
    } else {
        real no = med_o.n(h.p, inc_wave_normal, r.lambda, r.pol);
        real qo_im = sqrt(max(0.0, kt * kt - no * no));
        res.rs = (kn_in, qo_im) / (kn_in, -qo_im);
        res.ts = (0.0, 0.0);
        res.Rs = 1.0;
        res.To = 0.0;
    }

    // 2. Extraordinary channel (TM)
    res.d_e = trans_e.group_direction;
    res.wave_normal_e = trans_e.wave_normal;
    res.tir_e = trans_e.tir;
    res.medium_e = med_e;

    if (!res.tir_e) {
        real dn = abs(dot(trans_e.group_direction, N));
        real we_dot_de = abs(dot(trans_e.n * trans_e.wave_normal, trans_e.group_direction));
        real Ze = (we_dot_de > config.math.eps) ? (dn / we_dot_de) : (1.0 / trans_e.n);
        real Zi = kn_in / (n1 * n1);

        real denom_p = Zi + Ze;
        res.rp = (denom_p > config.math.eps) ? ((Zi - Ze) / denom_p, 0.0) : (-1.0, 0.0);
        res.tp = (denom_p > config.math.eps) ? ((2.0 * Zi) / denom_p, 0.0) : (0.0, 0.0);
        res.Rp = res.rp.x * res.rp.x;
        res.Te = 1.0 - res.Rp;
    } else {
        real ne = med_e.n(h.p, inc_wave_normal, r.lambda, r.pol);
        real sint_e = kt / ne;
        real sin2t = sint_e * sint_e;

        if (sin2t > 1.0) {
            real k_im = sqrt(sin2t - 1.0);
            pair cosi = (f.cosi, 0.0);
            pair cost = (0.0, -k_im);
            pair denom_p = ne * cosi + n1 * cost;
            res.rp = (abs(denom_p) > config.math.eps) ? (ne * cosi - n1 * cost) / denom_p : (1.0, 0.0);
        } else {
            res.rp = (-1.0, 0.0);
        }

        res.tp = (0.0, 0.0);
        res.Rp = 1.0;
        res.Te = 0.0;
    }

    res.R = weightedReflectance(r.pol, res.Rs, res.Rp);
    res.T = 1.0 - res.R;

    res.sync();
    return res;
}

pair boundaryCenter(path boundary) {
    if (length(boundary) == 0) return (0, 0);
    if (cyclic(boundary)) return 0.5 * (min(boundary) + max(boundary));
    return 0.5 * (point(boundary, 0) + point(boundary, length(boundary)));
}

pair boundaryAxis(path boundary) {
    if (length(boundary) == 0) return (1, 0);
    pair chord = point(boundary, length(boundary)) - point(boundary, 0);
    if (abs(chord) > config.geometry.eps) return perp(chord);
    pair tangent = dir(boundary, 0);
    if (abs(tangent) > config.geometry.eps) return perp(tangent);
    return (1, 0);
}

path linearBoundary(pair center, pair axis, real diameter) {
    pair A = unitOr(axis);
    pair q = (-A.y, A.x);
    return center - diameter / 2 * q -- center + diameter / 2 * q;
}

// STATELESS FLAT SURFACE PREPARATION
struct FlatSurfaceRecord {
    Surface s;
    transform inv_T; // Precomputed inverse!
    pair min_bound;
    pair max_bound;
}

private FlatSurfaceRecord[] prepareFlatSurfaces(Element[] elements) {
    FlatSurfaceRecord[] records;
    for (int i = 0; i < elements.length; ++i) {
        Surface[] flat = elements[i].flatten();
        for (int j = 0; j < flat.length; ++j) {
            FlatSurfaceRecord rec = new FlatSurfaceRecord;
            rec.s = flat[j];
            rec.inv_T = flat[j].inv_T; // Already inverted once by transformedSurface()
            if (length(flat[j].p) > 0) {
                path global_path = flat[j].T * flat[j].p;
                rec.min_bound = min(global_path);
                rec.max_bound = max(global_path);
            }
            records.push(rec);
        }
    }
    return records;
}

Element operator*(transform t, Element e) {
    Element out = copy(e);
    out.T = t * e.T;
    out.center = t * e.center;
    out.axis = transformedDirection(t, e.center, e.axis);
    return out;
}

// Propagation Models
// private
struct Segment {
    pair a;
    pair b;
    Ray before;
    Ray after;
    bool hit;
    path traj;
}

// private
Segment copy(Segment src) {
    if (src == null) return null;
    Segment s = new Segment;
    s.a = src.a;
    s.b = src.b;
    s.before = copy(src.before);
    s.after = copy(src.after);
    s.hit = src.hit;
    s.traj = src.traj;
    return s;
}

// 2D Ray-AABB (Slab Method) Intersection Test
private
bool rayIntersectsAABB(pair ro, pair rd, pair bmin, pair bmax, real max_dist) {
    real eps = config.geometry.eps_intersection;
    // X-axis slab
    real inv_dx = abs(rd.x) > 1e-12 ? 1.0 / rd.x : (rd.x >= 0 ? 1e12 : -1e12);
    real tx1 = (bmin.x - eps - ro.x) * inv_dx;
    real tx2 = (bmax.x + eps - ro.x) * inv_dx;
    real tmin = min(tx1, tx2);
    real tmax = max(tx1, tx2);
    // Y-axis slab
    real inv_dy = abs(rd.y) > 1e-12 ? 1.0 / rd.y : (rd.y >= 0 ? 1e12 : -1e12);
    real ty1 = (bmin.y - eps - ro.y) * inv_dy;
    real ty2 = (bmax.y + eps - ro.y) * inv_dy;
    tmin = max(tmin, min(ty1, ty2));
    tmax = min(tmax, max(ty1, ty2));
    // Must enter before exiting, exit in front of ray, and enter within reach
    return (tmax >= 0.0) && (tmin <= tmax) && (tmin <= max_dist + eps);
}

private
Hit nearestIntersection(Ray r, FlatSurfaceRecord[] flat_surfaces, real max_distance = infinity) {
    Hit best_hit = new Hit;
    best_hit.valid = false;
    best_hit.t = infinity;
    for (int i = 0; i < flat_surfaces.length; ++i) {
        FlatSurfaceRecord rec = flat_surfaces[i];
        // Fast Ray-Box Slab Test (With Distance Pruning)
        // If a hit was already found at distance best_hit.t, reject any surface further away!
        real current_reach = min(max_distance, best_hit.t);
        if (!rayIntersectsAABB(r.p, r.d, rec.min_bound, rec.max_bound, current_reach)) {
            continue; // Skipped instantly without touching Bezier math!
        }
        // Topological origin check (self-intersection prevention)
        bool ignore_origin = false;
        for (int idx : r.ignored_flat_indices) {
            if (i == idx) {
                ignore_origin = true;
                break;
            }
        }
        // Exact Bézier Spline Intersection (Only evaluated when ray actually hits the box)
        Hit h = rec.s.intersect(r, max_distance, ignore_origin);
        if (h.valid && h.t < best_hit.t) {
            best_hit = h;
            best_hit.flat_index = i;
        }
    }
    return best_hit;
}

// BEER-LAMBERT EXPONENTIAL DECAY
private
void applyAbsorption(Ray r, real dist) {
    if (r == null || !r.alive || dist <= 0 || r.medium == null) return;
    r.pol = r.medium.stepPolarization(r.pol, r.p, r.d, r.lambda, dist);
}

private
bool requiresNumericalIntegration(Medium m) {
    if (m == null) return false;
    for (VolumeProperty vp : m.props) {
        if (vp.spatial_index || vp.phase_gradient != null) return true;
    }
    return false;
}

private
PropagationResult propagate(Medium medium, Ray r, FlatSurfaceRecord[] flat_surfaces, real far, real ds) {
    if (requiresNumericalIntegration(medium)) {
        real t_total = 0;
        Ray curr_ray = copy(r);
        guide traj = curr_ray.p;
        Integrator solver = (config.physics.grin_solver == GRINSolverType.RK4) ? Integrator.RK4() : Integrator.RK45();
        real max_ds = 0.075, min_ds = 1e-5;
        real current_ds = max(min_ds, min(ds, max_ds));
        real tol = config.physics.eps_integration;
        Ray test_ray = new Ray;
        while (t_total < far && curr_ray.alive) {
            real n_center = curr_ray.medium.n(curr_ray.p, curr_ray.d, curr_ray.lambda, curr_ray.pol);
            if (n_center <= config.physics.min_refractive_index) {
                curr_ray.alive = false;
                break;
            }
            real step = min(current_ds, far - t_total);
            GRINStep step_result = solver.advance(curr_ray, step);
            if (solver.is_adaptive) {
                if (step_result.error > tol && step > min_ds) {
                    current_ds *= 0.5;
                    continue;
                }
                if (step_result.error > 0) {
                    real scale_factor = 0.9 * (tol / step_result.error) ^ 0.2;
                    current_ds = min(max_ds, current_ds * max(0.2, min(2.0, scale_factor)));
                } else {
                    current_ds = min(max_ds, current_ds * 2.0);
                }
            }
            Ray next_ray = step_result.ray;
            applyAbsorption(next_ray, step);
            if (next_ray.intensity() <= config.physics.min_intensity) {
                next_ray.alive = false;
                traj = traj -- next_ray.p;
                return propagationResult(next_ray, false, null, traj);
            }
            pair chord = next_ray.p - curr_ray.p;
            real chord_length = length(chord);
            Hit h;
            bool hit_valid = false;
            // Test every non-degenerate chord: skipping short chords (as an
            // integration-tolerance threshold did) lets tiny adaptive steps
            // tunnel through boundaries.
            if (chord_length > config.geometry.eps) {
                test_ray.p = curr_ray.p;
                test_ray.d = chord / chord_length;
                test_ray.ignored_flat_indices = copy(curr_ray.ignored_flat_indices);
                h = nearestIntersection(test_ray, flat_surfaces, chord_length);
                hit_valid = h.valid;
            }
            if (hit_valid) {
                real hit_step = step * (h.t / chord_length);
                Ray incident_ray = solver.advance(curr_ray, min(step, hit_step)).ray;
                incident_ray.p = h.p;
                // Midpoint-rule OPL over the partial step (same order as the
                // interior steps); the index at the step start biased it.
                real n_hit = curr_ray.medium.n(0.5 * (curr_ray.p + h.p),
                                               unitOr(curr_ray.d + incident_ray.d, curr_ray.d),
                                               curr_ray.lambda,
                                               curr_ray.pol);
                incident_ray.opl = curr_ray.opl + geometryToUm(hit_step) * n_hit;
                applyAbsorption(incident_ray, hit_step);
                traj = traj -- h.p;
                return propagationResult(incident_ray, true, h, traj);
            } else {
                pair mid_point = 0.5 * (curr_ray.p + next_ray.p);
                pair mid_direction = unitOr(curr_ray.d + next_ray.d, curr_ray.d);
                real n_mid = curr_ray.medium.n(mid_point, mid_direction, curr_ray.lambda, curr_ray.pol);
                curr_ray = next_ray;
                curr_ray.opl += geometryToUm(step) * n_mid;
                t_total += step;
                traj = traj -- curr_ray.p;
                curr_ray.ignored_flat_indices = new int[];
            }
        }
        return propagationResult(curr_ray, false, null, traj);
    } else {
        Hit h = nearestIntersection(r, flat_surfaces, far);
        // OPL accumulates k.dr along the group direction; for
        // extraordinary rays this is n_phase * cos(walk-off), not n_phase.
        real current_n = rayIndex(medium, r.p, r.wave_normal, r.d, r.lambda, r.pol);
        Ray next_ray = copy(r);
        if (!h.valid || h.t > far) {
            next_ray.opl += geometryToUm(far) * current_n;
            next_ray.p = r.p + far * r.d;
            applyAbsorption(next_ray, far);
            return propagationResult(next_ray, false, null, r.p -- next_ray.p);
        } else {
            next_ray.opl += geometryToUm(h.t) * current_n;
            next_ray.p = h.p;
            applyAbsorption(next_ray, h.t);
            return propagationResult(next_ray, true, h, r.p -- h.p);
        }
    }
}

// Geometric length a ray segment has travelled (arc length of its
// trajectory, falling back to the straight hit distance).
private
real travelledLength(Ray r) {
    if (r == null) return 0;
    if (length(r.trajectory) > 0) return arclength(r.trajectory);
    if (r.hit != null && r.hit.valid) return r.hit.t;
    return 0;
}

// The "Navigator" Trace Facade & Ray Tree
struct Trace {
    Ray[] ray_tree;
    int root_id = 0;
    int[] route;
    int[] root_ids;
    bool is_reflection = false;
    int[] activeBranchRayIds() {
        int[] branch;
        if (ray_tree.length == 0) return branch;
        int curr_id = root_id;
        if (curr_id < 0 || curr_id >= ray_tree.length) return branch;
        int route_step = 0;
        while (true) {
            branch.push(curr_id);
            Ray r = ray_tree[curr_id];
            if (r.children_ids.length == 0) break;
            int child_choice = r.default_child_index;
            if (route_step < route.length) {
                child_choice = route[route_step];
                ++route_step;
            }
            if (child_choice < 0 || child_choice >= r.children_ids.length) {
                child_choice = min(max(child_choice, 0), r.children_ids.length - 1);
            }
            curr_id = r.children_ids[child_choice];
        }
        return branch;
    }
    Segment[] segments() {
        Segment[] segs;
        int[] branch = activeBranchRayIds();
        for (int k = 0; k < branch.length; ++k) {
            Ray r = ray_tree[branch[k]];
            Segment seg = new Segment;
            seg.before = r;
            if (length(r.trajectory) > 0) {
                seg.traj = r.trajectory;
                seg.a = point(r.trajectory, 0);
                seg.b = point(r.trajectory, length(r.trajectory));
            } else {
                seg.a = r.p;
                seg.b = (r.hit != null && r.hit.valid) ? r.hit.p : r.p;
                seg.traj = seg.a -- seg.b;
            }
            seg.hit = (r.hit != null && r.hit.valid);
            if (k + 1 < branch.length) {
                seg.after = ray_tree[branch[k + 1]];
            } else {
                Ray after_ray = copy(r);
                after_ray.p = seg.b;
                if (r.terminated) after_ray.alive = false;
                seg.after = after_ray;
            }
            segs.push(seg);
        }
        return segs;
    }
    Interaction[] interactions() {
        Interaction[] ixs;
        int[] branch = activeBranchRayIds();
        for (int k = 0; k < branch.length; ++k) {
            Ray r = ray_tree[branch[k]];
            if (r.hit != null && r.hit.valid) {
                Interaction ix = new Interaction;
                ix.hit = r.hit;
                ix.incident = copy(r);
                ix.incident.p = r.hit.p;
                if (length(r.trajectory) > 0) {
                    ix.incident.d = dir(r.trajectory, length(r.trajectory));
                }
                // R and T are fractions of the power *arriving* at the
                // surface, i.e. after bulk absorption/gain along the segment.
                ix.incident.pol = resolvePolarizationAt(r, travelledLength(r));
                real parent_i = ix.incident.intensity();
                pair N = ix.normal();
                for (int c_id : r.children_ids) {
                    Ray c_ray = ray_tree[c_id];
                    if (dot(c_ray.d, N) > 0) {
                        ix.reflected.push(c_ray);
                        if (parent_i > config.math.eps) ix.R += c_ray.intensity() / parent_i;
                    } else {
                        ix.transmitted.push(c_ray);
                        if (parent_i > config.math.eps) ix.T += c_ray.intensity() / parent_i;
                    }
                }
                ixs.push(ix);
            }
        }
        return ixs;
    }
    int length() { return interactions().length; }
    Interaction interaction(int ix_index) {
        Interaction[] ixs = interactions();
        if (ixs.length == 0) return null;
        int idx = (ix_index < 0) ? ix_index + ixs.length : ix_index;
        if (idx < 0 || idx >= ixs.length) return null;
        return ixs[idx];
    }
    Ray incident(int ix_index) {
        Interaction ix = interaction(ix_index);
        return (ix != null) ? ix.incident : null;
    }
    Ray transmitted(int ix_index, int ray_index = 0) {
        Interaction ix = interaction(ix_index);
        if (ix == null || ix.transmitted.length == 0) return null;
        int idx = (ray_index < 0) ? ray_index + ix.transmitted.length : ray_index;
        if (idx < 0 || idx >= ix.transmitted.length) return null;
        return ix.transmitted[idx];
    }
    Ray[] Transmitted(int ix_index) {
        Interaction ix = interaction(ix_index);
        if (ix == null) return new Ray[];
        return ix.transmitted;
    }
    Ray reflected(int ix_index, int ray_index = 0) {
        Interaction ix = interaction(ix_index);
        if (ix == null || ix.reflected.length == 0) return null;
        int idx = (ray_index < 0) ? ray_index + ix.reflected.length : ray_index;
        if (idx < 0 || idx >= ix.reflected.length) return null;
        return ix.reflected[idx];
    }
    Ray[] Reflected(int ix_index) {
        Interaction ix = interaction(ix_index);
        if (ix == null) return new Ray[];
        return ix.reflected;
    }
    // GETTERS & OPERATORS
    pair hit(int i) {
        Interaction ix = interaction(i);
        return (ix != null) ? ix.hit.p : (0, 0);
    }
    pair X(int i) { return hit(i); } // Physical crossing alias
    pair normal(int i) {
        Interaction ix = interaction(i);
        return (ix != null) ? ix.normal() : (0, 0);
    }
    path normal(int i, real before, real after) {
        Interaction ix = interaction(i);
        if (ix == null) return nullpath;
        pair p = ix.hit.p;
        pair n = unit(ix.normal());
        return (p - n * after) -- (p + n * before);
    }
    Ray last() {
        int[] branch = activeBranchRayIds();
        if (branch.length == 0) return null;
        Ray r = copy(ray_tree[branch[branch.length - 1]]);
        if (r.hit != null && r.hit.valid) {
            r.p = r.hit.p;
            // hit.t is only the final chord for curved (GRIN) segments;
            // attenuation must use the full travelled arc length.
            r.pol = resolvePolarizationAt(r, travelledLength(r));
            // If the ray terminated at the collision (e.g. Absorber / Aperture screen):
            if (r.terminated || !r.alive) {
                r.alive = false;
                r.pol.scaleTo(0.0); // Post-absorption intensity is 0!
            }
        } else if (length(r.trajectory) > 0) {
            real end_t = length(r.trajectory);
            real total_l = arclength(r.trajectory);
            r.p = point(r.trajectory, end_t);
            r.d = dir(r.trajectory, end_t);
            r.wave_normal = r.medium.waveNormalFromGroup(r.d, r.p, r.lambda);
            r.pol = resolvePolarizationAt(r, total_l);
            if (r.terminated || !r.alive) {
                r.alive = false;
                r.pol.scaleTo(0.0);
            }
        }
        return r;
    }
    Ray first() {
        if (root_id >= 0 && root_id < ray_tree.length) {
            return copy(ray_tree[root_id]);
        }
        return null;
    }
    real opl() {
        Ray r = last();
        return (r != null) ? r.opl : 0.0;
    }
    Trace[] branches() {
        Trace[] b;
        if (root_ids.length > 1) {
            for (int r_id : root_ids) {
                Trace sub = new Trace;
                sub.ray_tree = ray_tree;
                sub.root_id = r_id;
                sub.route = copy(route);
                sub.root_ids = new int[];
                b.push(sub);
            }
        } else {
            int[] current_branch = activeBranchRayIds();
            int route_step = 0;
            for (int k = 0; k < current_branch.length; ++k) {
                Ray r = ray_tree[current_branch[k]];
                if (r.children_ids.length == 0) continue;
                int followed_index = r.default_child_index;
                bool was_explicit_step = (route_step < route.length);
                if (was_explicit_step) {
                    followed_index = route[route_step];
                    ++route_step;
                    continue; // Skip spawning side-branches for explicit follow() steps
                }
                for (int i = 0; i < r.children_ids.length; ++i) {
                    if (i == followed_index) continue;
                    Trace sub = new Trace;
                    sub.ray_tree = ray_tree;
                    sub.root_id = r.children_ids[i];
                    sub.route = new int[];
                    Ray child = ray_tree[r.children_ids[i]];
                    if (r.hit != null && r.hit.valid) {
                        pair N = r.hit.normal;
                        sub.is_reflection = (dot(child.d, N) * dot(r.d, N) < 0);
                    } else {
                        sub.is_reflection = false;
                    }
                    b.push(sub);
                }
            }
        }
        return b;
    }
    Trace follow( ... int[] newRoute) {
        Trace out = new Trace;
        out.ray_tree = ray_tree;
        out.root_id = root_id;
        out.route = copy(newRoute);
        out.root_ids = copy(root_ids);
        out.is_reflection = is_reflection;
        return out;
    }
    Trace subtrace(int i) {
        if (root_ids.length == 0) return this;
        int idx = (i < 0) ? i + root_ids.length : i;
        if (idx < 0 || idx >= root_ids.length) return null;
        Trace sub = new Trace;
        sub.ray_tree = ray_tree;
        sub.root_id = root_ids[idx];
        sub.route = copy(route);
        sub.root_ids = new int[];
        return sub;
    }
    if (VERSION < "3.11") {
        Trace get(int i) { return subtrace(i); }
    } else {
        // Native Subscript Operator: t[i] returns Trace view for i-th ray in a beam!
        Trace operator[](int i) { return subtrace(i); }
    }
} 

Trace copy(Trace src) {
    if (src == null) return null;
    Trace out = new Trace;
    out.ray_tree = src.ray_tree;
    out.root_id = src.root_id;
    out.route = copy(src.route);
    out.root_ids = copy(src.root_ids);
    out.is_reflection = src.is_reflection;
    return out;
}

// The Tracer Engine (Queue-Based BFS)
private
struct Tracer {
    FlatSurfaceRecord[] flat_surfaces;
    mediumResolver resolve;
    int max_depth = config.physics.max_depth;
    real min_intensity = config.physics.min_intensity;
    int branch_budget = config.physics.max_branches;
    real ds = config.physics.ds;
    real far = config.physics.far;
    Trace traceRays(Ray[] rays) {
        Trace out = new Trace;
        if (rays.length == 0) return out;
        Ray[] combined_ray_tree;
        int[] root_ids;
        for (int r_idx = 0; r_idx < rays.length; ++r_idx) {
            int start_tree_offset = combined_ray_tree.length;
            Ray root = copy(rays[r_idx]);
            root.initTree(start_tree_offset, -1, 0);
            if (root.medium == null) {
                root.medium = resolve(root.p, root, null, null, true);
            }
            combined_ray_tree.push(root);
            root_ids.push(start_tree_offset);
            int[] q = {start_tree_offset};
            int head = 0;
            int branches_count = 1;
            while (head < q.length) {
                int curr_id = q[head];
                ++head;
                Ray r = combined_ray_tree[curr_id];
                if (!r.alive || r.terminated || r.intensity() <= min_intensity || r.generation >= max_depth) continue;
                // 5 parameters (no sceneDomain)
                // 1. Propagate ray to boundary or far reach
                PropagationResult prop_res = propagate(r.medium, r, flat_surfaces, far, ds);
                r.trajectory = prop_res.trajectory;
                r.opl = prop_res.ray.opl;

                // 2. Handle Surface Interaction
                if (prop_res.hit_surface) {
                    r.hit = prop_res.hit;
                    Surface hit_surf = flat_surfaces[prop_res.hit.flat_index].s;

                    // Execute surface physics and receive InteractionResult:
                    InteractionResult ix_res = hit_surf.interact(prop_res.ray, prop_res.hit, resolve);
                    r.default_child_index = ix_res.default_index;

                    bool accepted_child = false;
                    if (!ix_res.absorbed && ix_res.rays.length > 0) {
                        for (int i = 0; i < ix_res.rays.length; ++i) {
                            Ray child = ix_res.rays[i];
                            if (child.intensity() > min_intensity && r.generation + 1 <= max_depth &&
                                (branch_budget <= 0 || branches_count < branch_budget)) {
                                accepted_child = true;
                                child.initTree(combined_ray_tree.length, curr_id, r.generation + 1);
                                r.children_ids.push(child.id);
                                combined_ray_tree.push(child);
                                q.push(child.id);
                                ++branches_count;
                            }
                        }
                    }
                    r.terminated = !accepted_child;
                }
                combined_ray_tree[curr_id] = r;
            }
        }
        out.ray_tree = combined_ray_tree;
        out.root_id = (rays.length == 1) ? root_ids[0] : -1;
        out.route = new int[];
        out.root_ids = root_ids;
        return out;
    }
}

private Tracer
makeTracer(FlatSurfaceRecord[] flat_surfaces, mediumResolver resolve, real ds, real far, int max_depth) {
    Tracer tracer = new Tracer;
    tracer.flat_surfaces = flat_surfaces;
    tracer.resolve = resolve;
    tracer.ds = ds;
    tracer.far = far;
    tracer.max_depth = max_depth;
    tracer.min_intensity = config.physics.min_intensity;
    tracer.branch_budget = config.physics.max_branches;
    return tracer;
}

// STATELESS OPTICAL SYSTEM
// ============================================================================
// OpticalSystem
// ============================================================================
// STATELESS OPTICAL SYSTEM
struct OpticalSystem {
    Element[] elements;
    Medium ambient = Materials.VACUUM;
    void add(Element e) { this.elements.push(e); }
    void add(Element[] eList) { this.elements.append(eList); }
    void add(OpticalSystem s) { this.add(s.elements); }
    int length() { return elements.length; }
    // STATELESS NATIVE TRACING
    Trace trace(Ray[] rays,
                real ds = config.physics.ds,
                real far = config.physics.far,
                int max_depth = config.physics.max_depth) {
        FlatSurfaceRecord[] flat_surfaces = prepareFlatSurfaces(elements);
        Element[] trace_elements = elements;
        Medium trace_ambient = ambient;
        mediumResolver resolve = new Medium(pair p, Ray r, Medium m1, Medium m2, bool from_m1) {
            return resolveMediumAt(trace_elements, trace_ambient, p, r, m1, m2, from_m1);
        };
        return makeTracer(flat_surfaces, resolve, ds, far, max_depth).traceRays(rays);
    }
    Trace trace(Ray ray,
                real ds = config.physics.ds,
                real far = config.physics.far,
                int max_depth = config.physics.max_depth) {
        return trace(new Ray[]{ray}, ds = ds, far = far, max_depth = max_depth);
    }
}

OpticalSystem
system(Element element, Medium ambient = Materials.VACUUM) {
    OpticalSystem s = new OpticalSystem;
    s.ambient = ambient;
    s.add(element);
    return s;
}

OpticalSystem system(Element[] elements, Medium ambient = Materials.VACUUM) {
    OpticalSystem s = new OpticalSystem;
    s.ambient = ambient;
    s.add(elements);
    return s;
}

OpticalSystem copy(OpticalSystem src) {
    if (src == null) return null;
    OpticalSystem out = new OpticalSystem;
    out.ambient = src.ambient;
    for (Element e : src.elements) {
        out.elements.push(copy(e));
    }
    return out;
}

OpticalSystem operator*(transform t, OpticalSystem sys) {
    OpticalSystem out = new OpticalSystem;
    out.ambient = sys.ambient;
    for (int i = 0; i < sys.elements.length; ++i) out.add(t * sys.elements[i]);
    return out;
}

// Segment rule shared by RayArrows and RayPolarization.
struct SegmentRule {
    int[] indices;
    // Common visibility override
    // true    = visible
    // false   = invisible
    // default = inherit global setting
    bool3 visible = default;
    // RayArrows
    // first/middle/last refer to the segment's position in the
    // complete ray.
    // RayPolarization
    // first/last refer to the first/last polarization marker
    // within this segment.
    bool3 first = default;
    bool3 middle = default;
    bool3 last = default;
    arrowhead head = null;
    real size = nan;
    real relative = nan;
    int count = -1;
    real spacing = nan;
    real before = nan;
    real after = nan;
    real scale = nan;
    bool3 show_s = default;
    bool3 show_p = default;
    pen p = nullpen;
}

SegmentRule
copy(SegmentRule src) {
    if (src == null) return null;
    SegmentRule out = new SegmentRule;
    out.indices = copy(src.indices);
    out.visible = src.visible;
    out.first = src.first;
    out.middle = src.middle;
    out.last = src.last;
    out.head = src.head;
    out.size = src.size;
    out.relative = src.relative;
    out.count = src.count;
    out.spacing = src.spacing;
    out.before = src.before;
    out.after = src.after;
    out.scale = src.scale;
    out.show_s = src.show_s;
    out.show_p = src.show_p;
    out.p = src.p;
    return out;
}

// segment()
// For RayArrows:
//   visible
//   first / middle / last
//   head / size / relative
// For RayPolarization:
//   visible
//   first / last
//   count / spacing / before / after
//   scale / show_s / show_p / p
//   head / size
// Unused fields are simply ignored by the corresponding resolver.
SegmentRule segment(bool3 visible = default,
                    explicit real size = nan,
                    explicit real relative = nan,
                    explicit real spacing = nan,
                    explicit real before = nan,
                    explicit real after = nan,
                    explicit real scale = nan,
                    bool3 keyword first = default,
                    bool3 keyword middle = default,
                    bool3 keyword last = default,
                    bool3 keyword show_s = default,
                    bool3 keyword show_p = default,
                    pen keyword p = nullpen,
                    int keyword count = -1,
                    arrowhead keyword head = null, ... int[] indices) {
    SegmentRule out = new SegmentRule;
    out.indices = copy(indices);
    // Common
    out.visible = visible;
    out.first = first;
    out.last = last;
    out.head = head;
    out.size = size;
    // Arrows
    out.middle = middle;
    out.relative = relative;
    // Polarization
    out.count = count;
    out.spacing = spacing;
    out.before = before;
    out.after = after;
    out.scale = scale;
    out.show_s = show_s;
    out.show_p = show_p;
    out.p = p;
    return out;
}

struct ResolvedArrow {
    bool visible = false;
    arrowhead head = config.rendering.arrow_head;
    real size = config.rendering.arrow_size;
    real relative = nan;
}

arrowbar
operator cast(ResolvedArrow a) {
    if (!a.visible) return None;
    real pos = isnan(a.relative) ? 0.5 : a.relative;
    arrowhead h = (a.head != null) ? a.head : config.rendering.arrow_head;
    real sz = isnan(a.size) ? config.rendering.arrow_size : a.size;
    return Arrow(h, sz, position = Relative(pos));
}

struct RayArrows {
    // Global visibility
    bool visible = true;
    // Global segment policy
    // These refer to segments of the ray.
    bool first = true;
    bool middle = true;
    bool last = true;
    // Global arrow appearance
    arrowhead head = config.rendering.arrow_head;
    real size = config.rendering.arrow_size;
    // Global arrow position
    real relative = nan;
    // Segment rules
    SegmentRule[] rules;
    // Resolve arrow configuration for one segment
    ResolvedArrow getArrow(int i, int total_segments) {
        ResolvedArrow out = new ResolvedArrow;
        bool is_first = (i == 0);
        bool is_last = (i == total_segments - 1);
        // Start with global visibility + positional policy.
        out.visible = this.visible && (is_last ? this.last : (is_first ? this.first : this.middle));
        // Global appearance.
        out.head = this.head;
        out.size = this.size;
        // Global position.
        // If no global relative position is specified:
        //   last segment -> end
        //   other        -> middle
        out.relative = isnan(this.relative) ? (is_last ? 1.0 : 0.5) : this.relative;
        // Apply matching segment rules.
        // Positional settings are applied first.
        // Explicit segment visibility is applied last.
        // Therefore:
        //   segment(visible = true, ...)
        // can force a segment on, while:
        //   segment(visible = false, ...)
        // can force it off.
        for (SegmentRule rule : rules) {
            for (int idx : rule.indices) {
                int resolved_idx = (idx < 0) ? idx + total_segments : idx;
                if (resolved_idx == i) {
                    // Segment positional policy.
                    if (is_first && rule.first != default) {
                        out.visible = (bool)rule.first;
                    }
                    if (!is_first && !is_last && rule.middle != default) {
                        out.visible = (bool)rule.middle;
                    }
                    if (is_last && rule.last != default) {
                        out.visible = (bool)rule.last;
                    }
                    // Explicit visibility is the final override.
                    if (rule.visible != default) {
                        out.visible = (bool)rule.visible;
                    }
                    // Appearance.
                    if (rule.head != null) out.head = rule.head;
                    if (!isnan(rule.size)) out.size = rule.size;
                    if (!isnan(rule.relative)) out.relative = rule.relative;
                }
            }
        }
        return out;
    }
}

RayArrows
copy(RayArrows src) {
    if (src == null) return null;
    RayArrows out = new RayArrows;
    out.visible = src.visible;
    out.first = src.first;
    out.middle = src.middle;
    out.last = src.last;
    out.head = src.head;
    out.size = src.size;
    out.relative = src.relative;
    for (SegmentRule r : src.rules) {
        out.rules.push(copy(r));
    }
    return out;
}

// RayArrows constructor
RayArrows rayArrows(bool visible = true,
                    bool3 first = default,
                    bool3 middle = default,
                    bool3 last = default,
                    arrowhead head = config.rendering.arrow_head,
                    real size = config.rendering.arrow_size,
                    real relative = nan, ... SegmentRule[] rules) {
    RayArrows out = new RayArrows;
    out.head = head;
    out.size = size;
    out.relative = relative;
    if (visible) {
        out.visible = true;
        out.first = (first != default) ? (bool)first : true;
        out.middle = (middle != default) ? (bool)middle : true;
        out.last = (last != default) ? (bool)last : true;
        for (SegmentRule r : rules) out.rules.push(r);
    } else {
        // When visible = false is requested:
        // Default all segments to false unless explicitly restored via first/middle/last or segment()
        bool any_explicit = (first == true || middle == true || last == true || rules.length > 0);
        out.visible = any_explicit;
        out.first = (first == true);
        out.middle = (middle == true);
        out.last = (last == true);
        for (SegmentRule r : rules) {
            SegmentRule rule_copy = copy(r);
            if (rule_copy.visible == default) {
                rule_copy.visible = true;
            }
            out.rules.push(rule_copy);
        }
    }
    return out;
}

struct RayPolarization {
    // Global visibility
    bool visible = config.rendering.draw_polarization;
    // Global segment selection
    int[] segments = {-1};
    // Global polarization settings
    int count = -1;
    real spacing = 0.5;
    real before = 0.2;
    real after = 0.2;
    // IMPORTANT:
    // For polarization these mean:
    // first  = first marker in each segment
    // middle = middle markers in each segment
    // last   = last marker in each segment
    bool first = true;
    bool middle = true;
    bool last = true;
    real scale = 0.5;
    bool show_s = true;
    bool show_p = true;
    pen p = nullpen;
    real min_amplitude = 0;
    real min_segment_length = config.rendering.min_polarization_segment_length;
    real min_separation = config.rendering.min_polarization_marker_separation;
    // Polarization arrow appearance.
    // The arrowbar itself is resolved only when needed.
    arrowhead head = config.rendering.arrow_pol_head;
    real size = config.rendering.arrow_pol_size;
    // Segment rules
    SegmentRule[] rules;
    // Resolve polarization settings for one segment
    RayPolarization getPolarization(int i, int total_segments) {
        // Create a fresh copy of current state
        RayPolarization out = new RayPolarization;
        out.visible = this.visible;
        out.count = this.count;
        out.spacing = this.spacing;
        out.before = this.before;
        out.after = this.after;
        out.first = this.first;
        out.middle = this.middle;
        out.last = this.last;
        out.scale = this.scale;
        out.show_s = this.show_s;
        out.show_p = this.show_p;
        out.p = this.p;
        out.min_amplitude = this.min_amplitude;
        out.head = this.head;
        out.size = this.size;
        out.min_segment_length = this.min_segment_length;
        out.min_separation = this.min_separation;
        // Apply rules only to this specific segment
        for (SegmentRule rule : rules) {
            for (int idx : rule.indices) {
                int resolved_idx = (idx < 0) ? idx + total_segments : idx;
                if (resolved_idx == i) {
                    if (rule.visible != default) out.visible = (bool)rule.visible;
                    if (rule.count >= 0) out.count = rule.count;
                    if (!isnan(rule.spacing)) out.spacing = rule.spacing;
                    if (!isnan(rule.before)) out.before = rule.before;
                    if (!isnan(rule.after)) out.after = rule.after;
                    if (rule.first != default) out.first = (bool)rule.first;
                    if (rule.middle != default) out.middle = (bool)rule.middle;
                    if (rule.last != default) out.last = (bool)rule.last;
                    if (!isnan(rule.scale)) out.scale = rule.scale;
                    if (rule.show_s != default) out.show_s = (bool)rule.show_s;
                    if (rule.show_p != default) out.show_p = (bool)rule.show_p;
                    if (rule.p != nullpen) out.p = rule.p;
                    if (rule.head != null) out.head = rule.head;
                    if (!isnan(rule.size)) out.size = rule.size;
                }
            }
        }
        return out;
    }
}

RayPolarization
copy(RayPolarization src) {
    if (src == null) return null;
    RayPolarization out = new RayPolarization;
    out.visible = src.visible;
    out.count = src.count;
    out.spacing = src.spacing;
    out.before = src.before;
    out.after = src.after;
    out.first = src.first;
    out.middle = src.middle;
    out.last = src.last;
    out.scale = src.scale;
    out.show_s = src.show_s;
    out.show_p = src.show_p;
    out.p = src.p;
    out.min_amplitude = src.min_amplitude;
    out.head = src.head;
    out.size = src.size;
    out.min_segment_length = src.min_segment_length;
    out.min_separation = src.min_separation;
    for (SegmentRule r : src.rules) {
        out.rules.push(copy(r));
    }
    return out;
}

// RayPolarization constructor
// When visible = false, any marker position (first/middle/last) explicitly
// passed as true is restored, so you can write:
//   rayPolarization(visible = false, first = true)
// to hide middle and last markers while keeping the first one visible.
RayPolarization rayPolarization(bool visible = config.rendering.draw_polarization,
                                int count = -1,
                                real spacing = 0.5,
                                real before = 0.2,
                                real after = 0.2,
                                bool3 first = default,
                                bool3 middle = default,
                                bool3 last = default,
                                real scale = 0.5,
                                bool show_s = true,
                                bool show_p = true,
                                pen p = nullpen,
                                arrowhead head = config.rendering.arrow_pol_head,
                                real size = config.rendering.arrow_pol_size, ... SegmentRule[] rules) {
    RayPolarization out = new RayPolarization;
    out.count = count;
    out.spacing = spacing;
    out.before = before;
    out.after = after;
    out.scale = scale;
    out.show_s = show_s;
    out.show_p = show_p;
    out.p = p;
    out.head = head;
    out.size = size;
    if (visible) {
        out.visible = true;
        out.first = (first != default) ? (bool)first : true;
        out.middle = (middle != default) ? (bool)middle : true;
        out.last = (last != default) ? (bool)last : true;
        for (SegmentRule r : rules) out.rules.push(r);
    } else {
        // visible = false: default all marker positions to hidden unless
        // explicitly restored via first/middle/last = true or segment() rules.
        bool any_explicit = (first == true || middle == true || last == true || rules.length > 0);
        out.visible = any_explicit;
        out.first = (first == true);
        out.middle = (middle == true);
        out.last = (last == true);
        for (SegmentRule r : rules) {
            SegmentRule rule_copy = copy(r);
            if (rule_copy.visible == default) rule_copy.visible = true;
            out.rules.push(rule_copy);
        }
    }
    return out;
}

private
path polarizationMarkerLine(pair pos, pair d, real Ap, RayPolarization pol, real scale = 1) {
    pair v_perp = perp(d) * (Ap * pol.scale * scale);
    return pos - v_perp -- pos + v_perp;
}

private
pen markerPen(RayPolarization pol, pen p) {
    if (pol.p == nullpen) return p;
    real a = opacity(p);
    if (a == 0) a = 1.0;
    return pol.p + opacity(a);
}

// Core implementation on frame
private
void drawPolarizationMarker(frame f, pair pos, pair d, real As, real Ap, real scale, pen p, RayPolarization pol) {
    pen q = markerPen(pol, p);
    if (pol.show_p && Ap > pol.min_amplitude) {
        path line = polarizationMarkerLine(pos, d, Ap, pol, scale);
        real line_len = arclength(line);
        // Resolve effective arrow size and head safely
        real arrow_sz = (pol.size > 0 && !isnan(pol.size)) ? pol.size : config.rendering.arrow_pol_size;
        arrowhead h = (pol.head != null) ? pol.head : config.rendering.arrow_pol_head;
        // Condition: Line must fit at least 2 arrowheads (one on each side)
        //    Otherwise draw a clean plain line (None) to avoid overlapping arrow blobs
        real min_len_needed = max(config.rendering.min_polarization_marker_with_arrow, 2.0 * arrow_sz);
        arrowbar ab = (line_len > min_len_needed) ? Arrows(h, arrow_sz) : None;
        draw(f, line, q, ab);
    }
    if (pol.show_s && As > pol.min_amplitude) {
        dot(f, pos, q);
    }
}

// Picture wrapper with correct user-to-frame scale factor
private
void drawPolarizationMarker(picture pic = currentpicture,
                            pair pos,
                            pair d,
                            real As,
                            real Ap,
                            real scale = 1,
                            pen p,
                            RayPolarization pol) {
    pic.add(new void(frame f, transform t) {
        // Multiply by length(t * (1, 0)) to convert user units -> canvas pt
        real user_to_pt = length(t * (1, 0)) * scale;
        pair pf = t * pos;
        pair df = t * (pos + d) - pf;
        drawPolarizationMarker(f, pf, df, As, Ap, user_to_pt, p, pol);
    });
}

private
void drawPolarization(picture pic = currentpicture,
                      Ray r,
                      path traj,
                      int segment_idx = 0,
                      int total_segments = 1,
                      pen p = currentpen,
                      RayPolarization pol = rayPolarization(),
                      bool apply_opacity = true) {
    if (r == null || !r.alive) return;
    RayPolarization eff_pol = pol.getPolarization(segment_idx, total_segments);
    if (!eff_pol.visible) return;

    real L = arclength(traj);
    if (L <= config.rendering.min_drawable_length) return;

    // Minimum Segment Length Gate: If the segment is too short to reasonably display a marker, skip entirely.
    real min_seg_l = !isnan(eff_pol.min_segment_length) ? eff_pol.min_segment_length
                                                        : config.rendering.min_polarization_segment_length;
    if (L < min_seg_l) return;

    // Margins and Spacing Limits
    real s_start = clamp(eff_pol.before, 0, L);
    real s_end = clamp(L - eff_pol.after, s_start, L);
    real range = s_end - s_start;
    real min_sep =
        !isnan(eff_pol.min_separation) ? eff_pol.min_separation : config.rendering.min_polarization_marker_separation;
    if (min_sep <= 0) min_sep = 0.1;

    real[] rawPositions;

    // Position Calculation with Overlap Prevention
    if (range < min_sep) {
        // Range is too tight for multiple markers: place exactly ONE centered marker
        if (eff_pol.middle || eff_pol.first || eff_pol.last) {
            rawPositions.push(L / 2.0);
        }
    } else if (eff_pol.count > 0) {
        if (eff_pol.count == 1) {
            rawPositions.push(s_start + range / 2.0);
        } else {
            // Clamp count to prevent stacking
            int max_allowed = max(1, floor(range / min_sep) + 1);
            int eff_count = min(eff_pol.count, max_allowed);
            if (eff_count == 1) {
                rawPositions.push(s_start + range / 2.0);
            } else {
                real step = range / (eff_count - 1);
                for (int i = 0; i < eff_count; ++i) {
                    if (i == 0 && !eff_pol.first) continue;
                    if (i == eff_count - 1 && !eff_pol.last) continue;
                    if (i > 0 && i < eff_count - 1 && !eff_pol.middle) continue;
                    rawPositions.push(s_start + i * step);
                }
            }
        }
    } else {
        // Spacing-based Placement
        real sp = max(min_sep, eff_pol.spacing);
        int num_intervals = floor(range / sp);
        if (num_intervals <= 0) {
            if (eff_pol.middle || eff_pol.first || eff_pol.last) {
                rawPositions.push(s_start + range / 2.0);
            }
        } else {
            real actual_start = s_start + (range - num_intervals * sp) / 2.0;
            for (int i = 0; i <= num_intervals; ++i) {
                if (i == 0 && !eff_pol.first) continue;
                if (i == num_intervals && !eff_pol.last) continue;
                if (i > 0 && i < num_intervals && !eff_pol.middle) continue;
                rawPositions.push(actual_start + i * sp);
            }
        }
    }

    // Strict Deduplication Filter (Guarantees no two markers are closer than min_sep)
    real[] positions;
    real last_s = -infinity;
    for (real s : rawPositions) {
        if (s - last_s >= min_sep - config.math.eps) {
            positions.push(s);
            last_s = s;
        }
    }
    if (positions.length == 0) return;

    // Render Filtered Markers
    for (real s : positions) {
        real t = arctime(traj, s);
        pair pos = point(traj, t);
        pair d = unitOr(dir(traj, t));

        // Evaluate local physical polarization at arc distance s
        Polarization local_pol = resolvePolarizationAt(r, s);

        // Derive S and P field amplitudes directly from Stokes parameters S0 and S1
        real local_as = sqrt(max(0.0, 0.5 * (local_pol.S0 + local_pol.S1)));
        real local_ap = sqrt(max(0.0, 0.5 * (local_pol.S0 - local_pol.S1)));
        real local_i = local_pol.intensity();

        // Skip marker if local intensity has decayed below threshold
        if (local_i < config.physics.min_intensity) continue;
        if (local_as <= eff_pol.min_amplitude && local_ap <= eff_pol.min_amplitude) continue;

        pen pen_local = p;
        if (apply_opacity && config.rendering.use_intensity_opacity) {
            real base_alpha = opacity(p);
            if (base_alpha == 0) base_alpha = 1.0;
            pen_local = p + opacity(base_alpha * max(0.0, min(1.0, local_i)));
        }

        drawPolarizationMarker(pic, pos, d, local_as, local_ap, 1, pen_local, eff_pol);
    }
}

// ============================================================================
// Huygens-Fresnel Wavefront System
// ============================================================================
struct RayWavefront {
    // Global visibility
    bool visible = config.rendering.draw_wavefronts;
    // Global segment selection
    bool first = true;
    bool middle = true;
    bool last = true;
    // Per-ray wavefront crests / Huygens wavelets
    int count = -1;     // Total wavefront crests/wavelets per segment
    real spacing = nan; // Vacuum wavelength lambda0 (auto-scaled by local n)
    real width = 14.0;  // Transverse bar half-width (crest mode)
    real before = 0.05; // Margin from segment start
    real after = 0.05;  // Margin from segment end
    pen p = nullpen;
    // Huygens Wavelet Settings
    bool wavelets = false;       // If true: draw circular Huygens secondary wavelet arcs
    real alpha = 30.0;           // Half-angle of arc (spans d +- alpha, default 30 deg, up to 90 deg)
    bool cut_at_boundary = true; // If true: arcs are clipped at physical surfaces
    bool dashed_behind = false;  // If true: portion behind boundary is dashed instead of omitted
    bool show_source = true;     // Draw center dot at emission origin P0
    // Bundle iso-OPL contours (draw(Trace[], ...) only)
    bool contours = false;      // Draw iso-OPL / iso-phase wavefront contours across the bundle
    bool phase_locked = false;  // Use common physical-OPL contours for bundle wavefronts
    pen p_contours = lightblue; // Pen for iso-OPL contour lines
    // Interference lines (draw(Trace[], ...) only, requires contours or wavelets)
    bool antinodes = false;            // Draw constructive interference channels (Delta phi = 2*pi*m)
    bool nodes = false;                // Draw destructive interference channels (Delta phi = pi*(2m+1))
    pen pen_line = gray(0.4) + dashed; // Pen for interference lines
    pen pen_dots = gray(0.4) + 2pt;    // Pen for intersection dots
    real max_dist = config.rendering.wavefront_max_dist;    // Max distance from source midpoint to trace interference
    int max_orders = config.rendering.wavefront_max_orders; // Maximum interference order |m| to draw
    bool show_dots = true;                                  // Show circle-crossing dots at antinodes
    bool show_labels = true;                                // Label each antinodal channel with m
    // Segment rules
    SegmentRule[] rules;

    // Resolve wavefront configuration for one segment
    RayWavefront getWavefront(int i, int total_segments) {
        bool is_first = (i == 0);
        bool is_last = (i == total_segments - 1);
        RayWavefront out = new RayWavefront;
        out.visible = this.visible && (is_last ? this.last : (is_first ? this.first : this.middle));
        out.count = this.count;
        out.spacing = this.spacing;
        out.width = this.width;
        out.before = this.before;
        out.after = this.after;
        out.first = this.first;
        out.middle = this.middle;
        out.last = this.last;
        out.wavelets = this.wavelets;
        out.alpha = this.alpha;
        out.cut_at_boundary = this.cut_at_boundary;
        out.dashed_behind = this.dashed_behind;
        out.show_source = this.show_source;
        out.p = this.p;
        // Bundle-level fields pass through unchanged
        out.contours = this.contours;
        out.phase_locked = this.phase_locked;
        out.p_contours = this.p_contours;
        out.antinodes = this.antinodes;
        out.nodes = this.nodes;
        out.pen_line = this.pen_line;
        out.pen_dots = this.pen_dots;
        out.max_dist = this.max_dist;
        out.max_orders = this.max_orders;
        out.show_dots = this.show_dots;
        out.show_labels = this.show_labels;

        for (SegmentRule rule : rules) {
            for (int idx : rule.indices) {
                int resolved_idx = (idx < 0) ? idx + total_segments : idx;
                if (resolved_idx == i) {
                    // Segment positional policy
                    if (is_first && rule.first != default) out.visible = (bool)rule.first;
                    if (!is_first && !is_last && rule.middle != default) out.visible = (bool)rule.middle;
                    if (is_last && rule.last != default) out.visible = (bool)rule.last;
                    // Explicit visibility is the final override
                    if (rule.visible != default) out.visible = (bool)rule.visible;
                    // Appearance
                    if (rule.count >= 0) out.count = rule.count;
                    if (!isnan(rule.spacing)) out.spacing = rule.spacing;
                    if (!isnan(rule.scale)) out.width = rule.scale * this.width;
                    if (!isnan(rule.before)) out.before = rule.before;
                    if (!isnan(rule.after)) out.after = rule.after;
                    if (rule.p != nullpen) out.p = rule.p;
                }
            }
        }
        return out;
    }
}

RayWavefront
copy(RayWavefront src) {
    if (src == null) return null;
    RayWavefront out = new RayWavefront;
    out.visible = src.visible;
    out.first = src.first;
    out.middle = src.middle;
    out.last = src.last;
    out.count = src.count;
    out.spacing = src.spacing;
    out.width = src.width;
    out.before = src.before;
    out.after = src.after;
    out.p = src.p;
    out.wavelets = src.wavelets;
    out.alpha = src.alpha;
    out.cut_at_boundary = src.cut_at_boundary;
    out.dashed_behind = src.dashed_behind;
    out.show_source = src.show_source;
    out.contours = src.contours;
    out.phase_locked = src.phase_locked;
    out.p_contours = src.p_contours;
    out.antinodes = src.antinodes;
    out.nodes = src.nodes;
    out.pen_line = src.pen_line;
    out.pen_dots = src.pen_dots;
    out.max_dist = src.max_dist;
    out.max_orders = src.max_orders;
    out.show_dots = src.show_dots;
    out.show_labels = src.show_labels;
    for (SegmentRule r : src.rules) out.rules.push(copy(r));
    return out;
}

// RayWavefront constructor
// Controls per-ray wavefront crests/wavelets AND bundle-level iso-OPL contours
// and interference lines when used with draw(Trace[], ...).
//
// When visible = false, any marker position (first/middle/last) explicitly
// passed as true is restored, so you can write:
//   rayWavefront(visible = false, last = true)
// to hide crests on all segments except the last one.
RayWavefront rayWavefront(bool visible = config.rendering.draw_wavefronts,
                          real spacing = nan,
                          int count = -1,
                          real width = 14.0,
                          real before = 0.05,
                          real after = 0.05,
                          bool wavelets = false,
                          real alpha = 30.0,
                          bool cut_at_boundary = true,
                          bool dashed_behind = false,
                          bool show_source = true,
                          pen p = nullpen,
                          // Bundle iso-OPL contours
                          bool contours = false,
                          bool phase_locked = false,
                          pen p_contours = lightblue,
                          // Interference lines
                          bool antinodes = false,
                          bool nodes = false,
                          pen pen_line = gray(0.4) + dashed,
                          pen pen_dots = gray(0.4) + 2.5pt,
                          real max_dist = config.rendering.wavefront_max_dist,
                          int max_orders = config.rendering.wavefront_max_orders,
                          bool show_dots = true,
                          bool show_labels = true,
                          bool3 first = default,
                          bool3 middle = default,
                          bool3 last = default, ... SegmentRule[] rules) {
    RayWavefront out = new RayWavefront;
    out.spacing = spacing;
    out.width = width;
    out.count = count;
    out.before = before;
    out.after = after;
    out.wavelets = wavelets;
    out.alpha = alpha;
    out.cut_at_boundary = cut_at_boundary;
    out.dashed_behind = dashed_behind;
    out.show_source = show_source;
    out.p = p;
    out.contours = contours;
    out.phase_locked = phase_locked;
    out.p_contours = p_contours;
    out.antinodes = antinodes;
    out.nodes = nodes;
    out.pen_line = pen_line;
    out.pen_dots = pen_dots;
    out.max_dist = max_dist;
    out.max_orders = max_orders;
    out.show_dots = show_dots;
    out.show_labels = show_labels;
    if (visible) {
        out.visible = true;
        out.first = (first != default) ? (bool)first : true;
        out.middle = (middle != default) ? (bool)middle : true;
        out.last = (last != default) ? (bool)last : true;
        for (SegmentRule r : rules) out.rules.push(r);
    } else {
        // visible = false: default all segment positions to hidden unless
        // explicitly restored via first/middle/last = true or segment() rules.
        bool any_explicit = (first == true || middle == true || last == true || rules.length > 0);
        out.visible = any_explicit;
        out.first = (first == true);
        out.middle = (middle == true);
        out.last = (last == true);
        for (SegmentRule r : rules) {
            SegmentRule rule_copy = copy(r);
            if (rule_copy.visible == default) rule_copy.visible = true;
            out.rules.push(rule_copy);
        }
    }
    return out;
}

// Standalone Wavefront Marker / Huygens Wavelet Renderer
void drawWavefrontMarker(picture pic = currentpicture,
                         pair center,        // Emission origin P0 (for wavelets) or bar center (for crests)
                         pair d,             // Wave propagation direction
                         real radius,        // Arc radius (for wavelets) or arc distance s
                         real local_n = 1.0, // Refractive index for wavelength scaling
                         Hit hit = null,     // Optional boundary for automatic clipping
                         pen p = currentpen,
                         RayWavefront wf = rayWavefront()) {
    pen q = (wf.p != nullpen) ? wf.p : p;
    pair d_hat = unitOr(d);

    if (wf.wavelets) {
        // =====================================================================
        // Mode B: Circular Huygens Secondary Wavelet Arc
        // =====================================================================
        real theta0 = degrees(atan2(d_hat.y, d_hat.x));
        real a_span = clamp(wf.alpha, 5.0, 90.0);
        path arc_full = Arc(center, radius, theta0 - a_span, theta0 + a_span, 160);

        // Boundary clipping: truncate the arc if it hits a physical surface
        if (wf.cut_at_boundary && hit != null && hit.valid) {
            pair hit_p = hit.p;
            pair N = unitOr(hit.normal);
            if (dot(d_hat, N) > 0) N = -N; // Opposes incident wave

            guide arc_visible, arc_behind;
            bool visible_run = false;
            bool hidden_run = false;
            int n_pts = 160;
            for (int j = 0; j <= n_pts; ++j) {
                pair pt = point(arc_full, reltime(arc_full, j / (real)n_pts));
                bool is_visible = dot(pt - hit_p, N) >= -config.math.eps;
                if (is_visible) {
                    if (wf.dashed_behind && hidden_run && length(arc_behind) > 1)
                        draw(pic, arc_behind, q + dashed);
                    arc_behind = nullpath;
                    hidden_run = false;
                    arc_visible = arc_visible -- pt;
                    visible_run = true;
                } else {
                    if (visible_run && length(arc_visible) > 1) draw(pic, arc_visible, q);
                    arc_visible = nullpath;
                    visible_run = false;
                    arc_behind = arc_behind -- pt;
                    hidden_run = true;
                }
            }
            if (length(arc_visible) > 1) draw(pic, arc_visible, q);
            if (wf.dashed_behind && length(arc_behind) > 1) draw(pic, arc_behind, q + dashed);
        } else {
            draw(pic, arc_full, q);
        }
    } else {
        // =====================================================================
        // Mode A: Planar / Orthogonal Wavefront Crest Bar
        // =====================================================================
        pair w = perp(d_hat); // Transverse vector perpendicular to ray
        real half_w = 0.5 * wf.width;
        draw(pic, (center - half_w * w) -- (center + half_w * w), q + 1.0bp);
    }
}

// Picture wrapper for canvas coordinate mapping
private
void drawWavefront(picture pic = currentpicture,
                   Ray r,
                   path traj,
                   int segment_idx = 0,
                   int total_segments = 1,
                   pen p = currentpen,
                   RayWavefront wf = rayWavefront()) {
    if (r == null || !r.alive || !wf.visible) return;
    RayWavefront eff_wf = wf.getWavefront(segment_idx, total_segments);
    if (!eff_wf.visible) return;

    real L = arclength(traj);
    if (L <= config.rendering.min_drawable_length) return;

    real local_n = (r.medium != null) ? r.medium.n(r.p, r.d, r.lambda, r.pol) : 1.0;
    if (local_n <= config.math.eps) local_n = 1.0;

    pen pen_draw = (eff_wf.p != nullpen) ? eff_wf.p : p;
    pair P0 = point(traj, 0);
    pair d0 = unitOr(dir(traj, 0));

    // Calculate expanding radii / distance steps:
    real[] radii;
    if (eff_wf.count > 0) {
        real step = L / (eff_wf.count + 1);
        for (int i = 1; i <= eff_wf.count; ++i) radii.push(i * step);
    } else {
        real physical_lambda = (r.lambda > 0) ? r.lambda : config.physics.lambda;
        real sp = (!isnan(eff_wf.spacing) && eff_wf.spacing > 0)
                      ? (eff_wf.spacing / local_n)
                      : resolveEffectiveLambda(nan, physical_lambda) / local_n;
        if (sp <= config.math.eps) sp = 12.0;
        int requested = max(0, ceil((L - config.math.eps) / sp));
        int marker_count = min(requested, max(1, config.rendering.wavefront_max_markers));
        if (requested <= marker_count) {
            for (real cur_r = sp; cur_r < L - config.math.eps; cur_r += sp) {
                radii.push(cur_r);
            }
        } else {
            // Extremely small rendered wavelengths can otherwise allocate
            // millions of paths. Preserve coverage while bounding rendering
            // work; this affects drawing only, never physical phase/OPL.
            real render_step = L / (marker_count + 1);
            for (int i = 1; i <= marker_count; ++i) {
                radii.push(i * render_step);
            }
        }
    }

    // Draw source dot once if wavelets are enabled
    if (eff_wf.wavelets && eff_wf.show_source) dot(pic, P0, pen_draw + 2.0pt);

    // Render each wavefront marker:
    for (real rad : radii) {
        if (eff_wf.wavelets) {
            // Wavelet mode: center is P0, radius expands outward
            drawWavefrontMarker(pic, P0, d0, rad, local_n, r.hit, pen_draw, eff_wf);
        } else {
            // Transverse bar mode: center is point along the trajectory
            real u = arctime(traj, rad);
            pair pos = point(traj, u);
            pair d_local = unitOr(dir(traj, u));
            drawWavefrontMarker(pic, pos, d_local, rad, local_n, r.hit, pen_draw, eff_wf);
        }
    }
}

// Helper packet for a point on a trace at a specific OPL
struct TracePoint {
    bool valid = false;
    pair p = (nan, nan);
    pair d = (1, 0);
    real intensity = 1.0;
    int seg_idx = -1;
    real seg_o_p_l_start = 0.0;
    real seg_o_p_l_end = 0.0;
}

// Robust OPL sampler along any ray trace
private TracePoint pointAtOPL(Segment[] segs, real target_o_p_l) {
    TracePoint tp;
    int n_segs = segs.length;
    if (n_segs == 0) return tp;

    real cur_opl = 0.0;
    for (int i = 0; i < n_segs; ++i) {
        Segment s = segs[i];
        real seg_len = arclength(s.traj);
        real n_med = rayIndex(s.before);
        real seg_opl = geometryToUm(seg_len) * n_med;
        real next_opl = cur_opl + seg_opl;

        // Strict range check: if target_o_p_l exceeds next_opl on the last segment, the ray has ENDED!
        if (target_o_p_l >= cur_opl - config.math.eps && target_o_p_l <= next_opl + config.math.eps) {
            real dist_along = (n_med > config.math.eps) ? (target_o_p_l - cur_opl) / n_med : 0.0;
            dist_along = clamp(dist_along, 0.0, seg_len);
            real u = (seg_len > 0) ? arctime(s.traj, dist_along) : 0.0;

            tp.valid = true;
            tp.p = point(s.traj, u);
            tp.d = unitOr(dir(s.traj, u), s.before.d);
            tp.intensity = s.before.intensity();
            tp.seg_idx = i;
            tp.seg_o_p_l_start = cur_opl;
            tp.seg_o_p_l_end = next_opl;
            return tp;
        }
        cur_opl = next_opl;
    }
    return tp; // Ray has ended or not reached yet -> invalid!
}

// Phase-locked Huygens wavelets from common physical OPL targets.
private
void drawPhaseLockedWavelets(picture pic = currentpicture, Trace[] traces, RayWavefront wf) {
    if (traces.length < 2) return;
    real source_lambda = (traces[0].first() != null && traces[0].first().lambda > 0)
                         ? traces[0].first().lambda : config.physics.lambda;
    real spacing_um = geometryToUm((!isnan(wf.spacing) && wf.spacing > 0)
                                   ? wf.spacing : resolveEffectiveLambda(nan, source_lambda));
    if (spacing_um <= config.math.eps) return;
    // For explicitly requested interference families, use the identical
    // paired-circle construction as drawBundleInterference. This guarantees
    // that visible wavelet intersections lie on the same node/antinode loci.
    if (wf.antinodes || wf.nodes) {
        Ray a = traces[0].first(), b = traces[1].first();
        real lambda_geom = spacing_um / metersToUm(config.physics.meters_per_unit);
        real separation = length(a.p - b.p);
        if (separation <= config.math.eps) return;
        int max_order = max(0, wf.max_orders);
        for (int m = -max_order; m <= max_order; ++m) {
            for (real offset = 0.0; offset <= 0.5; offset += 0.5) {
                if ((offset == 0.0 && !wf.antinodes) || (offset == 0.5 && !wf.nodes)) continue;
                int n_max = ceil((wf.max_dist + separation) / lambda_geom);
                for (int n = 1; n <= n_max; ++n) {
                    real r1 = n * lambda_geom;
                    real r2 = (n - m - offset) * lambda_geom;
                    if (r2 <= 0 || abs(r1 - r2) > separation || r1 + r2 < separation) continue;
                    Hit hit_a = (traces[0].last() != null) ? traces[0].last().hit : null;
                    Hit hit_b = (traces[1].last() != null) ? traces[1].last().hit : null;
                    drawWavefrontMarker(pic, a.p, a.d, r1, 1.0, hit_a,
                                        (wf.p != nullpen) ? wf.p : lightblue, wf);
                    drawWavefrontMarker(pic, b.p, b.d, r2, 1.0, hit_b,
                                        (wf.p != nullpen) ? wf.p : lightblue, wf);
                }
            }
        }
        return;
    }
    real max_opl = 0.0;
    for (Trace tr : traces) {
        if (tr == null || tr.first() == null) continue;
        real total = 0.0;
        for (Segment s : tr.segments()) total += geometryToUm(arclength(s.traj)) * rayIndex(s.before);
        max_opl = max(max_opl, total);
    }
    int count = min(max(1, ceil(max_opl / spacing_um)), max(1, config.rendering.wavefront_max_markers));
    for (int i = 1; i <= count; ++i) {
        real target = i * max_opl / (count + 1);
        for (int k = 0; k < traces.length; ++k) {
            Ray src = traces[k].first();
            if (src == null) continue;
            // Use the same effective wavelength as the interference-line
            // renderer. This keeps wavelet crossings aligned with nodes and
            // antinodes when spacing is a visual geometry wavelength.
            real k0 = 2.0 * pi / spacing_um;
            real phase_offset = (k0 > config.math.eps) ? src.phase / k0 : 0.0;
            real opl_from_source = target - phase_offset;
            if (wf.nodes && k == 1) opl_from_source -= 0.5 * spacing_um;
            if (opl_from_source <= config.math.eps) continue;
            real n_src = (src.medium != null) ? src.medium.n(src.p, src.d, src.lambda, src.pol) : 1.0;
            real radius = opl_from_source / max(config.math.eps, n_src) /
                          metersToUm(config.physics.meters_per_unit);
            drawWavefrontMarker(pic, src.p, src.d, radius, n_src, null,
                                (wf.p != nullpen) ? wf.p : lightblue, wf);
        }
    }
}

// Private: Iso-OPL / Iso-Phase Wavefront Contours across a ray bundle
private
void drawBundleContours(picture pic = currentpicture, Trace[] traces, RayWavefront wf) {
    int N = traces.length;
    if (N < 2) return;
    pen p = (wf.p_contours != nullpen) ? wf.p_contours : heavycyan;
    int count = (wf.count > 0) ? wf.count : 12;

    // 1. Resolve Incoming Direction and Reference Origin
    Trace t0 = traces[0];
    if (t0.first() == null) return;
    real source_lambda = (t0.first().lambda > 0) ? t0.first().lambda : config.physics.lambda;
    real spacing_geom = (!isnan(wf.spacing) && wf.spacing > 0) ? wf.spacing : resolveEffectiveLambda(nan, source_lambda);
    real spacing = geometryToUm(spacing_geom);
    pair p0_ref = t0.first().p;
    pair d_in = unitOr(t0.first().d);
    real n_amb = (t0.first().medium != null) ? t0.first().medium.n() : 1.0;

    // Segment lists are rebuilt on every Trace.segments() call; cache them
    // once per trace instead of once per (contour, ray) pair.
    Segment[][] segs_of = new Segment[N][];
    for (int k = 0; k < N; ++k) segs_of[k] = traces[k].segments();

    // 2. Compute Global OPL Range Across All Rays
    real[] launch_offsets = new real[N];
    real min_global_o_p_l = infinity;
    real max_global_o_p_l = -infinity;

    for (int k = 0; k < N; ++k) {
        Trace t = traces[k];
        if (t.first() == null || t.last() == null) continue;

        real k0 = vacuumWaveNumber(t.first().lambda);
        real phase_opl_offset = (k0 > config.math.eps) ? t.first().phase / k0 : 0.0;
        real offset = geometryToUm(dot(t.first().p - p0_ref, d_in) * n_amb) + phase_opl_offset;
        launch_offsets[k] = offset;

        real trace_total_o_p_l = 0.0;
        for (Segment s : segs_of[k]) {
            real n_med = rayIndex(s.before);
            trace_total_o_p_l += geometryToUm(arclength(s.traj)) * n_med;
        }

        min_global_o_p_l = min(min_global_o_p_l, offset);
        max_global_o_p_l = max(max_global_o_p_l, offset + trace_total_o_p_l);
    }

    if (min_global_o_p_l >= max_global_o_p_l - config.math.eps) return;

    // 3. Determine Target OPL values
    real[] targetOPLs;
    if (!isnan(spacing) && spacing > 0) {
        int requested = max(0, ceil((max_global_o_p_l - min_global_o_p_l) / spacing));
        int target_count = min(requested, max(1, config.rendering.wavefront_max_markers));
        if (requested <= target_count) {
            for (real L = min_global_o_p_l + spacing; L < max_global_o_p_l; L += spacing) {
                targetOPLs.push(L);
            }
        } else {
            real contour_step = (max_global_o_p_l - min_global_o_p_l) / (target_count + 1);
            for (int i = 1; i <= target_count; ++i) {
                targetOPLs.push(min_global_o_p_l + i * contour_step);
            }
        }
    } else {
        real d_l = (max_global_o_p_l - min_global_o_p_l) / (count + 1);
        for (int i = 1; i <= count; ++i) {
            targetOPLs.push(min_global_o_p_l + i * d_l);
        }
    }

    // 4. Trace Wavefront Contours at Each Target OPL
    for (real L_target : targetOPLs) {
        guide current_wf;
        int current_run = 0;
        TracePoint prev_t_p;
        int prev_k = -1;

        for (int k = 0; k < N; ++k) {
            real target_along_ray = L_target - launch_offsets[k];
            TracePoint tp = pointAtOPL(segs_of[k], target_along_ray);

            if (!tp.valid) {
                // Ray k has ended (or hasn't started):
                // Finish and draw the current contiguous run if it spans >= 2 rays:
                if (current_run >= 2) {
                    draw(pic, current_wf, p);
                }
                current_wf = nullpath;
                current_run = 0;
                prev_t_p.valid = false;
                prev_k = -1;
                continue;
            }

            // Ray k is active!
            if (prev_t_p.valid && k == prev_k + 1) {
                // Contiguous adjacent rays both active:
                // Check if crossing between air and glass:
                if (tp.seg_idx != prev_t_p.seg_idx) {
                    Segment[] s_curr = segs_of[k];
                    Segment[] s_prev = segs_of[prev_k];
                    int trans_seg = min(tp.seg_idx, prev_t_p.seg_idx);

                    if (trans_seg < s_curr.length && trans_seg < s_prev.length) {
                        pair S_hit_curr = s_curr[trans_seg].b;
                        pair S_hit_prev = s_prev[trans_seg].b;

                        real opl_hit_curr =
                            launch_offsets[k] + (tp.seg_idx > trans_seg ? tp.seg_o_p_l_start : tp.seg_o_p_l_end);
                        real opl_hit_prev =
                            launch_offsets[prev_k] +
                            (prev_t_p.seg_idx > trans_seg ? prev_t_p.seg_o_p_l_start : prev_t_p.seg_o_p_l_end);

                        real denom = opl_hit_curr - opl_hit_prev;
                        real tau =
                            (abs(denom) > config.math.eps) ? clamp((L_target - opl_hit_prev) / denom, 0.0, 1.0) : 0.5;
                        pair S_contact = S_hit_prev + tau * (S_hit_curr - S_hit_prev);

                        current_wf = current_wf -- S_contact; // Route through boundary contact point S
                    }
                }

                current_wf = current_wf -- tp.p;
                ++current_run;
            } else {
                // Start a new contiguous wavefront line at tp.p
                current_wf = tp.p;
                current_run = 1;
            }

            prev_t_p = tp;
            prev_k = k;
        }

        // Draw pending run at the end of the rays:
        if (current_run >= 2) {
            draw(pic, current_wf, p);
        }
    }
}

// Private: Interference antinodal/nodal channels between two point sources
// Draws loci of Delta-phi = 2*pi*m (antinodes) or 2*pi*(m+1/2) (nodes)
private
void drawBundleInterference(picture pic = currentpicture, Trace t1, Trace t2, RayWavefront wf) {
    if (t1 == null || t2 == null) return;
    if (t1.first() == null || t2.first() == null) return;

    pair S1 = t1.first().p;
    pair S2 = t2.first().p;
    pair mid_s = 0.5 * (S1 + S2);

    pair v_slit = S1 - S2;
    real d = length(v_slit);
    if (d <= config.math.eps) return;

    // Orthonormal coordinate system centred at mid_s
    pair u_slit = v_slit / d;
    pair u_axis = perp(u_slit);
    if (dot(u_axis, t1.first().d) < 0) u_axis = -u_axis;

    // Resolve wavelength from explicit spacing, global vis_lambda, or the
    // physical ray wavelength. All values here are geometric units.
    real source_lambda = (t1.first().lambda > 0) ? t1.first().lambda : config.physics.lambda;
    real lam = (!isnan(wf.spacing) && wf.spacing > 0) ? wf.spacing : resolveEffectiveLambda(nan, source_lambda);
    if (lam <= config.math.eps) return;

    int max_order_possible = floor(d / lam);
    int orders_to_draw = min(wf.max_orders, max_order_possible);

    // Collect half-integer offsets for nodes (0 = antinodes only, 0.5 = nodes only)
    real[] offsets;
    if (wf.antinodes) offsets.push(0.0);
    if (wf.nodes) offsets.push(0.5);

    for (real offset : offsets) {
        // Loop over each interference order m = 0, +-1, +-2, ...
        for (int m = -orders_to_draw; m <= orders_to_draw; ++m) {
            guide hyperbola_path;
            pair last_pt = (nan, nan);

            int n_max = ceil((wf.max_dist + d) / lam);
            for (int n1 = 1; n1 <= n_max; ++n1) {
                // Antinodes: r1 - r2 = m*lam
                // Nodes:     r1 - r2 = (m + 0.5)*lam
                real n2f = n1 - m - offset; // fractional n2
                if (n2f <= 0) continue;
                real r1 = n1 * lam;
                real r2 = n2f * lam;

                // Two circles intersect iff |r1 - r2| <= d <= r1 + r2
                if (abs(r1 - r2) <= d && (r1 + r2) >= d) {
                    real y_local = (r2 * r2 - r1 * r1) / (2.0 * d);
                    real x2 = r1 * r1 - (y_local - d / 2.0) ^ 2;
                    if (x2 >= 0) {
                        real x_local = sqrt(x2);
                        if (x_local <= wf.max_dist) {
                            pair P_cross = mid_s + x_local * u_axis + y_local * u_slit;
                            hyperbola_path = hyperbola_path -- P_cross;
                            if (wf.show_dots) dot(pic, P_cross, wf.pen_dots);
                            last_pt = P_cross;
                        }
                    }
                }
            }

            if (length(hyperbola_path) > 0) {
                draw(pic, hyperbola_path, wf.pen_line);
                if (wf.show_labels && !isnan(last_pt.x)) {
                    string lbl;
                    if (offset == 0) {
                        lbl = (m == 0) ? "$m=0$" : format("$m=%+d$", m);
                    } else {
                        lbl = (m >= 0) ? format("$m=%+d.5$", m) : format("$m=%d.5$", m + 1);
                    }
                    label(pic, Label("\tiny " + lbl), last_pt, E, wf.pen_line);
                }
            }
        }
    }
}
private
void drawRaySegment(picture pic = currentpicture,
                    Label L = Label(""),
                    Ray r,
                    path traj,
                    int segment_idx = 0,
                    int total_segments = 1,
                    pen p = currentpen,
                    ResolvedArrow arrow = new ResolvedArrow,
                    RayPolarization polarization = rayPolarization(),
                    RayWavefront wavefronts = rayWavefront(),
                    bool apply_opacity = true) {
    real L_traj = arclength(traj);
    if (L_traj < config.rendering.min_drawable_length) return;
    pen pen_line = p;
    if (apply_opacity) {
        real start_i = max(0.0, min(1.0, r.intensity()));
        real base_alpha = opacity(p);
        if (base_alpha == 0) base_alpha = 1.0;
        pen_line = p + opacity(base_alpha * start_i);
    }
    arrowbar ab =
        (arrow.visible && L_traj >= config.rendering.min_polarization_segment_length) ? (arrowbar)arrow : None;

    drawWavefront(pic, r, traj, segment_idx, total_segments, p, wavefronts);
    draw(pic, L, traj, pen_line, ab);
    drawPolarization(pic, r, traj, segment_idx, total_segments, p, polarization, apply_opacity);
}

private
void drawFadingSegment(picture pic = currentpicture,
                       Label L = Label(""),
                       Ray r,
                       path traj,
                       int segment_idx = 0,
                       int total_segments = 1,
                       pen p = currentpen,
                       ResolvedArrow arrow = new ResolvedArrow,
                       RayPolarization polarization = rayPolarization(),
                       RayWavefront wavefronts = rayWavefront(),
                       int steps = 16,
                       bool apply_opacity = true) {
    // Draw wavefronts (if requested)
    drawWavefront(pic, r, traj, segment_idx, total_segments, p, wavefronts);

    real L_traj = arclength(traj);
    real alpha = (r.medium != null) ? r.medium.absorption(r.p, r.d, r.lambda, r.pol) : 0.0;
    real start_i = r.intensity();
    if (L_traj <= config.rendering.min_drawable_length || steps <= 1 || alpha <= 0 || start_i <= config.math.eps) {
        drawRaySegment(pic, L, r, traj, segment_idx, total_segments, p, arrow, polarization, wavefronts, apply_opacity);
        return;
    }
    real base_alpha = opacity(p);
    if (base_alpha == 0) base_alpha = 1.0;
    real step_l = L_traj / steps;
    // Draw exponentially fading sub-segments without arrows
    for (int i = 0; i < steps; ++i) {
        real t1 = arctime(traj, i * step_l);
        real t2 = arctime(traj, (i + 1) * step_l);
        path sub = subpath(traj, t1, t2);
        real local_dist = (i + 0.5) * step_l;
        real local_i = start_i * exp(-alpha * geometryToUm(local_dist));
        pen pen_sub = p;
        if (apply_opacity) {
            pen_sub = p + opacity(base_alpha * max(0.0, min(1.0, local_i))) + squarecap;
        } else {
            pen_sub = p + squarecap;
        }
        draw(pic, (i == 0 ? L : Label("")), sub, pen_sub, None);
    }
    // Overlay a single full-sized arrowhead matching the local faded intensity
    if (arrow.visible) {
        arrowhead h = (arrow.head != null) ? arrow.head : config.rendering.arrow_head;
        real sz = isnan(arrow.size) ? config.rendering.arrow_size : arrow.size;
        real pos = isnan(arrow.relative) ? (segment_idx == total_segments - 1 ? 1.0 : 0.5) : arrow.relative;
        pen pen_arrow = p;
        if (apply_opacity) {
            real local_i = start_i * exp(-alpha * geometryToUm(pos * L_traj));
            pen_arrow = p + opacity(base_alpha * max(0.0, min(1.0, local_i)));
        }
        pic.add(arrow(h, traj, invisible, sz, FillDraw(pen_arrow, pen_arrow), position = Relative(pos)));
    }
    // Draw polarization markers
    drawPolarization(pic, r, traj, segment_idx, total_segments, p, polarization, apply_opacity);
}

// Fast 2D Ray-Viewport Slab Clipper (Device Space)
private
real[] clipRayToViewport(pair p, pair d, pair m, pair M) {
    real Tin = -infinity, Tout = infinity;
    if (abs(d.x) < config.rendering.clip_velocity_eps) {
        if (p.x < m.x || p.x > M.x) return new real[]{};
    } else {
        real inv_dx = 1.0 / d.x;
        real t1 = (m.x - p.x) * inv_dx;
        real t2 = (M.x - p.x) * inv_dx;
        Tin = max(Tin, min(t1, t2));
        Tout = min(Tout, max(t1, t2));
    }
    if (abs(d.y) < config.rendering.clip_velocity_eps) {
        if (p.y < m.y || p.y > M.y) return new real[]{};
    } else {
        real inv_dy = 1.0 / d.y;
        real t1 = (m.y - p.y) * inv_dy;
        real t2 = (M.y - p.y) * inv_dy;
        Tin = max(Tin, min(t1, t2));
        Tout = min(Tout, max(t1, t2));
    }
    if (Tin <= Tout) return new real[]{Tin, Tout};
    return new real[]{};
}

private
void drawClippedRay(picture pic,
                    pair pt_world,
                    pair dir_world,
                    real margin,
                    pen p,
                    ResolvedArrow arrow,
                    bool clip_backward,
                    bool clip_forward) {
    pic.add(new void(frame f, transform t, transform T_trans, pair m, pair M) {
        pair m_eff = m + (margin, margin);
        pair M_eff = M - (margin, margin);
        pair pf = t * pt_world;
        pair df = t * (pt_world + dir_world) - pf;
        real[] clip = clipRayToViewport(pf, df, m_eff, M_eff);
        if (clip.length == 2) {
            real t0 = clip_backward ? clip[0] : 0.0;
            real t1 = clip_forward ? clip[1] : 0.0;
            if (t0 < t1) {
                draw(f, (pf + t0 * df) -- (pf + t1 * df), p, (arrowbar)arrow);
            }
        }
    });
}

// Helper to modulate pen opacity by the ray's physical intensity (if enabled)
private
pen resolveRayPen(pen pen_base, Ray r, bool apply_op) {
    if (!apply_op || r == null) return pen_base;
    real base_alpha = opacity(pen_base);
    if (base_alpha == 0) base_alpha = 1.0;
    real I = clamp(r.intensity(), 0.0, 1.0);
    return pen_base + opacity(base_alpha * I);
}

// Complete Renderer Struct
// Default Surface Drawing Delegate
private
void drawSurfaces(picture pic = currentpicture, Element e) {
    for (Surface s : e.surfaces) {
        if (length(s.p) > 0) {
            transform total_T = e.T * s.T;
            if (s.fill != nullpen && cyclic(s.p)) {
                fill(pic, total_T * s.p, s.fill);
            }
            draw(pic, total_T * s.p, s.stroke);
        }
    }
}

struct Renderer {
    pen p = currentpen;
    RayArrows arrows = rayArrows();
    RayPolarization polarization = rayPolarization();
    RayWavefront wavefronts = rayWavefront();
    bool intensity_opacity = true;
    bool draw_absorption_fading = config.rendering.draw_absorption_fading;
    int fade_subpath_steps = config.rendering.fade_subpath_steps;
    real before = nan;
    real after = nan;
    real margin = 0;
    bool draw_reflections = config.rendering.draw_reflections;
    void drawTrace(picture pic = currentpicture, Trace T) {
        Segment[] segs = T.segments();
        Trace[] brs = T.branches();
        if (segs.length == 0 && brs.length == 0) return;
        for (int i = 0; i < segs.length; ++i) {
            Segment s = segs[i];
            path tr = s.traj;
            real L_tr = arclength(tr);
            bool is_first = (i == 0);
            bool is_last = (i == segs.length - 1);
            ResolvedArrow current_arr = arrows.getArrow(i, segs.length);
            bool apply_op = config.rendering.use_intensity_opacity && intensity_opacity;
            pen pen_eff = resolveRayPen(p, s.before, apply_op);
            // Special Case: Single unhit segment passing entirely through space
            if (is_first && is_last && (before == infinity || before == -infinity) && after == infinity) {
                pair p0 = (length(tr) > 0) ? point(tr, 0) : s.before.p;
                pair d0 = (length(tr) > 0) ? dir(tr, 0) : s.before.d;
                drawClippedRay(pic, p0, d0, this.margin, pen_eff, current_arr, true, true);
                continue;
            }
            // A. First Segment (Entrance Leg)
            if (is_first) {
                if (before == 0) {
                    tr = nullpath;
                } else if (before == infinity || before == -infinity) {
                    pair b = point(tr, length(tr));
                    pair d_in = (length(tr) > 0) ? dir(tr, length(tr)) : s.before.d;
                    drawClippedRay(pic, b, d_in, this.margin, pen_eff, current_arr, true, false);
                    if (!is_last) continue;
                    tr = nullpath;
                } else if (!isnan(before) && before > 0) {
                    if (before < L_tr) {
                        // Shorten: start closer to the surface
                        tr = subpath(tr, arctime(tr, L_tr - before), length(tr));
                    } else if (before > L_tr) {
                        // Extend backwards past the source position P0:
                        pair p0 = point(tr, 0);
                        pair d0 = (length(tr) > 0) ? dir(tr, 0) : s.before.d;
                        pair p_ext = p0 - (before - L_tr) * unitOr(d0);

                        if (length(tr) == 1 && straight(tr, 0)) {
                            // Keep as a single clean straight segment to the surface:
                            tr = p_ext -- point(tr, 1);
                        } else {
                            // Prepend backward extension for multi-segment paths:
                            tr = p_ext -- tr;
                        }
                    }
                }
            }
            // B. Last Segment (Exit Leg)
            if (is_last && (s.after == null || s.after.alive)) {
                if (after == 0) {
                    tr = nullpath;
                } else if (after == infinity) {
                    pair a = point(tr, 0);
                    pair d_out = (length(tr) > 0) ? dir(tr, 0) : s.before.d;
                    drawClippedRay(pic, a, d_out, this.margin, pen_eff, current_arr, false, true);
                    continue;
                } else if (!isnan(after) && after > 0 && after < L_tr) {
                    tr = subpath(tr, 0, arctime(tr, after));
                }
            }
            if (arclength(tr) < config.rendering.min_drawable_length) continue;
            if (draw_absorption_fading && s.before.medium != null &&
                s.before.medium.absorption(s.before.p, s.before.d, s.before.lambda, s.before.pol) > 0) {
                drawFadingSegment(pic,
                                  Label(""),
                                  s.before,
                                  tr,
                                  i,
                                  segs.length,
                                  p,
                                  current_arr,
                                  polarization,
                                  wavefronts,
                                  fade_subpath_steps,
                                  apply_op);
            } else {
                drawRaySegment(pic,
                               Label(""),
                               s.before,
                               tr,
                               i,
                               segs.length,
                               pen_eff,
                               current_arr,
                               polarization,
                               wavefronts,
                               false);
            }
        }
        // Branch rendering for reflections
        for (Trace br : brs) {
            if (br.is_reflection && !draw_reflections) continue;
            Renderer branch_renderer = new Renderer;
            branch_renderer.p = p;
            branch_renderer.arrows = copy(arrows);
            if (segs.length > 0) branch_renderer.arrows.first = false;
            branch_renderer.before = (segs.length == 0) ? before : nan;
            branch_renderer.after = after;
            branch_renderer.margin = margin;
            branch_renderer.draw_reflections = draw_reflections;
            branch_renderer.polarization = polarization;
            branch_renderer.drawTrace(pic, br);
        }
    }
    void drawElement(picture pic = currentpicture, Element e) {
        if (e.fill != nullpen && cyclic(e.boundary)) {
            fill(pic, e.T * e.boundary, e.fill);
        }
        if (e.custom_draw != null) {
            e.custom_draw(pic, e);
        } else {
            drawSurfaces(pic, e);
        }
        for (Element child : e.children) {
            drawElement(pic, e.T * child);
        }
    }
}

// Global Public Trace Drawing Wrapper
void draw(picture pic = currentpicture,
          Trace T,
          real before = nan,
          real after = nan,
          real margin = 0,
          pen p = currentpen,
          RayArrows arrows = rayArrows(),
          RayPolarization polarization = rayPolarization(),
          RayWavefront wavefronts = rayWavefront(),
          bool draw_reflections = config.rendering.draw_reflections) {
    Renderer r = new Renderer;
    r.p = p;
    r.before = before;
    r.after = after;
    r.margin = margin;
    r.arrows = arrows;
    r.polarization = polarization;
    r.wavefronts = wavefronts;
    r.draw_reflections = draw_reflections;
    r.drawTrace(pic, T);
}

// Global Public Trace Bundle Drawing Wrapper
void draw(picture pic = currentpicture,
          Trace[] traces,
          real before = nan,
          real after = nan,
          real margin = 0,
          pen p = currentpen,
          RayArrows arrows = rayArrows(),
          RayPolarization polarization = rayPolarization(),
          RayWavefront wavefronts = rayWavefront(),
          bool draw_reflections = config.rendering.draw_reflections) {
    // 1. Draw wavefronts before rays so rays and arrows remain legible on top.
    // Phase-locked bundles draw common OPL contours instead of crest bars.
    RayWavefront individual_wavefronts = wavefronts;
    // Wavelet arcs are intentionally retained as a visual Huygens overlay;
    // only ordinary per-ray crest bars are replaced by common OPL contours.
    if (wavefronts.phase_locked) individual_wavefronts.visible = false;
    if (individual_wavefronts.visible) {
        for (Trace t : traces) {
            draw(pic, t, before, after, margin, nullpen, rayArrows(visible = false),
                 rayPolarization(visible = false), individual_wavefronts, draw_reflections);
        }
    }

    if (wavefronts.phase_locked && wavefronts.wavelets) {
        drawPhaseLockedWavelets(pic, traces, wavefronts);
    }

    RayWavefront no_wavefronts = wavefronts;
    no_wavefronts.visible = false;
    for (Trace t : traces) {
        draw(pic, t, before, after, margin, p, arrows, polarization, no_wavefronts, draw_reflections);
    }
    // 2. Bundle iso-OPL contours
    if ((wavefronts.contours || wavefronts.phase_locked) && !wavefronts.wavelets && traces.length >= 2) {
        drawBundleContours(pic, traces, wavefronts);
    }
    // 3. Interference lines between the first two traces
    if ((wavefronts.antinodes || wavefronts.nodes) && traces.length >= 2) {
        drawBundleInterference(pic, traces[0], traces[1], wavefronts);
    }
}

void draw(picture pic = currentpicture, Element e) {
    Renderer r = new Renderer;
    r.drawElement(pic, e);
}

void draw(picture pic = currentpicture, Element[] elements) {
    for (int i = 0; i < elements.length; ++i) draw(pic, elements[i]);
}

void draw(picture pic = currentpicture, OpticalSystem sys) { draw(pic, sys.elements); }
void draw(picture pic = currentpicture, Segment s, pen p = currentpen) { draw(pic, s.a -- s.b, p); }
void draw(picture pic = currentpicture, Segment[] segments, pen p = currentpen) {
    for (int i = 0; i < segments.length; ++i) draw(pic, segments[i], p);
}

// Standalone Ray Drawer with Viewport Clipping for infinity
void draw(picture pic = currentpicture,
          Label L = Label(""),
          Ray r,
          real length = infinity,
          real margin = 0,
          pen p = currentpen,
          RayArrows arrows = rayArrows(),
          pair anchor = (infinity, infinity),
          RayPolarization polarization = rayPolarization(),
          RayWavefront wavefronts = rayWavefront()) {
    pair z = (abs(anchor.x) == infinity || abs(anchor.y) == infinity) ? r.p : projection(r, anchor);
    ResolvedArrow arrow = arrows.getArrow(0, 1);
    pen pen_eff = resolveRayPen(p, r, config.rendering.use_intensity_opacity);
    if (length == infinity) {
        drawClippedRay(pic, z, r.d, margin, pen_eff, arrow, false, true);
    } else if (length == -infinity) {
        drawClippedRay(pic, z, r.d, margin, pen_eff, arrow, true, false);
    } else if (length >= 0) {
        drawRaySegment(pic, L, r, z -- (z + length * r.d), 0, 1, pen_eff, arrow, polarization, wavefronts, false);
    } else {
        drawRaySegment(pic,
                       L,
                       r,
                       (z + length * r.d) -- z,
                       0,
                       1,
                       pen_eff,
                       arrow,
                       polarization,
                       wavefronts,
                       false);
    }
}

// ============================================================================
// PART II: FACTORIES, SOLVERS & EXTENSIONS
// ============================================================================
// ----------------------------------------------------------------------------
// Section 7: Optical Component & Source Factories
// ----------------------------------------------------------------------------
// Beam Generator Factories
struct Beam {
    void operator init() {}
    static Ray[] parallel(pair start,
                          pair step,
                          int count,
                          pair d,
                          Medium medium = null, // null: resolved at the start point by the tracer
                          real lambda = config.physics.lambda,
                          Polarization pol = Polarization.unpolarized(),
                          real phase = 0) {
        Ray[] rays;
        for (int i = 0; i < count; ++i) {
            rays.push(Ray(start + i * step, d, medium, lambda, pol, phase));
        }
        return rays;
    }
    static Ray[] diverging(pair source,
                           pair d,
                           real angle,
                           int count,
                           Medium medium = null, // null: resolved at the source by the tracer
                           real lambda = config.physics.lambda,
                           Polarization pol = Polarization.unpolarized(),
                           real phase = 0) {
        Ray[] rays;
        if (count <= 1) {
            rays.push(Ray(source, d, medium, lambda, pol, phase));
            return rays;
        }
        for (int i = 0; i < count; ++i) {
            real a = -angle / 2 + angle * i / (count - 1);
            rays.push(Ray(source, rotate(a) * d, medium, lambda, pol, phase));
        }
        return rays;
    }
}

// Surface Property: Refract (Updated with Complex Fresnel Phase Retardance)
Property Refract(Medium inside, Medium ambient = null, real reflectance = nan) {
    Property p = new Property;
    p.name = "Refract";

    p.generate = new Ray[](Ray r, Hit h, mediumResolver resolve) {
        FresnelResult res = fresnelIsotropic(inside, ambient, r, h, resolve, reflectance);
        Ray[] rays;

        // 1. Transmitted Ray
        if (!res.tir && r.intensity() * res.T > config.physics.min_intensity) {
            real intensity = r.intensity() * res.T;
            Ray trans_ray = spawnRay(r, h, res.d, intensity, res.medium, 0, res.wave_normal);

            if (isnan(reflectance)) {
                applyDiattenuation(trans_ray.pol, r.pol, res.Ts, res.Tp);
            }
            rays.push(trans_ray);
        }

        // 2. Reflected Ray (handles TIR and specular reflection)
        if (r.intensity() * res.R > config.physics.min_intensity &&
            (res.tir || config.physics.trace_reflections || config.rendering.draw_reflections)) {
            real intensity = r.intensity() * res.R;
            real refl_phase = atan2(res.rs.y, res.rs.x); // pi for external reflection (n1 < n2)
            Ray refl_ray = spawnRay(r,
                                    h,
                                    res.reflected_direction,
                                    intensity,
                                    res.source_medium,
                                    refl_phase,
                                    res.reflected_wave_normal);

            if (isnan(reflectance)) {
                // Mueller reflection transformation for Stokes vector
                applyDiattenuation(refl_ray.pol, r.pol, res.rs, res.rp);
            }
            rays.push(refl_ray);
        }

        return rays;
    };

    return p;
}

// Ideal 2D polarizing splitter: independent scalar ordinary/extraordinary
// modes. This is deliberately not presented as a general anisotropic solver.
Property DoubleRefraction(Medium med_o, Medium med_e, Medium ambient = null) {
    Property p = new Property;
    p.name = "DoubleRefraction";
    p.generate = new Ray[](Ray r, Hit h, mediumResolver resolve) {
        Ray[] rays;

        // Resolve interface frame
        InterfaceFrame f = interfaceFrame(med_o, ambient, r, h, resolve);

        // ---------------------------------------------------------------------
        // Definitive test: If the ray is NOT already in med_o or med_e,
        // it is ENTERING the crystal from the outside!
        // ---------------------------------------------------------------------
        bool is_entering = (r.medium != med_o && r.medium != med_e);

        if (is_entering) {
            // =================================================================
            // CASE 1: ENTERING CRYSTAL (External Ray splits into O and E modes)
            // =================================================================
            real Is_inc = max(0.0, 0.5 * (r.pol.S0 + r.pol.S1));
            real Ip_inc = max(0.0, 0.5 * (r.pol.S0 - r.pol.S1));

            FresnelResult b = fresnelUniaxial(r, h, f.source_medium, med_o, med_e, resolve);

            // 1. Ordinary Ray (Pure S-polarization in med_o)
            if (!b.tir_o) {
                real I_o = Is_inc * b.To;
                if (I_o > config.physics.min_intensity) {
                    Ray o_ray = spawnRay(r, h, b.d_o, I_o, med_o, 0, b.wave_normal_o);
                    o_ray.pol.S0 = I_o;
                    o_ray.pol.S1 = I_o;
                    o_ray.pol.S2 = 0.0;
                    o_ray.pol.S3 = 0.0;
                    rays.push(o_ray);
                }
            }

            // 2. Extraordinary Ray (Pure P-polarization in med_e with walk-off)
            if (!b.tir_e) {
                real I_e = Ip_inc * b.Te;
                if (I_e > config.physics.min_intensity) {
                    Ray e_ray = spawnRay(r, h, b.d_e, I_e, med_e, 0, b.wave_normal_e);
                    e_ray.pol.S0 = I_e;
                    e_ray.pol.S1 = -I_e;
                    e_ray.pol.S2 = 0.0;
                    e_ray.pol.S3 = 0.0;
                    rays.push(e_ray);
                }
            }

            // 3. Reflected Ray (Reflecting back into source medium)
            if (config.physics.trace_reflections || config.rendering.draw_reflections) {
                real I_refl = b.Rs * Is_inc + b.Rp * Ip_inc;
                if (I_refl > config.physics.min_intensity) {
                    Ray refl_ray =
                        spawnRay(r, h, b.reflected_direction, I_refl, f.source_medium, 0, b.reflected_wave_normal);
                    pair C = b.rs * conj(b.rp);

                    refl_ray.pol.S0 = I_refl;
                    refl_ray.pol.S1 = b.Rs * Is_inc - b.Rp * Ip_inc;
                    refl_ray.pol.S2 = C.x * r.pol.S2 - C.y * r.pol.S3;
                    refl_ray.pol.S3 = C.y * r.pol.S2 + C.x * r.pol.S3;
                    rays.push(refl_ray);
                }
            }
        } else {
            // =================================================================
            // CASE 2: EXITING CRYSTAL (Ray is already in med_o or med_e)
            // =================================================================
            Medium src_med = (r.medium != null) ? r.medium : med_o;
            FresnelResult res = fresnelIsotropic(src_med, ambient, r, h, resolve, nan);

            if (!res.tir) {
                real I_trans = r.intensity() * res.T;
                if (I_trans > config.physics.min_intensity) {
                    Ray trans_ray = spawnRay(r, h, res.d, I_trans, res.medium, 0, res.wave_normal);
                    trans_ray.pol = copy(r.pol);
                    trans_ray.pol.scaleTo(I_trans);
                    rays.push(trans_ray);
                }
            }

            // Internal reflection (e.g. TIR of ordinary ray in Nicol prism)
            if (config.physics.trace_reflections || config.rendering.draw_reflections || res.tir) {
                real I_refl = r.intensity() * res.R;
                if (I_refl > config.physics.min_intensity) {
                    Ray refl_ray = spawnRay(r,
                                            h,
                                            res.reflected_direction,
                                            I_refl,
                                            res.source_medium,
                                            0,
                                            res.reflected_wave_normal);
                    refl_ray.pol = copy(r.pol);
                    refl_ray.pol.scaleTo(I_refl);
                    rays.push(refl_ray);
                }
            }
        }

        return rays;
    };
    return p;
}

// Ideal (perfect-conductor) mirror in the Fresnel s/p convention used by
// fresnelIsotropic: rs = -1, rp = +1, i.e. Mueller diag(1, 1, -1, -1).
// This is the n2 -> infinity limit of dielectric reflection and reverses
// the handedness of circular light.
void applyMirrorPolarization(Polarization pol) {
    if (pol == null) return;
    pol.S2 = -pol.S2;
    pol.S3 = -pol.S3;
}

Property Reflect(real phase_shift = pi) {
    Property p = new Property;
    p.name = "Reflect";
    p.generate = new Ray[](Ray r, Hit h, mediumResolver resolve) {
        pair N = phaseGradientNormal(r, h);
        pair reflected_wave_normal = unitOr(r.wave_normal - 2 * dot(r.wave_normal, N) * N, r.d);
        pair reflected = r.medium.groupDirection(reflected_wave_normal);
        Ray out = spawnRay(r, h, reflected, r.intensity(), r.medium, phase_shift, reflected_wave_normal);
        applyMirrorPolarization(out.pol);
        return new Ray[]{out};
    };
    // Specular Reflection alignment: Bisector normal!
    p.align = new pair(pair in_dir, pair target_dir, real lambda, int order_index) { return target_dir - in_dir; };
    return p;
}

Property IdealReflect(real f, real phase_shift = pi) {
    Property p = new Property;
    p.name = "IdealReflect";
    p.generate = new Ray[](Ray r, Hit h, mediumResolver resolve) {
        if (abs(f) <= config.math.eps) return new Ray[];
        pair A = unitOr(h.axis); // Inward optical axis
        // Theoretical focus F
        pair F = h.center + f * A;
        pair u_F = unitOr(F - h.p, A);
        // Exact parabolic bisector normal at the curved hit point h.p
        pair N_ideal = unitOr(-A - u_F, A);
        // Exact specular reflection across the ideal normal
        pair v_out = unitOr(r.d - 2.0 * dot(r.d, N_ideal) * N_ideal, r.d);
        Ray out = spawnRay(r, h, v_out, r.intensity(), r.medium, phase_shift);
        applyMirrorPolarization(out.pol);
        return new Ray[]{out};
    };
    p.align = new pair(pair in_dir, pair target_dir, real lambda, int order_index) { return target_dir - in_dir; };
    return p;
}

Property Focus(real f) {
    Property p = new Property;
    p.name = "Focus";
    p.stage = PropertyStage.Deflect;
    p.generate = new Ray[](Ray r, Hit h, mediumResolver resolve) {
        pair v_out = thinLensRefractDir(r, h, f);
        Ray out = spawnRay(r, h, v_out, r.intensity(), r.medium);
        out.wave_normal = v_out;
        return new Ray[]{out};
    };
    return p;
}

Property IdealFocus(real f) {
    Property p = new Property;
    p.name = "IdealFocus";
    p.stage = PropertyStage.Deflect;
    p.generate = new Ray[](Ray r, Hit h, mediumResolver resolve) {
        pair A = unitOr(h.axis);
        pair F = h.center + f * A;
        pair v_out = unitOr(F - h.p, r.d);
        Ray out = spawnRay(r, h, v_out, r.intensity(), r.medium);
        out.wave_normal = v_out;
        return new Ray[]{out};
    };
    return p;
}

private
pair wavevectorShiftAxis(pair in_dir, pair target_dir, real delta_k, bool reflective) {
    pair v = target_dir - in_dir;
    real v_mag = length(v);
    if (v_mag < config.physics.wavevector_alignment_tolerance || abs(delta_k) > v_mag) {
        return (infinity, infinity);
    }
    real angle_t = aCos(clamp(delta_k / v_mag, -1, 1));
    pair v_hat = v / v_mag;
    pair v_rot1 = rotate(angle_t) * v_hat;
    pair N1 = (-v_rot1.y, v_rot1.x);
    pair v_rot2 = rotate(-angle_t) * v_hat;
    pair N2 = (-v_rot2.y, v_rot2.x);
    real side_product = dot(in_dir, N1) * dot(target_dir, N1);
    if (reflective) return side_product < 0 ? N1 : N2;
    return side_product > 0 ? N1 : N2;
}

Property
Diffract(real lines_per_mm, int[] orders = new int[]{1}, bool reflective = false, real[] efficiencies = new real[]{}) {
    Property p = new Property;
    p.name = "Diffract";
    real grating_spacing = 1000.0 / lines_per_mm;
    // Generation Closure
    p.generate = new Ray[](Ray r, Hit h, mediumResolver resolve) {
        Ray[] validRays;
        bool use_efficiencies = (efficiencies.length == orders.length);
        pair N = phaseGradientNormal(r, h);
        pair T = (-N.y, N.x);
        Medium m1 = (r.medium != null) ? r.medium : Materials.VACUUM;
        real n1 = m1.n(h.p, r.d, r.lambda, r.pol);
        Medium m2 = m1;
        if (!reflective) {
            pair target_point = h.p + config.geometry.eps_intersection * r.d;
            m2 = resolve(target_point, r, null, null, true);
            if (m2 == null) m2 = Materials.VACUUM;
        }
        real n2 = reflective ? n1 : m2.n(h.p, r.d, r.lambda, r.pol);
        real k_in = dot(r.d, T);
        real lambda_um = r.lambda * nm_to_um;
        real phase_shift = reflective ? pi : 0;
        for (int i = 0; i < orders.length; ++i) {
            real delta_k = orders[i] * lambda_um / grating_spacing;
            real k_out = (n1 * k_in + delta_k) / n2;
            if (abs(k_out) <= 1.0) {
                pair v_out = diffractionOutputDir(N, T, k_out, reflective);
                real eff = use_efficiencies ? max(0.0, efficiencies[i]) : 1.0;
                real intensity = r.intensity() * eff;
                Ray diff_ray = spawnRay(r, h, v_out, intensity, m2, phase_shift);
                validRays.push(diff_ray);
            } else {
                // Optional diagnostic:
                write(format("Notice: Order m = %d", orders[i]) + format(" is evanescent at lambda = %.1f", r.lambda) +
                      format(" (|k_out| = %.2f > 1.0)", abs(k_out)));
            }
        }
        if (validRays.length > 0 && !use_efficiencies) {
            real split_intensity = r.intensity() / validRays.length;
            for (int i = 0; i < validRays.length; ++i) {
                validRays[i].pol.scaleTo(split_intensity);
            }
        }
        return validRays;
    };
    // Alignment Closure (Decouples align() completely!)
    p.align = new pair(pair in_dir, pair target_dir, real lambda, int order_index) {
        int target_order = (orders.length > order_index) ? orders[order_index] : 0;
        real delta_k = target_order * (lambda * nm_to_um) / grating_spacing;
        return wavevectorShiftAxis(in_dir, target_dir, delta_k, reflective);
    };
    return p;
}

Property PhaseShift(real delta_k, bool reflective = false) {
    Property p = new Property;
    p.name = "PhaseShift";
    p.generate = new Ray[](Ray r, Hit h, mediumResolver resolve) {
        pair N = phaseGradientNormal(r, h);
        pair T = (-N.y, N.x);
        real k_in = dot(r.d, T);
        // real k_out = k_in + delta_k; // Initial version - without medium correction
        // real n = r.medium.n(r.p, r.d, r.lambda, r.pol); // Get the local index of the medium
        // real k_out = k_in + delta_k / n;                   // Scale the transverse wavevector s
        Medium m1 = (r.medium != null) ? r.medium : Materials.VACUUM;
        real n1 = m1.n(h.p, r.d, r.lambda, r.pol);
        Medium m2 = reflective ? m1 : resolve(h.p + config.geometry.eps_intersection * r.d, r, null, null, true);
        real n2 = m2.n(h.p, r.d, r.lambda, r.pol);
        real k_out = (n1 * dot(r.d, T) + delta_k) / n2;
        if (abs(k_out) > 1.0) return new Ray[];
        pair v_out = diffractionOutputDir(N, T, k_out, reflective);
        real p_shift = reflective ? pi : 0;
        return new Ray[]{spawnRay(r, h, v_out, r.intensity(), reflective ? r.medium : m2, p_shift)};
    };
    p.align = new pair(pair in_dir, pair target_dir, real lambda, int order_index) {
        return wavevectorShiftAxis(in_dir, target_dir, delta_k, reflective);
    };
    return p;
}

Property Absorb() {
    Property p = new Property;
    p.name = "Absorb";
    p.generate = new Ray[](Ray r, Hit h, mediumResolver resolve) { return new Ray[]; };
    return p;
}

// Universal Refractive Interface Primitive
Element Interface(path boundary,
                  Medium inside,
                  Medium ambient = null,
                  pen stroke = black,
                  pen fill = nullpen,
                  string name = "Interface") {
    Property[] props = {Refract(inside, ambient)};
    return element(boundary,
                   props,
                   stroke,
                   fill,
                   name,
                   boundaryCenter(boundary),
                   boundaryAxis(boundary),
                   inside); // <-- Passes inside directly
}

// Master Lens Factories
struct LensPath {
    // Helper: Calculate the maximum physically realizable diameter for given curvatures and thickness
    static real thickLensMaxDiameter(real R1, real R2, real thickness, real min_edge_thickness = 0.0) {
        if (thickness <= config.math.eps) return 0.0;
        // 1. Hemisphere bound: h cannot exceed sphere radius
        real h_max = infinity;
        if (abs(R1) > config.math.eps) h_max = min(h_max, config.geometry.hemisphere_margin * abs(R1));
        if (abs(R2) > config.math.eps) h_max = min(h_max, config.geometry.hemisphere_margin * abs(R2));
        // Helper to evaluate edge thickness at height y
        real edgeThick(real y) {
            real s1 = (abs(R1) > config.math.eps) ? (abs(R1) - sqrt(max(0.0, R1 * R1 - y * y))) : 0.0;
            real s2 = (abs(R2) > config.math.eps) ? (abs(R2) - sqrt(max(0.0, R2 * R2 - y * y))) : 0.0;
            real dx1 = (R1 > config.math.eps) ? s1 : ((R1 < -config.math.eps) ? -s1 : 0.0);
            real dx2 = (R2 < -config.math.eps) ? -s2 : ((R2 > config.math.eps) ? s2 : 0.0);
            return thickness + dx2 - dx1;
        }
        // 2. Crossover bound: edge thickness must be >= min_edge_thickness
        if (h_max < infinity && edgeThick(h_max) < min_edge_thickness) {
            // Fast bisection to find knife-edge height where edgeThick(y) == min_edge_thickness
            real y_low = 0.0;
            real y_high = h_max;
            for (int i = 0; i < config.geometry.max_bisection_steps; ++i) {
                real y_mid = 0.5 * (y_low + y_high);
                if (edgeThick(y_mid) < min_edge_thickness)
                    y_high = y_mid;
                else
                    y_low = y_mid;
            }
            h_max = y_low;
        }
        return (h_max < infinity) ? (2.0 * h_max) : infinity;
    }
    // Robust thick lens boundary path generator
    static path thickLens(real R1, real R2, real thickness, real diameter, real min_edge_thickness = 0.0) {
        real t_safe = max(thickness, config.geometry.eps);
        real max_d = thickLensMaxDiameter(R1, R2, t_safe, min_edge_thickness);
        // Auto-clamp requested diameter to physical maximum:
        real eff_diameter = min(diameter, max_d);
        real h = eff_diameter / 2.0;
        pair front_vertex = (-t_safe / 2.0, 0);
        pair back_vertex = (t_safe / 2.0, 0);
        // 1. FRONT SURFACE ARC (curves from y = -h to y = +h)
        path front_arc;
        if (R1 > config.math.eps) {
            pair C1 = front_vertex + (R1, 0);
            real a1 = aSin(clamp(h / R1, -1.0, 1.0));
            front_arc = reverse(Arc(C1, R1, 180 - a1, 180 + a1, config.geometry.arc_subdivisions));
        } else if (R1 < -config.math.eps) {
            pair C1 = front_vertex + (R1, 0);
            real a1 = aSin(clamp(h / abs(R1), -1.0, 1.0));
            front_arc = Arc(C1, abs(R1), -a1, a1, config.geometry.arc_subdivisions);
        } else {
            front_arc = (front_vertex.x, -h) -- (front_vertex.x, h);
        }
        // 2. REAR SURFACE ARC (curves from y = -h to y = +h)
        path back_arc;
        if (R2 < -config.math.eps) {
            pair C2 = back_vertex + (R2, 0);
            real a2 = aSin(clamp(h / abs(R2), -1.0, 1.0));
            back_arc = Arc(C2, abs(R2), -a2, a2, config.geometry.arc_subdivisions);
        } else if (R2 > config.math.eps) {
            pair C2 = back_vertex + (R2, 0);
            real a2 = aSin(clamp(h / R2, -1.0, 1.0));
            back_arc = reverse(Arc(C2, R2, 180 - a2, 180 + a2, config.geometry.arc_subdivisions));
        } else {
            back_arc = (back_vertex.x, -h) -- (back_vertex.x, h);
        }
        // 3. ASSEMBLE CLOSED NON-SELF-INTERSECTING BOUNDARY
        // front_arc goes (-h -> +h), reverse(back_arc) goes (+h -> -h)
        return front_arc -- reverse(back_arc) -- cycle;
    }
}

// Default Surface Drawing
// Native Asymptote Arrowhead for Thin Lenses (Convex & Concave)
arrowhead LensHead(bool concave = false, real angle = 20, real bow = 0.3, real barb = 0.5) {
    arrowhead ah = new arrowhead;
    ah.head = new path(path g, position position = EndPoint, pen p = currentpen, real size = 0, real angle = angle) {
        if (size == 0) size = config.rendering.min_lens_head > 0 ? config.rendering.min_lens_head : 5.0;
        // Resolve path endpoint P0 and outward unit tangent d
        bool relative = position.relative;
        real pos = position.position.x;
        real t = relative ? reltime(g, pos) : pos;
        pair P0 = point(g, t);
        pair d = unit(dir(g, t));
        // Orientation:
        pair arrow_dir = concave ? -d : d;
        pair n = (-arrow_dir.y, arrow_dir.x);
        real w = size * tan(angle * pi / 180.0);
        // Construct standard hook shape locally with tip at (0, 0)
        pair M = (0, 0);
        pair tip = size * barb * arrow_dir;
        pair p_l = size * (barb - 1) * arrow_dir + w * n;
        pair p_r = size * (barb - 1) * arrow_dir - w * n;
        pair c_l = size * (barb - 0.5) * arrow_dir + (w * bow) * n;
        pair c_r = size * (barb - 0.5) * arrow_dir - (w * bow) * n;
        path hook_path = (concave) ? M -- p_l .. c_l .. tip -- tip .. c_r .. p_r -- cycle
                                   : (-tip) -- M -- p_l .. c_l .. tip -- tip .. c_r .. p_r -- M -- cycle;
        // Attachment Point:
        pair attach_point = (concave) ? -(size * barb) * arrow_dir : (size * barb) * arrow_dir;
        return shift(P0 - attach_point) * hook_path;
    };
    return ah;
}

// Convenient pre-configured arrowheads:
arrowhead LensConvexHead = LensHead(concave = false);
arrowhead LensConcaveHead = LensHead(concave = true);

// Thin Lens Drawing Delegate (Arrowhead hooks & Focal points)
void drawThinLens(picture pic, Element e) {
    if (e.surfaces.length == 0) return;
    Surface lens_surface = transformedSurface(e.surfaces[0], e.T);
    pair center = lens_surface.center;
    pair axis = lens_surface.axis;
    pair q = (-axis.y, axis.x);
    real radius = 0;
    if (length(lens_surface.p) > 0) {
        path global_path = lens_surface.T * lens_surface.p;
        radius = abs(point(global_path, length(global_path)) - point(global_path, 0)) / 2.0;
    }
    pair a = center - radius * q;
    pair b = center + radius * q;
    real head_size = config.rendering.min_lens_head;
    real f_val = !isnan(e.focal_length) ? e.focal_length : 0;
    // Single declarative draw line using native Asymptote Arrows!
    draw(pic,
         a -- b,
         e.stroke,
         Arrows(LensHead(concave = (f_val < 0)), size = head_size, filltype = FillDraw));
    // Optional focal points F1, F2
    if (e.show_info) {
        pair F1 = center - f_val * axis;
        pair F2 = center + f_val * axis;
        dot(pic, F1, e.stroke);
        dot(pic, F2, e.stroke);
        label(pic, "$F_1$", F1, -q, e.stroke);
        label(pic, "$F_2$", F2, -q, e.stroke);
    }
}

// Grating Drawing Delegate (Tick marks)
void drawGrating(picture pic = currentpicture, Element e) {
    if (e.surfaces.length == 0) return;
    Surface grating_surface = transformedSurface(e.surfaces[0], e.T);
    path global_path = grating_surface.T * grating_surface.p;
    real L_total = arclength(global_path);
    if (L_total <= config.rendering.min_drawable_length) return;
    int nlines = config.rendering.grating_tick_count;
    real tick_length = config.rendering.grating_tick_length;
    real alpha = config.rendering.grating_tick_angle;
    pair tangent = unitOr(dir(global_path, 0));
    pair nv = -unitOr(grating_surface.axis);
    pair tick_dir = rotate(alpha) * nv;
    pair tick_vector = -tick_length * tick_dir;
    real proj = dot(tick_vector, tangent);
    real margin = abs(proj);
    real s_start = margin;
    real s_end = L_total - margin;
    real range = s_end - s_start;
    if (range > 0 && nlines > 0) {
        real step_l = range / nlines;
        for (int i = 0; i <= nlines; ++i) {
            real s = s_start + i * step_l;
            real t = arctime(global_path, s);
            pair pt = point(global_path, t);
            draw(pic, pt -- pt + tick_vector, e.stroke + 0.5 * linewidth(e.stroke));
        }
    }
    draw(pic, global_path, e.stroke);
}

struct Lens {
    static Element thin(pair center,
                        real angle,
                        real f,
                        real diameter,
                        pen stroke = config.rendering.pen_lens,
                        string name = "Thin Lens") {
        pair A = dir(angle);
        path boundary = linearBoundary(center, A, diameter);
        Property[] props = {Focus(f)};
        Element e = element(boundary, props, stroke, nullpen, name, center, A);
        e.focal_length = f;
        e.custom_draw = drawThinLens;
        return e;
    }

    static Element thick(pair center,
                         real angle,
                         real R1,
                         real R2,
                         real thickness,
                         real diameter,
                         Medium glass,
                         Medium ambient = null,
                         pen stroke = config.rendering.pen_lens,
                         pen fill = nullpen,
                         string name = "Thick Lens") {
        path local_p = LensPath.thickLens(R1, R2, thickness, diameter);
        transform T = shift(center) * rotate(angle);
        Element e = Interface(T * local_p, glass, ambient, stroke, fill, name);

        // Explicitly assign the real optical axis and center (do not let Interface guess!):
        e.center = center;
        e.axis = dir(angle);
        if (e.surfaces.length > 0) {
            e.surfaces[0].center = center;
            e.surfaces[0].axis = dir(angle);
        }
        return e;
    }

    static Element biconvex(pair center,
                            real angle,
                            real R1,
                            real R2,
                            real thickness,
                            real diameter,
                            Medium glass,
                            Medium ambient = null,
                            pen stroke = config.rendering.pen_lens,
                            pen fill = nullpen,
                            string name = "Biconvex Lens") {
        return thick(center, angle, R1, -R2, thickness, diameter, glass, ambient, stroke, fill, name);
    }

    static Element biconcave(pair center,
                             real angle,
                             real R1,
                             real R2,
                             real thickness,
                             real diameter,
                             Medium glass,
                             Medium ambient = null,
                             pen stroke = config.rendering.pen_lens,
                             pen fill = nullpen,
                             string name = "Biconcave Lens") {
        return thick(center, angle, -R1, R2, thickness, diameter, glass, ambient, stroke, fill, name);
    }

    static Element custom(path boundary,
                          Medium glass,
                          Medium ambient = null,
                          pen stroke = config.rendering.pen_lens,
                          pen fill = nullpen,
                          string name = "Custom Lens") {
        return Interface(boundary, glass, ambient, stroke, fill, name);
    }
}

// Master Mirror Factories
struct Mirror {
    static Element flat(pair center,
                        real angle,
                        real diameter,
                        pen stroke = config.rendering.pen_mirror,
                        string name = "Flat Mirror") {
        pair A = dir(angle);
        path boundary = linearBoundary(center, A, diameter);
        Property[] props = {Reflect()};
        return element(boundary, props, stroke, nullpen, name, center, A);
    }

    // TRUE Physical Spherical Mirror (Has real spherical aberration!)
    static Element spherical(pair vertex,
                             real angle,
                             real f,
                             real diameter,
                             pen stroke = config.rendering.pen_mirror,
                             string name = "Spherical Mirror") {
        real R = 2 * abs(f);
        pair N = dir(angle);
        pair C = vertex + sgn(f) * R * N;
        real half = aSin(min(1.0, diameter / (2 * R)));
        real theta = degrees(vertex - C);
        path boundary = Arc(C, R, theta - half, theta + half, config.geometry.arc_subdivisions);
        Property[] props = {Reflect()};
        return element(boundary, props, stroke, nullpen, name, vertex, N);
    }

    // IDEAL Paraxial Spherical Mirror (Zero spherical aberration)
    static Element idealSpherical(pair vertex,
                                  real angle,
                                  real f,
                                  real diameter,
                                  pen stroke = config.rendering.pen_mirror,
                                  string name = "Ideal Spherical Mirror") {
        real R = 2 * abs(f);
        pair N = dir(angle);
        pair C = vertex + sgn(f) * R * N;
        real half = aSin(min(1.0, diameter / (2 * R)));
        real theta = degrees(vertex - C);
        path boundary = Arc(C, R, theta - half, theta + half, config.geometry.arc_subdivisions);
        Property[] props = {IdealReflect(f)};
        return element(boundary, props, stroke, nullpen, name, vertex, N);
    }

    static Element custom(path boundary, pen stroke = config.rendering.pen_mirror, string name = "Custom Mirror") {
        Property[] props = {Reflect()};
        return element(boundary, props, stroke, nullpen, name, boundaryCenter(boundary), boundaryAxis(boundary));
    }
}

// Beam Splitter Factory Namespace
struct BeamSplitter {
    // Plate Beam Splitter (Single flat partially-reflective surface)
    static Element plate(pair center,
                         real angle,
                         real diameter,
                         real R = 0.5,
                         Medium left = null,  // medium.ambient = null; default
                         Medium right = null, // medium.ambient = null; default
                         pen stroke = gray(0.3),
                         string name = "Plate Beam Splitter") {
        // angle = tilt angle (e.g. 45 deg tilts the plate to reflect horizontal rays UPWARD at 90 deg!)
        pair A = dir(-angle);
        path boundary = linearBoundary(center, A, diameter);
        // Attach a Refract property with a custom 'R' override!
        Property[] props = {Refract(left, right, R)};
        return element(boundary, props, stroke, nullpen, name, center, dir(-angle));
    }

    // Cube Beam Splitter (Glass cube with internal beam-splitting coating)
    static Element cube(pair center,
                        real angle,
                        real size,
                        real R = 0.5,
                        Medium glass = Materials.BK7,
                        Medium ambient = null,
                        pen stroke = gray(0.3),
                        pen fill = 0.08paleblue + white,
                        string name = "Cube Beam Splitter") {
        real h = size / 2.0;
        // Local boundaries
        path local_box = (-h, -h) -- (h, -h) -- (h, h) -- (-h, h) -- cycle;
        path local_diag = (-h, -h) -- (h, h); // Bottom-left to top-right (Reflects +X rays UPWARD)
        transform T = shift(center) * rotate(angle);
        Element e = new Element;
        e.name = name;
        e.stroke = stroke;
        e.fill = fill;
        e.center = center;
        e.axis = dir(angle);
        e.boundary = T * local_box;
        e.medium = glass;
        // Surface 1: Outer Glass Cube Boundary (Glass <-> Ambient)
        Surface s_outer = new Surface;
        s_outer.name = name + " Outer";
        s_outer.p = ensureCCW(T * local_box);
        s_outer.stroke = stroke;
        s_outer.fill = fill;
        s_outer.center = center;
        s_outer.axis = e.axis;
        s_outer.add(Refract(glass, ambient)); // Standard Refraction
        e.add(s_outer);
        // Surface 2: Internal Coating (Glass <-> Glass with specified R ratio!)
        Surface s_diag = new Surface;
        s_diag.name = name + " Coating";
        s_diag.p = T * local_diag;
        s_diag.stroke = stroke;
        s_diag.fill = nullpen;
        s_diag.center = center;
        s_diag.axis = e.axis;
        s_diag.add(Refract(glass, glass, R)); // Custom Ratio Refraction!
        e.add(s_diag);
        return e;
    }
}

// Thin polarization optics.  `axis` / `fast_axis` are angles in the
// Stokes (s/p) frame of the transmitted light; `angle` is the direction of
// the plate's optical axis (its surface normal) in the drawing plane.  The
// two are independent: rotating a waveplate's fast axis does not tilt the
// plate.
// General linear retarder: retardance delta (rad), fast axis at fast_axis.
Property Retard(real delta, real fast_axis = 0.0) {
    real ax = radians(fast_axis);
    return property(
        "Retarder",
        new Ray[](Ray r, Hit h, mediumResolver res) {
            r.pol.applyRotation(-ax);
            r.pol.applyRetardance(delta);
            r.pol.applyRotation(ax);
            return new Ray[]{r};
        },
        PropertyStage.Interact);
}

struct Polarizer {
    // 1. Ideal Linear Polarizer Plate (transmits linear light at axis)
    static Element LP(pair center,
                      real axis = 0.0,
                      real diameter = 30.0,
                      real angle = 0.0,
                      real thickness = 2.0,
                      real extinction = 1000.0) {
        // Index-free dichroic medium: no Fresnel loss, only diattenuation.
        Medium pol_med = medium("Polarizer").add(Dichroic(0.0, extinction, axis));
        path plate = box((-thickness / 2, -diameter / 2), (thickness / 2, diameter / 2));
        Element e = Interface(shift(center) * rotate(angle) * plate, pol_med);
        e.name = "Linear Polarizer";
        e.center = center;
        e.axis = dir(angle);
        e.surfaces[0].center = center;
        e.surfaces[0].axis = dir(angle);
        return e;
    }

    // Ideal zero-thickness retarder plate (no Fresnel loss, no medium change).
    static Element retarder(pair center,
                            real delta,
                            real fast_axis = 0.0,
                            real diameter = 30.0,
                            real angle = 0.0,
                            pen stroke = black,
                            string name = "Retarder") {
        pair A = dir(angle);
        Property[] props = {Retard(delta, fast_axis)};
        return element(linearBoundary(center, A, diameter), props, stroke, nullpen, name, center, A);
    }

    // 2. Quarter-Wave Plate (QWP) - Retardance = pi/2
    static Element QWP(pair center, real fast_axis = 0.0, real diameter = 30.0, real angle = 0.0) {
        return retarder(center, pi / 2.0, fast_axis, diameter, angle, name = "QWP");
    }

    // 3. Half-Wave Plate (HWP) - Retardance = pi (Rotates linear polarization by 2*theta)
    static Element HWP(pair center, real fast_axis = 0.0, real diameter = 30.0, real angle = 0.0) {
        return retarder(center, pi, fast_axis, diameter, angle, name = "HWP");
    }
}

// Absorber, Aperture, Grating & Phase Gradient Factories
Element Absorber(path boundary, pen stroke = config.rendering.pen_absorber, string name = "Absorber") {
    Property[] props = {Absorb()};
    return element(boundary, props, stroke, nullpen, name, boundaryCenter(boundary), boundaryAxis(boundary));
}

Element Aperture(pair center, real diameter, real hole, real angle = 0, pen stroke = black, string name = "Aperture") {
    pair ax = dir(angle);
    pair q = (-ax.y, ax.x);
    pair top1 = center + hole / 2 * q;
    pair top2 = center + diameter / 2 * q;
    pair bot1 = center - hole / 2 * q;
    pair bot2 = center - diameter / 2 * q;
    Element e = new Element;
    e.name = name;
    e.center = center;
    e.axis = ax;
    e.add(Absorber(top1 -- top2, stroke, name + " top"));
    e.add(Absorber(bot1 -- bot2, stroke, name + " bottom"));
    return e;
}

Element PhaseGradient(pair center,
                      real diameter,
                      real angle,
                      real delta_k,
                      bool reflective = false,
                      pen stroke = lightblue,
                      string name = "Phase Gradient") {
    pair ax = dir(angle);
    path boundary = linearBoundary(center, ax, diameter);
    Property[] props = {PhaseShift(delta_k, reflective)};
    return element(boundary, props, stroke, nullpen, name, center, ax);
}

Element Grating(pair center,
                real diameter,
                real angle,
                real lines_per_mm,
                bool reflective = false,
                int[] orders = new int[]{1},
                real[] efficiencies = new real[]{},
                pen stroke = lightblue,
                string name = "Grating") {
    pair ax = dir(angle);
    path boundary = linearBoundary(center, ax, diameter);
    Property[] props = {Diffract(lines_per_mm, orders, reflective, efficiencies)};
    Element e = element(boundary, props, stroke, nullpen, name, center, ax);
    e.custom_draw = drawGrating;
    return e;
}

// Generic Two-Mode Birefringent Crystal (Accepts any ordinary + extraordinary media)
Element Crystal(path boundary,
                Medium med_o,
                Medium med_e,
                Medium ambient = null,
                pen stroke = black,
                pen fill = nullpen,
                string name = "Crystal") {
    Property[] props = {DoubleRefraction(med_o, med_e, ambient)};
    return element(boundary,
                   props,
                   stroke,
                   fill,
                   name,
                   boundaryCenter(boundary),
                   boundaryAxis(boundary),
                   med_o,  // Ordinary mode
                   med_e); // Extraordinary mode
}

// Convenience Uniaxial Crystal Factory (Builds med_o and med_e in 1 line)
Element UniaxialCrystal(path boundary,
                        real no,
                        real ne,
                        pair optic_axis = (1, 0),
                        Medium ambient = null,
                        pen stroke = black,
                        pen fill = nullpen,
                        string name = "Uniaxial Crystal") {
    Medium med_o = medium(name + " Ordinary", no);
    Medium med_e = uniaxialExtraordinary(name + " Extraordinary", no, ne, optic_axis);
    return Crystal(boundary, med_o, med_e, ambient, stroke, fill, name);
}

// ----------------------------------------------------------------------------
// Section 8: Geometric Alignment, Solvers & Optical Analysis
// ----------------------------------------------------------------------------
pair lineIntersection(Ray r1, Ray r2) { return lineIntersection(r1.p, r1.d, r2.p, r2.d); }
pair lineIntersection(Ray r, pair line_origin, pair line_dir) {
    return lineIntersection(r.p, r.d, line_origin, line_dir);
}
pair intersection(Ray r1, Ray r2) { return lineIntersection(r1, r2); }

// Dynamic Geometric Alignment API
private
Element alignAxis(Element e, pair target_axis) {
    pair a = unitOr(e.axis);
    pair n = unitOr(target_axis, a);
    real angle = aCos(clamp(dot(a, n), -1, 1));
    if (cross(a, n) < 0) angle = -angle;
    return rotate(angle, e.center) * e;
}

Element align(Element e, pair from_pos, pair to_pos, real lambda = config.physics.lambda, int order_index = 0) {
    pair in_dir = unitOr(e.center - from_pos, e.axis);
    pair target_dir = unitOr(to_pos - e.center, e.axis);
    // Check if any attached surface property defines a custom alignment rule:
    //    - Mirrors -> use specular bisector (target_dir - in_dir)
    //    - Gratings -> use wavevectorShiftAxis()
    //    - Phase Gradients -> use delta_k shift
    for (Surface s : e.surfaces) {
        for (Property prop : s.props) {
            if (prop.align != null) {
                pair target_axis = prop.align(in_dir, target_dir, lambda, order_index);
                if (abs(xpart(target_axis)) != infinity && abs(ypart(target_axis)) != infinity) {
                    return alignAxis(e, target_axis);
                }
            }
        }
    }
    // Default: In-line transmission alignment (Lenses, Apertures, Interfaces)
    return alignAxis(e, in_dir);
}

Element align(Element e, pair from_pos, Element to_element, real lambda = config.physics.lambda, int order_index = 0) {
    return align(e, from_pos, to_element.center, lambda, order_index);
}

Element align(Element e, Element from_element, pair to_pos, real lambda = config.physics.lambda, int order_index = 0) {
    return align(e, from_element.center, to_pos, lambda, order_index);
}

Element
align(Element e, Element from_element, Element to_element, real lambda = config.physics.lambda, int order_index = 0) {
    return align(e, from_element.center, to_element.center, lambda, order_index);
}

path operator cast(Element e) { return e.T * e.boundary; }
path operator cast(Ray r) {
    return (r.p - config.geometry.ray_to_path_extension * unit(r.d)) -- (r.p + config.geometry.ray_to_path_extension * unit(r.d));
}

private
pair visibleDirection(Ray r, pair O) {
    if (r == null) return (0, 0);
    return unitOr(dot(r.p - O, r.d) >= -1e-9 ? r.d : -r.d);
}

real markAcute(Label L = "",
               int n = 1,
               real radius = 20,
               real space = 0,
               pair A,
               pair O,
               pair B,
               pen p = currentpen,
               arrowbar arrow = None,
               bool opposite = false) {
    pair v1 = unitOr(A - O), v2 = unitOr(B - O);
    if (dot(v1, v2) < 0) v2 = -v2;
    if (opposite) {
        v1 = -v1;
        v2 = -v2;
    }
    if (v1.x * v2.y - v1.y * v2.x < 0) {
        pair temp = v1;
        v1 = v2;
        v2 = temp;
    }
    markangle(L, n, radius, space, O + v1, O, O + v2, arrow, p);
    return aCos(clamp(dot(v1, v2), -1.0, 1.0));
}

real markAcute(Label L = "",
               int n = 1,
               real radius = 20,
               real space = 0,
               pair O,
               pair normal,
               Ray r,
               pen p = currentpen,
               arrowbar arrow = None,
               bool opposite = false) {
    return markAcute(L, n, radius, space, O + visibleDirection(r, O), O, O + normal, p, arrow, opposite);
}

real markAcute(Label L = "",
               int n = 1,
               real radius = 20,
               real space = 0,
               Interaction inter,
               int ray_type = Ray.Incident,
               int ray_index = 0,
               pen p = currentpen,
               arrowbar arrow = None,
               bool opposite = false) {
    if (inter == null || inter.hit == null) return 0;
    Ray r = (ray_type == Ray.Incident) ? inter.incident
                                       : ((ray_type == Ray.Reflected && ray_index < inter.reflected.length)
                                              ? inter.reflected[ray_index]
                                              : ((ray_type == Ray.Transmitted && ray_index < inter.transmitted.length)
                                                     ? inter.transmitted[ray_index]
                                                     : null));
    if (r == null) return 0;
    pair O = inter.hit.p;
    pair v_ray = (ray_type == Ray.Incident) ? -unitOr(r.d) : unitOr(r.d);
    return markAcute(L, n, radius, space, O + v_ray, O, O + inter.normal(), p, arrow, opposite);
}

real markAcute(Label L = "",
               int n = 1,
               real radius = 20,
               real space = 0,
               Ray r1,
               Ray r2,
               pen p = currentpen,
               arrowbar arrow = None,
               bool opposite = false) {
    pair O = lineIntersection(r1, r2);
    if (abs(O.x) == infinity) return 0;
    return markAcute(L,
                     n,
                     radius,
                     space,
                     O + visibleDirection(r1, O),
                     O,
                     O + visibleDirection(r2, O),
                     p,
                     arrow,
                     opposite);
}

real markAcute(Label L = "",
               int n = 1,
               real radius = 20,
               real space = 0,
               path a,
               path b,
               int index = 0,
               pen p = currentpen,
               arrowbar arrow = None,
               bool opposite = false) {
    if (length(a) == 0 || length(b) == 0) return 0.0;

    pair O = (nan, nan);
    pair v1 = (0, 0), v2 = (0, 0);

    // 1. Find intersection times directly using Asymptote's intersections()
    real[][] t = intersections(a, b);
    if (t.length > 0) {
        int idx = max(0, min(t.length - 1, index)); // Integer clamp
        real ta = t[idx][0];                        // Time parameter on path a
        real tb = t[idx][1];                        // Time parameter on path b
        O = point(a, ta);
        v1 = unitOr(dir(a, ta)); // Exact local tangent on path a
        v2 = unitOr(dir(b, tb)); // Exact local tangent on path b
    }
    // 2. If segments do not physically touch, compute virtual line intersection:
    else if (length(a) == 1 && length(b) == 1 && straight(a, 0) && straight(b, 0)) {
        O = lineIntersection(point(a, 0), dir(a, 0), point(b, 0), dir(b, 0));
        if (abs(O.x) == infinity) return 0.0; // Parallel lines
        v1 = unitOr(dir(a, 0));
        v2 = unitOr(dir(b, 0));
    } else {
        return 0.0;
    }

    return markAcute(L, n, radius, space, O + v1, O, O + v2, p, arrow, opposite);
}

real getAngle(Interaction inter, int ray_type = Ray.Incident, int ray_index = 0) {
    if (inter == null || inter.hit == null) return 0.0;
    Ray r;
    if (ray_type == Ray.Incident)
        r = inter.incident;
    else if (ray_type == Ray.Reflected) {
        if (ray_index < 0 || inter.reflected.length <= ray_index) return 0.0;
        r = inter.reflected[ray_index];
    } else if (ray_type == Ray.Transmitted) {
        if (ray_index < 0 || inter.transmitted.length <= ray_index) return 0.0;
        r = inter.transmitted[ray_index];
    } else {
        return 0.0;
    }
    if (r == null) return 0.0;
    pair vr = (length(r.p - inter.hit.p) > 1e-4) ? -unit(r.d) : unit(r.d);
    pair vn = unit(inter.normal());
    if (dot(vr, vn) < 0) vn = -vn;
    real cosine_val = max(min(dot(vr, vn), 1.0), -1.0);
    return aCos(cosine_val);
}

// Evaluates signed perpendicular miss distance from target B to the ray exiting the system
private
real rayMissDistance(pair A, real angle, OpticalSystem sys, pair B, real lambda = 0) {
    real lam = (lambda == 0) ? config.physics.lambda : lambda;
    Ray r = Ray.angle(A, angle, lambda = lam);
    Trace t = sys.trace(r);
    Segment[] segs = t.segments();
    if (segs.length == 0) return nan;
    real best_miss = nan;
    real min_distance = infinity;
    // Search ALL segments for the closest approach to target point B
    for (int i = 0; i < segs.length; ++i) {
        Segment s = segs[i];
        if (s.before == null || !s.before.alive) continue;
        pair p1 = s.a;
        pair p2 = s.b;
        pair d_seg = unitOr(p2 - p1, s.before.d);
        real seg_len = length(p2 - p1);
        pair to_b = B - p1;
        real proj = dot(to_b, d_seg);
        bool is_last = (i == segs.length - 1);
        if (proj >= -config.math.eps && (proj <= seg_len + config.math.eps || is_last)) {
            real perp_dist = abs(cross(to_b, d_seg));
            if (perp_dist < min_distance) {
                min_distance = perp_dist;
                best_miss = cross(to_b, d_seg);
            }
        }
    }
    return best_miss;
}

// Robust Brent's Method Ray Aimer
real aimAngle(pair A,
              pair B,
              OpticalSystem sys,
              real angle_min = -45.0,
              real angle_max = 45.0,
              real lambda = 0,
              real tol = 1e-6,
              int max_iter = 40) {
    real a = min(angle_min, angle_max);
    real b = max(angle_min, angle_max);
    real fa = rayMissDistance(A, a, sys, B, lambda);
    real fb = rayMissDistance(A, b, sys, B, lambda);
    // Auto-Bracketing with NaN Safeguards
    if (isnan(fa) || isnan(fb) || fa * fb > 0) {
        int n_samples = 20;
        real step = (b - a) / (real)n_samples;
        real best_ang = a;
        real min_miss = infinity;
        bool found_bracket = false;
        real prev_ang = a;
        real prev_f = fa;
        for (int i = 1; i <= n_samples; ++i) {
            real cur_ang = a + i * step;
            real cur_f = rayMissDistance(A, cur_ang, sys, B, lambda);
            if (!isnan(cur_f)) {
                if (abs(cur_f) < min_miss) {
                    min_miss = abs(cur_f);
                    best_ang = cur_ang;
                }
                if (!isnan(prev_f) && prev_f * cur_f <= 0) {
                    a = prev_ang;
                    fa = prev_f;
                    b = cur_ang;
                    fb = cur_f;
                    found_bracket = true;
                    break;
                }
            }
            prev_ang = cur_ang;
            prev_f = cur_f;
        }
        if (!found_bracket) return best_ang;
    }
    if (isnan(fa) || isnan(fb)) return 0.5 * (a + b);
    // Safe Brent's Algorithm
    if (abs(fa) < abs(fb)) {
        real temp = a;
        a = b;
        b = temp;
        real ftemp = fa;
        fa = fb;
        fb = ftemp;
    }
    real c = a, fc = fa;
    bool mflag = true;
    real d_step = 0;
    for (int iter = 0; iter < max_iter; ++iter) {
        if (abs(fb) < tol || abs(b - a) < tol) return b;
        real s = 0;
        real denom1 = (fa - fb) * (fa - fc);
        real denom2 = (fb - fa) * (fb - fc);
        real denom3 = (fc - fa) * (fc - fb);
        if (fa != fc && fb != fc && abs(denom1) > config.math.solver_denominator_eps &&
            abs(denom2) > config.math.solver_denominator_eps && abs(denom3) > config.math.solver_denominator_eps) {
            // Inverse Quadratic Interpolation
            s = (a * fb * fc) / denom1 + (b * fa * fc) / denom2 + (c * fa * fb) / denom3;
        } else {
            // Secant Method
            real denom_sec = fb - fa;
            if (abs(denom_sec) > config.math.solver_denominator_eps)
                s = b - fb * (b - a) / denom_sec;
            else
                s = 0.5 * (a + b);
        }
        bool cond1 = (s < min((3.0 * a + b) / 4.0, b)) || (s > max((3.0 * a + b) / 4.0, b));
        bool cond2 = mflag && (abs(s - b) >= 0.5 * abs(b - c));
        bool cond3 = !mflag && (abs(s - b) >= 0.5 * abs(c - d_step));
        bool cond4 = mflag && (abs(b - c) < tol);
        bool cond5 = !mflag && (abs(c - d_step) < tol);
        if (cond1 || cond2 || cond3 || cond4 || cond5 || isnan(s)) {
            s = 0.5 * (a + b);
            mflag = true;
        } else {
            mflag = false;
        }
        real fs = rayMissDistance(A, s, sys, B, lambda);
        if (isnan(fs)) {
            s = 0.5 * (a + b);
            fs = rayMissDistance(A, s, sys, B, lambda);
            if (isnan(fs)) return b; // Safe termination
        }
        d_step = c;
        c = b;
        fc = fb;
        if (fa * fs < 0) {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }
        if (abs(fa) < abs(fb)) {
            real temp = a;
            a = b;
            b = temp;
            real ftemp = fa;
            fa = fb;
            fb = ftemp;
        }
    }
    return b;
}

// Global standalone helper (Defined after Ray and OpticalSystem)
Ray aim(pair A, pair B, OpticalSystem sys, real angle_min = -45.0, real angle_max = 45.0, real lambda = 0) {
    real ang = aimAngle(A, B, sys, angle_min, angle_max, lambda);
    real lam = (lambda == 0) ? config.physics.lambda : lambda;
    return Ray.angle(A, ang, lambda = lam);
}

// Element[] array overload
Ray aim(pair A, pair B, Element[] elements, real angle_min = -45.0, real angle_max = 45.0, real lambda = 0) {
    return aim(A, B, system(elements), angle_min, angle_max, lambda);
}

// ============================================================================
// Unified Optical Imaging & Object Path Transformation Engine
// ============================================================================
// Conjugate Image Point of a Single Point P through System/Subsystem
pair imagePoint(OpticalSystem sys, pair P, real delta_angle = config.physics.differential_angle, real lambda = 0) {
    if (sys.elements.length == 0) return P;
    real lam = (lambda == 0) ? config.physics.lambda : lambda;
    pair target_optic = sys.elements[0].center;
    pair d_center = unitOr(target_optic - P, sys.elements[0].axis);
    real theta_center = degrees(d_center);
    Ray r1 = Ray.angle(P, theta_center - delta_angle, lambda = lam);
    Ray r2 = Ray.angle(P, theta_center + delta_angle, lambda = lam);
    Trace t1 = sys.trace(r1);
    Trace t2 = sys.trace(r2);
    Ray exit1 = t1.last();
    Ray exit2 = t2.last();
    if (exit1 == null || exit2 == null || !exit1.alive || !exit2.alive) return (infinity, infinity);
    return lineIntersection(exit1.p, exit1.d, exit2.p, exit2.d);
}

// ============================================================================
// Unified Optical Imaging & Gaussian Ray Engine (Replaces lines 1358-1488)
// ============================================================================
// Cardinal points
// Struct storing the 6 Gaussian cardinal points, vertices, and focal lengths
struct CardinalPoints {
    pair F1, F2; // Front & Rear Focal points
    pair H1, H2; // Front & Rear Principal points
    pair N1, N2; // Front & Rear Nodal points
    pair V1, V2; // Front & Rear Surface Vertices
    real f1;     // Front focal length (H1 -> F1)
    real f2;     // Rear focal length / EFL (H2 -> F2)
    real BFL;    // Back focal length (V2 -> F2)
    real FFL;    // Front focal distance (V1 -> F1)
    pair axis_origin;
    pair axis_dir;
    pair axis_norm;
    bool is_thick;
}

// Numerical Paraxial Solver: Traces forward and backward paraxial rays
// Numerical Paraxial Solver: Traces forward and backward paraxial rays
CardinalPoints cardinalPoints(OpticalSystem sys,
                              pair axis_origin = (0, 0),
                              pair axis_dir = (1, 0),
                              real paraxial_h = config.physics.paraxial_height) {
    CardinalPoints cp;
    cp.axis_origin = axis_origin;
    cp.axis_dir = unit(axis_dir);
    cp.axis_norm = (-cp.axis_dir.y, cp.axis_dir.x); // Transverse axis
    if (sys.elements.length == 0) return cp;

    real launch_distance = 0.5 * config.physics.far;
    real min_proj = infinity;
    real max_proj = -infinity;
    FlatSurfaceRecord[] recs = prepareFlatSurfaces(sys.elements);
    for (FlatSurfaceRecord rec : recs) {
        if (rec.min_bound.x < infinity && rec.max_bound.x > -infinity) {
            pair[] pts = {rec.min_bound,
                          (rec.min_bound.x, rec.max_bound.y),
                          rec.max_bound,
                          (rec.max_bound.x, rec.min_bound.y)};
            for (pair pt : pts) {
                real pr = dot(pt - axis_origin, cp.axis_dir);
                if (pr < min_proj) min_proj = pr;
                if (pr > max_proj) max_proj = pr;
            }
        }
    }
    if (min_proj < infinity && max_proj > -infinity) {
        real extent = max(max_proj - min_proj, 1.0);
        real margin = max(min(20.0, 2.0 * extent), 0.2 * extent);
        real safe_dist = max(abs(min_proj), abs(max_proj)) + margin;
        launch_distance = min(safe_dist, 0.5 * config.physics.far);
    }

    // -------------------------------------------------------------------------
    // 1. RESOLVE VERTICES V1 & V2 (Axial intercepts of front & rear surfaces)
    // -------------------------------------------------------------------------
    Ray r_axis_fwd = Ray(p = axis_origin - launch_distance * cp.axis_dir, d = cp.axis_dir);
    Trace t_axis_fwd = sys.trace(r_axis_fwd, far = config.physics.far);
    if (t_axis_fwd.segments().length > 1) {
        cp.V1 = t_axis_fwd.segments()[1].a; // First surface hit point
    } else {
        cp.V1 = sys.elements[0].center;
    }

    Ray r_axis_bwd = Ray(p = axis_origin + launch_distance * cp.axis_dir, d = -cp.axis_dir);
    Trace t_axis_bwd = sys.trace(r_axis_bwd, far = config.physics.far);
    if (t_axis_bwd.segments().length > 1) {
        cp.V2 = t_axis_bwd.segments()[1].a; // Rear surface hit point
    } else {
        cp.V2 = sys.elements[sys.elements.length - 1].center;
    }

    // -------------------------------------------------------------------------
    // 2. FORWARD TRACE: Rear Focal Point F2 & Rear Principal Point H2
    // -------------------------------------------------------------------------
    Ray r_fwd = Ray(p = axis_origin - launch_distance * cp.axis_dir + paraxial_h * cp.axis_norm, d = cp.axis_dir);
    Trace t_fwd = sys.trace(r_fwd, far = config.physics.far);
    Ray exit_fwd = t_fwd.last();

    // Directly intersect the exit ray with the optical axis line:
    cp.F2 = lineIntersection(exit_fwd, axis_origin, cp.axis_dir);

    if (abs(cp.F2.x) != infinity) {
        // Principal plane = where the extended input and output rays meet,
        // so f = h / tan(u) (not h / sin(u)).
        real slope_exit = dot(exit_fwd.d, cp.axis_norm) / dot(exit_fwd.d, cp.axis_dir);
        cp.f2 = -paraxial_h / slope_exit;
        cp.H2 = cp.F2 - cp.f2 * cp.axis_dir;
    } else {
        // Afocal / Flat window system
        cp.f2 = infinity;
        cp.F2 = (infinity, infinity);
        cp.H2 = cp.V2;
    }

    // -------------------------------------------------------------------------
    // 3. BACKWARD TRACE: Front Focal Point F1 & Front Principal Point H1
    // -------------------------------------------------------------------------
    Ray r_bwd = Ray(p = axis_origin + launch_distance * cp.axis_dir + paraxial_h * cp.axis_norm, d = -cp.axis_dir);
    Trace t_bwd = sys.trace(r_bwd, far = config.physics.far);
    Ray exit_bwd = t_bwd.last();

    // Directly intersect the backward exit ray with the optical axis line:
    cp.F1 = lineIntersection(exit_bwd, axis_origin, -cp.axis_dir);

    if (abs(cp.F1.x) != infinity) {
        real slope_bwd = dot(exit_bwd.d, cp.axis_norm) / abs(dot(exit_bwd.d, cp.axis_dir));
        cp.f1 = -paraxial_h / slope_bwd;
        cp.H1 = cp.F1 + cp.f1 * cp.axis_dir;
    } else {
        cp.f1 = infinity;
        cp.F1 = (-infinity, -infinity);
        cp.H1 = cp.V1;
    }

    // -------------------------------------------------------------------------
    // 4. NODAL POINTS N1 & N2
    // -------------------------------------------------------------------------
    // In Gaussian optics: N1 = H1 + (1 - n1/n2) * f2 * axis_dir
    // When immersed in the same medium on both sides (n1 == n2): N1 == H1, N2 == H2
    // Object-space medium as resolved by the tracer (the launch ray itself
    // carries medium == null until traced).
    Ray launched = t_fwd.first();
    real n_obj = (launched != null && launched.medium != null) ? launched.medium.n() : 1.0;
    real n_img = (exit_fwd.medium != null) ? exit_fwd.medium.n() : 1.0;
    real nodal_shift = (abs(n_img) > config.math.eps) ? (1.0 - n_obj / n_img) * cp.f2 : 0.0;

    cp.N1 = cp.H1 + nodal_shift * cp.axis_dir;
    cp.N2 = cp.H2 + nodal_shift * cp.axis_dir;

    // -------------------------------------------------------------------------
    // 5. DISTANCES (BFL, FFL) & THICKNESS FLAG
    // -------------------------------------------------------------------------
    cp.BFL = dot(cp.F2 - cp.V2, cp.axis_dir);
    cp.FFL = dot(cp.V1 - cp.F1, cp.axis_dir);
    cp.is_thick = (length(cp.H2 - cp.H1) > config.geometry.eps);

    return cp;
}

// 1. One-Line 3-Ray Gaussian Generator (replaces ad-hoc ray construction)
Ray[] gaussianRays(OpticalSystem sys, pair obj_tip, real axis_angle = nan) {
    if (sys.elements.length == 0) return new Ray[];
    Element first_optic = sys.elements[0];
    pair axis = !isnan(axis_angle) ? dir(axis_angle) : unitOr(first_optic.axis, (1, 0));
    pair q = (-axis.y, axis.x);

    CardinalPoints cp;
    if (!isnan(first_optic.focal_length) && abs(first_optic.focal_length) > config.math.eps) {
        cp.H1 = first_optic.center;
        cp.F1 = first_optic.center - first_optic.focal_length * axis;
    } else {
        cp = cardinalPoints(system(first_optic, sys.ambient), first_optic.center, axis);
    }

    Ray r_par = Ray(obj_tip, axis);
    Ray r_ctr = Ray.between(obj_tip, cp.H1);

    // Focal ray with aperture-clipping safeguard:
    Ray r_foc;
    if (abs(cp.F1.x) != infinity) {
        r_foc = Ray.between(obj_tip, cp.F1);
    } else {
        r_foc = Ray(obj_tip, axis);
    }

    return new Ray[]{r_par, r_ctr, r_foc};
}

// 2. Rich ImageResult Struct (replaces raw path return)
struct ImageResult {
    path image_path;
    real magnification = 1.0;
    pair base = (nan, nan);
    pair tip = (nan, nan);
    bool is_real = true;
}

// Implicit cast: 'path p = image(...)' continues to work 100% identically!
path operator cast(ImageResult res) {
    return res.image_path;
}

// 3. Upgraded image() (returns ImageResult with .magnification, .base, .tip)
ImageResult image(OpticalSystem sys,
                  path object_path,
                  int samples = 20,
                  real delta_angle = config.physics.differential_angle,
                  real lambda = 0) {
    ImageResult res;
    int n_segs = length(object_path);
    if (sys.elements.length == 0 || n_segs == 0) return res;

    guide img_guide;
    real total_l = arclength(object_path);
    int n_points = max(samples, n_segs * 4);

    for (int i = 0; i <= n_points; ++i) {
        real s = i * total_l / (real)n_points;
        pair obj_pt = point(object_path, arctime(object_path, s));
        pair img_pt = imagePoint(sys, obj_pt, delta_angle, lambda);
        if (abs(img_pt.x) != infinity && abs(img_pt.y) != infinity) {
            img_guide = img_guide -- img_pt;
        }
    }
    if (cyclic(object_path)) img_guide = img_guide -- cycle;
    res.image_path = img_guide;

    if (length(res.image_path) > 0) {
        pair obj_base = point(object_path, 0);
        pair obj_tip = point(object_path, length(object_path));
        res.base = point(res.image_path, 0);
        res.tip = point(res.image_path, length(res.image_path));

        real h_obj = length(obj_tip - obj_base);
        real h_img = length(res.tip - res.base);
        res.magnification = (h_obj > config.math.eps) ? (h_img / h_obj) : 1.0;
        if (dot(res.tip - res.base, obj_tip - obj_base) < 0) res.magnification = -res.magnification;
    }
    return res;
}

// Renders the folded optical axis through any system with before/after bounds
void drawAxis(picture pic = currentpicture,
              OpticalSystem sys,
              real before = infinity,
              real after = infinity,
              pen p = gray(0.5) + dashed,
              real margin = 0) {
    if (sys.elements.length == 0) return;

    Element first_optic = sys.elements[0];
    pair axis0 = unitOr(first_optic.axis, (1, 0));

    // 1. Resolve front vertex V1
    real search_dist = 0.5 * config.physics.far;
    Ray r_find_v1 = Ray(p = first_optic.center - search_dist * axis0, d = axis0);
    Trace t_find = sys.trace(r_find_v1, far = config.physics.far);
    pair V1 = (t_find.segments().length > 1) ? t_find.segments()[1].a : first_optic.center;

    // 2. Trace optical axis from front vertex V1 through the entire system
    Ray axis_ray = Ray(V1, axis0);
    Trace t_axis = sys.trace(axis_ray, far = config.physics.far);
    Segment[] segs = t_axis.segments();
    if (segs.length == 0) return;

    int n_segs = segs.length;

    // -------------------------------------------------------------------------
    // A. ENTRANCE EXTENSION (Before front vertex V1)
    // -------------------------------------------------------------------------
    if (before == infinity || before == -infinity) {
        // Extend to viewport edge before the first element
        drawClippedRay(pic, V1, -axis0, margin, p, new ResolvedArrow, false, true);
    } else if (!isnan(before) && before > 0) {
        // Extend by exact user distance before V1
        draw(pic, (V1 - before * axis0) -- V1, p);
    }
    // If isnan(before): do nothing! Starts strictly at V1.

    // -------------------------------------------------------------------------
    // B. INTERNAL PROPAGATION (Between V1 and V2 through all lenses & mirrors)
    // -------------------------------------------------------------------------
    // Segments 0 to n_segs - 2 are internal to the system:
    for (int i = 0; i < n_segs - 1; ++i) {
        draw(pic, segs[i].traj, p);
    }

    // -------------------------------------------------------------------------
    // C. EXIT EXTENSION (After rear vertex V2)
    // -------------------------------------------------------------------------
    Segment last_seg = segs[n_segs - 1];
    pair V2 = last_seg.a; // Rear vertex of the final element
    pair d_exit = unitOr(last_seg.before.d, axis0);

    if (after == infinity) {
        // Extend to viewport edge after the last element
        drawClippedRay(pic, V2, d_exit, margin, p, new ResolvedArrow, false, true);
    } else if (!isnan(after) && after > 0) {
        // Extend by exact user distance after V2
        draw(pic, V2 -- (V2 + after * d_exit), p);
    }
    // If isnan(after): do nothing! Ends strictly at V2.
}

// 4. Compact draw(sys, object) (shrunk from 110 lines to 45 lines)

void draw(picture pic = currentpicture,
          OpticalSystem sys,
          path object,
          pen p_ray = black,
          pen p_virtual = p_ray + dashed,
          pen p_object = black,
          pen p_image = p_object,
          pen p_axis = gray(0.5) + dashed,
          bool intermediate = true,
          bool show_axis = true,
          bool scale_arrows = true,
          real arrow_size = 5pt,
          real before = nan,
          real after = infinity,
          real axis_before = infinity, // Default: extends to canvas edge
          real axis_after = infinity,  // Default: extends to canvas edge
          real far = config.physics.far) {
    if (sys.elements.length == 0 || length(object) == 0) return;

    pair obj_base = point(object, 0);
    pair obj_tip = point(object, length(object));
    real h_obj = length(obj_tip - obj_base);

    // Axis and Ray tracing
    pair first_axis = unitOr(sys.elements[0].axis, (1, 0));
    Trace t_axis = sys.trace(Ray(obj_base, first_axis), far = far);
    Segment[] axis_segs = t_axis.segments();

    Ray[] rays = gaussianRays(sys, obj_tip);
    Trace[] traces;
    for (Ray r : rays) traces.push(sys.trace(r, far = far));

    // Resolve true optical stages (skipping internal glass propagation)
    int num_segs = traces[0].segments().length;
    pair[] stageTips, stageBases;
    bool[] stageIsReal;

    for (int s = 1; s < num_segs; ++s) {
        Segment seg0 = traces[0].segments()[s];
        Segment seg1 = traces[1].segments()[s];
        if (seg0.before.medium != sys.ambient && s < num_segs - 1) continue; // Skip glass

        pair tip = lineIntersection(seg0.a, seg0.before.d, seg1.a, seg1.before.d);
        if (abs(tip.x) == infinity || abs(tip.y) == infinity) continue;

        bool is_real = (dot(tip - seg0.a, seg0.before.d) > config.geometry.eps);
        int ax_idx = min(s, axis_segs.length - 1);
        pair ax_p1 = axis_segs[ax_idx].a;
        pair ax_d = unitOr(axis_segs[ax_idx].b - ax_p1, axis_segs[ax_idx].before.d);
        pair base = ax_p1 + dot(tip - ax_p1, ax_d) * ax_d;

        stageTips.push(tip);
        stageBases.push(base);
        stageIsReal.push(is_real);
    }

    // Render Axis, Rays, Extensions & Arrows
    // Render Folded Optical Axis (with clean start and end bounds)
    if (show_axis) {
        drawAxis(pic, sys, before = axis_before, after = axis_after, p = p_axis);
    }

    for (Trace t : traces)
        draw(pic,
             t,
             before = before,
             after = after,
             p = p_ray,
             arrows = rayArrows(first = true, middle = false, last = true));

    for (int k = 0; k < stageTips.length; ++k) {
        if (!stageIsReal[k]) {
            for (Trace t : traces)
                if (k + 1 < t.segments().length) draw(pic, t.segments()[k + 1].a -- stageTips[k], p_virtual);
        }
    }

    draw(pic, obj_base -- obj_tip, p_object, EndArrow(DefaultHead, size = arrow_size, filltype = Draw));

    for (int k = 0; k < stageTips.length; ++k) {
        if (!intermediate && k < stageTips.length - 1) continue;
        real h_img = length(stageTips[k] - stageBases[k]);
        if (h_img < config.math.eps) continue;
        real eff_arrow = scale_arrows ? (h_img / h_obj) * arrow_size : arrow_size;
        pen pen_arr = (k == stageTips.length - 1) ? p_image : (0.7 * p_image + 0.3 * white);
        if (!stageIsReal[k]) pen_arr = pen_arr + dashed;
        draw(pic,
             stageBases[k] -- stageTips[k],
             pen_arr,
             EndArrow(DefaultHead, size = eff_arrow, filltype = Draw));
    }
    draw(pic, sys);
}

void draw(picture pic = currentpicture,
          CardinalPoints cp,
          real plane_height = 35.0,
          pen p_axis = gray(0.5) + dashed + 0.6bp,
          pen p_planes = rgb(0.1, 0.4, 0.8) + dashed + 0.8bp,
          pen p_points = heavyred) {
    if (abs(cp.f2) == infinity) return;
    // 1. Draw Optical Axis
    real x_left = min(cp.F1.x, cp.V1.x) - 30.0;
    real x_right = max(cp.F2.x, cp.V2.x) + 30.0;
    pair ax_a = cp.axis_origin + x_left * cp.axis_dir;
    pair ax_b = cp.axis_origin + x_right * cp.axis_dir;
    draw(pic, ax_a -- ax_b, p_axis);
    // 2. Draw Principal Planes (H1 & H2)
    pair H1_top = cp.H1 + (plane_height / 2.0) * cp.axis_norm;
    pair H1_bot = cp.H1 - (plane_height / 2.0) * cp.axis_norm;
    pair H2_top = cp.H2 + (plane_height / 2.0) * cp.axis_norm;
    pair H2_bot = cp.H2 - (plane_height / 2.0) * cp.axis_norm;
    draw(pic, H1_bot -- H1_top, p_planes);
    draw(pic, H2_bot -- H2_top, p_planes);
    label(pic, "\small $\mathbf{H}_1$", H1_top, N, rgb(0.1, 0.4, 0.8));
    label(pic, "\small $\mathbf{H}_2$", H2_top, N, rgb(0.1, 0.4, 0.8));
    // 3. Mark Cardinal Points on Axis
    dot(pic, cp.F1, p_points + 3pt);
    dot(pic, cp.F2, p_points + 3pt);
    dot(pic, cp.H1, blue + 3pt);
    dot(pic, cp.H2, blue + 3pt);
    label(pic, "\small $\mathbf{F}_1$", cp.F1, 2S, p_points);
    label(pic, "\small $\mathbf{F}_2$", cp.F2, 2S, p_points);
    label(pic, "\small $H_1$", cp.H1, 2S, blue);
    label(pic, "\small $H_2$", cp.H2, 2S, blue);
    // 4. Mark Vertices
    dot(pic, cp.V1, black + 2pt);
    dot(pic, cp.V2, black + 2pt);
    label(pic, "\tiny $V_1$", cp.V1, NW, gray(0.3));
    label(pic, "\tiny $V_2$", cp.V2, NE, gray(0.3));
    // 5. Dimension Line for Effective Focal Length (EFL)
    real dim_y = (plane_height / 2.0) + 12.0;
    pair dim_a = cp.H2 + dim_y * cp.axis_norm;
    pair dim_b = cp.F2 + dim_y * cp.axis_norm;
    draw(pic, dim_a -- dim_b, p_points + 0.6bp, Arrows(HookHead, 2.5pt));
    label(pic, format("\small $f = %#.2f$", cp.f2), midpoint(dim_a -- dim_b), N, p_points);
}

// ----------------------------------------------------------------------------
// Section 9: Visual Insets, Legend Panels & Declarative DSL
// ----------------------------------------------------------------------------
real resolveSegmentOffset(Segment seg,
                          int idx,
                          int total_segments,
                          real relative = nan,
                          real absolute = nan,
                          string caller = "resolveSegmentOffset") {
    real L = arclength(seg.traj);
    // Single Unbounded Ray (No optical collisions) -> Allow standard relative or absolute
    if (total_segments == 1) {
        if (!isnan(relative)) return clamp(relative, 0.0, 1.0) * L;
        if (!isnan(absolute)) return clamp(abs(absolute) < config.math.eps ? 0.5 * L : abs(absolute), 0.0, L);
        return 0.5 * L;
    }
    // Open Entrance Ray (First Segment) -> Anchored at first surface hit
    if (idx == 0 && total_segments > 1) {
        if (isnan(absolute)) {
            abort("[" + caller +
                  "] Error: Please specify 'absolute' distance from the first surface on the entrance ray.");
        }
        real s;
        if (abs(absolute) < config.math.eps) {
            // special case: absolute = 0 means exact midpoint (relative = 0.5)
            s = 0.5 * L;
        } else {
            s = L - abs(absolute);
        }
        return clamp(s, 0.0, L);
    }
    // Open Exit Ray (Last Segment) -> Anchored at exit surface
    if (idx == total_segments - 1) {
        if (isnan(absolute)) {
            abort("[" + caller + "] Error: Please specify 'absolute' distance from the exit surface on the exit ray.");
        }
        real s;
        if (abs(absolute) < config.math.eps) {
            // special case: absolute = 0 means exact midpoint (relative = 0.5)
            s = 0.5 * L;
        } else {
            s = abs(absolute);
        }
        return clamp(s, 0.0, L);
    }
    // Intermediate Bounded Segments (Between optics) -> Both ends known
    if (!isnan(relative)) {
        return clamp(relative, 0.0, 1.0) * L;
    }
    if (!isnan(absolute)) {
        // real s = (absolute >= 0) ? absolute : L + absolute;
        real s;
        if (absolute > 0) {
            // Positive absolute: distance from beginning of segment (along ray)
            s = absolute;
        } else if (absolute < 0) {
            // Negative absolute: distance from end of segment (opposite ray)
            s = L + absolute;
        } else {
            // special case: absolute = 0 means exact midpoint (relative = 0.5)
            s = 0.5 * L;
        }
        return clamp(s, 0.0, L);
    }
    // Default for bounded segments: exact geometric midpoint (relative = 0.5)
    return 0.5 * L;
}

// ============================================================================
// Clean General-Purpose Inset Component
// ============================================================================
struct Inset {
    picture pic;              // Content canvas
    pair pos = (nan, nan);    // Local badge/inset position on diagram
    pair target = (nan, nan); // Physical target point on ray
    string title = "";        // Optional title
    string tag = "";          // Badge tag (e.g. "A", "1")
    real size = 20.0;         // Content scale/radius
    real title_content_gap = 5.0;
    real box_margin = 1.0; // User-configurable padding inside the inset box
    // Color Palette
    pen border = rgb(0.0, 0.2, 0.6); // Box outline, connector line, & target dot color
    pen fill = white;                // Inset box background fill
    pen connector = nullpen;         // Connector line pen (defaults to border + dashed)
    pen pen_badge = nullpen;         // Badge fill color (defaults to border)
    pen badge_text = white;          // Badge text label color
    void add(Label L) { label(this.pic, L); }
    void addGraphic(string filename) { label(this.pic, graphic(filename)); }
    // Helpers to resolve fallback pens
    pen getConnectorPen() { return (this.connector != nullpen) ? this.connector : (this.border + dashed); }
    pen getBadgePen() { return (this.pen_badge != nullpen) ? this.pen_badge : this.border; }
    // Standalone Drawing
    void draw(picture dest = currentpicture, pair position = (nan, nan), pen p = nullpen) {
        pair draw_pos = (!isnan(position.x)) ? position : this.pos;
        dest.add(new void(frame f, transform t) {
            pair pos_pt = t * draw_pos;
            frame content = this.pic.fit();
            frame box_frame;
            // Wrap content in box and shift to target canvas position
            filltype ft = (this.fill != nullpen) ? FillDraw(this.fill, this.border) : Fill(this.border);
            box(box_frame, content, this.box_margin, this.box_margin, this.border, ft, above = false);
            add(box_frame, content);
            box_frame = shift(pos_pt) * box_frame;
            // Draw connector line and dot at target point
            if (!isnan(this.target.x) && !isnan(this.target.y)) {
                pair target_pt = t * this.target;
                if (target_pt != pos_pt) {
                    pen pen_line = (p != nullpen) ? p : this.getConnectorPen();
                    pen pen_dot = (this.border != nullpen) ? this.border : pen_line;
                    pair intersection_pt = point(box_frame, unitOr(target_pt - pos_pt));
                    draw(f, target_pt -- intersection_pt, pen_line);
                    dot(f, target_pt, pen_dot);
                }
            }
            // Stamp boxed inset on canvas
            add(f, box_frame);
        });
    }
}

// Shared inset header: every inset uses the same tag pill and title baseline.
private void drawInsetHeader(picture pic, real y, string title, string tag, pen border) {
    frame header;
    real title_x = 0.0;
    if (tag != "") {
        frame tag_frame;
        label(tag_frame, Label("\bfseries{\tiny " + tag + "}", white));
        roundbox(tag_frame, tag_frame, 0.0, 0.0, border, Fill(border), above = false);
        add(header, tag_frame, N);
        title_x = max(tag_frame).x;
    }
    if (title != "") {
        // Keep the title vertically centered with or without a badge.  N
        // alignment makes a badge-less title extend downward from its anchor,
        // creating an unnecessarily large title/content gap.
        pair title_anchor = (tag != "") ? (title_x, 0) : (0, 0);
        label(header, Label("\bfseries{\tiny\strut " + title + "}", border), title_anchor, (tag != "") ? E : W);
    }
    if (title != "" || tag != "") {
        pair header_center = ((min(header).x + max(header).x) / 2.0, 0);
        add(pic, shift((-header_center.x, y)) * header);
    }
}

// Polarization Inset Factory (From raw Polarization state)
Inset PolarizationInset(Polarization pol,
                        pair pos = (nan, nan),
                        pair target = (nan, nan),
                        real size = 20.0,
                        real ref_amplitude = 1.0,
                        string title = "",
                        string tag = "",
                        pen p = currentpen,
                        pen border = rgb(0.0, 0.2, 0.6),
                        pen fill = white,
                        pen badge_text = white) {
    Inset ins;
    ins.pos = pos;
    ins.target = target;
    ins.size = size;
    ins.title = title;
    ins.tag = tag;
    ins.border = border;
    ins.fill = fill;
    ins.badge_text = badge_text;
    real r = size;

    // Reference Axes
    draw(ins.pic, (-r, 0) -- (r, 0), gray(0.8), Arrow(HookHead, 2pt));
    draw(ins.pic, (0, -r) -- (0, r), gray(0.8), Arrow(HookHead, 2pt));
    label(ins.pic, "\tiny $S$", (r, 0), 2S + W, gray(0.4));
    label(ins.pic, "\tiny $P$", (0, r), S + 2W, gray(0.4));

    // Render tag and title together so their baseline is independent of the
    // surrounding panel layout.
    real header_y = r + ins.title_content_gap;
    drawInsetHeader(ins.pic, header_y, title, tag, border);

    real local_i = pol.intensity();
    real orbit_r = 0.75 * r;

    // Guaranteed auto-scale: if ref_amplitude is nan, scale by sqrt(local_i)
    real scale_norm =
        (!isnan(ref_amplitude) && ref_amplitude > 0) ? ref_amplitude : sqrt(max(config.math.eps, local_i));
    real display_factor = orbit_r / scale_norm;

    real dop = pol.dop();
    real I_pol = dop * local_i;
    real min_arrow_size = config.rendering.polarization_inset_min_arrow_fraction * orbit_r;

    // =========================================================================
    // 1. STARBURST SPOKES (Unpolarized & Partially Polarized)
    // =========================================================================
    // Only drawn when there is a significant unpolarized fraction (DoP < 95%)
    if (dop < 0.95 && local_i > config.math.eps) {
        real unpol_r = display_factor * sqrt(local_i);
        int spokes = config.rendering.polarization_inset_unpolarized_spokes;
        if (spokes <= 0) spokes = 8;

        real psi = pol.psi();
        real min_spoke_alpha = 0.05; // Faint floor for orthogonal/unpolarized directions
        real max_spoke_alpha = 0.95; // Bold ceiling for dominant polarization direction

        for (int k = 0; k < spokes; ++k) {
            real ang = radians(k * 180.0 / spokes);
            pair dir = (cos(ang), sin(ang)) * unpol_r * 0.85;

            // Alignment factor along the spoke angle: 1.0 = aligned with psi, 0.0 = perpendicular
            real cos_diff = cos(ang - psi);
            real align_weight = cos_diff * cos_diff;

            // Physical opacity modulation:
            // - Pure unpolarized (dop = 0): all spokes have equal alpha (0.70)
            // - Partially polarized (dop > 0): aligned spoke stays bold, orthogonal fades to min_spoke_alpha
            real spoke_factor = (1.0 - dop) * 0.75 + dop * align_weight;
            real spoke_alpha = min_spoke_alpha + (max_spoke_alpha - min_spoke_alpha) * spoke_factor;
            spoke_alpha = clamp(spoke_alpha, min_spoke_alpha, 1.0);

            // Give dominant spokes slightly more line weight for clarity
            // real spoke_width = (dop > 0.2 && align_weight > 0.7) ? 0.6bp : 0.4bp;
            // pen pen_spoke = p + opacity(spoke_alpha);
            real s = 1.0 - spoke_alpha;
            pen pen_spoke = p + gray(s);

            if (length(dir) > min_arrow_size) {
                draw(ins.pic, (0, 0) -- - dir, pen_spoke, Arrow(HookHead, 2.0pt), BeginPenMargin);
                draw(ins.pic, (0, 0) -- dir, pen_spoke, Arrow(HookHead, 2.0pt), BeginPenMargin);
            } else
                draw(ins.pic, -dir -- dir, pen_spoke);
        }
    }

    // =========================================================================
    // 2. POLARIZED SHAPE (Linear Line or Ellipse)
    // =========================================================================
    // Condition to avoid double-drawing:
    // - Pure polarized states (dop >= 0.95): always draw the clean line/ellipse.
    // - Partially polarized states (dop < 0.95): ONLY draw if it has non-zero ellipticity
    //   (so the viewer can see the ellipse and rotation arrow). If linear, the modulated
    //   spokes already represent the state without a duplicate line.
    bool has_ellipticity = (abs(pol.chi()) >= radians(1.5));
    bool draw_polarized_shape = (dop >= 0.95) || (has_ellipticity && I_pol > config.math.eps);

    if (draw_polarized_shape && I_pol > config.math.eps) {
        real psi = pol.psi();
        real chi = pol.chi();
        real a = display_factor * sqrt(I_pol) * cos(chi);
        real b = display_factor * sqrt(I_pol) * abs(sin(chi));

        // Linear polarization (only runs for pure polarized states):
        if (!has_ellipticity) {
            pair tip = a * (cos(psi), sin(psi));
            if (a > min_arrow_size)
                draw(ins.pic, -tip -- tip, p, Arrows(HookHead, 2.5pt));
            else
                draw(ins.pic, -tip -- tip, p);
        } else {
            // Elliptical / Circular:
            guide ellipse_path;
            int steps = config.rendering.polarization_inset_ellipse_steps;
            for (int k = 0; k <= steps; ++k) {
                real tau = k * 2.0 * pi / steps;
                real x_loc = a * cos(tau);
                real y_loc = (pol.S3 >= 0 ? 1.0 : -1.0) * b * sin(tau);
                real S_t = x_loc * cos(psi) - y_loc * sin(psi);
                real P_t = x_loc * sin(psi) + y_loc * cos(psi);
                ellipse_path = ellipse_path .. (S_t, P_t);
            }
            if (a > min_arrow_size)
                draw(ins.pic, ellipse_path, p, Arrow(HookHead(), 2.5pt, Relative(0.85)));
            else
                draw(ins.pic, ellipse_path, p);
        }
    }

    // =========================================================================
    // 3. READOUT LABEL
    // =========================================================================
    object[] lines;
    lines.push(object(Label(format("\tiny $I = %#.3f$", local_i))));
    if (dop >= config.rendering.polarization_inset_unpolarized_threshold && dop < 0.999) {
        lines.push(object(Label(format("\tiny $\mathrm{DoP} = %.0f\%$", dop * 100.0))));
    } else if (dop < config.rendering.polarization_inset_unpolarized_threshold) {
        lines.push(object(Label("\tiny (Unpolarized)")));
    }

    frame f_readout = pack(align = S, ... lines);
    add(ins.pic, f_readout, (0, -r - 2.0), S);
    return ins;
}

// Trace Overload: Automatically extracts target point from the ray segment
Inset PolarizationInset(Trace tr,
                        int segment_index = 0,
                        real relative = nan,
                        real absolute = nan,
                        pair pos = (nan, nan),
                        real size = 20.0,
                        real ref_amplitude = nan,
                        string title = "",
                        string tag = "",
                        pen p = currentpen,
                        pen border = rgb(0.0, 0.2, 0.6),
                        pen fill = white,
                        pen badge_text = white) {
    Segment[] segs = tr.segments();
    if (segs.length == 0) return new Inset;
    int idx = (segment_index < 0) ? segment_index + segs.length : segment_index;
    if (idx < 0 || idx >= segs.length) return new Inset;
    Segment seg = segs[idx];
    if (seg.before == null || !seg.before.alive) return new Inset;
    real L = arclength(seg.traj);
    real s = resolveSegmentOffset(seg,
                                  idx,
                                  segs.length,
                                  relative = relative,
                                  absolute = absolute,
                                  caller = "PolarizationInset");
    pair target_point = (L > 0) ? point(seg.traj, arctime(seg.traj, s)) : seg.a;
    Polarization local_pol = resolvePolarizationAt(seg.before, s);
    real ref_amp = ref_amplitude;
    if (isnan(ref_amp) && tr.first() != null) {
        real start_i = tr.first().intensity();
        ref_amp = (start_i > config.math.eps) ? sqrt(start_i) : 1.0;
    }
    return PolarizationInset(local_pol,
                             pos = pos,
                             target = target_point,
                             size = size,
                             ref_amplitude = ref_amp,
                             title = title,
                             tag = tag,
                             p = p,
                             border = border,
                             fill = fill,
                             badge_text = badge_text);
}

// ============================================================================
// InterferenceInset Factory (Type Inset, using interfere/interfereState)
// ============================================================================
Inset InterferenceInset(Ray[] sources,
                        pair screen_center = (0, 0),
                        pair screen_axis = (0, 1),
                        real span = nan,
                        real vis_lambda = nan,
                        int periods = 5,
                        int samples = 400,
                        pair pos = (nan, nan),
                        pair target = (nan, nan),
                        real size = 28.0,
                        string title = "Screen Interference",
                        string tag = "",
                        pen pen_laser = Red(1.0),
                        pen border = rgb(0.0, 0.2, 0.6),
                        pen fill = white,
                        pen badge_text = white,
                        real rotation = 0.0) {
    if (samples < 1) abort("InterferenceInset: samples must be positive.");
    if (periods < 1) abort("InterferenceInset: periods must be positive.");
    if (sources.length < 2) abort("InterferenceInset: at least two sources are required.");
    if (!isnan(vis_lambda) && vis_lambda <= 0) abort("InterferenceInset: vis_lambda must be positive or nan.");

    Inset ins;
    picture content;
    picture title_content;
    ins.pos = pos;
    ins.target = target;
    ins.size = size;
    ins.title = title;
    ins.tag = tag;
    ins.border = border;
    ins.fill = fill;
    ins.badge_text = badge_text;

    real r = size;
    real w = r * 1.5;
    real h = r;

    // ========================================================================
    // Header Title & Tag Badge (Sharing the exact same typographical baseline)
    // ========================================================================
    real y_header = h * 0.85 + ins.title_content_gap;
    drawInsetHeader(title_content, y_header, title, tag, border);

    // ========================================================================
    // 1. Effective Wavelength & Universal N-Slit Auto-Span
    // ========================================================================
    pair probe_center = (!isnan(target.x) && !isnan(target.y)) ? target : screen_center;
    pair u_screen = unitOr(screen_axis, (0, 1));
    pair n_screen = perp(u_screen);

    real ref_lambda_nm = (sources.length > 0 && sources[0] != null) ? sources[0].lambda : config.physics.lambda;
    bool has_visual_lambda = !isnan(vis_lambda) && vis_lambda > 0;
    real lambda_geom_units = resolveEffectiveLambda(vis_lambda, ref_lambda_nm);
    // Keep the physical source wavelength in the public nm convention. Only
    // an explicitly supplied rendering wavelength needs conversion from
    // geometry units before material queries are made.
    real lambda_eff_nm =
        has_visual_lambda ? lambda_geom_units * metersToUm(config.physics.meters_per_unit) * um_to_nm : ref_lambda_nm;

    real actual_span = span;
    if (isnan(actual_span) || actual_span <= 0) {
        real min_y = infinity, max_y = -infinity;
        real avg_D = 0.0;
        int count = 0;

        for (Ray src : sources) {
            if (src == null) continue;
            real y_proj = dot(src.p - probe_center, u_screen);
            real x_proj = abs(dot(src.p - probe_center, n_screen));
            if (y_proj < min_y) min_y = y_proj;
            if (y_proj > max_y) max_y = y_proj;
            avg_D += x_proj;
            ++count;
        }

        // Universal pitch for N >= 2 collinear sources: d = Delta_y / (N - 1)
        real d_char = (count > 1 && max_y > min_y) ? ((max_y - min_y) / (count - 1)) : 10.0;
        real D_char = (count > 0) ? (avg_D / count) : 100.0;

        // Exact fringe period: Delta_y = lambda * D / d
        real delta_y = (d_char > config.math.eps) ? (lambda_geom_units * D_char / d_char) : 10.0;
        actual_span = periods * delta_y;
    }

    // ========================================================================
    // 2. Physical Wave Interference Sampling
    // ========================================================================
    real max_i = 0.0;
    real[] intensities = new real[samples + 1];

    for (int i = 0; i <= samples; ++i) {
        real t = -1.0 + 2.0 * i / samples; // [-1, 1]
        pair P_y = probe_center + (t * actual_span / 2.0) * u_screen;

        Ray[] raysAtP;
        for (Ray src : sources) {
            if (src == null) continue;
            pair d_k = unitOr(P_y - src.p);
            real dist = length(P_y - src.p);
            real n_med = rayIndex(src.medium, src.p, d_k, d_k, lambda_eff_nm, src.pol);

            Ray r_k = Ray(src.p, d_k, medium = src.medium, lambda = lambda_eff_nm, pol = src.pol, phase = src.phase);
            r_k.coherence_id = src.coherence_id;
            r_k.opl = src.opl + geometryToUm(dist) * n_med;
            r_k.p = P_y;

            raysAtP.push(r_k);
        }

        Polarization pol_P = interfereState(raysAtP);
        real I_P = pol_P.intensity();

        intensities[i] = I_P;
        if (I_P > max_i) max_i = I_P;
    }

    if (max_i <= 0) max_i = 1.0;

    // ========================================================================
    // 3. Render Inset: Horizontal vs Vertical Layout
    // ========================================================================
    // --------------------------------------------------------------------
    // Standard Position (I(y) on top, Fringes directly UNDER)
    // --------------------------------------------------------------------
    real plot_half_w = w * 0.75;
    real plot_h = h * 0.70;
    real y_base = h * -0.15; // baseline between graph and fringes

    // 3A. Draw 1D Intensity Plot Curve I(y) on Top
    guide curve;
    for (int i = 0; i <= samples; ++i) {
        real t = -1.0 + 2.0 * i / samples;
        real norm_i = clamp(intensities[i] / max_i, 0.0, 1.0);
        pair pt = (t * plot_half_w, y_base + norm_i * plot_h);
        if (i == 0)
            curve = pt;
        else
            curve = curve -- pt;
    }

    // Fill area under curve and stroke line
    fill(content,
         (-plot_half_w, y_base) -- curve -- (plot_half_w, y_base) -- cycle,
                                                               interp(white, pen_laser, 0.15));
    draw(content, curve, pen_laser);

    // Axes: Horizontal position axis (y) & Vertical Intensity axis (I)
    draw(content,
         (-plot_half_w - 2, y_base) -- (plot_half_w + 6, y_base), gray(0.5), Arrow(HookHead, size = 3pt)); // y-axis
    draw(content,
         (-plot_half_w,
          y_base) -- (-plot_half_w, y_base + plot_h + 4), gray(0.5), Arrow(HookHead, size = 3pt)); // I-axis

    // label(content, "\tiny $y$", (plot_half_w + 6, y_base), E, gray(0.4));
    // label(content, "\tiny $I$", (-plot_half_w, y_base + plot_h + 4), N, gray(0.4));

    // Center tick (y = 0)
    draw(content, (0, y_base - 1.2) -- (0, y_base + 1.2), gray(0.5));
    label(content, "\tiny $0$", (0, y_base + 1.5), N, gray(0.4));

    // 3B. Draw 2D Fringe Camera Strip directly UNDER the graph
    real strip_top_y = y_base - 4.0;
    real strip_h = h * 0.45;
    real strip_bot_y = strip_top_y - strip_h;

    fill(content, box((-plot_half_w, strip_bot_y), (plot_half_w, strip_top_y)), black);

    for (int i = 0; i < samples; ++i) {
        real t1 = -1.0 + 2.0 * i / samples;
        real t2 = -1.0 + 2.0 * (i + 1) / samples;
        real x1 = t1 * plot_half_w;
        real x2 = t2 * plot_half_w;
        real norm_i = clamp(0.5 * (intensities[i] + intensities[i + 1]) / max_i, 0.0, 1.0);

        pen pen_slice = interp(black, pen_laser, norm_i);
        fill(content, box((x1, strip_bot_y), (x2, strip_top_y)), pen_slice);
    }

    draw(content, box((-plot_half_w, strip_bot_y), (plot_half_w, strip_top_y)), gray(0.4) + 0.6pt);
    // label(content, "\tiny Fringes", (0, strip_bot_y - 1.5), S, gray(0.3));

    // Dotted vertical projection line aligning central maximum down to fringe strip
    draw(content, (0, y_base + plot_h) -- (0, strip_bot_y), gray(0.5) + dotted);

    // Keep the header clear when the plot is rotated into a vertical layout.
    real rotated_content_shift = (abs(sin(radians(rotation))) > 0.5) ? -0.55 * h : 0.0;
    ins.pic.add(shift((0, rotated_content_shift)) * rotate(rotation) * content);
    ins.pic.add(title_content);

    return ins;
}

// Trace[] overload for InterferenceInset
Inset InterferenceInset(Trace[] traces,
                        pair screen_center = (0, 0),
                        pair screen_axis = (0, 1),
                        real span = nan,
                        real vis_lambda = nan,
                        int periods = 5,
                        int samples = 400,
                        pair pos = (nan, nan),
                        pair target = (nan, nan),
                        real size = 28.0,
                        string title = "Screen Interference",
                        string tag = "",
                        pen pen_laser = Red(1.0),
                        pen border = rgb(0.0, 0.2, 0.6),
                        pen fill = white,
                        pen badge_text = white,
                        real rotation = 0.0) {
    Ray[] sources;
    for (Trace tr : traces) {
        if (tr != null && tr.last() != null) sources.push(tr.last());
    }
    return InterferenceInset(sources,
                             screen_center,
                             screen_axis,
                             span,
                             vis_lambda,
                             periods,
                             samples,
                             pos,
                             target,
                             size,
                             title,
                             tag,
                             pen_laser,
                             border,
                             fill,
                             badge_text,
                             rotation);
}

// ============================================================================
// Clean Auto-Grid Inset Panel
// ============================================================================
struct InsetPanel {
    Inset[] insets;
    int perrow = 1;     // Grid columns (1 = vertical list, 2 = 2-col table)
    real margin = 2.0;  // Panel padding in pt
    real spacing = 8.0; // Spacing between cells in pt
    // Panel Global Styling
    pen border = nullpen;    // Panel outline (nullpen = defaults to insets[0].border)
    pen fill = nullpen;      // Panel background fill (nullpen = defaults to insets[0].fill)
    pen badge_text = white;  // Default badge text color
    bool draw_badges = true; // Draw badges and connectors on diagram
    InsetPanel add(Inset ins, string tag = "") {
        if (tag != "")
            ins.tag = tag;
        else if (ins.tag == "")
            ins.tag = string(this.insets.length + 1);
        this.insets.push(ins);
        return this;
    }
    // Direct Trace overload
    InsetPanel add(Trace tr,
                   int segment_index = 0,
                   real relative = nan,
                   real absolute = nan,
                   pair pos = (nan, nan),
                   real size = 20.0,
                   real ref_amplitude = 1.0,
                   string tag = "",
                   string title = "",
                   pen p = currentpen,
                   pen border = rgb(0.0, 0.2, 0.6),
                   pen fill = white,
                   pen badge_text = white) {
        Inset ins = PolarizationInset(tr,
                                      segment_index,
                                      relative,
                                      absolute,
                                      pos,
                                      size,
                                      ref_amplitude,
                                      title,
                                      tag,
                                      p,
                                      border,
                                      fill,
                                      badge_text);
        return this.add(ins, tag = tag);
    }
    // Draws the complete panel and diagram badges
    void draw(picture dest = currentpicture,
              pair position = (nan, nan),
              pair rel = (nan, nan),
              pair align = NW,
              pen p = nullpen) {
        if (this.insets.length == 0) return;
        dest.add(new void(frame f, transform t, transform T, pair m, pair M) {
            int n = this.insets.length;
            int cols = min(n, max(1, this.perrow));
            int rows = (int)ceil(n / (real)cols);
            real[] col_w = array(cols, 0.0);
            real[] row_h = array(rows, 0.0);
            frame[] cells = new frame[n];
            // [Section 1] Build Individual Legend Cells (Badge + Inset Content)
            for (int i = 0; i < n; ++i) {
                frame cf;
                Inset ins = this.insets[i];
                // Inset content next to badge
                frame content = ins.pic.fit();
                add(cf, shift(-min(content).x, min(content).y) * content);
                cells[i] = cf;
                int r = floor(i / cols), c = i % cols;
                pair sz = max(cf) - min(cf);
                if (sz.x > col_w[c]) col_w[c] = sz.x;
                if (sz.y > row_h[r]) row_h[r] = sz.y;
            }
            // [Section 2] Assemble the Multi-Column Grid
            frame grid_frame;
            real y_curr = 0;
            for (int r = 0; r < rows; ++r) {
                real x_curr = 0;
                for (int c = 0; c < cols; ++c) {
                    int idx = r * cols + c;
                    if (idx < n) {
                        pair offset = (x_curr - min(cells[idx]).x, y_curr - max(cells[idx]).y);
                        add(grid_frame, shift(offset) * cells[idx]);
                    }
                    x_curr += col_w[c] + this.spacing;
                }
                y_curr -= (row_h[r] + this.spacing);
            }
            // [Section 3] Wrap Grid in Panel Background Container Box
            frame panel_box;
            // Resolve panel outline & fill defaults from the first inset if unset
            pen p_border = (this.border != nullpen) ? this.border : (n > 0 ? this.insets[0].border : gray(0.5));
            pen p_fill = (this.fill != nullpen) ? this.fill : (n > 0 ? this.insets[0].fill : white);
            filltype ft = (p_fill != nullpen) ? FillDraw(p_fill, p_border) : Fill(p_border);
            box(panel_box, grid_frame, this.margin, p_border, ft, above = false);
            add(panel_box, grid_frame);
            // [Section 4] Resolve Canvas / User Position & 1-Line Alignment
            pair pos_pt = (!isnan(rel.x) && !isnan(rel.y)) ? (m.x + rel.x * (M.x - m.x), m.y + rel.y * (M.y - m.y))
                                                           : t * position;
            pair shift_vec = pos_pt - point(panel_box, align);
            // [Section 5] Draw Badges, Target Dots & Bridges on the Diagram
            if (this.draw_badges) {
                for (Inset ins : this.insets) {
                    if (!isnan(ins.target.x) && !isnan(ins.target.y)) {
                        pair target_pt = t * ins.target;
                        pair badge_pt = (!isnan(ins.pos.x) && ins.pos != ins.target) ? t * ins.pos : target_pt;
                        pen pen_badge = ins.getBadgePen();
                        pen b_txt = (ins.badge_text != nullpen) ? ins.badge_text : this.badge_text;
                        // Draw on-diagram badge box
                        frame b_frame;
                        label(b_frame, Label("\tiny\bfseries " + ins.tag, b_txt));
                        roundbox(b_frame, b_frame, 0.0, 0.0, pen_badge, Fill(pen_badge), above = false);
                        b_frame = shift(badge_pt) * b_frame;
                        // If badge is offset from target, draw connector line and target dot
                        if (badge_pt != target_pt) {
                            pen pen_line = (p != nullpen) ? p : ins.getConnectorPen();
                            pen pen_dot = (ins.border != nullpen) ? ins.border : pen_line;
                            pair intersection_pt = point(b_frame, unitOr(target_pt - badge_pt));
                            draw(f, target_pt -- intersection_pt, pen_line);
                            dot(f, target_pt, pen_dot);
                        }
                        add(f, b_frame);
                    }
                }
            }
            // [Section 6] Stamp the Legend Panel onto the Canvas
            add(f, shift(shift_vec) * panel_box);
        });
    }
}

// ============================================================================
// Implicit String-to-Medium Casting (Case-Insensitive & Type-Safe)
// ============================================================================
Medium operator cast(string name) {
    string n = UpperCase(name);
    if (n == "VACUUM" || n == "VAC") return Materials.VACUUM;
    if (n == "AIR") return Materials.AIR;
    if (n == "WATER") return Materials.WATER;
    if (n == "BK7") return Materials.BK7;
    if (n == "FLINT") return Materials.FLINT;
    if (n == "SILICA") return Materials.SILICA;
    abort("Unknown optical material string: '" + name + "'. Use medium(name, n) for custom media.");
    return Materials.VACUUM;
}

// ============================================================================
// SYSTEM ALGEBRA & DECLARATIVE DSL OPERATORS
// ============================================================================

// 1. Ray Tracing via Multiplication Operator (*)
// Usage: Trace t = sys * r;  or  Trace t = sys * rays;
Trace operator*(OpticalSystem sys, Ray r) { return sys.trace(r); }
Trace operator*(OpticalSystem sys, Ray[] rays) { return sys.trace(rays); }

// 2. System Assembly via Addition Operator (+)

// OpticalSystem + OpticalSystem -> Combines elements of both systems
OpticalSystem operator+(OpticalSystem sys1, OpticalSystem sys2) {
    OpticalSystem out = copy(sys1);
    out.add(sys2.elements);
    return out;
}

// OpticalSystem + Element -> Appends a single element
OpticalSystem operator+(OpticalSystem sys, Element e) {
    OpticalSystem out = copy(sys);
    out.add(e);
    return out;
}

// OpticalSystem + Element[] -> Appends an array of elements
OpticalSystem operator+(OpticalSystem sys, Element[] es) {
    OpticalSystem out = copy(sys);
    out.add(es);
    return out;
}

// Element + Element -> Creates a new system from two elements
OpticalSystem operator+(Element e1, Element e2) {
    OpticalSystem out = system(e1);
    out.add(e2);
    return out;
}

// Element + OpticalSystem -> Prepends an element to a system
OpticalSystem operator+(Element e, OpticalSystem sys) {
    OpticalSystem out = system(e, sys.ambient);
    out.add(sys.elements);
    return out;
}

// Medium + Element / Element[] -> Initializes system with specific ambient medium
OpticalSystem operator+(Medium ambient, Element e) { return system(e, ambient); }
OpticalSystem operator+(Medium ambient, Element[] es) { return system(es, ambient); }

// 3. Transformation Operator (*) is defined next to OpticalSystem:
// Usage: OpticalSystem shifted_sys = shift(10, 0) * rotate(45) * sys;

// 4. Angle-to-Polarization Cast Operator
// Usage: Ray r = Ray(p, d, pol = 45); // Automatically creates linear polarization at 45 degrees
Polarization operator cast(real angle) { return Polarization.linear(angle); }

// Draw dimension lines and labels
void distance(picture pic = currentpicture,
              string L = "",
              pair A,
              pair B,
              // 1. Structure & Geometry (How it's built)
              real offset = 0bp,
              bool rotated = true,
              bool extensions = true,
              // 2. Label Styling (How the text looks)
              align align = Center,
              filltype fill = UnFill(0pt),
              pen textpen = currentpen,
              // 3. Line Styling (How the line looks)
              pen linepen = grey + dashed,
              arrowbar arrow = Arrows(HookHead, 5pt),
              arrowbar bar = Bars,
              margin margin = PenMargins) {
    real D = abs(B - A) / 1mm;
    add(
        pic,
        new void(frame f, transform t) {
            pair A = t * A;
            pair B = t * B;

            string s_l = L;
            if (s_l == "") s_l = format("%#.2f", D);

            pair d = B - A;
            real dist = abs(d);
            if (dist == 0) return;

            pair u = unit(d);
            pair n = rotate(90) * u;

            // True positions of the dimension line in the frame
            pair A_dim = A + offset * n;
            pair B_dim = B + offset * n;

            // Extension lines
            if (extensions) draw(f, A -- A_dim ^^ B -- B_dim, linepen);

            // Local coordinates for the temporary picture
            pair p_a = rotated ? (0, 0) : A_dim;
            pair p_b = rotated ? (dist, 0) : B_dim;

            picture tmp;

            // Build the Label using only the function's explicit arguments
            Label lab = Label(s_l, align = align, p = textpen, filltype = fill);

            draw(tmp, lab, p_a -- p_b, linepen + solid, arrow = arrow, bar = bar, margin = margin);

            // Transform: shift by true start point if rotated, otherwise do nothing
            transform trans = rotated ? shift(A_dim) * rotate(degrees(u)) : identity();
            add(f, trans * tmp.fit());
        });
}

// Second overload wrapper for paths
void distance(picture pic = currentpicture,
              string L = "",
              path g,
              real offset = 0bp,
              bool rotated = true,
              bool extensions = true,
              align align = Center,
              filltype fill = UnFill(0pt),
              pen textpen = currentpen,
              pen linepen = grey + dashed,
              arrowbar arrow = Arrows(HookHead, 5pt),
              arrowbar bar = Bars,
              margin margin = PenMargins) {
    pair A = point(g, 0);
    pair B = point(g, length(g));

    // Explicitly forward everything using named arguments to ensure no mismatches
    distance(pic = pic,
             L = L,
             A = A,
             B = B,
             offset = offset,
             rotated = rotated,
             extensions = extensions,
             align = align,
             fill = fill,
             textpen = textpen,
             linepen = linepen,
             arrow = arrow,
             bar = bar,
             margin = margin);
}
