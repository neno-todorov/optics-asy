// optics.asy
//
// An attempt at a 2D optical ray tracing engine
//
// TODO:
//   1. make it work
//   2. make documentation
//   3. make more tests

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
    real stokesIntensityEps = 1e-12;
    real solverDenominatorEps = 1e-12;
}

struct PhysicsConfig {
    real lambda = 633;
    real minIntensity = 1e-3;
    real minRefractiveIndex = 0.1;
    real epsIntegration = 1e-4;
    real wavevectorAlignmentTolerance = 1e-6;
    real aimingDirectionEps = 1e-7;
    real aimingHeightTolerance = 1e-5;
    real far = 1000;
    real ds = 0.01;
    real paraxialHeight = 0.05;
    real differentialAngle = 0.5;
    int grinSolver = GRINSolverType.RK45;
    int maxDepth = 20;
    int maxBranches = 128;
    bool traceReflections = true;
    bool traceDiffraction = true;
}

struct GeometryConfig {
    real eps = 1e-8;
    real epsIntersection = 1e-8;
    real selfIntersectionFactor = 100;
    real nearIntersectionFactor = 10;
    real rayToPathExtension = 1000;
    real hemisphereMargin = 0.999;
    int maxBisectionSteps = 25;
    int arcSubdivisions = 64;
}

struct RenderingConfig {
    bool useIntensityOpacity = true;
    bool drawPolarization = false;
    bool drawWavefronts = false;
    bool drawReflections = false;
    bool drawAbsorptionFading = true;
    pen penRay = black;
    pen penBeam = black;
    pen penMirror = black;
    pen penAbsorber = black;
    pen penLens = black;
    real clipVelocityEps = 1e-9;
    real minLensHead = 10pt;
    real minDrawableLength = 1e-10;
    int fadeSubpathSteps = 32;
    int gratingTickCount = 12;
    real gratingTickAngle = 20;
    real gratingTickLength = 0.12;
    real polarizationSpacing = 0.5;
    real minPolarizationMarkerWithArrow = 0.5;
    real minPolarizationSegmentLength = 0.40;    // Segments shorter than this show 0 markers
    real minPolarizationMarkerSeparation = 0.25; // Minimum distance between consecutive markers
    real arrowSize = 5pt;
    real arrowPolSize = arrowSize / 3;
    arrowhead arrowHead = HookHead;
    arrowhead arrowPolHead = arrowHead;
    real polarizationInsetMinArrowFraction = 0.18;
    real polarizationInsetUnpolarizedThreshold = 1e-4;
    int polarizationInsetUnpolarizedSpokes = 4;
    int polarizationInsetEllipseSteps = 36;
    real wavefrontMaxDist = 120.0;
    int wavefrontMaxOrders = 3;
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
real wavelengthInUnits(real lambda_nm) { return lambda_nm / 1000.0; }

real vacuumWaveNumber(real lambda_nm) {
    real lambda = wavelengthInUnits(lambda_nm);
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

private
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

    real dop() {
        if (S0 <= config.math.eps) return 0.0;
        return clamp(sqrt(S1 * S1 + S2 * S2 + S3 * S3) / S0, 0.0, 1.0);
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
    void scaleTo(real newIntensity) {
        if (this.S0 > config.math.eps) {
            real factor = newIntensity / this.S0;
            this.S0 = newIntensity;
            this.S1 *= factor;
            this.S2 *= factor;
            this.S3 *= factor;
        } else {
            this.S0 = newIntensity;
        }
    }

    // Phase retardance between p and s components (Fresnel Rhomb / Waveplates)
    void applyRetardance(real delta) {
        real c = cos(delta);
        real s = sin(delta);
        real newS2 = S2 * c - S3 * s;
        real newS3 = S2 * s + S3 * c;
        S2 = newS2;
        S3 = newS3;
    }

    // Rotation of the polarization plane by theta (Optical Activity / Faraday)
    void applyRotation(real theta) {
        real c2 = cos(2.0 * theta);
        real s2 = sin(2.0 * theta);
        real newS1 = S1 * c2 - S2 * s2;
        real newS2 = S1 * s2 + S2 * c2;
        S1 = newS1;
        S2 = newS2;
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

    // 2. Linear polarization at angle theta_deg
    static Polarization linear(real theta_deg, real intensity = 1.0) {
        Polarization p = new Polarization;
        real th = radians(theta_deg);
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
    static Polarization fromAngles(real psi_deg, real chi_deg, real intensity = 1.0) {
        Polarization p = new Polarization;
        real psi = radians(psi_deg);
        real chi = radians(chi_deg);
        p.S0 = intensity;
        p.S1 = intensity * cos(2.0 * chi) * cos(2.0 * psi);
        p.S2 = intensity * cos(2.0 * chi) * sin(2.0 * psi);
        p.S3 = intensity * sin(2.0 * chi);
        return p;
    }

    // 5. Elliptical polarization via semi-axes (a, b), phase shift, and tilt
    static Polarization elliptical(real a = nan,
                                   real b = nan,
                                   real phase_deg = 90.0,
                                   real theta_deg = 0.0,
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
        real delta = radians(phase_deg);
        real theta = radians(theta_deg);

        real s0 = a_val * a_val + b_val * b_val;
        real s1 = a_val * a_val - b_val * b_val;
        real s2 = 2.0 * a_val * b_val * cos(delta);
        real s3 = 2.0 * a_val * b_val * sin(delta);

        real c2 = cos(2.0 * theta);
        real s2_th = sin(2.0 * theta);

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
};

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

// Compatibility wrapper for polarization rotation
void applyPolarizationRotation(Polarization pol, real theta) {
    if (pol != null) pol.applyRotation(theta);
}

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

// Applies a physical polarization angle rotation 'theta' (in radians)
// to both the coherent Jones field amplitudes and the Mueller Stokes vector.
void applyPolarizationRotation(Polarization pol, real theta) {
    if (pol == null) return;

    // Rotate Stokes vector (S1, S2) by 2 * theta (Mueller Rotator Matrix)
    real c = cos(2.0 * theta);
    real s = sin(2.0 * theta);

    real newS1 = pol.S1 * c - pol.S2 * s;
    real newS2 = pol.S1 * s + pol.S2 * c;

    pol.S1 = newS1;
    pol.S2 = newS2;
}

// Volume Properties (Medium Components)
// Anisotropic Boundary Packet & Function Types
struct AnisotropicTransmission {
    bool valid = false;
    pair waveNormal;
    pair groupDirection;
    real phaseIndex;
    real cosNormal;
}

struct PhaseMatch {
    bool tir = false;
    pair waveNormal;
    pair groupDirection;
    real n;
    real cosTheta;
}

// Function signatures
typedef real indexProfile(pair p, real lambda);
typedef real volumeIndexFunc(pair p, pair d, real lambda, Polarization pol);
typedef real volumeAbsorptionFunc(pair p, pair d, real lambda, Polarization pol);
typedef pair gradientProfile(pair p, real lambda);
typedef pair volumeGradientFunc(pair p, pair d, real lambda);
typedef pair groupDirectionFunc(pair waveNormal, pair p, real lambda);
typedef pair waveNormalFunc(pair groupDir, pair p, real lambda);
typedef AnisotropicTransmission anisotropicSolverFunc(pair N, pair T, real kt, pair forward, pair p, real lambda);
typedef Polarization volumeStepFunc(Polarization pol, pair p, pair d, real lambda, real ds);

struct VolumeProperty {
    string name = "";
    // Lower values run first; equal priorities retain add() order.
    int priority = 0;
    bool spatialIndex = false;
    volumeIndexFunc phaseIndex = null;
    volumeGradientFunc phaseGradient = null;
    volumeAbsorptionFunc attenuation = null;
    volumeStepFunc polarizationStep = null;
    groupDirectionFunc groupDir = null;
    waveNormalFunc waveNormal = null;
    anisotropicSolverFunc anisotropicSolver = null;
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
RefractiveIndex(real A = 1.0, real B = 0.0, indexProfile nProfile = null, gradientProfile gradProfile = null) {
    VolumeProperty vp = new VolumeProperty;
    vp.name = "RefractiveIndex";
    vp.spatialIndex = (nProfile != null);
    if (nProfile != null) {
        vp.phaseIndex = new real(pair p, pair d, real lambda, Polarization pol) {
            real l_um = wavelengthInUnits(lambda);
            real baseN = A + (l_um > 0 ? B / (l_um * l_um) : 0);
            return baseN * nProfile(p, lambda);
        };
    } else if (A != 1.0 || B != 0.0) {
        vp.phaseIndex = new real(pair p, pair d, real lambda, Polarization pol) {
            real l_um = wavelengthInUnits(lambda);
            return A + (l_um > 0 ? B / (l_um * l_um) : 0);
        };
    }
    if (gradProfile != null) {
        vp.phaseGradient = new pair(pair p, pair d, real lambda) {
            real l_um = wavelengthInUnits(lambda);
            real baseN = A + (l_um > 0 ? B / (l_um * l_um) : 0);
            return baseN * gradProfile(p, lambda);
        };
    }
    return vp;
}

VolumeProperty GRIN(indexProfile nProfile, gradientProfile gradProfile = null) {
    VolumeProperty vp = new VolumeProperty;
    vp.name = "GRIN";
    vp.spatialIndex = true;
    vp.phaseIndex = new real(pair p, pair d, real lambda, Polarization pol) { return nProfile(p, lambda); };
    if (gradProfile != null) {
        vp.phaseGradient = new pair(pair p, pair d, real lambda) { return gradProfile(p, lambda); };
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

VolumeProperty Dichroic(real alpha_s, real alpha_p) {
    VolumeProperty vp = new VolumeProperty;
    vp.name = "Dichroic";

    vp.polarizationStep = new Polarization(Polarization inPol, pair p, pair d, real lambda, real ds) {
        if (inPol == null || ds <= 0) return inPol;

        Polarization out = copy(inPol);
        real Ts = exp(-alpha_s * ds);
        real Tp = exp(-alpha_p * ds);
        // Transform Stokes parameters using the Mueller diattenuation matrix
        applyDiattenuation(out, inPol, Ts, Tp);
        return out;
    };
    return vp;
}

VolumeProperty UniaxialExtraordinary(real no, real ne, pair opticAxis = (1, 0)) {
    if (no <= 0 || ne <= 0) abort("Uniaxial refractive indices must be positive.");
    VolumeProperty vp = new VolumeProperty;
    vp.name = "Uniaxial Extraordinary";
    pair c_axis = unitOr(opticAxis);
    pair q_axis = perp(c_axis);
    vp.phaseIndex = new real(pair p, pair waveNormal, real lambda, Polarization pol) {
        pair u = unitOr(waveNormal);
        real invN2 = dot(u, c_axis) ^ 2 / (no * no) + dot(u, q_axis) ^ 2 / (ne * ne);
        return 1.0 / sqrt(max(config.math.eps, invN2));
    };
    vp.groupDir = new pair(pair waveNormal, pair p, real lambda) {
        pair k = unitOr(waveNormal);
        return unitOr(dot(k, c_axis) * c_axis / (no * no) + dot(k, q_axis) * q_axis / (ne * ne), k);
    };
    vp.waveNormal = new pair(pair groupDir, pair p, real lambda) {
        pair d = unitOr(groupDir);
        return unitOr((dot(d, c_axis) * no * no) * c_axis + (dot(d, q_axis) * ne * ne) * q_axis, d);
    };
    vp.anisotropicSolver = new AnisotropicTransmission(pair N, pair T, real kt, pair forward, pair p, real lambda) {
        AnisotropicTransmission out = new AnisotropicTransmission;
        real a = dot(N, c_axis) ^ 2 / (no * no) + dot(N, q_axis) ^ 2 / (ne * ne);
        real b = 2 * kt * (dot(N, c_axis) * dot(T, c_axis) / (no * no) + dot(N, q_axis) * dot(T, q_axis) / (ne * ne));
        real cc = kt * kt * (dot(T, c_axis) ^ 2 / (no * no) + dot(T, q_axis) ^ 2 / (ne * ne)) - 1;
        real discriminant = b * b - 4 * a * cc;
        if (a <= config.math.eps || discriminant < -config.math.eps) return out;
        real root = sqrt(max(0.0, discriminant));
        real candidates[] = {(-b + root) / (2 * a), (-b - root) / (2 * a)};
        real bestForward = -infinity;
        for (real kn : candidates) {
            pair k = kn * N + kt * T;
            pair waveNormal = unitOr(k);
            pair group =
                unitOr(dot(waveNormal, c_axis) * c_axis / (no * no) + dot(waveNormal, q_axis) * q_axis / (ne * ne),
                       waveNormal);
            real score = dot(group, forward);
            if (score > bestForward) {
                bestForward = score;
                out.valid = (score > config.math.eps);
                out.waveNormal = waveNormal;
                out.groupDirection = group;
                out.phaseIndex = length(k);
                out.cosNormal = abs(dot(waveNormal, N));
            }
        }
        return out;
    };
    return vp;
}

// Optical Activity / Natural Polarization Rotation
VolumeProperty OpticalActivity(real rot_deg_per_unit) {
    VolumeProperty vp = new VolumeProperty;
    vp.name = "OpticalActivity";

    vp.polarizationStep = new Polarization(Polarization inPol, pair p, pair d, real lambda, real ds) {
        if (inPol == null || ds <= 0) return inPol;
        Polarization out = copy(inPol);

        real theta = radians(rot_deg_per_unit * ds);
        applyPolarizationRotation(out, theta);

        return out;
    };
    return vp;
}

// Faraday Rotation: theta = V * (B . d) * ds (Non-reciprocal magnetic rotation)
VolumeProperty FaradayRotator(real V, pair B = (1, 0)) {
    VolumeProperty vp = new VolumeProperty;
    vp.name = "FaradayRotator";

    vp.polarizationStep = new Polarization(Polarization inPol, pair p, pair d, real lambda, real ds) {
        if (inPol == null || ds <= 0) return inPol;
        Polarization out = copy(inPol);

        real theta = radians(V * dot(B, d) * ds);
        applyPolarizationRotation(out, theta);

        return out;
    };
    return vp;
}

// Modular Medium Struct
// Bulk Volume Physics:
// Phase Velocity / Eikonal Equation   -->  n(p, d, lambda)
// Attenuation Coefficient             -->  absorption(p, d, lambda)
// Continuous Field Evolution (ds)     -->  stepPolarization(pol, ds)
// Constitutive Relation (k <-> S)     -->  groupDirection(), waveNormalFromGroup()
struct Medium {
    string name = "";
    VolumeProperty[] props;
    Medium add(VolumeProperty prop) {
        props.push(prop);
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
    // Refractive-index profiles multiply in stage/priority/declaration order.
    // This models a base material with dimensionless spatial perturbations.
    real n(pair p = (0, 0), pair d = (1, 0), real lambda = config.physics.lambda, Polarization pol = null) {
        real val = 1.0;
        bool hasCustom = false;
        for (VolumeProperty vp : orderedProperties()) {
            if (vp.phaseIndex != null) {
                real curN = vp.phaseIndex(p, d, lambda, pol);
                val = hasCustom ? (val * curN) : curN;
                hasCustom = true;
            }
        }
        return val;
    }
    // Cumulative Absorption/Gain
    // The Beer–Lambert absorption coefficient alpha.
    real absorption(pair p = (0, 0), pair d = (1, 0), real lambda = config.physics.lambda, Polarization pol = null) {
        real totalAlpha = 0.0;
        for (VolumeProperty vp : orderedProperties()) {
            if (vp.attenuation != null) totalAlpha += vp.attenuation(p, d, lambda, pol);
        }
        return totalAlpha;
    }
    // Polarization & Loss Stepping (Chains all steps + Beer-Lambert decay)
    // Evolves the (Es, Ep) Jones vector over a small distance ds (applying dichroism, optical activity/rotation, and
    // continuous absorption decay).
    Polarization stepPolarization(Polarization inPol, pair p, pair d, real lambda, real ds) {
        if (inPol == null || ds <= 0) return inPol;
        Polarization currPol = copy(inPol);
        for (VolumeProperty vp : orderedProperties()) {
            if (vp.polarizationStep != null) {
                currPol = vp.polarizationStep(currPol, p, d, lambda, ds);
            }
        }
        real a = absorption(p, d, lambda, inPol);
        if (a != 0) {
            // Field amplitudes scale as exp(-alpha*ds/2),
            // so intensity scales as exp(-alpha*ds).
            real decay = exp(-a * ds / 2.0);

            // Keep Stokes intensity consistent with the field amplitudes.
            // Unpolarized intensity() reads S0 directly.
            currPol.S0 *= decay * decay;
        }

        return currPol;
    }
    // Anisotropic Direction & Normal Conversion (Pure delegation)
    // Maps wavevector k to Poynting vector S (and vice versa).
    pair groupDirection(pair k, pair p = (0, 0), real lambda = config.physics.lambda) {
        pair out = unitOr(k);
        for (VolumeProperty vp : orderedProperties()) {
            if (vp.groupDir != null) out = vp.groupDir(k, p, lambda);
        }
        return out;
    }
    pair waveNormalFromGroup(pair d, pair p = (0, 0), real lambda = config.physics.lambda) {
        pair out = unitOr(d);
        for (VolumeProperty vp : orderedProperties()) {
            if (vp.waveNormal != null) out = vp.waveNormal(d, p, lambda);
        }
        return out;
    }
}

Medium
medium(string name = "",
       real A = 1,
       real B = 0,
       real alpha = 0.0,
       indexProfile nProfile = null,
       gradientProfile gradProfile = null) {
    Medium m = new Medium;
    m.name = name;
    if (A != 1.0 || B != 0.0 || nProfile != null) m.add(RefractiveIndex(A, B, nProfile, gradProfile));
    if (alpha != 0) m.add(Absorption(alpha));
    return m;
}

Medium uniaxialExtraordinary(string name, real no, real ne, pair opticAxis = (1, 0)) {
    Medium m = new Medium;
    m.name = name;
    m.add(UniaxialExtraordinary(no, ne, opticAxis));
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
    int flatIndex = -1;
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
    out.flatIndex = src.flatIndex;
    return out;
}

struct Ray {
    int id = -1;
    int parentId = -1;
    int generation = 0;
    int[] childrenIds = new int[];
    int defaultChildIndex = 0;
    pair p;
    // Wave-normal direction; `d` is the geometric/group direction.  They are
    // identical in isotropic media and differ for extraordinary rays.
    pair waveNormal;
    pair d;
    real lambda;
    real opl;
    real phase;
    Medium medium;
    bool alive = true;
    bool terminated = false;
    int[] ignoredFlatIndices = new int[];
    Hit hit = null;
    path trajectory = nullpath;
    // REPLACED: Es, Ep, and unpolarized are now encapsulated!
    Polarization pol;
    static int Incident = 1;
    static int Reflected = 2;
    static int Transmitted = 3;
    pair point(real t) { return p + t * d; }
    real intensity() { return pol.intensity(); } // Delegates to pol!
    void initTree(int newId, int newParentId, int newGeneration) {
        this.id = newId;
        this.parentId = newParentId;
        this.generation = newGeneration;
        this.childrenIds = new int[];
        this.defaultChildIndex = 0;
    }
    void operator init(pair p = (0, 0),
                       pair d = (1, 0),
                       Medium medium = null,
                       real lambda = config.physics.lambda,
                       Polarization pol = Polarization.unpolarized(),
                       real phase = 0) {
        this.id = -1;
        this.parentId = -1;
        this.generation = 0;
        this.childrenIds = new int[];
        this.defaultChildIndex = 0;
        this.p = p;
        this.d = unitOr(d);
        this.waveNormal = this.d;
        this.medium = medium;
        this.lambda = lambda;
        this.phase = phase;
        // Directly assign the incoming polarization state
        this.pol = copy(pol);
        this.opl = 0;
        this.alive = true;
        this.terminated = false;
        this.ignoredFlatIndices = new int[];
        this.hit = null;
        this.trajectory = nullpath;
    }
    // STATIC FACTORY NAMESPACES
    static Ray angle(pair p,
                     real angle,
                     Medium medium = null,
                     real lambda = config.physics.lambda,
                     Polarization pol = Polarization.unpolarized(),
                     real phase = 0) {
        return Ray(p, dir(angle), medium, lambda, pol, phase);
    }
    static Ray between(pair a,
                       pair b,
                       Medium medium = null,
                       real lambda = config.physics.lambda,
                       Polarization pol = Polarization.unpolarized(),
                       real phase = 0) {
        return Ray(a, b - a, medium, lambda, pol, phase);
    }
    static Ray toPoint(pair target,
                       pair d,
                       real distance,
                       Medium medium = null,
                       real lambda = config.physics.lambda,
                       Polarization pol = Polarization.unpolarized(),
                       real phase = 0) {
        pair dHat = unitOr(d);
        return Ray(target - distance * dHat, dHat, medium, lambda, pol, phase);
    }
    static Ray toPoint(pair target,
                       real angle,
                       real distance,
                       Medium medium = null,
                       real lambda = config.physics.lambda,
                       Polarization pol = Polarization.unpolarized(),
                       real phase = 0) {
        return toPoint(target, dir(angle), distance, medium, lambda, pol, phase);
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
    r.parentId = src.parentId;
    r.generation = src.generation;
    r.childrenIds = copy(src.childrenIds);
    r.defaultChildIndex = src.defaultChildIndex;
    r.p = src.p;
    r.waveNormal = src.waveNormal;
    r.d = src.d;
    r.medium = src.medium;
    r.lambda = src.lambda;
    r.phase = src.phase;
    r.opl = src.opl;
    r.pol = copy(src.pol);
    r.alive = src.alive;
    r.terminated = src.terminated;
    r.ignoredFlatIndices = copy(src.ignoredFlatIndices);
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
    out.waveNormal = -r.waveNormal;
    return out;
}

Ray spawnRay(Ray src,
             Hit h,
             pair direction,
             real intensity,
             Medium medium,
             real phaseShift = 0.0,
             pair waveNormal = (0, 0)) {
    Ray out = copy(src);
    out.initTree(-1, src.id, src.generation + 1);
    out.d = unitOr(direction, src.d);
    out.waveNormal = unitOr(waveNormal, out.d);
    out.p = h.p;
    out.medium = medium;
    // Scale polarization state
    out.pol.scaleTo(intensity);
    // Boundary phase shift
    out.phase += phaseShift;
    out.terminated = false;
    // KEEP EXISTING flatIndex LOGIC EXACTLY AS IT IS:
    if (h.t > config.geometry.epsIntersection) {
        if (h.flatIndex >= 0)
            out.ignoredFlatIndices = new int[]{h.flatIndex};
        else
            out.ignoredFlatIndices = new int[];
    } else {
        out.ignoredFlatIndices = copy(src.ignoredFlatIndices);
        if (h.flatIndex >= 0) out.ignoredFlatIndices.push(h.flatIndex);
    }
    out.hit = null;
    out.trajectory = nullpath;
    return out;
}

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
    if (m == null) {
        res.n = 1.0;
        res.grad = (0, 0);
        return res;
    }
    real indexProduct = 1.0;
    pair gradient = (0, 0);
    bool hasIndex = false;
    bool hasCompleteAnalyticGradient = true;
    for (VolumeProperty vp : m.orderedProperties()) {
        if (vp.phaseIndex == null) continue;
        real localIndex = vp.phaseIndex(p, d, lambda, pol);
        if (vp.phaseGradient == null) {
            hasCompleteAnalyticGradient = false;
        } else {
            // Product rule for Medium.n(), whose index profiles compose
            // multiplicatively: d(a*b) = b da + a db.
            gradient = localIndex * gradient + indexProduct * vp.phaseGradient(p, d, lambda);
        }
        indexProduct *= localIndex;
        hasIndex = true;
    }
    res.n = hasIndex ? indexProduct : m.n(p, d, lambda, pol);
    if (hasIndex && hasCompleteAnalyticGradient) {
        res.grad = gradient;
        return res;
    }
    real eps = config.physics.epsIntegration;
    real dx = (m.n(p + (eps, 0), d, lambda, pol) - m.n(p - (eps, 0), d, lambda, pol)) / (2 * eps);
    real dy = (m.n(p + (0, eps), d, lambda, pol) - m.n(p - (0, eps), d, lambda, pol)) / (2 * eps);
    res.grad = (dx, dy);
    return res;
}

private
pair mediumGradient(Medium m, pair p, pair d, real lambda) { return mediumIndexAndGradient(m, p, d, lambda).grad; }

private
pair grinDirectionDerivative(Ray r, pair p, pair d) {
    pair u = (length(d) > config.math.eps) ? unit(d) : r.d;
    IndexAndGradient ig = mediumIndexAndGradient(r.medium, p, u, r.lambda, r.pol);
    if (ig.n <= config.physics.minRefractiveIndex) return (0, 0);
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
    res.ray.waveNormal = r.medium.waveNormalFromGroup(res.ray.d, res.ray.p, res.ray.lambda);
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
    res.ray.waveNormal = r.medium.waveNormalFromGroup(res.ray.d, res.ray.p, res.ray.lambda);
    res.error = length(p_5th - p_4th) + length(d_5th - d_4th);
    return res;
}

private
struct Integrator {
    string name;
    bool isAdaptive;
    GRINStepFunction advance;
    static Integrator RK4() {
        Integrator ig = new Integrator;
        ig.name = "RK4";
        ig.isAdaptive = false;
        ig.advance = grinAdvanceRK4;
        return ig;
    }
    static Integrator RK45() {
        Integrator ig = new Integrator;
        ig.name = "RK45";
        ig.isAdaptive = true;
        ig.advance = grinAdvanceRK45;
        return ig;
    }
}

// Collision Packet & Plugin Results
// Collision Packet & Plugin Results
private struct InteractionResult {
    Ray[] rays;
    int defaultIndex;
    bool absorbed;
}

private InteractionResult
interactionResult(Ray[] rays = new Ray[], int defaultIndex = 0, bool absorbed = false) {
    InteractionResult res = new InteractionResult;
    res.rays = rays;
    res.defaultIndex = defaultIndex;
    res.absorbed = absorbed;
    return res;
}

private
struct PropagationResult {
    Ray ray;
    bool hitSurface;
    Hit hit;
    path trajectory;
}

private PropagationResult
propagationResult(Ray ray, bool hitSurface = false, Hit hit = null, path trajectory = nullpath) {
    PropagationResult res = new PropagationResult;
    res.ray = ray;
    res.hitSurface = hitSurface;
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
    real cosForward;
    real cosi;
    bool fromMedium1;
    Medium sourceMedium;
    Medium targetMedium;
}

typedef Medium
mediumResolver(pair, Ray, Medium, Medium, bool);

private
InterfaceFrame interfaceFrame(Medium m1, Medium m2, Ray r, Hit h, mediumResolver resolve);

// private
struct FresnelResult {
    // -------------------------------------------------------------------------
    // Backward-Compatible Primary Fields (Used by standard lenses, mirrors, tests)
    // -------------------------------------------------------------------------
    pair d;
    pair waveNormal;
    Medium medium;
    real T = 1.0;
    real R = 0.0;
    bool totalInternalReflection = false;

    // Boundary Context
    pair normal;
    Medium sourceMedium;
    real n1 = 1.0;
    real n2 = 1.0;
    real cosi = 1.0;
    real cost = 1.0;

    // Reflection
    pair reflectedDirection;
    pair reflectedWaveNormal;
    pair rs = (0, 0);
    pair rp = (0, 0);
    real Rs = 0.0;
    real Rp = 0.0;

    // -------------------------------------------------------------------------
    // Mode-Specific Fields (For Anisotropic Crystals / Birefringence)
    // -------------------------------------------------------------------------
    bool isBirefringent = false;

    // Channel 1: Ordinary / TE / S
    pair d_o;
    pair waveNormal_o;
    pair ts = (0, 0); // or 'to'
    real To = 0.0;
    bool tir_o = false;
    Medium medium_o;

    // Channel 2: Extraordinary / TM / P
    pair d_e;
    pair waveNormal_e;
    pair tp = (0, 0); // or 'te'
    real Te = 0.0;
    bool tir_e = false;
    Medium medium_e;

    real Ts() { return 1.0 - Rs; }
    real Tp() { return 1.0 - Rp; }

    // -------------------------------------------------------------------------
    // Synchronizer: Call at the end of solver to guarantee all aliases match
    // -------------------------------------------------------------------------
    void sync() {
        if (!this.isBirefringent) {
            // In isotropic media: both modes share the primary path
            this.d_o = this.d;
            this.d_e = this.d;
            this.waveNormal_o = this.waveNormal;
            this.waveNormal_e = this.waveNormal;
            this.medium_o = this.medium;
            this.medium_e = this.medium;
            this.To = this.Ts();
            this.Te = this.Tp();
            this.tir_o = this.totalInternalReflection;
            this.tir_e = this.totalInternalReflection;
        } else {
            // In crystals: primary aliases default to the ordinary mode
            this.d = this.d_o;
            this.waveNormal = this.waveNormal_o;
            this.medium = this.medium_o;
            this.totalInternalReflection = (this.tir_o && this.tir_e);
            this.T = this.To + this.Te;
        }
    }
};

// Refraction
// fresnelIsotropic(Medium m1, Medium m2, Ray r, Hit h, mediumResolver resolve, real reflectance = nan);

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
pair diffractionOutputDir(pair N, pair T, real kOut, bool reflective) {
    real kNormal = sqrt(max(0.0, 1.0 - kOut * kOut));
    return kOut * T + (reflective ? kNormal * N : -kNormal * N);
}

typedef Ray[] rayGenerator(Ray r, Hit h, mediumResolver resolve);
typedef pair propertyAlignFunc(pair inDir, pair targetDir, real lambda, int orderIndex);

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
    transform invT = identity();
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
    Hit intersect(Ray r, real maxDistance = infinity, bool ignoreOrigin = false) {
        Hit h = new Hit;
        h.valid = false;
        h.t = infinity;
        if (length(this.p) == 0) return h;
        // Uses the precomputed inverse! No inverse() computation during tracing!
        real rayReach = (maxDistance < infinity) ? maxDistance : config.physics.far;
        pair local_p = this.invT * r.p;
        pair local_end = this.invT * (r.p + rayReach * r.d);
        real[] hits = intersections(this.p, local_p, local_end);
        real bestTime = -1;
        real bestDistance = infinity;
        for (real u : hits) {
            pair pt = point(this.p, u);
            pair global_pt = this.T * pt;
            pair global_q = global_pt - r.p;
            real d = length(global_q);
            if (ignoreOrigin && d <= config.geometry.selfIntersectionFactor * config.geometry.epsIntersection) {
                continue;
            }
            if (d <= maxDistance + config.geometry.epsIntersection && d < bestDistance) {
                if (d <= config.geometry.nearIntersectionFactor * config.geometry.epsIntersection ||
                    dot(global_q, r.d) > config.geometry.eps) {
                    bestDistance = d;
                    bestTime = u;
                }
            }
        }
        if (bestTime < 0) return h;
        h.valid = true;
        pair local_hit_pt = point(this.p, bestTime);
        pair local_tangent = dir(this.p, bestTime);
        pair local_normal = (local_tangent.y, -local_tangent.x);
        h.p = this.T * local_hit_pt;
        h.tangent = unitOr(this.T * (local_hit_pt + local_tangent) - h.p, local_tangent);
        h.normal = unitOr(this.T * (local_hit_pt + local_normal) - h.p, local_normal);
        h.center = this.center;
        h.axis = this.axis;
        h.t = bestDistance;
        h.entering = dot(r.d, h.normal) >= 0;
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

        // Seed the ray at the boundary:
        Ray[] active = {spawnRay(r, h, r.d, r.intensity(), r.medium)};
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
        // - defaultIndex: 0 (default ray to follow in tree)
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
Surface transformedSurface(Surface src, transform parentT) {
    Surface out = copy(src);
    out.T = parentT * src.T;
    out.invT = inverse(out.T); // Precalculated ONCE here!
    out.center = parentT * src.center;
    out.axis = transformedDirection(parentT, src.center, src.axis);
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
    real focalLength = nan; // Optional visual metadata for lenses/mirrors
    pen stroke;
    pen fill = nullpen;
    bool showInfo;
    void customDraw(picture pic, Element e) = null;
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
    Surface[] flatten(transform parentT = identity()) {
        Surface[] result;
        transform totalT = parentT * this.T;
        // Transform and add own surfaces
        for (Surface s : this.surfaces) {
            result.push(transformedSurface(s, totalT));
        }
        // Recursively flatten children and append in 1 line
        for (Element child : this.children) {
            result.append(child.flatten(totalT)); // <-- Native append()
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
    e.showInfo = src.showInfo;
    e.boundary = src.boundary;
    e.medium = src.medium;
    e.med_e = src.med_e;
    e.customDraw = src.customDraw;
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
searchEnclosing(Element[] list, pair p, transform currentT, ElementMatch best) {
    for (int i = 0; i < list.length; ++i) {
        Element e = list[i];
        transform totalT = currentT * e.T;
        // Only cyclic (closed) boundaries define a volumetric medium
        if (length(e.boundary) > 0 && cyclic(e.boundary)) {
            path b = totalT * e.boundary;
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
            searchEnclosing(e.children, p, totalT, best);
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
Medium resolveMediumAt(Element[] elements, Medium ambient, pair p, Ray r, Medium m1, Medium m2, bool fromM1) {
    // Explicit target override:
    // If the surface explicitly hardcodes an outside medium (m2 != null), use it.
    // Otherwise, null means "dynamically resolve the medium at point p".
    Medium explicitTarget = fromM1 ? m2 : m1;
    if (explicitTarget != null) return explicitTarget;
    // System Ambient Fallback
    Medium sysAmbient = (ambient != null) ? ambient : Materials.VACUUM;
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
    return sysAmbient;
}

interfaceFrame = new InterfaceFrame(Medium m1, Medium m2, Ray r, Hit h, mediumResolver resolve) {
    InterfaceFrame f = new InterfaceFrame;
    pair N = unitOr(h.normal);
    real cosForward = dot(r.d, N);
    bool fromM1 = cosForward >= 0;
    f.normal = fromM1 ? -N : N;
    f.cosForward = cosForward;
    f.cosi = abs(dot(unitOr(r.waveNormal, r.d), N));
    f.fromMedium1 = fromM1;
    // Source medium is ALWAYS the ray's current physical medium
    f.sourceMedium = (r.medium != null) ? r.medium : Materials.VACUUM;
    // Target medium is resolved spatially in transmission direction
    pair targetPoint = h.p + config.geometry.epsIntersection * r.d;
    f.targetMedium = resolve(targetPoint, r, m1, m2, fromM1);
    // Failsafe
    if (f.sourceMedium == null) f.sourceMedium = Materials.VACUUM;
    if (f.targetMedium == null) f.targetMedium = Materials.VACUUM;
    return f;
};

private
PhaseMatch solvePhaseMatch(Medium m, pair N, pair T, real kt, pair forward, pair p, real lambda) {
    PhaseMatch pm = new PhaseMatch;
    real localN = m.n(p, forward, lambda);
    pm.n = localN;
    // A constitutive anisotropic solver defines the dispersion relation; two
    // independent solvers cannot be composed without an explicit model.
    anisotropicSolverFunc solver = null;
    for (VolumeProperty vp : m.orderedProperties()) {
        if (vp.anisotropicSolver != null) {
            if (solver != null) abort("A Medium may have one anisotropic constitutive solver.");
            solver = vp.anisotropicSolver;
        }
    }
    if (solver != null) {
        AnisotropicTransmission at = solver(N, T, kt, forward, p, lambda);
        pm.tir = !at.valid;
        pm.waveNormal = at.waveNormal;
        pm.groupDirection = at.groupDirection;
        pm.n = (at.phaseIndex > 0) ? at.phaseIndex : localN;
        pm.cosTheta = at.cosNormal;
        return pm;
    }
    // Isotropic dispersion sphere
    real kn2 = localN * localN - kt * kt;
    if (kn2 < 0) {
        pm.tir = true;
        pm.cosTheta = 0;
        return pm;
    }
    real kn = sqrt(max(0.0, kn2));
    real sign = (dot(N, forward) >= 0) ? 1.0 : -1.0;
    pair kVec = kt * T + (sign * kn) * N;
    pm.tir = false;
    pm.waveNormal = unitOr(kVec, forward);
    pm.groupDirection = pm.waveNormal;
    pm.cosTheta = (localN > config.math.eps) ? (kn / localN) : 0.0;
    return pm;
}

// ============================================================================
// Maxwell Boundary Solver for Isotropic <-> Uniaxial Planar Interfaces
// ============================================================================
private
FresnelResult fresnelUniaxial(Ray r, Hit h, Medium med_iso, Medium med_o, Medium med_e, mediumResolver resolve) {
    FresnelResult res;
    res.isBirefringent = true;

    // 1. Establish orthonormal interface frame
    InterfaceFrame f = interfaceFrame(med_iso, med_o, r, h, resolve);
    pair N = f.normal;
    pair T = (-N.y, N.x); // Unit tangent along surface

    pair incWaveNormal = unitOr(r.waveNormal, r.d);
    real n1 = f.sourceMedium.n(h.p, incWaveNormal, r.lambda, r.pol);
    real kt = n1 * dot(incWaveNormal, T);
    real kn_in = abs(n1 * dot(incWaveNormal, N));

    // Boundary metadata
    res.normal = f.normal;
    res.sourceMedium = f.sourceMedium;
    res.n1 = n1;
    res.cosi = f.cosi;

    // Reflected wave kinematics
    PhaseMatch reflMatch = solvePhaseMatch(f.sourceMedium, N, T, kt, N, h.p, r.lambda);
    res.reflectedDirection = reflMatch.groupDirection;
    res.reflectedWaveNormal = reflMatch.waveNormal;

    // 2. Transmitted Kinematics using existing functional boundary matcher
    pair intoMed = -N; // Normal pointing into crystal
    PhaseMatch trans_o = solvePhaseMatch(med_o, N, T, kt, intoMed, h.p, r.lambda);
    PhaseMatch trans_e = solvePhaseMatch(med_e, N, T, kt, intoMed, h.p, r.lambda);

    res.d_o = trans_o.groupDirection;
    res.waveNormal_o = trans_o.waveNormal;
    res.tir_o = trans_o.tir;
    res.medium_o = med_o;

    res.d_e = trans_e.groupDirection; // Contains exact Poynting walk-off
    res.waveNormal_e = trans_e.waveNormal;
    res.tir_e = trans_e.tir;
    res.medium_e = med_e;

    // ========================================================================
    // 3. Ordinary Mode (TE / S-Polarization)
    // ========================================================================
    if (!res.tir_o) {
        real qo = abs(trans_o.n * dot(trans_o.waveNormal, N));
        real denom = kn_in + qo;
        res.rs = (denom > config.math.eps) ? ((kn_in - qo) / denom, 0.0) : (-1.0, 0.0);
        res.ts = (denom > config.math.eps) ? ((2.0 * kn_in) / denom, 0.0) : (0.0, 0.0);
        // res.to = res.ts;
        res.Rs = abs2(res.rs);
        res.To = 1.0 - res.Rs;
    } else {
        real no = med_o.n(h.p, incWaveNormal, r.lambda, r.pol);
        real qo2 = no * no - kt * kt;
        real qo_im = sqrt(max(0.0, -qo2));
        pair denom = (kn_in, qo_im);
        res.rs = (kn_in, -qo_im) / denom;
        res.ts = (0.0, 0.0);
        // res.to = (0.0, 0.0);
        res.Rs = 1.0;
        res.To = 0.0;
    }

    // ========================================================================
    // 4. Extraordinary Mode (TM / P-Polarization via Kinematic Impedance)
    // ========================================================================
    if (!res.tir_e) {
        pair we = trans_e.n * trans_e.waveNormal;
        pair de = trans_e.groupDirection;

        real dn = abs(dot(de, N));
        real we_dot_de = abs(dot(we, de));

        // Exact Maxwell boundary surface impedance: Ze = dn / (we . de)
        real Ze = (we_dot_de > config.math.eps) ? (dn / we_dot_de) : (1.0 / trans_e.n);
        real Zi = kn_in / (n1 * n1); // Incident TM wave impedance

        real denom_p = Zi + Ze;
        res.rp = (abs(denom_p) > config.math.eps) ? ((Zi - Ze) / denom_p, 0.0) : (-1.0, 0.0);
        res.tp = (abs(denom_p) > config.math.eps) ? ((2.0 * Zi) / denom_p, 0.0) : (0.0, 0.0);
        // res.te = res.tp;
        res.Rp = abs2(res.rp);
        res.Te = 1.0 - res.Rp;
    } else {
        res.rp = (-1.0, 0.0);
        res.tp = (0.0, 0.0);
        // res.te = (0.0, 0.0);
        res.Rp = 1.0;
        res.Te = 0.0;
    }

    // Total incident-weighted reflectance
    real I_in = r.intensity();
    if (I_in > config.math.eps) {
        real Is = max(0.0, 0.5 * (r.pol.S0 + r.pol.S1));
        real Ip = max(0.0, 0.5 * (r.pol.S0 - r.pol.S1));
        res.R = (res.Rs * Is + res.Rp * Ip) / I_in;
    } else {
        res.R = 0.5 * (res.Rs + res.Rp);
    }
    res.T = 1.0 - res.R;

    res.sync();
    return res;
}

FresnelResult fresnelIsotropic(Medium m1, Medium m2, Ray r, Hit h, mediumResolver resolve, real reflectance = nan) {
    FresnelResult res;
    res.isBirefringent = false;

    InterfaceFrame f = interfaceFrame(m1, m2, r, h, resolve);
    pair N = f.normal;
    pair T = (-N.y, N.x);
    pair incWaveNormal = unitOr(r.waveNormal, r.d);
    real n1 = f.sourceMedium.n(h.p, incWaveNormal, r.lambda, r.pol);
    real kt = n1 * dot(incWaveNormal, T);

    // Unified Transmitted & Reflected Modes
    PhaseMatch trans = solvePhaseMatch(f.targetMedium, N, T, kt, -N, h.p, r.lambda);
    PhaseMatch refl = solvePhaseMatch(f.sourceMedium, N, T, kt, N, h.p, r.lambda);
    real n2 = max(config.physics.minRefractiveIndex, trans.n);

    res.normal = f.normal;
    res.sourceMedium = f.sourceMedium;
    res.n1 = n1;
    res.n2 = n2;
    res.cosi = f.cosi;
    res.reflectedWaveNormal = refl.waveNormal;
    res.reflectedDirection = refl.groupDirection;

    pair cosi = (f.cosi, 0);

    if (trans.tir) {
        // Total Internal Reflection (TIR)
        res.totalInternalReflection = true;
        res.d = refl.groupDirection;
        res.waveNormal = refl.waveNormal;
        res.medium = f.sourceMedium;

        // Exact complex evanescent wave normal: cost = -i * sqrt((kt/n2)^2 - 1)
        real k_im = (n2 > config.math.eps) ? sqrt(max(0.0, (kt / n2) ^ 2 - 1.0)) : 0.0;
        pair cost = (0, -k_im);
        pair denom_s = n1 * cosi + n2 * cost;
        pair denom_p = n2 * cosi + n1 * cost;

        res.rs = (abs(denom_s) > config.math.eps) ? (n1 * cosi - n2 * cost) / denom_s : (1, 0);
        res.rp = (abs(denom_p) > config.math.eps) ? (n2 * cosi - n1 * cost) / denom_p : (1, 0);
        res.ts = (0, 0);
        res.tp = (0, 0);
        // res.to = (0, 0);
        // res.te = (0, 0);
        res.Rs = 1.0;
        res.Rp = 1.0;
        res.R = 1.0;
        res.T = 0.0;
        res.sync();
        return res;
    }

    // Standard Isotropic Refraction
    res.totalInternalReflection = false;
    res.cost = trans.cosTheta;
    res.waveNormal = trans.waveNormal;
    res.d = trans.groupDirection;
    res.medium = f.targetMedium;

    pair cost = (trans.cosTheta, 0);
    pair denom_s = n1 * cosi + n2 * cost;
    pair denom_p = n2 * cosi + n1 * cost;

    res.rs = (abs(denom_s) > config.math.eps) ? (n1 * cosi - n2 * cost) / denom_s : (0, 0);
    res.rp = (abs(denom_p) > config.math.eps) ? (n2 * cosi - n1 * cost) / denom_p : (0, 0);
    res.ts = (abs(denom_s) > config.math.eps) ? (2 * n1 * cosi) / denom_s : (0, 0);
    res.tp = (abs(denom_p) > config.math.eps) ? (2 * n1 * cosi) / denom_p : (0, 0);
    // res.to = res.ts;
    // res.te = res.tp;

    res.Rs = abs2(res.rs);
    res.Rp = abs2(res.rp);

    real I_in = r.intensity();
    if (I_in > config.math.eps) {
        real Is = max(0.0, 0.5 * (r.pol.S0 + r.pol.S1));
        real Ip = max(0.0, 0.5 * (r.pol.S0 - r.pol.S1));
        res.R = (res.Rs * Is + res.Rp * Ip) / I_in;
    } else {
        res.R = 0.5 * (res.Rs + res.Rp);
    }

    res.R = isnan(reflectance) ? clamp(res.R, 0.0, 1.0) : clamp(reflectance, 0.0, 1.0);
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
    transform invT; // Precomputed inverse!
    pair minBound;
    pair maxBound;
}

private FlatSurfaceRecord[] prepareFlatSurfaces(Element[] elements) {
    FlatSurfaceRecord[] records;
    for (int i = 0; i < elements.length; ++i) {
        Surface[] flat = elements[i].flatten();
        for (int j = 0; j < flat.length; ++j) {
            FlatSurfaceRecord rec = new FlatSurfaceRecord;
            rec.s = flat[j];
            rec.invT = inverse(flat[j].T); // Inverted ONCE at setup!
            if (length(flat[j].p) > 0) {
                path globalPath = flat[j].T * flat[j].p;
                rec.minBound = min(globalPath);
                rec.maxBound = max(globalPath);
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
bool rayIntersectsAABB(pair ro, pair rd, pair bmin, pair bmax, real maxDist) {
    real eps = config.geometry.epsIntersection;
    // X-axis slab
    real invDx = abs(rd.x) > 1e-12 ? 1.0 / rd.x : (rd.x >= 0 ? 1e12 : -1e12);
    real tx1 = (bmin.x - eps - ro.x) * invDx;
    real tx2 = (bmax.x + eps - ro.x) * invDx;
    real tmin = min(tx1, tx2);
    real tmax = max(tx1, tx2);
    // Y-axis slab
    real invDy = abs(rd.y) > 1e-12 ? 1.0 / rd.y : (rd.y >= 0 ? 1e12 : -1e12);
    real ty1 = (bmin.y - eps - ro.y) * invDy;
    real ty2 = (bmax.y + eps - ro.y) * invDy;
    tmin = max(tmin, min(ty1, ty2));
    tmax = min(tmax, max(ty1, ty2));
    // Must enter before exiting, exit in front of ray, and enter within reach
    return (tmax >= 0.0) && (tmin <= tmax) && (tmin <= maxDist + eps);
}

private
Hit nearestIntersection(Ray r, FlatSurfaceRecord[] flatSurfaces, real maxDistance = infinity) {
    Hit bestHit = new Hit;
    bestHit.valid = false;
    bestHit.t = infinity;
    for (int i = 0; i < flatSurfaces.length; ++i) {
        FlatSurfaceRecord rec = flatSurfaces[i];
        // Fast Ray-Box Slab Test (With Distance Pruning)
        // If a hit was already found at distance bestHit.t, reject any surface further away!
        real currentReach = min(maxDistance, bestHit.t);
        if (!rayIntersectsAABB(r.p, r.d, rec.minBound, rec.maxBound, currentReach)) {
            continue; // Skipped instantly without touching Bezier math!
        }
        // Topological origin check (self-intersection prevention)
        bool ignoreOrigin = false;
        for (int idx : r.ignoredFlatIndices) {
            if (i == idx) {
                ignoreOrigin = true;
                break;
            }
        }
        // Exact Bézier Spline Intersection (Only evaluated when ray actually hits the box)
        Hit h = rec.s.intersect(r, maxDistance, ignoreOrigin);
        if (h.valid && h.t < bestHit.t) {
            bestHit = h;
            bestHit.flatIndex = i;
        }
    }
    return bestHit;
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
        if (vp.spatialIndex || vp.phaseGradient != null) return true;
    }
    return false;
}

private
PropagationResult propagate(Medium medium, Ray r, FlatSurfaceRecord[] flatSurfaces, real far, real ds) {
    if (requiresNumericalIntegration(medium)) {
        real tTotal = 0;
        Ray currRay = copy(r);
        guide traj = currRay.p;
        Integrator solver = (config.physics.grinSolver == GRINSolverType.RK4) ? Integrator.RK4() : Integrator.RK45();
        real max_ds = 0.075, min_ds = 1e-5;
        real current_ds = max(min_ds, min(ds, max_ds));
        real tol = config.physics.epsIntegration;
        Ray testRay = new Ray;
        while (tTotal < far && currRay.alive) {
            real nCenter = currRay.medium.n(currRay.p, currRay.d, currRay.lambda, currRay.pol);
            if (nCenter <= config.physics.minRefractiveIndex) {
                currRay.alive = false;
                break;
            }
            real step = min(current_ds, far - tTotal);
            GRINStep stepResult = solver.advance(currRay, step);
            if (solver.isAdaptive) {
                if (stepResult.error > tol && step > min_ds) {
                    current_ds *= 0.5;
                    continue;
                }
                if (stepResult.error > 0) {
                    real scaleFactor = 0.9 * (tol / stepResult.error) ^ 0.2;
                    current_ds = min(max_ds, current_ds * max(0.2, min(2.0, scaleFactor)));
                } else {
                    current_ds = min(max_ds, current_ds * 2.0);
                }
            }
            Ray nextRay = stepResult.ray;
            applyAbsorption(nextRay, step);
            if (nextRay.intensity() <= config.physics.minIntensity) {
                nextRay.alive = false;
                traj = traj -- nextRay.p;
                return propagationResult(nextRay, false, null, traj);
            }
            pair chord = nextRay.p - currRay.p;
            real chordLength = length(chord);
            Hit h;
            bool hitValid = false;
            if (chordLength > config.physics.epsIntegration) {
                testRay.p = currRay.p;
                testRay.d = chord / chordLength;
                testRay.ignoredFlatIndices = copy(currRay.ignoredFlatIndices);
                h = nearestIntersection(testRay, flatSurfaces, chordLength);
                hitValid = h.valid;
            }
            if (hitValid) {
                real hitStep = step * (h.t / chordLength);
                Ray incidentRay = solver.advance(currRay, min(step, hitStep)).ray;
                incidentRay.opl = currRay.opl + hitStep * nCenter;
                applyAbsorption(incidentRay, hitStep);
                traj = traj -- h.p;
                return propagationResult(incidentRay, true, h, traj);
            } else {
                pair midPoint = 0.5 * (currRay.p + nextRay.p);
                pair midDirection = unitOr(currRay.d + nextRay.d, currRay.d);
                real nMid = currRay.medium.n(midPoint, midDirection, currRay.lambda, currRay.pol);
                currRay = nextRay;
                currRay.opl += step * nMid;
                tTotal += step;
                traj = traj -- currRay.p;
                currRay.ignoredFlatIndices = new int[];
            }
        }
        return propagationResult(currRay, false, null, traj);
    } else {
        Hit h = nearestIntersection(r, flatSurfaces, far);
        real current_n = medium.n(r.p, r.d, r.lambda, r.pol);
        Ray nextRay = copy(r);
        if (!h.valid || h.t > far) {
            nextRay.opl += far * current_n;
            nextRay.p = r.p + far * r.d;
            applyAbsorption(nextRay, far);
            return propagationResult(nextRay, false, null, r.p -- nextRay.p);
        } else {
            nextRay.opl += h.t * current_n;
            nextRay.p = h.p;
            applyAbsorption(nextRay, h.t);
            return propagationResult(nextRay, true, h, r.p -- h.p);
        }
    }
}

// The "Navigator" Trace Facade & Ray Tree
struct Trace {
    Ray[] rayTree;
    int rootId = 0;
    int[] route;
    int[] rootIds;
    bool isReflection = false;
    int[] activeBranchRayIds() {
        int[] branch;
        if (rayTree.length == 0) return branch;
        int currId = rootId;
        if (currId < 0 || currId >= rayTree.length) return branch;
        int routeStep = 0;
        while (true) {
            branch.push(currId);
            Ray r = rayTree[currId];
            if (r.childrenIds.length == 0) break;
            int childChoice = r.defaultChildIndex;
            if (routeStep < route.length) {
                childChoice = route[routeStep];
                ++routeStep;
            }
            if (childChoice < 0 || childChoice >= r.childrenIds.length) {
                childChoice = min(max(childChoice, 0), r.childrenIds.length - 1);
            }
            currId = r.childrenIds[childChoice];
        }
        return branch;
    }
    Segment[] segments() {
        Segment[] segs;
        int[] branch = activeBranchRayIds();
        for (int k = 0; k < branch.length; ++k) {
            Ray r = rayTree[branch[k]];
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
                seg.after = rayTree[branch[k + 1]];
            } else {
                Ray afterRay = copy(r);
                afterRay.p = seg.b;
                if (r.terminated) afterRay.alive = false;
                seg.after = afterRay;
            }
            segs.push(seg);
        }
        return segs;
    }
    Interaction[] interactions() {
        Interaction[] ixs;
        int[] branch = activeBranchRayIds();
        for (int k = 0; k < branch.length; ++k) {
            Ray r = rayTree[branch[k]];
            if (r.hit != null && r.hit.valid) {
                Interaction ix = new Interaction;
                ix.hit = r.hit;
                ix.incident = copy(r);
                ix.incident.p = r.hit.p;
                if (length(r.trajectory) > 0) {
                    ix.incident.d = dir(r.trajectory, length(r.trajectory));
                }
                real parentI = r.intensity();
                pair N = ix.normal();
                for (int cId : r.childrenIds) {
                    Ray cRay = rayTree[cId];
                    if (dot(cRay.d, N) > 0) {
                        ix.reflected.push(cRay);
                        if (parentI > config.math.eps) ix.R += cRay.intensity() / parentI;
                    } else {
                        ix.transmitted.push(cRay);
                        if (parentI > config.math.eps) ix.T += cRay.intensity() / parentI;
                    }
                }
                ixs.push(ix);
            }
        }
        return ixs;
    }
    int length() { return interactions().length; }
    Interaction interaction(int ixIndex) {
        Interaction[] ixs = interactions();
        if (ixs.length == 0) return null;
        int idx = (ixIndex < 0) ? ixIndex + ixs.length : ixIndex;
        if (idx < 0 || idx >= ixs.length) return null;
        return ixs[idx];
    }
    Ray incident(int ixIndex) {
        Interaction ix = interaction(ixIndex);
        return (ix != null) ? ix.incident : null;
    }
    Ray transmitted(int ixIndex, int rayIndex = 0) {
        Interaction ix = interaction(ixIndex);
        if (ix == null || ix.transmitted.length == 0) return null;
        int idx = (rayIndex < 0) ? rayIndex + ix.transmitted.length : rayIndex;
        if (idx < 0 || idx >= ix.transmitted.length) return null;
        return ix.transmitted[idx];
    }
    Ray[] Transmitted(int ixIndex) {
        Interaction ix = interaction(ixIndex);
        if (ix == null) return new Ray[];
        return ix.transmitted;
    }
    Ray reflected(int ixIndex, int rayIndex = 0) {
        Interaction ix = interaction(ixIndex);
        if (ix == null || ix.reflected.length == 0) return null;
        int idx = (rayIndex < 0) ? rayIndex + ix.reflected.length : rayIndex;
        if (idx < 0 || idx >= ix.reflected.length) return null;
        return ix.reflected[idx];
    }
    Ray[] Reflected(int ixIndex) {
        Interaction ix = interaction(ixIndex);
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
        Ray r = copy(rayTree[branch[branch.length - 1]]);
        if (r.hit != null && r.hit.valid) {
            r.p = r.hit.p;
            r.pol = resolvePolarizationAt(r, r.hit.t);
            // If the ray terminated at the collision (e.g. Absorber / Aperture screen):
            if (r.terminated || !r.alive) {
                r.alive = false;
                r.pol.scaleTo(0.0); // Post-absorption intensity is 0!
            }
        } else if (length(r.trajectory) > 0) {
            real endT = length(r.trajectory);
            real totalL = arclength(r.trajectory);
            r.p = point(r.trajectory, endT);
            r.d = dir(r.trajectory, endT);
            r.waveNormal = r.medium.waveNormalFromGroup(r.d, r.p, r.lambda);
            r.pol = resolvePolarizationAt(r, totalL);
            if (r.terminated || !r.alive) {
                r.alive = false;
                r.pol.scaleTo(0.0);
            }
        }
        return r;
    }
    Ray first() {
        if (rootId >= 0 && rootId < rayTree.length) {
            return copy(rayTree[rootId]);
        }
        return null;
    }
    real opl() {
        Ray r = last();
        return (r != null) ? r.opl : 0.0;
    }
    Trace[] branches() {
        Trace[] b;
        if (rootIds.length > 1) {
            for (int rId : rootIds) {
                Trace sub = new Trace;
                sub.rayTree = rayTree;
                sub.rootId = rId;
                sub.route = copy(route);
                sub.rootIds = new int[];
                b.push(sub);
            }
        } else {
            int[] currentBranch = activeBranchRayIds();
            int routeStep = 0;
            for (int k = 0; k < currentBranch.length; ++k) {
                Ray r = rayTree[currentBranch[k]];
                if (r.childrenIds.length == 0) continue;
                int followedIndex = r.defaultChildIndex;
                bool wasExplicitStep = (routeStep < route.length);
                if (wasExplicitStep) {
                    followedIndex = route[routeStep];
                    ++routeStep;
                    continue; // Skip spawning side-branches for explicit follow() steps
                }
                for (int i = 0; i < r.childrenIds.length; ++i) {
                    if (i == followedIndex) continue;
                    Trace sub = new Trace;
                    sub.rayTree = rayTree;
                    sub.rootId = r.childrenIds[i];
                    sub.route = new int[];
                    Ray child = rayTree[r.childrenIds[i]];
                    if (r.hit != null && r.hit.valid) {
                        pair N = r.hit.normal;
                        sub.isReflection = (dot(child.d, N) * dot(r.d, N) < 0);
                    } else {
                        sub.isReflection = false;
                    }
                    b.push(sub);
                }
            }
        }
        return b;
    }
    Trace follow( ... int[] newRoute) {
        Trace out = new Trace;
        out.rayTree = rayTree;
        out.rootId = rootId;
        out.route = copy(newRoute);
        out.rootIds = copy(rootIds);
        out.isReflection = isReflection;
        return out;
    }
    Trace subtrace(int i) {
        if (rootIds.length == 0) return this;
        int idx = (i < 0) ? i + rootIds.length : i;
        if (idx < 0 || idx >= rootIds.length) return null;
        Trace sub = new Trace;
        sub.rayTree = rayTree;
        sub.rootId = rootIds[idx];
        sub.route = copy(route);
        sub.rootIds = new int[];
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
    out.rayTree = src.rayTree;
    out.rootId = src.rootId;
    out.route = copy(src.route);
    out.rootIds = copy(src.rootIds);
    out.isReflection = src.isReflection;
    return out;
}

// The Tracer Engine (Queue-Based BFS)
private
struct Tracer {
    FlatSurfaceRecord[] flatSurfaces;
    mediumResolver resolve;
    int maxDepth = config.physics.maxDepth;
    real minIntensity = config.physics.minIntensity;
    int branchBudget = config.physics.maxBranches;
    real ds = config.physics.ds;
    real far = config.physics.far;
    Trace traceRays(Ray[] rays) {
        Trace out = new Trace;
        if (rays.length == 0) return out;
        Ray[] combinedRayTree;
        int[] rootIds;
        for (int rIdx = 0; rIdx < rays.length; ++rIdx) {
            int startTreeOffset = combinedRayTree.length;
            Ray root = copy(rays[rIdx]);
            root.initTree(startTreeOffset, -1, 0);
            if (root.medium == null) {
                root.medium = resolve(root.p, root, null, null, true);
            }
            combinedRayTree.push(root);
            rootIds.push(startTreeOffset);
            int[] q = {startTreeOffset};
            int head = 0;
            int branchesCount = 1;
            while (head < q.length) {
                int currId = q[head];
                ++head;
                Ray r = combinedRayTree[currId];
                if (!r.alive || r.terminated || r.intensity() <= minIntensity || r.generation >= maxDepth) continue;
                // 5 parameters (no sceneDomain)
                // 1. Propagate ray to boundary or far reach
                PropagationResult propRes = propagate(r.medium, r, flatSurfaces, far, ds);
                r.trajectory = propRes.trajectory;
                r.opl = propRes.ray.opl;

                // 2. Handle Surface Interaction
                if (propRes.hitSurface) {
                    r.hit = propRes.hit;
                    Surface hitSurf = flatSurfaces[propRes.hit.flatIndex].s;

                    // Execute surface physics and receive InteractionResult:
                    InteractionResult ixRes = hitSurf.interact(propRes.ray, propRes.hit, resolve);
                    r.defaultChildIndex = ixRes.defaultIndex;

                    bool acceptedChild = false;
                    if (!ixRes.absorbed && ixRes.rays.length > 0) {
                        for (int i = 0; i < ixRes.rays.length; ++i) {
                            Ray child = ixRes.rays[i];
                            if (child.intensity() > minIntensity && r.generation + 1 <= maxDepth &&
                                (branchBudget <= 0 || branchesCount < branchBudget)) {
                                acceptedChild = true;
                                child.initTree(combinedRayTree.length, currId, r.generation + 1);
                                r.childrenIds.push(child.id);
                                combinedRayTree.push(child);
                                q.push(child.id);
                                ++branchesCount;
                            }
                        }
                    }
                    r.terminated = !acceptedChild;
                }
                combinedRayTree[currId] = r;
            }
        }
        out.rayTree = combinedRayTree;
        out.rootId = (rays.length == 1) ? rootIds[0] : -1;
        out.route = new int[];
        out.rootIds = rootIds;
        return out;
    }
}

private Tracer
makeTracer(FlatSurfaceRecord[] flatSurfaces, mediumResolver resolve, real ds, real far, int maxDepth) {
    Tracer tracer = new Tracer;
    tracer.flatSurfaces = flatSurfaces;
    tracer.resolve = resolve;
    tracer.ds = ds;
    tracer.far = far;
    tracer.maxDepth = maxDepth;
    tracer.minIntensity = config.physics.minIntensity;
    tracer.branchBudget = config.physics.maxBranches;
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
                int maxDepth = config.physics.maxDepth) {
        FlatSurfaceRecord[] flatSurfaces = prepareFlatSurfaces(elements);
        Element[] traceElements = elements;
        Medium traceAmbient = ambient;
        mediumResolver resolve = new Medium(pair p, Ray r, Medium m1, Medium m2, bool fromM1) {
            return resolveMediumAt(traceElements, traceAmbient, p, r, m1, m2, fromM1);
        };
        return makeTracer(flatSurfaces, resolve, ds, far, maxDepth).traceRays(rays);
    }
    Trace
    trace(Ray ray, real ds = config.physics.ds, real far = config.physics.far, int maxDepth = config.physics.maxDepth) {
        return trace(new Ray[]{ray}, ds = ds, far = far, maxDepth = maxDepth);
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
    bool3 showS = default;
    bool3 showP = default;
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
    out.showS = src.showS;
    out.showP = src.showP;
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
//   scale / showS / showP / p
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
                    bool3 keyword showS = default,
                    bool3 keyword showP = default,
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
    out.showS = showS;
    out.showP = showP;
    out.p = p;
    return out;
}

struct ResolvedArrow {
    bool visible = false;
    arrowhead head = config.rendering.arrowHead;
    real size = config.rendering.arrowSize;
    real relative = nan;
}

arrowbar
operator cast(ResolvedArrow a) {
    if (!a.visible) return None;
    real pos = isnan(a.relative) ? 0.5 : a.relative;
    arrowhead h = (a.head != null) ? a.head : config.rendering.arrowHead;
    real sz = isnan(a.size) ? config.rendering.arrowSize : a.size;
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
    arrowhead head = config.rendering.arrowHead;
    real size = config.rendering.arrowSize;
    // Global arrow position
    real relative = nan;
    // Segment rules
    SegmentRule[] rules;
    // Resolve arrow configuration for one segment
    ResolvedArrow getArrow(int i, int totalSegments) {
        ResolvedArrow out = new ResolvedArrow;
        bool isFirst = (i == 0);
        bool isLast = (i == totalSegments - 1);
        // Start with global visibility + positional policy.
        out.visible = this.visible && (isLast ? this.last : (isFirst ? this.first : this.middle));
        // Global appearance.
        out.head = this.head;
        out.size = this.size;
        // Global position.
        // If no global relative position is specified:
        //   last segment -> end
        //   other        -> middle
        out.relative = isnan(this.relative) ? (isLast ? 1.0 : 0.5) : this.relative;
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
                int resolvedIdx = (idx < 0) ? idx + totalSegments : idx;
                if (resolvedIdx == i) {
                    // Segment positional policy.
                    if (isFirst && rule.first != default) {
                        out.visible = (bool)rule.first;
                    }
                    if (!isFirst && !isLast && rule.middle != default) {
                        out.visible = (bool)rule.middle;
                    }
                    if (isLast && rule.last != default) {
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
                    arrowhead head = config.rendering.arrowHead,
                    real size = config.rendering.arrowSize,
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
        bool anyExplicit = (first == true || middle == true || last == true || rules.length > 0);
        out.visible = anyExplicit;
        out.first = (first == true);
        out.middle = (middle == true);
        out.last = (last == true);
        for (SegmentRule r : rules) {
            SegmentRule ruleCopy = copy(r);
            if (ruleCopy.visible == default) {
                ruleCopy.visible = true;
            }
            out.rules.push(ruleCopy);
        }
    }
    return out;
}

struct RayPolarization {
    // Global visibility
    bool visible = config.rendering.drawPolarization;
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
    bool showS = true;
    bool showP = true;
    pen p = nullpen;
    real minAmplitude = 0;
    real minSegmentLength = config.rendering.minPolarizationSegmentLength;
    real minSeparation = config.rendering.minPolarizationMarkerSeparation;
    // Polarization arrow appearance.
    // The arrowbar itself is resolved only when needed.
    arrowhead head = config.rendering.arrowPolHead;
    real size = config.rendering.arrowPolSize;
    // Segment rules
    SegmentRule[] rules;
    // Resolve polarization settings for one segment
    RayPolarization getPolarization(int i, int totalSegments) {
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
        out.showS = this.showS;
        out.showP = this.showP;
        out.p = this.p;
        out.minAmplitude = this.minAmplitude;
        out.head = this.head;
        out.size = this.size;
        out.minSegmentLength = this.minSegmentLength;
        out.minSeparation = this.minSeparation;
        // Apply rules only to this specific segment
        for (SegmentRule rule : rules) {
            for (int idx : rule.indices) {
                int resolvedIdx = (idx < 0) ? idx + totalSegments : idx;
                if (resolvedIdx == i) {
                    if (rule.visible != default) out.visible = (bool)rule.visible;
                    if (rule.count >= 0) out.count = rule.count;
                    if (!isnan(rule.spacing)) out.spacing = rule.spacing;
                    if (!isnan(rule.before)) out.before = rule.before;
                    if (!isnan(rule.after)) out.after = rule.after;
                    if (rule.first != default) out.first = (bool)rule.first;
                    if (rule.middle != default) out.middle = (bool)rule.middle;
                    if (rule.last != default) out.last = (bool)rule.last;
                    if (!isnan(rule.scale)) out.scale = rule.scale;
                    if (rule.showS != default) out.showS = (bool)rule.showS;
                    if (rule.showP != default) out.showP = (bool)rule.showP;
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
    out.showS = src.showS;
    out.showP = src.showP;
    out.p = src.p;
    out.minAmplitude = src.minAmplitude;
    out.head = src.head;
    out.size = src.size;
    out.minSegmentLength = src.minSegmentLength;
    out.minSeparation = src.minSeparation;
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
RayPolarization rayPolarization(bool visible = config.rendering.drawPolarization,
                                int count = -1,
                                real spacing = 0.5,
                                real before = 0.2,
                                real after = 0.2,
                                bool3 first = default,
                                bool3 middle = default,
                                bool3 last = default,
                                real scale = 0.5,
                                bool showS = true,
                                bool showP = true,
                                pen p = nullpen,
                                arrowhead head = config.rendering.arrowPolHead,
                                real size = config.rendering.arrowPolSize, ... SegmentRule[] rules) {
    RayPolarization out = new RayPolarization;
    out.count = count;
    out.spacing = spacing;
    out.before = before;
    out.after = after;
    out.scale = scale;
    out.showS = showS;
    out.showP = showP;
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
        bool anyExplicit = (first == true || middle == true || last == true || rules.length > 0);
        out.visible = anyExplicit;
        out.first = (first == true);
        out.middle = (middle == true);
        out.last = (last == true);
        for (SegmentRule r : rules) {
            SegmentRule ruleCopy = copy(r);
            if (ruleCopy.visible == default) ruleCopy.visible = true;
            out.rules.push(ruleCopy);
        }
    }
    return out;
}

private
path polarizationMarkerLine(pair pos, pair d, real Ap, RayPolarization pol, real scale = 1) {
    pair vPerp = perp(d) * (Ap * pol.scale * scale);
    return pos - vPerp -- pos + vPerp;
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
    if (pol.showP && Ap > pol.minAmplitude) {
        path line = polarizationMarkerLine(pos, d, Ap, pol, scale);
        real lineLen = arclength(line);
        // Resolve effective arrow size and head safely
        real arrowSz = (pol.size > 0 && !isnan(pol.size)) ? pol.size : config.rendering.arrowPolSize;
        arrowhead h = (pol.head != null) ? pol.head : config.rendering.arrowPolHead;
        // Condition: Line must fit at least 2 arrowheads (one on each side)
        //    Otherwise draw a clean plain line (None) to avoid overlapping arrow blobs
        real minLenNeeded = max(config.rendering.minPolarizationMarkerWithArrow, 2.0 * arrowSz);
        arrowbar ab = (lineLen > minLenNeeded) ? Arrows(h, arrowSz) : None;
        draw(f, line, q, ab);
    }
    if (pol.showS && As > pol.minAmplitude) {
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
        real userToPt = length(t * (1, 0)) * scale;
        pair pf = t * pos;
        pair df = t * (pos + d) - pf;
        drawPolarizationMarker(f, pf, df, As, Ap, userToPt, p, pol);
    });
}

private
void drawPolarization(picture pic = currentpicture,
                      Ray r,
                      path traj,
                      int segmentIdx = 0,
                      int totalSegments = 1,
                      pen p = currentpen,
                      RayPolarization pol = rayPolarization(),
                      bool applyOpacity = true) {
    if (r == null || !r.alive) return;
    RayPolarization effPol = pol.getPolarization(segmentIdx, totalSegments);
    if (!effPol.visible) return;

    real L = arclength(traj);
    if (L <= config.rendering.minDrawableLength) return;

    // Minimum Segment Length Gate: If the segment is too short to reasonably display a marker, skip entirely.
    real minSegL =
        !isnan(effPol.minSegmentLength) ? effPol.minSegmentLength : config.rendering.minPolarizationSegmentLength;
    if (L < minSegL) return;

    // Margins and Spacing Limits
    real s_start = clamp(effPol.before, 0, L);
    real s_end = clamp(L - effPol.after, s_start, L);
    real range = s_end - s_start;
    real minSep =
        !isnan(effPol.minSeparation) ? effPol.minSeparation : config.rendering.minPolarizationMarkerSeparation;
    if (minSep <= 0) minSep = 0.1;

    real[] rawPositions;

    // Position Calculation with Overlap Prevention
    if (range < minSep) {
        // Range is too tight for multiple markers: place exactly ONE centered marker
        if (effPol.middle || effPol.first || effPol.last) {
            rawPositions.push(L / 2.0);
        }
    } else if (effPol.count > 0) {
        if (effPol.count == 1) {
            rawPositions.push(s_start + range / 2.0);
        } else {
            // Clamp count to prevent stacking
            int maxAllowed = max(1, floor(range / minSep) + 1);
            int effCount = min(effPol.count, maxAllowed);
            if (effCount == 1) {
                rawPositions.push(s_start + range / 2.0);
            } else {
                real step = range / (effCount - 1);
                for (int i = 0; i < effCount; ++i) {
                    if (i == 0 && !effPol.first) continue;
                    if (i == effCount - 1 && !effPol.last) continue;
                    if (i > 0 && i < effCount - 1 && !effPol.middle) continue;
                    rawPositions.push(s_start + i * step);
                }
            }
        }
    } else {
        // Spacing-based Placement
        real sp = max(minSep, effPol.spacing);
        int numIntervals = floor(range / sp);
        if (numIntervals <= 0) {
            if (effPol.middle || effPol.first || effPol.last) {
                rawPositions.push(s_start + range / 2.0);
            }
        } else {
            real actualStart = s_start + (range - numIntervals * sp) / 2.0;
            for (int i = 0; i <= numIntervals; ++i) {
                if (i == 0 && !effPol.first) continue;
                if (i == numIntervals && !effPol.last) continue;
                if (i > 0 && i < numIntervals && !effPol.middle) continue;
                rawPositions.push(actualStart + i * sp);
            }
        }
    }

    // Strict Deduplication Filter (Guarantees no two markers are closer than minSep)
    real[] positions;
    real last_s = -infinity;
    for (real s : rawPositions) {
        if (s - last_s >= minSep - config.math.eps) {
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
        Polarization localPol = resolvePolarizationAt(r, s);

        // Derive S and P field amplitudes directly from Stokes parameters S0 and S1
        real localAs = sqrt(max(0.0, 0.5 * (localPol.S0 + localPol.S1)));
        real localAp = sqrt(max(0.0, 0.5 * (localPol.S0 - localPol.S1)));
        real localI = localPol.intensity();

        // Skip marker if local intensity has decayed below threshold
        if (localI < config.physics.minIntensity) continue;
        if (localAs <= effPol.minAmplitude && localAp <= effPol.minAmplitude) continue;

        pen localPen = p;
        if (applyOpacity && config.rendering.useIntensityOpacity) {
            real baseAlpha = opacity(p);
            if (baseAlpha == 0) baseAlpha = 1.0;
            localPen = p + opacity(baseAlpha * max(0.0, min(1.0, localI)));
        }

        drawPolarizationMarker(pic, pos, d, localAs, localAp, 1, localPen, effPol);
    }
}

// ============================================================================
// Huygens-Fresnel Wavefront System
// ============================================================================
struct RayWavefront {
    // Global visibility
    bool visible = config.rendering.drawWavefronts;
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
    bool wavelets = false;     // If true: draw circular Huygens secondary wavelet arcs
    real alpha = 30.0;         // Half-angle of arc (spans d +- alpha, default 30 deg, up to 90 deg)
    bool cutAtBoundary = true; // If true: arcs are clipped at physical surfaces
    bool dashedBehind = false; // If true: portion behind boundary is dashed instead of omitted
    bool showSource = true;    // Draw center dot at emission origin P0
    // Bundle iso-OPL contours (draw(Trace[], ...) only)
    bool contours = false;     // Draw iso-OPL / iso-phase wavefront contours across the bundle
    pen pContours = lightblue; // Pen for iso-OPL contour lines
    // Interference lines (draw(Trace[], ...) only, requires contours or wavelets)
    bool antinodes = false;                           // Draw constructive interference channels (Delta phi = 2*pi*m)
    bool nodes = false;                               // Draw destructive interference channels (Delta phi = pi*(2m+1))
    pen pLine = gray(0.4) + dashed;                   // Pen for interference lines
    pen pDots = gray(0.4) + 2pt;                      // Pen for intersection dots
    real maxDist = config.rendering.wavefrontMaxDist; // Max distance from source midpoint to trace interference
    int maxOrders = config.rendering.wavefrontMaxOrders; // Maximum interference order |m| to draw
    bool showDots = true;                                // Show circle-crossing dots at antinodes
    bool showLabels = true;                              // Label each antinodal channel with m
    // Segment rules
    SegmentRule[] rules;

    // Resolve wavefront configuration for one segment
    RayWavefront getWavefront(int i, int totalSegments) {
        bool isFirst = (i == 0);
        bool isLast = (i == totalSegments - 1);
        RayWavefront out = new RayWavefront;
        out.visible = this.visible && (isLast ? this.last : (isFirst ? this.first : this.middle));
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
        out.cutAtBoundary = this.cutAtBoundary;
        out.dashedBehind = this.dashedBehind;
        out.showSource = this.showSource;
        out.p = this.p;
        // Bundle-level fields pass through unchanged
        out.contours = this.contours;
        out.pContours = this.pContours;
        out.antinodes = this.antinodes;
        out.nodes = this.nodes;
        out.pLine = this.pLine;
        out.pDots = this.pDots;
        out.maxDist = this.maxDist;
        out.maxOrders = this.maxOrders;
        out.showDots = this.showDots;
        out.showLabels = this.showLabels;

        for (SegmentRule rule : rules) {
            for (int idx : rule.indices) {
                int resolvedIdx = (idx < 0) ? idx + totalSegments : idx;
                if (resolvedIdx == i) {
                    // Segment positional policy
                    if (isFirst && rule.first != default) out.visible = (bool)rule.first;
                    if (!isFirst && !isLast && rule.middle != default) out.visible = (bool)rule.middle;
                    if (isLast && rule.last != default) out.visible = (bool)rule.last;
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
    out.cutAtBoundary = src.cutAtBoundary;
    out.dashedBehind = src.dashedBehind;
    out.showSource = src.showSource;
    out.contours = src.contours;
    out.pContours = src.pContours;
    out.antinodes = src.antinodes;
    out.nodes = src.nodes;
    out.pLine = src.pLine;
    out.pDots = src.pDots;
    out.maxDist = src.maxDist;
    out.maxOrders = src.maxOrders;
    out.showDots = src.showDots;
    out.showLabels = src.showLabels;
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
RayWavefront rayWavefront(bool visible = config.rendering.drawWavefronts,
                          real spacing = 12.0,
                          int count = -1,
                          real width = 14.0,
                          real before = 0.05,
                          real after = 0.05,
                          bool wavelets = false,
                          real alpha = 30.0,
                          bool cutAtBoundary = true,
                          bool dashedBehind = false,
                          bool showSource = true,
                          pen p = nullpen,
                          // Bundle iso-OPL contours
                          bool contours = false,
                          pen pContours = lightblue,
                          // Interference lines
                          bool antinodes = false,
                          bool nodes = false,
                          pen pLine = gray(0.4) + dashed,
                          pen pDots = gray(0.4) + 2.5pt,
                          real maxDist = config.rendering.wavefrontMaxDist,
                          int maxOrders = config.rendering.wavefrontMaxOrders,
                          bool showDots = true,
                          bool showLabels = true,
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
    out.cutAtBoundary = cutAtBoundary;
    out.dashedBehind = dashedBehind;
    out.showSource = showSource;
    out.p = p;
    out.contours = contours;
    out.pContours = pContours;
    out.antinodes = antinodes;
    out.nodes = nodes;
    out.pLine = pLine;
    out.pDots = pDots;
    out.maxDist = maxDist;
    out.maxOrders = maxOrders;
    out.showDots = showDots;
    out.showLabels = showLabels;
    if (visible) {
        out.visible = true;
        out.first = (first != default) ? (bool)first : true;
        out.middle = (middle != default) ? (bool)middle : true;
        out.last = (last != default) ? (bool)last : true;
        for (SegmentRule r : rules) out.rules.push(r);
    } else {
        // visible = false: default all segment positions to hidden unless
        // explicitly restored via first/middle/last = true or segment() rules.
        bool anyExplicit = (first == true || middle == true || last == true || rules.length > 0);
        out.visible = anyExplicit;
        out.first = (first == true);
        out.middle = (middle == true);
        out.last = (last == true);
        for (SegmentRule r : rules) {
            SegmentRule ruleCopy = copy(r);
            if (ruleCopy.visible == default) ruleCopy.visible = true;
            out.rules.push(ruleCopy);
        }
    }
    return out;
}

// Standalone Wavefront Marker / Huygens Wavelet Renderer
void drawWavefrontMarker(picture pic = currentpicture,
                         pair center,       // Emission origin P0 (for wavelets) or bar center (for crests)
                         pair d,            // Wave propagation direction
                         real radius,       // Arc radius (for wavelets) or arc distance s
                         real localN = 1.0, // Refractive index for wavelength scaling
                         Hit hit = null,    // Optional boundary for automatic clipping
                         pen p = currentpen,
                         RayWavefront wf = rayWavefront()) {
    pen q = (wf.p != nullpen) ? wf.p : p;
    pair dHat = unitOr(d);

    if (wf.wavelets) {
        // =====================================================================
        // Mode B: Circular Huygens Secondary Wavelet Arc
        // =====================================================================
        real theta0 = degrees(atan2(dHat.y, dHat.x));
        real aSpan = clamp(wf.alpha, 5.0, 90.0);
        path arcFull = Arc(center, radius, theta0 - aSpan, theta0 + aSpan, 100);

        // Boundary clipping: truncate the arc if it hits a physical surface
        if (wf.cutAtBoundary && hit != null && hit.valid) {
            pair hitP = hit.p;
            pair N = unitOr(hit.normal);
            if (dot(dHat, N) > 0) N = -N; // Opposes incident wave

            guide arcVisible, arcBehind;
            int nPts = 100;
            for (int j = 0; j <= nPts; ++j) {
                pair pt = point(arcFull, reltime(arcFull, j / (real)nPts));
                if (dot(pt - hitP, N) >= -config.math.eps)
                    arcVisible = arcVisible -- pt;
                else
                    arcBehind = arcBehind -- pt;
            }
            if (length(arcVisible) > 0) draw(pic, arcVisible, q);
            if (wf.dashedBehind && length(arcBehind) > 0) draw(pic, arcBehind, q + dashed);
        } else {
            draw(pic, arcFull, q);
        }
    } else {
        // =====================================================================
        // Mode A: Planar / Orthogonal Wavefront Crest Bar
        // =====================================================================
        pair w = perp(dHat); // Transverse vector perpendicular to ray
        real halfW = 0.5 * wf.width;
        draw(pic, (center - halfW * w) -- (center + halfW * w), q + 1.0bp);
    }
}

// Picture wrapper for canvas coordinate mapping
private
void drawWavefront(picture pic = currentpicture,
                   Ray r,
                   path traj,
                   int segmentIdx = 0,
                   int totalSegments = 1,
                   pen p = currentpen,
                   RayWavefront wf = rayWavefront()) {
    if (r == null || !r.alive || !wf.visible) return;
    RayWavefront effWf = wf.getWavefront(segmentIdx, totalSegments);
    if (!effWf.visible) return;

    real L = arclength(traj);
    if (L <= config.rendering.minDrawableLength) return;

    real localN = (r.medium != null) ? r.medium.n(r.p, r.d, r.lambda, r.pol) : 1.0;
    if (localN <= config.math.eps) localN = 1.0;

    pen drawPen = (effWf.p != nullpen) ? effWf.p : p;
    pair P0 = point(traj, 0);
    pair d0 = unitOr(dir(traj, 0));

    // Calculate expanding radii / distance steps:
    real[] radii;
    if (effWf.count > 0) {
        real step = L / (effWf.count + 1);
        for (int i = 1; i <= effWf.count; ++i) radii.push(i * step);
    } else {
        real sp = (!isnan(effWf.spacing) && effWf.spacing > 0) ? (effWf.spacing / localN) : 12.0;
        for (real cur_r = sp; cur_r < L - config.math.eps; cur_r += sp) {
            radii.push(cur_r);
        }
    }

    // Draw source dot once if wavelets are enabled
    if (effWf.wavelets && effWf.showSource) dot(pic, P0, drawPen + 2.0pt);

    // Render each wavefront marker:
    for (real rad : radii) {
        if (effWf.wavelets) {
            // Wavelet mode: center is P0, radius expands outward
            drawWavefrontMarker(pic, P0, d0, rad, localN, r.hit, drawPen, effWf);
        } else {
            // Transverse bar mode: center is point along the trajectory
            real u = arctime(traj, rad);
            pair pos = point(traj, u);
            pair d_local = unitOr(dir(traj, u));
            drawWavefrontMarker(pic, pos, d_local, rad, localN, r.hit, drawPen, effWf);
        }
    }
}

// Helper packet for a point on a trace at a specific OPL
struct TracePoint {
    bool valid = false;
    pair p = (nan, nan);
    pair d = (1, 0);
    real intensity = 1.0;
    int segIdx = -1;
    real segOPLStart = 0.0;
    real segOPLEnd = 0.0;
}

// Robust OPL sampler along any ray trace
private TracePoint pointAtOPL(Trace t, real targetOPL) {
    TracePoint tp;
    Segment[] segs = t.segments();
    int nSegs = segs.length;
    if (nSegs == 0) return tp;

    real cur_opl = 0.0;
    for (int i = 0; i < nSegs; ++i) {
        Segment s = segs[i];
        real segLen = arclength(s.traj);
        real n_med = (s.before.medium != null) ? s.before.medium.n(s.before.p, s.before.d, s.before.lambda) : 1.0;
        real seg_opl = segLen * n_med;
        real next_opl = cur_opl + seg_opl;

        // Strict range check: if targetOPL exceeds next_opl on the last segment, the ray has ENDED!
        if (targetOPL >= cur_opl - config.math.eps && targetOPL <= next_opl + config.math.eps) {
            real distAlong = (n_med > config.math.eps) ? (targetOPL - cur_opl) / n_med : 0.0;
            distAlong = clamp(distAlong, 0.0, segLen);
            real u = (segLen > 0) ? arctime(s.traj, distAlong) : 0.0;

            tp.valid = true;
            tp.p = point(s.traj, u);
            tp.d = unitOr(dir(s.traj, u), s.before.d);
            tp.intensity = s.before.intensity();
            tp.segIdx = i;
            tp.segOPLStart = cur_opl;
            tp.segOPLEnd = next_opl;
            return tp;
        }
        cur_opl = next_opl;
    }
    return tp; // Ray has ended or not reached yet -> invalid!
}

// Private: Iso-OPL / Iso-Phase Wavefront Contours across a ray bundle
private
void drawBundleContours(picture pic = currentpicture, Trace[] traces, RayWavefront wf) {
    int N = traces.length;
    if (N < 2) return;
    pen p = (wf.pContours != nullpen) ? wf.pContours : heavycyan;
    real spacing = wf.spacing;
    int count = (wf.count > 0) ? wf.count : 12;

    // 1. Resolve Incoming Direction and Reference Origin
    Trace t0 = traces[0];
    if (t0.first() == null) return;
    pair p0_ref = t0.first().p;
    pair d_in = unitOr(t0.first().d);
    real n_amb = (t0.first().medium != null) ? t0.first().medium.n() : 1.0;

    // 2. Compute Global OPL Range Across All Rays
    real[] launchOffsets = new real[N];
    real minGlobalOPL = infinity;
    real maxGlobalOPL = -infinity;

    for (int k = 0; k < N; ++k) {
        Trace t = traces[k];
        if (t.first() == null || t.last() == null) continue;

        real offset = dot(t.first().p - p0_ref, d_in) * n_amb;
        launchOffsets[k] = offset;

        real traceTotalOPL = 0.0;
        for (Segment s : t.segments()) {
            real n_med = (s.before.medium != null) ? s.before.medium.n(s.before.p, s.before.d, s.before.lambda) : 1.0;
            traceTotalOPL += arclength(s.traj) * n_med;
        }

        minGlobalOPL = min(minGlobalOPL, offset);
        maxGlobalOPL = max(maxGlobalOPL, offset + traceTotalOPL);
    }

    if (minGlobalOPL >= maxGlobalOPL - config.math.eps) return;

    // 3. Determine Target OPL values
    real[] targetOPLs;
    if (!isnan(spacing) && spacing > 0) {
        for (real L = minGlobalOPL + spacing; L < maxGlobalOPL; L += spacing) {
            targetOPLs.push(L);
        }
    } else {
        real dL = (maxGlobalOPL - minGlobalOPL) / (count + 1);
        for (int i = 1; i <= count; ++i) {
            targetOPLs.push(minGlobalOPL + i * dL);
        }
    }

    // 4. Trace Wavefront Contours at Each Target OPL
    for (real L_target : targetOPLs) {
        guide currentWf;
        int currentRun = 0;
        TracePoint prevTP;
        int prevK = -1;

        for (int k = 0; k < N; ++k) {
            real targetAlongRay = L_target - launchOffsets[k];
            TracePoint tp = pointAtOPL(traces[k], targetAlongRay);

            if (!tp.valid) {
                // Ray k has ended (or hasn't started):
                // Finish and draw the current contiguous run if it spans >= 2 rays:
                if (currentRun >= 2) {
                    draw(pic, currentWf, p);
                }
                currentWf = nullpath;
                currentRun = 0;
                prevTP.valid = false;
                prevK = -1;
                continue;
            }

            // Ray k is active!
            if (prevTP.valid && k == prevK + 1) {
                // Contiguous adjacent rays both active:
                // Check if crossing between air and glass:
                if (tp.segIdx != prevTP.segIdx) {
                    Trace t_curr = traces[k];
                    Trace t_prev = traces[prevK];
                    int transSeg = min(tp.segIdx, prevTP.segIdx);

                    if (transSeg < t_curr.segments().length && transSeg < t_prev.segments().length) {
                        pair S_hit_curr = t_curr.segments()[transSeg].b;
                        pair S_hit_prev = t_prev.segments()[transSeg].b;

                        real opl_hit_curr = launchOffsets[k] + (tp.segIdx > transSeg ? tp.segOPLStart : tp.segOPLEnd);
                        real opl_hit_prev =
                            launchOffsets[prevK] + (prevTP.segIdx > transSeg ? prevTP.segOPLStart : prevTP.segOPLEnd);

                        real denom = opl_hit_curr - opl_hit_prev;
                        real tau =
                            (abs(denom) > config.math.eps) ? clamp((L_target - opl_hit_prev) / denom, 0.0, 1.0) : 0.5;
                        pair S_contact = S_hit_prev + tau * (S_hit_curr - S_hit_prev);

                        currentWf = currentWf -- S_contact; // Route through boundary contact point S
                    }
                }

                currentWf = currentWf -- tp.p;
                ++currentRun;
            } else {
                // Start a new contiguous wavefront line at tp.p
                currentWf = tp.p;
                currentRun = 1;
            }

            prevTP = tp;
            prevK = k;
        }

        // Draw pending run at the end of the rays:
        if (currentRun >= 2) {
            draw(pic, currentWf, p);
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
    pair midS = 0.5 * (S1 + S2);

    pair v_slit = S1 - S2;
    real d = length(v_slit);
    if (d <= config.math.eps) return;

    // Orthonormal coordinate system centred at midS
    pair u_slit = v_slit / d;
    pair u_axis = perp(u_slit);
    if (dot(u_axis, t1.first().d) < 0) u_axis = -u_axis;

    // Resolve wavelength from wf.spacing or ray lambda
    real lam = (!isnan(wf.spacing) && wf.spacing > 0) ? wf.spacing : (t1.first().lambda > 0 ? t1.first().lambda : 8.0);

    int maxOrderPossible = floor(d / lam);
    int ordersToDraw = min(wf.maxOrders, maxOrderPossible);

    // Collect half-integer offsets for nodes (0 = antinodes only, 0.5 = nodes only)
    real[] offsets;
    if (wf.antinodes) offsets.push(0.0);
    if (wf.nodes) offsets.push(0.5);

    for (real offset : offsets) {
        // Loop over each interference order m = 0, +-1, +-2, ...
        for (int m = -ordersToDraw; m <= ordersToDraw; ++m) {
            guide hyperbolaPath;
            pair lastPt = (nan, nan);

            int n_max = ceil((wf.maxDist + d) / lam);
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
                        if (x_local <= wf.maxDist) {
                            pair P_cross = midS + x_local * u_axis + y_local * u_slit;
                            hyperbolaPath = hyperbolaPath -- P_cross;
                            if (wf.showDots) dot(pic, P_cross, wf.pDots);
                            lastPt = P_cross;
                        }
                    }
                }
            }

            if (length(hyperbolaPath) > 0) {
                draw(pic, hyperbolaPath, wf.pLine);
                if (wf.showLabels && !isnan(lastPt.x)) {
                    string lbl;
                    if (offset == 0) {
                        lbl = (m == 0) ? "$m=0$" : format("$m=%+d$", m);
                    } else {
                        lbl = (m >= 0) ? format("$m=%+d.5$", m) : format("$m=%d.5$", m + 1);
                    }
                    label(pic, Label("\tiny " + lbl), lastPt, E, wf.pLine);
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
                    int segmentIdx = 0,
                    int totalSegments = 1,
                    pen p = currentpen,
                    ResolvedArrow arrow = new ResolvedArrow,
                    RayPolarization polarization = rayPolarization(),
                    RayWavefront wavefronts = rayWavefront(),
                    bool applyOpacity = true) {
    real L_traj = arclength(traj);
    if (L_traj < config.rendering.minDrawableLength) return;
    pen linePen = p;
    if (applyOpacity) {
        real startI = max(0.0, min(1.0, r.intensity()));
        real baseAlpha = opacity(p);
        if (baseAlpha == 0) baseAlpha = 1.0;
        linePen = p + opacity(baseAlpha * startI);
    }
    real arrowSz = isnan(arrow.size) ? config.rendering.arrowSize : arrow.size;
    arrowbar ab = (arrow.visible && L_traj >= config.rendering.minPolarizationSegmentLength) ? (arrowbar)arrow : None;

    drawWavefront(pic, r, traj, segmentIdx, totalSegments, p, wavefronts);
    draw(pic, L, traj, linePen, ab);
    drawPolarization(pic, r, traj, segmentIdx, totalSegments, p, polarization, applyOpacity);
}

private
void drawFadingSegment(picture pic = currentpicture,
                       Label L = Label(""),
                       Ray r,
                       path traj,
                       int segmentIdx = 0,
                       int totalSegments = 1,
                       pen p = currentpen,
                       ResolvedArrow arrow = new ResolvedArrow,
                       RayPolarization polarization = rayPolarization(),
                       RayWavefront wavefronts = rayWavefront(),
                       int steps = 16,
                       bool applyOpacity = true) {
    // Draw wavefronts (if requested)
    drawWavefront(pic, r, traj, segmentIdx, totalSegments, p, wavefronts);

    real L_traj = arclength(traj);
    real alpha = (r.medium != null) ? r.medium.absorption(r.p, r.d, r.lambda, r.pol) : 0.0;
    real startI = r.intensity();
    if (L_traj <= config.rendering.minDrawableLength || steps <= 1 || alpha <= 0 || startI <= config.math.eps) {
        drawRaySegment(pic, L, r, traj, segmentIdx, totalSegments, p, arrow, polarization, wavefronts, applyOpacity);
        return;
    }
    real baseAlpha = opacity(p);
    if (baseAlpha == 0) baseAlpha = 1.0;
    real stepL = L_traj / steps;
    // Draw exponentially fading sub-segments without arrows
    for (int i = 0; i < steps; ++i) {
        real t1 = arctime(traj, i * stepL);
        real t2 = arctime(traj, (i + 1) * stepL);
        path sub = subpath(traj, t1, t2);
        real localDist = (i + 0.5) * stepL;
        real localI = startI * exp(-alpha * localDist);
        pen subPen = p;
        if (applyOpacity) {
            subPen = p + opacity(baseAlpha * max(0.0, min(1.0, localI))) + squarecap;
        } else {
            subPen = p + squarecap;
        }
        draw(pic, (i == 0 ? L : Label("")), sub, subPen, None);
    }
    // Overlay a single full-sized arrowhead matching the local faded intensity
    if (arrow.visible) {
        arrowhead h = (arrow.head != null) ? arrow.head : config.rendering.arrowHead;
        real sz = isnan(arrow.size) ? config.rendering.arrowSize : arrow.size;
        real pos = isnan(arrow.relative) ? (segmentIdx == totalSegments - 1 ? 1.0 : 0.5) : arrow.relative;
        pen arrowPen = p;
        if (applyOpacity) {
            real localI = startI * exp(-alpha * pos * L_traj);
            arrowPen = p + opacity(baseAlpha * max(0.0, min(1.0, localI)));
        }
        pic.add(arrow(h, traj, invisible, sz, FillDraw(arrowPen, arrowPen), position = Relative(pos)));
    }
    // Draw polarization markers
    drawPolarization(pic, r, traj, segmentIdx, totalSegments, p, polarization, applyOpacity);
}

// Fast 2D Ray-Viewport Slab Clipper (Device Space)
private
real[] clipRayToViewport(pair p, pair d, pair m, pair M) {
    real Tin = -infinity, Tout = infinity;
    if (abs(d.x) < config.rendering.clipVelocityEps) {
        if (p.x < m.x || p.x > M.x) return new real[]{};
    } else {
        real invDx = 1.0 / d.x;
        real t1 = (m.x - p.x) * invDx;
        real t2 = (M.x - p.x) * invDx;
        Tin = max(Tin, min(t1, t2));
        Tout = min(Tout, max(t1, t2));
    }
    if (abs(d.y) < config.rendering.clipVelocityEps) {
        if (p.y < m.y || p.y > M.y) return new real[]{};
    } else {
        real invDy = 1.0 / d.y;
        real t1 = (m.y - p.y) * invDy;
        real t2 = (M.y - p.y) * invDy;
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
                    bool clipBackward,
                    bool clipForward) {
    pic.add(new void(frame f, transform t, transform T_trans, pair m, pair M) {
        pair m_eff = m + (margin, margin);
        pair M_eff = M - (margin, margin);
        pair pf = t * pt_world;
        pair df = t * (pt_world + dir_world) - pf;
        real[] clip = clipRayToViewport(pf, df, m_eff, M_eff);
        if (clip.length == 2) {
            real t0 = clipBackward ? clip[0] : 0.0;
            real t1 = clipForward ? clip[1] : 0.0;
            if (t0 < t1) {
                draw(f, (pf + t0 * df) -- (pf + t1 * df), p, (arrowbar)arrow);
            }
        }
    });
}

// Helper to modulate pen opacity by the ray's physical intensity (if enabled)
private
pen resolveRayPen(pen basePen, Ray r, bool applyOp) {
    if (!applyOp || r == null) return basePen;
    real baseAlpha = opacity(basePen);
    if (baseAlpha == 0) baseAlpha = 1.0;
    real I = clamp(r.intensity(), 0.0, 1.0);
    return basePen + opacity(baseAlpha * I);
}

// Complete Renderer Struct
// Default Surface Drawing Delegate
private
void drawSurfaces(picture pic = currentpicture, Element e) {
    for (Surface s : e.surfaces) {
        if (length(s.p) > 0) {
            transform totalT = e.T * s.T;
            if (s.fill != nullpen && cyclic(s.p)) {
                fill(pic, totalT * s.p, s.fill);
            }
            draw(pic, totalT * s.p, s.stroke);
        }
    }
}

struct Renderer {
    pen p = currentpen;
    RayArrows arrows = rayArrows();
    RayPolarization polarization = rayPolarization();
    RayWavefront wavefronts = rayWavefront();
    bool intensityOpacity = true;
    bool drawAbsorptionFading = config.rendering.drawAbsorptionFading;
    int fadeSubpathSteps = config.rendering.fadeSubpathSteps;
    real before = nan;
    real after = nan;
    real margin = 0;
    bool drawReflections = config.rendering.drawReflections;
    void drawTrace(picture pic = currentpicture, Trace T) {
        Segment[] segs = T.segments();
        Trace[] brs = T.branches();
        if (segs.length == 0 && brs.length == 0) return;
        for (int i = 0; i < segs.length; ++i) {
            Segment s = segs[i];
            path tr = s.traj;
            real L_tr = arclength(tr);
            bool isFirst = (i == 0);
            bool isLast = (i == segs.length - 1);
            ResolvedArrow currentArr = arrows.getArrow(i, segs.length);
            bool applyOp = config.rendering.useIntensityOpacity && intensityOpacity;
            pen effPen = resolveRayPen(p, s.before, applyOp);
            // Special Case: Single unhit segment passing entirely through space
            if (isFirst && isLast && (before == infinity || before == -infinity) && after == infinity) {
                pair p0 = (length(tr) > 0) ? point(tr, 0) : s.before.p;
                pair d0 = (length(tr) > 0) ? dir(tr, 0) : s.before.d;
                drawClippedRay(pic, p0, d0, this.margin, effPen, currentArr, true, true);
                continue;
            }
            // A. First Segment (Entrance Leg)
            if (isFirst) {
                if (before == 0) {
                    tr = nullpath;
                } else if (before == infinity || before == -infinity) {
                    pair b = point(tr, length(tr));
                    pair d_in = (length(tr) > 0) ? dir(tr, length(tr)) : s.before.d;
                    drawClippedRay(pic, b, d_in, this.margin, effPen, currentArr, true, false);
                    if (!isLast) continue;
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
            if (isLast && (s.after == null || s.after.alive)) {
                if (after == 0) {
                    tr = nullpath;
                } else if (after == infinity) {
                    pair a = point(tr, 0);
                    pair d_out = (length(tr) > 0) ? dir(tr, 0) : s.before.d;
                    drawClippedRay(pic, a, d_out, this.margin, effPen, currentArr, false, true);
                    continue;
                } else if (!isnan(after) && after > 0 && after < L_tr) {
                    tr = subpath(tr, 0, arctime(tr, after));
                }
            }
            if (arclength(tr) < config.rendering.minDrawableLength) continue;
            if (drawAbsorptionFading && s.before.medium != null &&
                s.before.medium.absorption(s.before.p, s.before.d, s.before.lambda, s.before.pol) > 0) {
                drawFadingSegment(pic,
                                  Label(""),
                                  s.before,
                                  tr,
                                  i,
                                  segs.length,
                                  p,
                                  currentArr,
                                  polarization,
                                  wavefronts,
                                  fadeSubpathSteps,
                                  applyOp);
            } else {
                drawRaySegment(pic,
                               Label(""),
                               s.before,
                               tr,
                               i,
                               segs.length,
                               effPen,
                               currentArr,
                               polarization,
                               wavefronts,
                               false);
            }
        }
        // Branch rendering for reflections
        for (Trace br : brs) {
            if (br.isReflection && !drawReflections) continue;
            Renderer branchRenderer = new Renderer;
            branchRenderer.p = p;
            branchRenderer.arrows = copy(arrows);
            if (segs.length > 0) branchRenderer.arrows.first = false;
            branchRenderer.before = (segs.length == 0) ? before : nan;
            branchRenderer.after = after;
            branchRenderer.margin = margin;
            branchRenderer.drawReflections = drawReflections;
            branchRenderer.polarization = polarization;
            branchRenderer.drawTrace(pic, br);
        }
    }
    void drawElement(picture pic = currentpicture, Element e) {
        if (e.fill != nullpen && cyclic(e.boundary)) {
            fill(pic, e.T * e.boundary, e.fill);
        }
        if (e.customDraw != null) {
            e.customDraw(pic, e);
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
          bool drawReflections = config.rendering.drawReflections) {
    Renderer r = new Renderer;
    r.p = p;
    r.before = before;
    r.after = after;
    r.margin = margin;
    r.arrows = arrows;
    r.polarization = polarization;
    r.wavefronts = wavefronts;
    r.drawReflections = drawReflections;
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
          bool drawReflections = config.rendering.drawReflections) {
    // 1. Draw each individual trace (rays + per-ray crests/wavelets)
    for (Trace t : traces) {
        draw(pic, t, before, after, margin, p, arrows, polarization, wavefronts, drawReflections);
    }
    // 2. Bundle iso-OPL contours
    if (wavefronts.contours && traces.length >= 2) {
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
    pen effPen = resolveRayPen(p, r, config.rendering.useIntensityOpacity);
    if (length == infinity) {
        drawClippedRay(pic, z, r.d, margin, effPen, arrow, false, true);
    } else if (length == -infinity) {
        drawClippedRay(pic, z, r.d, margin, effPen, arrow, true, false);
    } else if (length >= 0) {
        drawRaySegment(pic, L, r, z -- (z + length * r.d), 0, 1, effPen, arrow, polarization, wavefronts, false);
    } else {
        drawRaySegment(pic,
                       L,
                       r,
                       (z + length * r.d) -- z,
                       0,
                       1,
                       effPen,
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
                          Medium medium = Materials.VACUUM,
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
                           Medium medium = Materials.VACUUM,
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
        if (!res.totalInternalReflection && r.intensity() * res.T > config.physics.minIntensity) {
            real intensity = r.intensity() * res.T;
            Ray transRay = spawnRay(r, h, res.d, intensity, res.medium, 0, res.waveNormal);

            if (isnan(reflectance)) {
                real sqrtFactor = res.cosi > config.math.eps ? sqrt((res.n2 * res.cost) / (res.n1 * res.cosi)) : 0.0;

                // Mueller transmission transformation for Stokes vector
                real Ts = 1.0 - abs2(res.rs);
                real Tp = 1.0 - abs2(res.rp);
                pair C_T = res.ts * conj(res.tp) * sqrtFactor;
                applyDiattenuation(transRay.pol, r.pol, Ts, Tp);
            }
            rays.push(transRay);
        }

        // 2. Reflected Ray (handles TIR and specular reflection)
        if (r.intensity() * res.R > config.physics.minIntensity &&
            (res.totalInternalReflection || config.physics.traceReflections || config.rendering.drawReflections)) {
            real intensity = r.intensity() * res.R;
            real reflPhase = atan2(res.rs.y, res.rs.x); // pi for external reflection (n1 < n2)
            Ray reflRay =
                spawnRay(r, h, res.reflectedDirection, intensity, res.sourceMedium, reflPhase, res.reflectedWaveNormal);

            if (isnan(reflectance)) {
                // Mueller reflection transformation for Stokes vector
                applyDiattenuation(reflRay.pol, r.pol, res.rs, res.rp);
            }
            rays.push(reflRay);
        }

        return rays;
    };

    return p;
}

Property Reflect(real phaseShift = pi) {
    Property p = new Property;
    p.name = "Reflect";
    p.generate = new Ray[](Ray r, Hit h, mediumResolver resolve) {
        pair N = phaseGradientNormal(r, h);
        pair reflectedWaveNormal = unitOr(r.waveNormal - 2 * dot(r.waveNormal, N) * N, r.d);
        pair reflected = r.medium.groupDirection(reflectedWaveNormal);
        return new Ray[]{spawnRay(r, h, reflected, r.intensity(), r.medium, phaseShift, reflectedWaveNormal)};
    };
    // Specular Reflection alignment: Bisector normal!
    p.align = new pair(pair inDir, pair targetDir, real lambda, int orderIndex) { return targetDir - inDir; };
    return p;
}

Property IdealReflect(real f, real phaseShift = pi) {
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
        pair vOut = unitOr(r.d - 2.0 * dot(r.d, N_ideal) * N_ideal, r.d);
        return new Ray[]{spawnRay(r, h, vOut, r.intensity(), r.medium, phaseShift)};
    };
    p.align = new pair(pair inDir, pair targetDir, real lambda, int orderIndex) { return targetDir - inDir; };
    return p;
}

Property Focus(real f) {
    Property p = new Property;
    p.name = "Focus";
    p.stage = PropertyStage.Deflect;
    p.generate = new Ray[](Ray r, Hit h, mediumResolver resolve) {
        pair vOut = thinLensRefractDir(r, h, f);
        Ray out = spawnRay(r, h, vOut, r.intensity(), r.medium);
        out.waveNormal = vOut;
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
        pair vOut = unitOr(F - h.p, r.d);
        Ray out = spawnRay(r, h, vOut, r.intensity(), r.medium);
        out.waveNormal = vOut;
        return new Ray[]{out};
    };
    return p;
}

private
pair wavevectorShiftAxis(pair inDir, pair targetDir, real deltaK, bool reflective) {
    pair v = targetDir - inDir;
    real vMag = length(v);
    if (vMag < config.physics.wavevectorAlignmentTolerance || abs(deltaK) > vMag) {
        return (infinity, infinity);
    }
    real angleT = aCos(clamp(deltaK / vMag, -1, 1));
    pair vHat = v / vMag;
    pair vRot1 = rotate(angleT) * vHat;
    pair N1 = (-vRot1.y, vRot1.x);
    pair vRot2 = rotate(-angleT) * vHat;
    pair N2 = (-vRot2.y, vRot2.x);
    real sideProduct = dot(inDir, N1) * dot(targetDir, N1);
    if (reflective) return sideProduct < 0 ? N1 : N2;
    return sideProduct > 0 ? N1 : N2;
}

Property
Diffract(real lines_per_mm, int[] orders = new int[]{1}, bool reflective = false, real[] efficiencies = new real[]{}) {
    Property p = new Property;
    p.name = "Diffract";
    real gratingSpacing = 1000.0 / lines_per_mm;
    // Generation Closure
    p.generate = new Ray[](Ray r, Hit h, mediumResolver resolve) {
        Ray[] validRays;
        bool useEfficiencies = (efficiencies.length == orders.length);
        pair N = phaseGradientNormal(r, h);
        pair T = (-N.y, N.x);
        Medium m1 = (r.medium != null) ? r.medium : Materials.VACUUM;
        real n1 = m1.n(h.p, r.d, r.lambda, r.pol);
        Medium m2 = m1;
        if (!reflective) {
            pair targetPoint = h.p + config.geometry.epsIntersection * r.d;
            m2 = resolve(targetPoint, r, null, null, true);
            if (m2 == null) m2 = Materials.VACUUM;
        }
        real n2 = reflective ? n1 : m2.n(h.p, r.d, r.lambda, r.pol);
        real kIn = dot(r.d, T);
        real lambda_um = wavelengthInUnits(r.lambda);
        real phaseShift = reflective ? pi : 0;
        for (int i = 0; i < orders.length; ++i) {
            real deltaK = orders[i] * lambda_um / gratingSpacing;
            real kOut = (n1 * kIn + deltaK) / n2;
            if (abs(kOut) <= 1.0) {
                pair vOut = diffractionOutputDir(N, T, kOut, reflective);
                real eff = useEfficiencies ? max(0.0, efficiencies[i]) : 1.0;
                real intensity = r.intensity() * eff;
                Ray diffRay = spawnRay(r, h, vOut, intensity, m2, phaseShift);
                validRays.push(diffRay);
            } else {
                // Optional diagnostic:
                write(format("Notice: Order m = %d", orders[i]) + format(" is evanescent at lambda = %.1f", r.lambda) +
                      format(" (|kOut| = %.2f > 1.0)", abs(kOut)));
            }
        }
        if (validRays.length > 0 && !useEfficiencies) {
            real splitIntensity = r.intensity() / validRays.length;
            for (int i = 0; i < validRays.length; ++i) {
                validRays[i].pol.scaleTo(splitIntensity);
            }
        }
        return validRays;
    };
    // Alignment Closure (Decouples align() completely!)
    p.align = new pair(pair inDir, pair targetDir, real lambda, int orderIndex) {
        int targetOrder = (orders.length > orderIndex) ? orders[orderIndex] : 0;
        real deltaK = targetOrder * wavelengthInUnits(lambda) / gratingSpacing;
        return wavevectorShiftAxis(inDir, targetDir, deltaK, reflective);
    };
    return p;
}

Property PhaseShift(real deltaK, bool reflective = false) {
    Property p = new Property;
    p.name = "PhaseShift";
    p.generate = new Ray[](Ray r, Hit h, mediumResolver resolve) {
        pair N = phaseGradientNormal(r, h);
        pair T = (-N.y, N.x);
        real kIn = dot(r.d, T);
        // real kOut = kIn + deltaK; // Initial version - without medium correction
        // real n = r.medium.n(r.p, r.d, r.lambda, r.pol); // Get the local index of the medium
        // real kOut = kIn + deltaK / n;                   // Scale the transverse wavevector s
        Medium m1 = (r.medium != null) ? r.medium : Materials.VACUUM;
        real n1 = m1.n(h.p, r.d, r.lambda, r.pol);
        Medium m2 = reflective ? m1 : resolve(h.p + config.geometry.epsIntersection * r.d, r, null, null, true);
        real n2 = m2.n(h.p, r.d, r.lambda, r.pol);
        real kOut = (n1 * dot(r.d, T) + deltaK) / n2;
        if (abs(kOut) > 1.0) return new Ray[];
        pair vOut = diffractionOutputDir(N, T, kOut, reflective);
        real pShift = reflective ? pi : 0;
        return new Ray[]{spawnRay(r, h, vOut, r.intensity(), reflective ? r.medium : m2, pShift)};
    };
    p.align = new pair(pair inDir, pair targetDir, real lambda, int orderIndex) {
        return wavevectorShiftAxis(inDir, targetDir, deltaK, reflective);
    };
    return p;
}

Property Absorb() {
    Property p = new Property;
    p.name = "Absorb";
    p.generate = new Ray[](Ray r, Hit h, mediumResolver resolve) { return new Ray[]; };
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
        bool isEntering = (r.medium != med_o && r.medium != med_e);

        if (isEntering) {
            // =================================================================
            // CASE 1: ENTERING CRYSTAL (External Ray splits into O and E modes)
            // =================================================================
            real Is_inc = max(0.0, 0.5 * (r.pol.S0 + r.pol.S1));
            real Ip_inc = max(0.0, 0.5 * (r.pol.S0 - r.pol.S1));

            FresnelResult b = fresnelUniaxial(r, h, f.sourceMedium, med_o, med_e, resolve);

            // 1. Ordinary Ray (Pure S-polarization in med_o)
            if (!b.tir_o) {
                real I_o = Is_inc * b.To;
                if (I_o > config.physics.minIntensity) {
                    Ray oRay = spawnRay(r, h, b.d_o, I_o, med_o, 0, b.waveNormal_o);
                    oRay.pol.S0 = I_o;
                    oRay.pol.S1 = I_o;
                    oRay.pol.S2 = 0.0;
                    oRay.pol.S3 = 0.0;
                    rays.push(oRay);
                }
            }

            // 2. Extraordinary Ray (Pure P-polarization in med_e with walk-off)
            if (!b.tir_e) {
                real I_e = Ip_inc * b.Te;
                if (I_e > config.physics.minIntensity) {
                    Ray eRay = spawnRay(r, h, b.d_e, I_e, med_e, 0, b.waveNormal_e);
                    eRay.pol.S0 = I_e;
                    eRay.pol.S1 = -I_e;
                    eRay.pol.S2 = 0.0;
                    eRay.pol.S3 = 0.0;
                    rays.push(eRay);
                }
            }

            // 3. Reflected Ray (Reflecting back into source medium)
            if (config.physics.traceReflections || config.rendering.drawReflections) {
                real I_refl = b.Rs * Is_inc + b.Rp * Ip_inc;
                if (I_refl > config.physics.minIntensity) {
                    Ray reflRay =
                        spawnRay(r, h, b.reflectedDirection, I_refl, f.sourceMedium, 0, b.reflectedWaveNormal);
                    pair C = b.rs * conj(b.rp);

                    reflRay.pol.S0 = I_refl;
                    reflRay.pol.S1 = b.Rs * Is_inc - b.Rp * Ip_inc;
                    reflRay.pol.S2 = C.x * r.pol.S2 - C.y * r.pol.S3;
                    reflRay.pol.S3 = C.y * r.pol.S2 + C.x * r.pol.S3;
                    rays.push(reflRay);
                }
            }
        } else {
            // =================================================================
            // CASE 2: EXITING CRYSTAL (Ray is already in med_o or med_e)
            // =================================================================
            Medium srcMed = (r.medium != null) ? r.medium : med_o;
            FresnelResult res = fresnelIsotropic(srcMed, ambient, r, h, resolve, nan);

            if (!res.totalInternalReflection) {
                real I_trans = r.intensity() * res.T;
                if (I_trans > config.physics.minIntensity) {
                    Ray transRay = spawnRay(r, h, res.d, I_trans, res.medium, 0, res.waveNormal);
                    transRay.pol = copy(r.pol);
                    transRay.pol.scaleTo(I_trans);
                    rays.push(transRay);
                }
            }

            // Internal reflection (e.g. TIR of ordinary ray in Nicol prism)
            if (config.physics.traceReflections || config.rendering.drawReflections || res.totalInternalReflection) {
                real I_refl = r.intensity() * res.R;
                if (I_refl > config.physics.minIntensity) {
                    Ray reflRay =
                        spawnRay(r, h, res.reflectedDirection, I_refl, res.sourceMedium, 0, res.reflectedWaveNormal);
                    reflRay.pol = copy(r.pol);
                    reflRay.pol.scaleTo(I_refl);
                    rays.push(reflRay);
                }
            }
        }

        return rays;
    };
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
    static real thickLensMaxDiameter(real R1, real R2, real thickness, real minEdgeThickness = 0.0) {
        if (thickness <= config.math.eps) return 0.0;
        // 1. Hemisphere bound: h cannot exceed sphere radius
        real h_max = infinity;
        if (abs(R1) > config.math.eps) h_max = min(h_max, config.geometry.hemisphereMargin * abs(R1));
        if (abs(R2) > config.math.eps) h_max = min(h_max, config.geometry.hemisphereMargin * abs(R2));
        // Helper to evaluate edge thickness at height y
        real edgeThick(real y) {
            real s1 = (abs(R1) > config.math.eps) ? (abs(R1) - sqrt(max(0.0, R1 * R1 - y * y))) : 0.0;
            real s2 = (abs(R2) > config.math.eps) ? (abs(R2) - sqrt(max(0.0, R2 * R2 - y * y))) : 0.0;
            real dx1 = (R1 > config.math.eps) ? s1 : ((R1 < -config.math.eps) ? -s1 : 0.0);
            real dx2 = (R2 < -config.math.eps) ? -s2 : ((R2 > config.math.eps) ? s2 : 0.0);
            return thickness + dx2 - dx1;
        }
        // 2. Crossover bound: edge thickness must be >= minEdgeThickness
        if (h_max < infinity && edgeThick(h_max) < minEdgeThickness) {
            // Fast bisection to find knife-edge height where edgeThick(y) == minEdgeThickness
            real yLow = 0.0;
            real yHigh = h_max;
            for (int i = 0; i < config.geometry.maxBisectionSteps; ++i) {
                real yMid = 0.5 * (yLow + yHigh);
                if (edgeThick(yMid) < minEdgeThickness)
                    yHigh = yMid;
                else
                    yLow = yMid;
            }
            h_max = yLow;
        }
        return (h_max < infinity) ? (2.0 * h_max) : infinity;
    }
    // Robust thick lens boundary path generator
    static path thickLens(real R1, real R2, real thickness, real diameter, real minEdgeThickness = 0.0) {
        real t_safe = max(thickness, config.geometry.eps);
        real maxD = thickLensMaxDiameter(R1, R2, t_safe, minEdgeThickness);
        // Auto-clamp requested diameter to physical maximum:
        real effDiameter = min(diameter, maxD);
        real h = effDiameter / 2.0;
        pair frontVertex = (-t_safe / 2.0, 0);
        pair backVertex = (t_safe / 2.0, 0);
        // 1. FRONT SURFACE ARC (curves from y = -h to y = +h)
        path frontArc;
        if (R1 > config.math.eps) {
            pair C1 = frontVertex + (R1, 0);
            real a1 = aSin(clamp(h / R1, -1.0, 1.0));
            frontArc = reverse(Arc(C1, R1, 180 - a1, 180 + a1, config.geometry.arcSubdivisions));
        } else if (R1 < -config.math.eps) {
            pair C1 = frontVertex + (R1, 0);
            real a1 = aSin(clamp(h / abs(R1), -1.0, 1.0));
            frontArc = Arc(C1, abs(R1), -a1, a1, config.geometry.arcSubdivisions);
        } else {
            frontArc = (frontVertex.x, -h) -- (frontVertex.x, h);
        }
        // 2. REAR SURFACE ARC (curves from y = -h to y = +h)
        path backArc;
        if (R2 < -config.math.eps) {
            pair C2 = backVertex + (R2, 0);
            real a2 = aSin(clamp(h / abs(R2), -1.0, 1.0));
            backArc = Arc(C2, abs(R2), -a2, a2, config.geometry.arcSubdivisions);
        } else if (R2 > config.math.eps) {
            pair C2 = backVertex + (R2, 0);
            real a2 = aSin(clamp(h / R2, -1.0, 1.0));
            backArc = reverse(Arc(C2, R2, 180 - a2, 180 + a2, config.geometry.arcSubdivisions));
        } else {
            backArc = (backVertex.x, -h) -- (backVertex.x, h);
        }
        // 3. ASSEMBLE CLOSED NON-SELF-INTERSECTING BOUNDARY
        // frontArc goes (-h -> +h), reverse(backArc) goes (+h -> -h)
        return frontArc -- reverse(backArc) -- cycle;
    }
}

// Default Surface Drawing
// Native Asymptote Arrowhead for Thin Lenses (Convex & Concave)
arrowhead LensHead(bool concave = false, real angle = 20, real bow = 0.3, real barb = 0.5) {
    arrowhead ah = new arrowhead;
    ah.head = new path(path g, position position = EndPoint, pen p = currentpen, real size = 0, real angle = angle) {
        if (size == 0) size = config.rendering.minLensHead > 0 ? config.rendering.minLensHead : 5.0;
        // Resolve path endpoint P0 and outward unit tangent d
        bool relative = position.relative;
        real pos = position.position.x;
        real t = relative ? reltime(g, pos) : pos;
        pair P0 = point(g, t);
        pair d = unit(dir(g, t));
        // Orientation:
        pair arrowDir = concave ? -d : d;
        pair n = (-arrowDir.y, arrowDir.x);
        real w = size * tan(angle * pi / 180.0);
        // Construct standard hook shape locally with tip at (0, 0)
        pair M = (0, 0);
        pair tip = size * barb * arrowDir;
        pair pL = size * (barb - 1) * arrowDir + w * n;
        pair pR = size * (barb - 1) * arrowDir - w * n;
        pair cL = size * (barb - 0.5) * arrowDir + (w * bow) * n;
        pair cR = size * (barb - 0.5) * arrowDir - (w * bow) * n;
        path hookPath = (concave) ? M -- pL .. cL .. tip -- tip .. cR .. pR -- cycle
                                  : (-tip) -- M -- pL .. cL .. tip -- tip .. cR .. pR -- M -- cycle;
        // Attachment Point:
        pair attachPoint = (concave) ? -(size * barb) * arrowDir : (size * barb) * arrowDir;
        return shift(P0 - attachPoint) * hookPath;
    };
    return ah;
}

// Convenient pre-configured arrowheads:
arrowhead LensConvexHead = LensHead(concave = false);
arrowhead LensConcaveHead = LensHead(concave = true);

// Thin Lens Drawing Delegate (Arrowhead hooks & Focal points)
void drawThinLens(picture pic, Element e) {
    if (e.surfaces.length == 0) return;
    Surface lensSurface = transformedSurface(e.surfaces[0], e.T);
    pair center = lensSurface.center;
    pair axis = lensSurface.axis;
    pair q = (-axis.y, axis.x);
    real radius = 0;
    if (length(lensSurface.p) > 0) {
        path globalPath = lensSurface.T * lensSurface.p;
        radius = abs(point(globalPath, length(globalPath)) - point(globalPath, 0)) / 2.0;
    }
    pair a = center - radius * q;
    pair b = center + radius * q;
    real headSize = config.rendering.minLensHead;
    real fVal = !isnan(e.focalLength) ? e.focalLength : 0;
    // Single declarative draw line using native Asymptote Arrows!
    draw(pic, a -- b, e.stroke, Arrows(LensHead(concave = (fVal < 0)), size = headSize, filltype = FillDraw));
    // Optional focal points F1, F2
    if (e.showInfo) {
        pair F1 = center - fVal * axis;
        pair F2 = center + fVal * axis;
        dot(pic, F1, e.stroke);
        dot(pic, F2, e.stroke);
        label(pic, "$F_1$", F1, -q, e.stroke);
        label(pic, "$F_2$", F2, -q, e.stroke);
    }
}

// Grating Drawing Delegate (Tick marks)
void drawGrating(picture pic = currentpicture, Element e) {
    if (e.surfaces.length == 0) return;
    Surface gratingSurface = transformedSurface(e.surfaces[0], e.T);
    path globalPath = gratingSurface.T * gratingSurface.p;
    real L_total = arclength(globalPath);
    if (L_total <= config.rendering.minDrawableLength) return;
    int nlines = config.rendering.gratingTickCount;
    real tickLength = config.rendering.gratingTickLength;
    real alpha = config.rendering.gratingTickAngle;
    pair tangent = unitOr(dir(globalPath, 0));
    pair nv = -unitOr(gratingSurface.axis);
    pair tickDir = rotate(alpha) * nv;
    pair tickVector = -tickLength * tickDir;
    real proj = dot(tickVector, tangent);
    real margin = abs(proj);
    real s_start = margin;
    real s_end = L_total - margin;
    real range = s_end - s_start;
    if (range > 0 && nlines > 0) {
        real stepL = range / nlines;
        for (int i = 0; i <= nlines; ++i) {
            real s = s_start + i * stepL;
            real t = arctime(globalPath, s);
            pair pt = point(globalPath, t);
            draw(pic, pt -- pt + tickVector, e.stroke + 0.5 * linewidth(e.stroke));
        }
    }
    draw(pic, globalPath, e.stroke);
}

struct Lens {
    static Element thin(pair center,
                        real angle,
                        real f,
                        real diameter,
                        pen stroke = config.rendering.penLens,
                        string name = "Thin Lens") {
        pair A = dir(angle);
        path boundary = linearBoundary(center, A, diameter);
        Property[] props = {Focus(f)};
        Element e = element(boundary, props, stroke, nullpen, name, center, A);
        e.focalLength = f;
        e.customDraw = drawThinLens;
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
                         pen stroke = config.rendering.penLens,
                         pen fill = nullpen,
                         string name = "Thick Lens") {
        path localP = LensPath.thickLens(R1, R2, thickness, diameter);
        transform T = shift(center) * rotate(angle);
        Element e = Interface(T * localP, glass, ambient, stroke, fill, name);

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
                            pen stroke = config.rendering.penLens,
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
                             pen stroke = config.rendering.penLens,
                             pen fill = nullpen,
                             string name = "Biconcave Lens") {
        return thick(center, angle, -R1, R2, thickness, diameter, glass, ambient, stroke, fill, name);
    }

    static Element custom(path boundary,
                          Medium glass,
                          Medium ambient = null,
                          pen stroke = config.rendering.penLens,
                          pen fill = nullpen,
                          string name = "Custom Lens") {
        return Interface(boundary, glass, ambient, stroke, fill, name);
    }
}

// Master Mirror Factories
struct Mirror {
    static Element
    flat(pair center, real angle, real diameter, pen stroke = config.rendering.penMirror, string name = "Flat Mirror") {
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
                             pen stroke = config.rendering.penMirror,
                             string name = "Spherical Mirror") {
        real R = 2 * abs(f);
        pair N = dir(angle);
        pair C = vertex + sgn(f) * R * N;
        real half = aSin(min(1.0, diameter / (2 * R)));
        real theta = degrees(vertex - C);
        path boundary = Arc(C, R, theta - half, theta + half, config.geometry.arcSubdivisions);
        Property[] props = {Reflect()};
        return element(boundary, props, stroke, nullpen, name, vertex, N);
    }

    // IDEAL Paraxial Spherical Mirror (Zero spherical aberration)
    static Element idealSpherical(pair vertex,
                                  real angle,
                                  real f,
                                  real diameter,
                                  pen stroke = config.rendering.penMirror,
                                  string name = "Ideal Spherical Mirror") {
        real R = 2 * abs(f);
        pair N = dir(angle);
        pair C = vertex + sgn(f) * R * N;
        real half = aSin(min(1.0, diameter / (2 * R)));
        real theta = degrees(vertex - C);
        path boundary = Arc(C, R, theta - half, theta + half, config.geometry.arcSubdivisions);
        Property[] props = {IdealReflect(f)};
        return element(boundary, props, stroke, nullpen, name, vertex, N);
    }

    static Element custom(path boundary, pen stroke = config.rendering.penMirror, string name = "Custom Mirror") {
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
        path localBox = (-h, -h) -- (h, -h) -- (h, h) -- (-h, h) -- cycle;
        path localDiag = (-h, -h) -- (h, h); // Bottom-left to top-right (Reflects +X rays UPWARD)
        transform T = shift(center) * rotate(angle);
        Element e = new Element;
        e.name = name;
        e.stroke = stroke;
        e.fill = fill;
        e.center = center;
        e.axis = dir(angle);
        e.boundary = T * localBox;
        e.medium = glass;
        // Surface 1: Outer Glass Cube Boundary (Glass <-> Ambient)
        Surface sOuter = new Surface;
        sOuter.name = name + " Outer";
        sOuter.p = ensureCCW(T * localBox);
        sOuter.stroke = stroke;
        sOuter.fill = fill;
        sOuter.center = center;
        sOuter.axis = e.axis;
        sOuter.add(Refract(glass, ambient)); // Standard Refraction
        e.add(sOuter);
        // Surface 2: Internal Coating (Glass <-> Glass with specified R ratio!)
        Surface sDiag = new Surface;
        sDiag.name = name + " Coating";
        sDiag.p = T * localDiag;
        sDiag.stroke = stroke;
        sDiag.fill = nullpen;
        sDiag.center = center;
        sDiag.axis = e.axis;
        sDiag.add(Refract(glass, glass, R)); // Custom Ratio Refraction!
        e.add(sDiag);
        return e;
    }
}

struct Polarizer {
    // 1. Ideal Linear Polarizer Plate
    static Element LP(pair center, real angle_deg = 0.0, real diameter = 30.0) {
        Medium polMed = medium("Polarizer").add(Dichroic(0.0, 1000.0));
        // Rotated to angle_deg
        return Interface(box(center - (1.0, diameter / 2), center + (1.0, diameter / 2)), polMed, Materials.VACUUM);
    }

    // 2. Quarter-Wave Plate (QWP) - Retardance = pi/2
    static Element QWP(pair center, real fastAxis_deg = 0.0, real diameter = 30.0) {
        // Retardance property delta = pi/2 rotated by fastAxis_deg
        Property qwp = property(
            "QWP",
            new Ray[](Ray r, Hit h, mediumResolver res) {
                r.pol.applyRotation(radians(-fastAxis_deg));
                r.pol.applyRetardance(pi / 2.0);
                r.pol.applyRotation(radians(fastAxis_deg));
                return new Ray[]{r};
            },
            PropertyStage.Interact);

        Element plate = Interface(
            linearBoundary(center, (cos(radians(fastAxis_deg + 90)), sin(radians(fastAxis_deg + 90))), diameter),
            Materials.VACUUM);
        plate.surfaces[0].add(qwp);
        return plate;
    }

    // 3. Half-Wave Plate (HWP) - Retardance = pi (Rotates linear polarization by 2*theta)
    static Element HWP(pair center, real fastAxis_deg = 0.0, real diameter = 30.0) {
        Property hwp = property(
            "HWP",
            new Ray[](Ray r, Hit h, mediumResolver res) {
                r.pol.applyRotation(radians(-fastAxis_deg));
                r.pol.applyRetardance(pi);
                r.pol.applyRotation(radians(fastAxis_deg));
                return new Ray[]{r};
            },
            PropertyStage.Interact);

        Element plate = Interface(
            linearBoundary(center, (cos(radians(fastAxis_deg + 90)), sin(radians(fastAxis_deg + 90))), diameter),
            Materials.VACUUM);
        plate.surfaces[0].add(hwp);
        return plate;
    }
}

// Absorber, Aperture, Grating & Phase Gradient Factories
Element Absorber(path boundary, pen stroke = config.rendering.penAbsorber, string name = "Absorber") {
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
                      real deltaK,
                      bool reflective = false,
                      pen stroke = lightblue,
                      string name = "Phase Gradient") {
    pair ax = dir(angle);
    path boundary = linearBoundary(center, ax, diameter);
    Property[] props = {PhaseShift(deltaK, reflective)};
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
    e.customDraw = drawGrating;
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
                        pair opticAxis = (1, 0),
                        Medium ambient = null,
                        pen stroke = black,
                        pen fill = nullpen,
                        string name = "Uniaxial Crystal") {
    Medium med_o = medium(name + " Ordinary", no);
    Medium med_e = uniaxialExtraordinary(name + " Extraordinary", no, ne, opticAxis);
    return Crystal(boundary, med_o, med_e, ambient, stroke, fill, name);
}

// ----------------------------------------------------------------------------
// Section 8: Geometric Alignment, Solvers & Optical Analysis
// ----------------------------------------------------------------------------
pair lineIntersection(Ray r1, Ray r2) { return lineIntersection(r1.p, r1.d, r2.p, r2.d); }
pair lineIntersection(Ray r, pair lineOrigin, pair lineDir) { return lineIntersection(r.p, r.d, lineOrigin, lineDir); }
pair intersection(Ray r1, Ray r2) { return lineIntersection(r1, r2); }

// Dynamic Geometric Alignment API
private
Element alignAxis(Element e, pair targetAxis) {
    pair a = unitOr(e.axis);
    pair n = unitOr(targetAxis, a);
    real angle = aCos(clamp(dot(a, n), -1, 1));
    if (cross(a, n) < 0) angle = -angle;
    return rotate(angle, e.center) * e;
}

Element align(Element e, pair fromPos, pair toPos, real lambda = config.physics.lambda, int orderIndex = 0) {
    pair inDir = unitOr(e.center - fromPos, e.axis);
    pair targetDir = unitOr(toPos - e.center, e.axis);
    // Check if any attached surface property defines a custom alignment rule:
    //    - Mirrors -> use specular bisector (targetDir - inDir)
    //    - Gratings -> use wavevectorShiftAxis()
    //    - Phase Gradients -> use deltaK shift
    for (Surface s : e.surfaces) {
        for (Property prop : s.props) {
            if (prop.align != null) {
                pair targetAxis = prop.align(inDir, targetDir, lambda, orderIndex);
                if (abs(xpart(targetAxis)) != infinity && abs(ypart(targetAxis)) != infinity) {
                    return alignAxis(e, targetAxis);
                }
            }
        }
    }
    // Default: In-line transmission alignment (Lenses, Apertures, Interfaces)
    return alignAxis(e, inDir);
}

Element align(Element e, pair fromPos, Element toElement, real lambda = config.physics.lambda, int orderIndex = 0) {
    return align(e, fromPos, toElement.center, lambda, orderIndex);
}

Element align(Element e, Element fromElement, pair toPos, real lambda = config.physics.lambda, int orderIndex = 0) {
    return align(e, fromElement.center, toPos, lambda, orderIndex);
}

Element
align(Element e, Element fromElement, Element toElement, real lambda = config.physics.lambda, int orderIndex = 0) {
    return align(e, fromElement.center, toElement.center, lambda, orderIndex);
}

path operator cast(Element e) { return e.T * e.boundary; }
path operator cast(Ray r) {
    return (r.p - config.geometry.rayToPathExtension * unit(r.d)) -- (r.p + config.geometry.rayToPathExtension * unit(r.d));
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
               int rayType = Ray.Incident,
               int rayIndex = 0,
               pen p = currentpen,
               arrowbar arrow = None,
               bool opposite = false) {
    if (inter == null || inter.hit == null) return 0;
    Ray r =
        (rayType == Ray.Incident)
            ? inter.incident
            : ((rayType == Ray.Reflected && rayIndex < inter.reflected.length)
                   ? inter.reflected[rayIndex]
                   : ((rayType == Ray.Transmitted && rayIndex < inter.transmitted.length) ? inter.transmitted[rayIndex]
                                                                                          : null));
    if (r == null) return 0;
    pair O = inter.hit.p;
    pair v_ray = (rayType == Ray.Incident) ? -unitOr(r.d) : unitOr(r.d);
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

real getAngle(Interaction inter, int rayType = Ray.Incident, int rayIndex = 0) {
    if (inter == null || inter.hit == null) return 0.0;
    Ray r;
    if (rayType == Ray.Incident)
        r = inter.incident;
    else if (rayType == Ray.Reflected) {
        if (rayIndex < 0 || inter.reflected.length <= rayIndex) return 0.0;
        r = inter.reflected[rayIndex];
    } else if (rayType == Ray.Transmitted) {
        if (rayIndex < 0 || inter.transmitted.length <= rayIndex) return 0.0;
        r = inter.transmitted[rayIndex];
    } else {
        return 0.0;
    }
    if (r == null) return 0.0;
    pair vr = (length(r.p - inter.hit.p) > 1e-4) ? -unit(r.d) : unit(r.d);
    pair vn = unit(inter.normal());
    if (dot(vr, vn) < 0) vn = -vn;
    real cosineVal = max(min(dot(vr, vn), 1.0), -1.0);
    return aCos(cosineVal);
}

// Evaluates signed perpendicular miss distance from target B to the ray exiting the system
private
real rayMissDistance(pair A, real angle, OpticalSystem sys, pair B, real lambda = 0) {
    real lam = (lambda == 0) ? config.physics.lambda : lambda;
    Ray r = Ray.angle(A, angle, lambda = lam);
    Trace t = sys.trace(r);
    Segment[] segs = t.segments();
    if (segs.length == 0) return nan;
    real bestMiss = nan;
    real minDistance = infinity;
    // Search ALL segments for the closest approach to target point B
    for (int i = 0; i < segs.length; ++i) {
        Segment s = segs[i];
        if (s.before == null || !s.before.alive) continue;
        pair p1 = s.a;
        pair p2 = s.b;
        pair d_seg = unitOr(p2 - p1, s.before.d);
        real segLen = length(p2 - p1);
        pair toB = B - p1;
        real proj = dot(toB, d_seg);
        bool isLast = (i == segs.length - 1);
        if (proj >= -config.math.eps && (proj <= segLen + config.math.eps || isLast)) {
            real perpDist = abs(cross(toB, d_seg));
            if (perpDist < minDistance) {
                minDistance = perpDist;
                bestMiss = cross(toB, d_seg);
            }
        }
    }
    return bestMiss;
}

// Robust Brent's Method Ray Aimer
real aimAngle(pair A,
              pair B,
              OpticalSystem sys,
              real angleMin = -45.0,
              real angleMax = 45.0,
              real lambda = 0,
              real tol = 1e-6,
              int maxIter = 40) {
    real a = min(angleMin, angleMax);
    real b = max(angleMin, angleMax);
    real fa = rayMissDistance(A, a, sys, B, lambda);
    real fb = rayMissDistance(A, b, sys, B, lambda);
    // Auto-Bracketing with NaN Safeguards
    if (isnan(fa) || isnan(fb) || fa * fb > 0) {
        int nSamples = 20;
        real step = (b - a) / (real)nSamples;
        real best_ang = a;
        real min_miss = infinity;
        bool foundBracket = false;
        real prev_ang = a;
        real prev_f = fa;
        for (int i = 1; i <= nSamples; ++i) {
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
                    foundBracket = true;
                    break;
                }
            }
            prev_ang = cur_ang;
            prev_f = cur_f;
        }
        if (!foundBracket) return best_ang;
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
    for (int iter = 0; iter < maxIter; ++iter) {
        if (abs(fb) < tol || abs(b - a) < tol) return b;
        real s = 0;
        real denom1 = (fa - fb) * (fa - fc);
        real denom2 = (fb - fa) * (fb - fc);
        real denom3 = (fc - fa) * (fc - fb);
        if (fa != fc && fb != fc && abs(denom1) > config.math.solverDenominatorEps &&
            abs(denom2) > config.math.solverDenominatorEps && abs(denom3) > config.math.solverDenominatorEps) {
            // Inverse Quadratic Interpolation
            s = (a * fb * fc) / denom1 + (b * fa * fc) / denom2 + (c * fa * fb) / denom3;
        } else {
            // Secant Method
            real denomSec = fb - fa;
            if (abs(denomSec) > config.math.solverDenominatorEps)
                s = b - fb * (b - a) / denomSec;
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
Ray aim(pair A, pair B, OpticalSystem sys, real angleMin = -45.0, real angleMax = 45.0, real lambda = 0) {
    real ang = aimAngle(A, B, sys, angleMin, angleMax, lambda);
    real lam = (lambda == 0) ? config.physics.lambda : lambda;
    return Ray.angle(A, ang, lambda = lam);
}

// Element[] array overload
Ray aim(pair A, pair B, Element[] elements, real angleMin = -45.0, real angleMax = 45.0, real lambda = 0) {
    return aim(A, B, system(elements), angleMin, angleMax, lambda);
}

// ============================================================================
// Unified Optical Imaging & Object Path Transformation Engine
// ============================================================================
// Conjugate Image Point of a Single Point P through System/Subsystem
pair imagePoint(OpticalSystem sys, pair P, real deltaAngle = config.physics.differentialAngle, real lambda = 0) {
    if (sys.elements.length == 0) return P;
    real lam = (lambda == 0) ? config.physics.lambda : lambda;
    pair targetOptic = sys.elements[0].center;
    pair d_center = unitOr(targetOptic - P, sys.elements[0].axis);
    real theta_center = degrees(d_center);
    Ray r1 = Ray.angle(P, theta_center - deltaAngle, lambda = lam);
    Ray r2 = Ray.angle(P, theta_center + deltaAngle, lambda = lam);
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
    pair axisOrigin;
    pair axisDir;
    pair axisNorm;
    bool isThick;
}

// Numerical Paraxial Solver: Traces forward and backward paraxial rays
// Numerical Paraxial Solver: Traces forward and backward paraxial rays
CardinalPoints cardinalPoints(OpticalSystem sys,
                              pair axisOrigin = (0, 0),
                              pair axisDir = (1, 0),
                              real paraxialH = config.physics.paraxialHeight) {
    CardinalPoints cp;
    cp.axisOrigin = axisOrigin;
    cp.axisDir = unit(axisDir);
    cp.axisNorm = (-cp.axisDir.y, cp.axisDir.x); // Transverse axis
    if (sys.elements.length == 0) return cp;

    real launchDistance = 0.5 * config.physics.far;
    real minProj = infinity;
    real maxProj = -infinity;
    FlatSurfaceRecord[] recs = prepareFlatSurfaces(sys.elements);
    for (FlatSurfaceRecord rec : recs) {
        if (rec.minBound.x < infinity && rec.maxBound.x > -infinity) {
            pair[] pts = {rec.minBound,
                          (rec.minBound.x, rec.maxBound.y),
                          rec.maxBound,
                          (rec.maxBound.x, rec.minBound.y)};
            for (pair pt : pts) {
                real pr = dot(pt - axisOrigin, cp.axisDir);
                if (pr < minProj) minProj = pr;
                if (pr > maxProj) maxProj = pr;
            }
        }
    }
    if (minProj < infinity && maxProj > -infinity) {
        real extent = max(maxProj - minProj, 1.0);
        real margin = max(min(20.0, 2.0 * extent), 0.2 * extent);
        real safeDist = max(abs(minProj), abs(maxProj)) + margin;
        launchDistance = min(safeDist, 0.5 * config.physics.far);
    }

    // -------------------------------------------------------------------------
    // 1. RESOLVE VERTICES V1 & V2 (Axial intercepts of front & rear surfaces)
    // -------------------------------------------------------------------------
    Ray r_axis_fwd = Ray(p = axisOrigin - launchDistance * cp.axisDir, d = cp.axisDir);
    Trace t_axis_fwd = sys.trace(r_axis_fwd, far = config.physics.far);
    if (t_axis_fwd.segments().length > 1) {
        cp.V1 = t_axis_fwd.segments()[1].a; // First surface hit point
    } else {
        cp.V1 = sys.elements[0].center;
    }

    Ray r_axis_bwd = Ray(p = axisOrigin + launchDistance * cp.axisDir, d = -cp.axisDir);
    Trace t_axis_bwd = sys.trace(r_axis_bwd, far = config.physics.far);
    if (t_axis_bwd.segments().length > 1) {
        cp.V2 = t_axis_bwd.segments()[1].a; // Rear surface hit point
    } else {
        cp.V2 = sys.elements[sys.elements.length - 1].center;
    }

    // -------------------------------------------------------------------------
    // 2. FORWARD TRACE: Rear Focal Point F2 & Rear Principal Point H2
    // -------------------------------------------------------------------------
    Ray r_fwd = Ray(p = axisOrigin - launchDistance * cp.axisDir + paraxialH * cp.axisNorm, d = cp.axisDir);
    Trace t_fwd = sys.trace(r_fwd, far = config.physics.far);
    Ray exit_fwd = t_fwd.last();

    // Directly intersect the exit ray with the optical axis line:
    cp.F2 = lineIntersection(exit_fwd, axisOrigin, cp.axisDir);

    if (abs(cp.F2.x) != infinity) {
        real slope_exit = dot(exit_fwd.d, cp.axisNorm);
        cp.f2 = -paraxialH / slope_exit;
        cp.H2 = cp.F2 - cp.f2 * cp.axisDir;
    } else {
        // Afocal / Flat window system
        cp.f2 = infinity;
        cp.F2 = (infinity, infinity);
        cp.H2 = cp.V2;
    }

    // -------------------------------------------------------------------------
    // 3. BACKWARD TRACE: Front Focal Point F1 & Front Principal Point H1
    // -------------------------------------------------------------------------
    Ray r_bwd = Ray(p = axisOrigin + launchDistance * cp.axisDir + paraxialH * cp.axisNorm, d = -cp.axisDir);
    Trace t_bwd = sys.trace(r_bwd, far = config.physics.far);
    Ray exit_bwd = t_bwd.last();

    // Directly intersect the backward exit ray with the optical axis line:
    cp.F1 = lineIntersection(exit_bwd, axisOrigin, -cp.axisDir);

    if (abs(cp.F1.x) != infinity) {
        real slope_bwd = dot(exit_bwd.d, cp.axisNorm);
        cp.f1 = -paraxialH / slope_bwd;
        cp.H1 = cp.F1 + cp.f1 * cp.axisDir;
    } else {
        cp.f1 = infinity;
        cp.F1 = (-infinity, -infinity);
        cp.H1 = cp.V1;
    }

    // -------------------------------------------------------------------------
    // 4. NODAL POINTS N1 & N2
    // -------------------------------------------------------------------------
    // In Gaussian optics: N1 = H1 + (1 - n1/n2) * f2 * axisDir
    // When immersed in the same medium on both sides (n1 == n2): N1 == H1, N2 == H2
    real n_obj = (r_fwd.medium != null) ? r_fwd.medium.n() : 1.0;
    real n_img = (exit_fwd.medium != null) ? exit_fwd.medium.n() : 1.0;
    real nodalShift = (abs(n_img) > config.math.eps) ? (1.0 - n_obj / n_img) * cp.f2 : 0.0;

    cp.N1 = cp.H1 + nodalShift * cp.axisDir;
    cp.N2 = cp.H2 + nodalShift * cp.axisDir;

    // -------------------------------------------------------------------------
    // 5. DISTANCES (BFL, FFL) & THICKNESS FLAG
    // -------------------------------------------------------------------------
    cp.BFL = dot(cp.F2 - cp.V2, cp.axisDir);
    cp.FFL = dot(cp.V1 - cp.F1, cp.axisDir);
    cp.isThick = (length(cp.H2 - cp.H1) > config.geometry.eps);

    return cp;
}

// 1. One-Line 3-Ray Gaussian Generator (replaces ad-hoc ray construction)
Ray[] gaussianRays(OpticalSystem sys, pair objTip, real axisAngle = nan) {
    if (sys.elements.length == 0) return new Ray[];
    Element firstOptic = sys.elements[0];
    pair axis = !isnan(axisAngle) ? dir(axisAngle) : unitOr(firstOptic.axis, (1, 0));
    pair q = (-axis.y, axis.x);

    CardinalPoints cp;
    if (!isnan(firstOptic.focalLength) && abs(firstOptic.focalLength) > config.math.eps) {
        cp.H1 = firstOptic.center;
        cp.F1 = firstOptic.center - firstOptic.focalLength * axis;
    } else {
        cp = cardinalPoints(system(firstOptic, sys.ambient), firstOptic.center, axis);
    }

    Ray r_par = Ray(objTip, axis);
    Ray r_ctr = Ray.between(objTip, cp.H1);

    // Focal ray with aperture-clipping safeguard:
    Ray r_foc;
    if (abs(cp.F1.x) != infinity) {
        r_foc = Ray.between(objTip, cp.F1);
    } else {
        r_foc = Ray(objTip, axis);
    }

    return new Ray[]{r_par, r_ctr, r_foc};
}

// 2. Rich ImageResult Struct (replaces raw path return)
struct ImageResult {
    path imagePath;
    real magnification = 1.0;
    pair base = (nan, nan);
    pair tip = (nan, nan);
    bool isReal = true;
}

// Implicit cast: 'path p = image(...)' continues to work 100% identically!
path operator cast(ImageResult res) {
    return res.imagePath;
}

// 3. Upgraded image() (returns ImageResult with .magnification, .base, .tip)
ImageResult image(OpticalSystem sys,
                  path objectPath,
                  int samples = 20,
                  real deltaAngle = config.physics.differentialAngle,
                  real lambda = 0) {
    ImageResult res;
    int nSegs = length(objectPath);
    if (sys.elements.length == 0 || nSegs == 0) return res;

    guide imgGuide;
    real totalL = arclength(objectPath);
    int nPoints = max(samples, nSegs * 4);

    for (int i = 0; i <= nPoints; ++i) {
        real s = i * totalL / (real)nPoints;
        pair objPt = point(objectPath, arctime(objectPath, s));
        pair imgPt = imagePoint(sys, objPt, deltaAngle, lambda);
        if (abs(imgPt.x) != infinity && abs(imgPt.y) != infinity) {
            imgGuide = imgGuide -- imgPt;
        }
    }
    if (cyclic(objectPath)) imgGuide = imgGuide -- cycle;
    res.imagePath = imgGuide;

    if (length(res.imagePath) > 0) {
        pair objBase = point(objectPath, 0);
        pair objTip = point(objectPath, length(objectPath));
        res.base = point(res.imagePath, 0);
        res.tip = point(res.imagePath, length(res.imagePath));

        real h_obj = length(objTip - objBase);
        real h_img = length(res.tip - res.base);
        res.magnification = (h_obj > config.math.eps) ? (h_img / h_obj) : 1.0;
        if (dot(res.tip - res.base, objTip - objBase) < 0) res.magnification = -res.magnification;
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

    Element firstOptic = sys.elements[0];
    pair axis0 = unitOr(firstOptic.axis, (1, 0));

    // 1. Resolve front vertex V1
    real searchDist = 0.5 * config.physics.far;
    Ray r_find_v1 = Ray(p = firstOptic.center - searchDist * axis0, d = axis0);
    Trace t_find = sys.trace(r_find_v1, far = config.physics.far);
    pair V1 = (t_find.segments().length > 1) ? t_find.segments()[1].a : firstOptic.center;

    // 2. Trace optical axis from front vertex V1 through the entire system
    Ray axisRay = Ray(V1, axis0);
    Trace t_axis = sys.trace(axisRay, far = config.physics.far);
    Segment[] segs = t_axis.segments();
    if (segs.length == 0) return;

    int nSegs = segs.length;

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
    // Segments 0 to nSegs - 2 are internal to the system:
    for (int i = 0; i < nSegs - 1; ++i) {
        draw(pic, segs[i].traj, p);
    }

    // -------------------------------------------------------------------------
    // C. EXIT EXTENSION (After rear vertex V2)
    // -------------------------------------------------------------------------
    Segment lastSeg = segs[nSegs - 1];
    pair V2 = lastSeg.a; // Rear vertex of the final element
    pair d_exit = unitOr(lastSeg.before.d, axis0);

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
          pen pRay = black,
          pen pVirtual = pRay + dashed,
          pen pObject = black,
          pen pImage = pObject,
          pen pAxis = gray(0.5) + dashed,
          bool intermediate = true,
          bool showAxis = true,
          bool scaleArrows = true,
          real arrowSize = 5pt,
          real before = nan,
          real after = infinity,
          real axisBefore = infinity, // Default: extends to canvas edge
          real axisAfter = infinity,  // Default: extends to canvas edge
          real far = config.physics.far) {
    if (sys.elements.length == 0 || length(object) == 0) return;

    pair objBase = point(object, 0);
    pair objTip = point(object, length(object));
    real h_obj = length(objTip - objBase);

    // Axis and Ray tracing
    pair firstAxis = unitOr(sys.elements[0].axis, (1, 0));
    Trace t_axis = sys.trace(Ray(objBase, firstAxis), far = far);
    Segment[] axisSegs = t_axis.segments();

    Ray[] rays = gaussianRays(sys, objTip);
    Trace[] traces;
    for (Ray r : rays) traces.push(sys.trace(r, far = far));

    // Resolve true optical stages (skipping internal glass propagation)
    int numSegs = traces[0].segments().length;
    pair[] stageTips, stageBases;
    bool[] stageIsReal;

    for (int s = 1; s < numSegs; ++s) {
        Segment seg0 = traces[0].segments()[s];
        Segment seg1 = traces[1].segments()[s];
        if (seg0.before.medium != sys.ambient && s < numSegs - 1) continue; // Skip glass

        pair tip = lineIntersection(seg0.a, seg0.before.d, seg1.a, seg1.before.d);
        if (abs(tip.x) == infinity || abs(tip.y) == infinity) continue;

        bool isReal = (dot(tip - seg0.a, seg0.before.d) > config.geometry.eps);
        int axIdx = min(s, axisSegs.length - 1);
        pair axP1 = axisSegs[axIdx].a;
        pair axD = unitOr(axisSegs[axIdx].b - axP1, axisSegs[axIdx].before.d);
        pair base = axP1 + dot(tip - axP1, axD) * axD;

        stageTips.push(tip);
        stageBases.push(base);
        stageIsReal.push(isReal);
    }

    // Render Axis, Rays, Extensions & Arrows
    // Render Folded Optical Axis (with clean start and end bounds)
    if (showAxis) {
        drawAxis(pic, sys, before = axisBefore, after = axisAfter, p = pAxis);
    }

    for (Trace t : traces)
        draw(pic,
             t,
             before = before,
             after = after,
             p = pRay,
             arrows = rayArrows(first = true, middle = false, last = true));

    for (int k = 0; k < stageTips.length; ++k) {
        if (!stageIsReal[k]) {
            for (Trace t : traces)
                if (k + 1 < t.segments().length) draw(pic, t.segments()[k + 1].a -- stageTips[k], pVirtual);
        }
    }

    draw(pic, objBase -- objTip, pObject, EndArrow(DefaultHead, size = arrowSize, filltype = Draw));

    for (int k = 0; k < stageTips.length; ++k) {
        if (!intermediate && k < stageTips.length - 1) continue;
        real h_img = length(stageTips[k] - stageBases[k]);
        if (h_img < config.math.eps) continue;
        real effArrow = scaleArrows ? (h_img / h_obj) * arrowSize : arrowSize;
        pen arrPen = (k == stageTips.length - 1) ? pImage : (0.7 * pImage + 0.3 * white);
        if (!stageIsReal[k]) arrPen = arrPen + dashed;
        draw(pic,
             stageBases[k] -- stageTips[k],
             arrPen,
             EndArrow(DefaultHead, size = effArrow, filltype = Draw));
    }
    draw(pic, sys);
}

void draw(picture pic = currentpicture,
          CardinalPoints cp,
          real planeHeight = 35.0,
          pen pAxis = gray(0.5) + dashed + 0.6bp,
          pen pPlanes = rgb(0.1, 0.4, 0.8) + dashed + 0.8bp,
          pen pPoints = heavyred) {
    if (abs(cp.f2) == infinity) return;
    // 1. Draw Optical Axis
    real x_left = min(cp.F1.x, cp.V1.x) - 30.0;
    real x_right = max(cp.F2.x, cp.V2.x) + 30.0;
    pair axA = cp.axisOrigin + x_left * cp.axisDir;
    pair axB = cp.axisOrigin + x_right * cp.axisDir;
    draw(pic, axA -- axB, pAxis);
    // 2. Draw Principal Planes (H1 & H2)
    pair H1_top = cp.H1 + (planeHeight / 2.0) * cp.axisNorm;
    pair H1_bot = cp.H1 - (planeHeight / 2.0) * cp.axisNorm;
    pair H2_top = cp.H2 + (planeHeight / 2.0) * cp.axisNorm;
    pair H2_bot = cp.H2 - (planeHeight / 2.0) * cp.axisNorm;
    draw(pic, H1_bot -- H1_top, pPlanes);
    draw(pic, H2_bot -- H2_top, pPlanes);
    label(pic, "\small $\mathbf{H}_1$", H1_top, N, rgb(0.1, 0.4, 0.8));
    label(pic, "\small $\mathbf{H}_2$", H2_top, N, rgb(0.1, 0.4, 0.8));
    // 3. Mark Cardinal Points on Axis
    dot(pic, cp.F1, pPoints + 3pt);
    dot(pic, cp.F2, pPoints + 3pt);
    dot(pic, cp.H1, blue + 3pt);
    dot(pic, cp.H2, blue + 3pt);
    label(pic, "\small $\mathbf{F}_1$", cp.F1, 2S, pPoints);
    label(pic, "\small $\mathbf{F}_2$", cp.F2, 2S, pPoints);
    label(pic, "\small $H_1$", cp.H1, 2S, blue);
    label(pic, "\small $H_2$", cp.H2, 2S, blue);
    // 4. Mark Vertices
    dot(pic, cp.V1, black + 2pt);
    dot(pic, cp.V2, black + 2pt);
    label(pic, "\tiny $V_1$", cp.V1, NW, gray(0.3));
    label(pic, "\tiny $V_2$", cp.V2, NE, gray(0.3));
    // 5. Dimension Line for Effective Focal Length (EFL)
    real dimY = (planeHeight / 2.0) + 12.0;
    pair dimA = cp.H2 + dimY * cp.axisNorm;
    pair dimB = cp.F2 + dimY * cp.axisNorm;
    draw(pic, dimA -- dimB, pPoints + 0.6bp, Arrows(HookHead, 2.5pt));
    label(pic, format("\small $f = %#.2f$", cp.f2), midpoint(dimA -- dimB), N, pPoints);
}

// ----------------------------------------------------------------------------
// Section 9: Visual Insets, Legend Panels & Declarative DSL
// ----------------------------------------------------------------------------
real resolveSegmentOffset(Segment seg,
                          int idx,
                          int totalSegments,
                          real relative = nan,
                          real absolute = nan,
                          string caller = "resolveSegmentOffset") {
    real L = arclength(seg.traj);
    // Single Unbounded Ray (No optical collisions) -> Allow standard relative or absolute
    if (totalSegments == 1) {
        if (!isnan(relative)) return clamp(relative, 0.0, 1.0) * L;
        if (!isnan(absolute)) return clamp(abs(absolute) < config.math.eps ? 0.5 * L : abs(absolute), 0.0, L);
        return 0.5 * L;
    }
    // Open Entrance Ray (First Segment) -> Anchored at first surface hit
    if (idx == 0 && totalSegments > 1) {
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
    if (idx == totalSegments - 1) {
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
    real margin = 2.0;        // Inner padding for box (pt)
    // Color Palette
    pen border = rgb(0.0, 0.2, 0.6); // Box outline, connector line, & target dot color
    pen fill = white;                // Inset box background fill
    pen connector = nullpen;         // Connector line pen (defaults to border + dashed)
    pen badgePen = nullpen;          // Badge fill color (defaults to border)
    pen badgeText = white;           // Badge text label color
    void add(Label L) { label(this.pic, L); }
    void addGraphic(string filename) { label(this.pic, graphic(filename)); }
    // Helpers to resolve fallback pens
    pen getConnectorPen() { return (this.connector != nullpen) ? this.connector : (this.border + dashed); }
    pen getBadgePen() { return (this.badgePen != nullpen) ? this.badgePen : this.border; }
    // Standalone Drawing
    void draw(picture dest = currentpicture, pair position = (nan, nan), pen p = nullpen) {
        pair draw_pos = (!isnan(position.x)) ? position : this.pos;
        dest.add(new void(frame f, transform t) {
            pair pos_pt = t * draw_pos;
            frame content = this.pic.fit();
            frame box_frame;
            // Wrap content in box and shift to target canvas position
            filltype ft = (this.fill != nullpen) ? FillDraw(this.fill, this.border) : Fill(this.border);
            box(box_frame, content, this.margin, this.margin, this.border, ft, above = false);
            add(box_frame, content);
            box_frame = shift(pos_pt) * box_frame;
            // Draw connector line and dot at target point
            if (!isnan(this.target.x) && !isnan(this.target.y)) {
                pair target_pt = t * this.target;
                if (target_pt != pos_pt) {
                    pen linePen = (p != nullpen) ? p : this.getConnectorPen();
                    pen dotPen = (this.border != nullpen) ? this.border : linePen;
                    pair intersection_pt = point(box_frame, unitOr(target_pt - pos_pt));
                    draw(f, target_pt -- intersection_pt, linePen);
                    dot(f, target_pt, dotPen);
                }
            }
            // Stamp boxed inset on canvas
            add(f, box_frame);
        });
    }
}

// Polarization Inset Factory (From raw Polarization state)
Inset PolarizationInset(Polarization pol,
                        pair pos = (nan, nan),
                        pair target = (nan, nan),
                        real size = 20.0,
                        real refAmplitude = 1.0,
                        string title = "",
                        string tag = "",
                        pen p = currentpen,
                        pen border = rgb(0.0, 0.2, 0.6),
                        pen fill = white,
                        pen badgeText = white) {
    Inset ins;
    ins.pos = pos;
    ins.target = target;
    ins.size = size;
    ins.title = title;
    ins.tag = tag;
    ins.border = border;
    ins.fill = fill;
    ins.badgeText = badgeText;
    real r = size;

    // Reference Axes
    draw(ins.pic, (-r, 0) -- (r, 0), gray(0.8), Arrow(HookHead, 2pt));
    draw(ins.pic, (0, -r) -- (0, r), gray(0.8), Arrow(HookHead, 2pt));
    label(ins.pic, "\tiny $S$", (r, 0), 2S + W, gray(0.4));
    label(ins.pic, "\tiny $P$", (0, r), S + 2W, gray(0.4));

    // Title
    if (title != "") label(ins.pic, Label("\bfseries{\strut\tiny " + title + "}"), (0, r + 2), N);

    real localI = pol.intensity();
    real orbitR = 0.75 * r;

    // Guaranteed auto-scale: if refAmplitude is nan, scale by sqrt(localI)
    real scaleNorm = (!isnan(refAmplitude) && refAmplitude > 0) ? refAmplitude : sqrt(max(config.math.eps, localI));
    real displayFactor = orbitR / scaleNorm;

    real dop = pol.dop();
    real I_pol = dop * localI;
    real I_unpol = (1.0 - dop) * localI;

    real minArrowSize = config.rendering.polarizationInsetMinArrowFraction * orbitR;

    // =========================================================================
    // 1. STARBURST SPOKES (Unpolarized & Partially Polarized)
    // =========================================================================
    // Only drawn when there is a significant unpolarized fraction (DoP < 95%)
    if (dop < 0.95 && localI > config.math.eps) {
        real unpolR = displayFactor * sqrt(localI);
        int spokes = config.rendering.polarizationInsetUnpolarizedSpokes;
        if (spokes <= 0) spokes = 8;

        real psi = pol.psi();
        real minSpokeAlpha = 0.05; // Faint floor for orthogonal/unpolarized directions
        real maxSpokeAlpha = 0.95; // Bold ceiling for dominant polarization direction

        for (int k = 0; k < spokes; ++k) {
            real ang = radians(k * 180.0 / spokes);
            pair dir = (cos(ang), sin(ang)) * unpolR * 0.85;

            // Alignment factor along the spoke angle: 1.0 = aligned with psi, 0.0 = perpendicular
            real cosDiff = cos(ang - psi);
            real alignWeight = cosDiff * cosDiff;

            // Physical opacity modulation:
            // - Pure unpolarized (dop = 0): all spokes have equal alpha (0.70)
            // - Partially polarized (dop > 0): aligned spoke stays bold, orthogonal fades to minSpokeAlpha
            real spokeFactor = (1.0 - dop) * 0.75 + dop * alignWeight;
            real spokeAlpha = minSpokeAlpha + (maxSpokeAlpha - minSpokeAlpha) * spokeFactor;
            spokeAlpha = clamp(spokeAlpha, minSpokeAlpha, 1.0);

            // Give dominant spokes slightly more line weight for clarity
            // real spokeWidth = (dop > 0.2 && alignWeight > 0.7) ? 0.6bp : 0.4bp;
            // pen spokePen = p + opacity(spokeAlpha);
            real s = 1.0 - spokeAlpha;
            pen spokePen = p + gray(s);

            if (length(dir) > minArrowSize) {
                draw(ins.pic, (0, 0) -- - dir, spokePen, Arrow(HookHead, 2.0pt), BeginPenMargin);
                draw(ins.pic, (0, 0) -- dir, spokePen, Arrow(HookHead, 2.0pt), BeginPenMargin);
            } else
                draw(ins.pic, -dir -- dir, spokePen);
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
    bool hasEllipticity = (abs(pol.chi()) >= radians(1.5));
    bool drawPolarizedShape = (dop >= 0.95) || (hasEllipticity && I_pol > config.math.eps);

    if (drawPolarizedShape && I_pol > config.math.eps) {
        real psi = pol.psi();
        real chi = pol.chi();
        real a = displayFactor * sqrt(I_pol) * cos(chi);
        real b = displayFactor * sqrt(I_pol) * abs(sin(chi));

        // Linear polarization (only runs for pure polarized states):
        if (!hasEllipticity) {
            pair tip = a * (cos(psi), sin(psi));
            if (a > minArrowSize)
                draw(ins.pic, -tip -- tip, p, Arrows(HookHead, 2.5pt));
            else
                draw(ins.pic, -tip -- tip, p);
        } else {
            // Elliptical / Circular:
            guide ellipsePath;
            int steps = config.rendering.polarizationInsetEllipseSteps;
            for (int k = 0; k <= steps; ++k) {
                real tau = k * 2.0 * pi / steps;
                real x_loc = a * cos(tau);
                real y_loc = (pol.S3 >= 0 ? 1.0 : -1.0) * b * sin(tau);
                real S_t = x_loc * cos(psi) - y_loc * sin(psi);
                real P_t = x_loc * sin(psi) + y_loc * cos(psi);
                ellipsePath = ellipsePath .. (S_t, P_t);
            }
            if (a > minArrowSize)
                draw(ins.pic, ellipsePath, p, Arrow(HookHead(), 2.5pt, Relative(0.85)));
            else
                draw(ins.pic, ellipsePath, p);
        }
    }

    // =========================================================================
    // 3. READOUT LABEL
    // =========================================================================
    object[] lines;
    lines.push(object(Label(format("\tiny $I = %#.3f$", localI))));
    if (dop >= config.rendering.polarizationInsetUnpolarizedThreshold && dop < 0.999) {
        lines.push(object(Label(format("\tiny $\mathrm{DoP} = %.0f\%$", dop * 100.0))));
    } else if (dop < config.rendering.polarizationInsetUnpolarizedThreshold) {
        lines.push(object(Label("\tiny (Unpolarized)")));
    }

    frame fReadout = pack(align = S, ... lines);
    add(ins.pic, fReadout, (0, -r - 2), S);
    return ins;
}
// Trace Overload: Automatically extracts target point from the ray segment
Inset PolarizationInset(Trace tr,
                        int segmentIndex = 0,
                        real relative = nan,
                        real absolute = nan,
                        pair pos = (nan, nan),
                        real size = 20.0,
                        real refAmplitude = nan,
                        string title = "",
                        string tag = "",
                        pen p = currentpen,
                        pen border = rgb(0.0, 0.2, 0.6),
                        pen fill = white,
                        pen badgeText = white) {
    Segment[] segs = tr.segments();
    if (segs.length == 0) return new Inset;
    int idx = (segmentIndex < 0) ? segmentIndex + segs.length : segmentIndex;
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
    pair targetPoint = (L > 0) ? point(seg.traj, arctime(seg.traj, s)) : seg.a;
    Polarization localPol = resolvePolarizationAt(seg.before, s);
    real refAmp = refAmplitude;
    if (isnan(refAmp) && tr.first() != null) {
        real startI = tr.first().intensity();
        refAmp = (startI > config.math.eps) ? sqrt(startI) : 1.0;
    }
    return PolarizationInset(localPol,
                             pos = pos,
                             target = targetPoint,
                             size = size,
                             refAmplitude = refAmp,
                             title = title,
                             tag = tag,
                             p = p,
                             border = border,
                             fill = fill,
                             badgeText = badgeText);
}

// ============================================================================
// Clean Auto-Grid Inset Panel
// ============================================================================
struct InsetPanel {
    Inset[] insets;
    int perrow = 1;     // Grid columns (1 = vertical list, 2 = 2-col table)
    real margin = 4.0;  // Panel padding in pt
    real spacing = 8.0; // Spacing between cells in pt
    // Panel Global Styling
    pen border = nullpen;   // Panel outline (nullpen = defaults to insets[0].border)
    pen fill = nullpen;     // Panel background fill (nullpen = defaults to insets[0].fill)
    pen badgeText = white;  // Default badge text color
    bool drawBadges = true; // Draw badges and connectors on diagram
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
                   int segmentIndex = 0,
                   real relative = nan,
                   real absolute = nan,
                   pair pos = (nan, nan),
                   real size = 20.0,
                   real refAmplitude = 1.0,
                   string tag = "",
                   string title = "",
                   pen p = currentpen,
                   pen border = rgb(0.0, 0.2, 0.6),
                   pen fill = white,
                   pen badgeText = white) {
        Inset ins = PolarizationInset(tr,
                                      segmentIndex,
                                      relative,
                                      absolute,
                                      pos,
                                      size,
                                      refAmplitude,
                                      title,
                                      tag,
                                      p,
                                      border,
                                      fill,
                                      badgeText);
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
            real[] colW = array(cols, 0.0);
            real[] rowH = array(rows, 0.0);
            frame[] cells = new frame[n];
            // [Section 1] Build Individual Legend Cells (Badge + Inset Content)
            for (int i = 0; i < n; ++i) {
                frame cf;
                Inset ins = this.insets[i];
                // Badge inside the panel: uses inset's badgePen/border & badgeText
                pen bPen = ins.getBadgePen();
                pen bTxt = (ins.badgeText != nullpen) ? ins.badgeText : this.badgeText;
                frame bf;
                label(bf, Label("\tiny\bfseries " + ins.tag, bTxt));
                roundbox(bf, bf, 0.0, 0.0, bPen, Fill(bPen), above = false);
                add(cf, bf);
                // Inset content next to badge
                frame content = ins.pic.fit();
                add(cf, shift(max(bf).x + 4.0 - min(content).x, min(content).y) * content);
                cells[i] = cf;
                int r = floor(i / cols), c = i % cols;
                pair sz = max(cf) - min(cf);
                if (sz.x > colW[c]) colW[c] = sz.x;
                if (sz.y > rowH[r]) rowH[r] = sz.y;
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
                    x_curr += colW[c] + this.spacing;
                }
                y_curr -= (rowH[r] + this.spacing);
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
            if (this.drawBadges) {
                for (Inset ins : this.insets) {
                    if (!isnan(ins.target.x) && !isnan(ins.target.y)) {
                        pair target_pt = t * ins.target;
                        pair badge_pt = (!isnan(ins.pos.x) && ins.pos != ins.target) ? t * ins.pos : target_pt;
                        pen bPen = ins.getBadgePen();
                        pen bTxt = (ins.badgeText != nullpen) ? ins.badgeText : this.badgeText;
                        // Draw on-diagram badge box
                        frame b_frame;
                        label(b_frame, Label("\tiny\bfseries " + ins.tag, bTxt));
                        roundbox(b_frame, b_frame, 0.0, 0.0, bPen, Fill(bPen), above = false);
                        b_frame = shift(badge_pt) * b_frame;
                        // If badge is offset from target, draw connector line and target dot
                        if (badge_pt != target_pt) {
                            pen linePen = (p != nullpen) ? p : ins.getConnectorPen();
                            pen dotPen = (ins.border != nullpen) ? ins.border : linePen;
                            pair intersection_pt = point(b_frame, unitOr(target_pt - badge_pt));
                            draw(f, target_pt -- intersection_pt, linePen);
                            dot(f, target_pt, dotPen);
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

// 3. Transformation Operator (*)
// Usage: OpticalSystem shiftedSys = shift(10, 0) * rotate(45) * sys;
OpticalSystem operator*(transform t, OpticalSystem sys) {
    OpticalSystem out = new OpticalSystem;
    out.ambient = sys.ambient;
    for (int i = 0; i < sys.elements.length; ++i) {
        out.add(t * sys.elements[i]);
    }
    return out;
}

// 4. Angle-to-Polarization Cast Operator
// Usage: Ray r = Ray(p, d, pol = 45); // Automatically creates linear polarization at 45 degrees
Polarization operator cast(real angle_deg) { return Polarization.linear(angle_deg); }

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

            string sL = L;
            if (sL == "") sL = format("%#.2f", D);

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
            pair pA = rotated ? (0, 0) : A_dim;
            pair pB = rotated ? (dist, 0) : B_dim;

            picture tmp;

            // Build the Label using only the function's explicit arguments
            Label lab = Label(sL, align = align, p = textpen, filltype = fill);

            draw(tmp, lab, pA -- pB, linepen + solid, arrow = arrow, bar = bar, margin = margin);

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
