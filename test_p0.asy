import optics;

// 1. Polarization.elliptical degree/radian fix
Polarization pe = Polarization.elliptical(1.0, 0.5, theta = 30.0);
assert(pe.isPhysical(), "elliptical Stokes vector must be physical");
assert(abs(pe.dop() - 1.0) < 1e-6, "elliptical DoP must be 1");

// 2. applyDiattenuation cross-term sign matches Polarization.jones
Polarization p45 = Polarization.linear(45.0);
Polarization qwp = copy(p45);
applyDiattenuation(qwp, p45, (1.0, 0.0), (0.0, 1.0));
Polarization j_ref = Polarization.jones((1.0 / sqrt(2), 0.0), (0.0, 1.0 / sqrt(2)));
assert(abs(qwp.S3 - j_ref.S3) < 1e-6, "applyDiattenuation S3 sign mismatch");

// 3. Oblique metal energy conservation (Rp + Tp == 1)
FresnelAmplitudes fa_au = fresnelAmplitudesKernel(1.0, (0.197, -3.09), cos(radians(45)));
assert(abs((fa_au.Rs + fa_au.Ts) - 1.0) < 1e-6, "Gold Rs+Ts != 1");
assert(abs((fa_au.Rp + fa_au.Tp) - 1.0) < 1e-6, "Gold Rp+Tp != 1");

// 4. ThinFilmStack: 0-layer matches Fresnel, thick Gold has R < 1, glass->air 45deg hits TIR
FresnelAmplitudes fa0 =
    thinFilmAmplitudes(new ThinFilmLayer[]{}, Materials.AIR, Materials.BK7, cos(radians(35)), 633.0, (0, 0));
FresnelAmplitudes fref =
    fresnelAmplitudesKernel(Materials.AIR.n(lambda = 633), Materials.BK7.n(lambda = 633), cos(radians(35)));
assert(abs(fa0.rs - fref.rs) < 1e-6 && abs(fa0.rp - fref.rp) < 1e-6, "0-layer TMM != Fresnel r");
assert(abs(fa0.ts - fref.ts) < 1e-6 && abs(fa0.tp - fref.tp) < 1e-6, "0-layer TMM != Fresnel t");

ThinFilmLayer[] thick_au = {ThinFilmLayer(Materials.GOLD, 5000.0)};
FresnelAmplitudes fa_thick =
    thinFilmAmplitudes(thick_au, Materials.AIR, Materials.BK7, cos(radians(45)), 633.0, (0, 0));
assert(fa_thick.Rs < 1.0 && fa_thick.Rp < 1.0, "Thick Gold film must have R < 1");

FresnelAmplitudes fa_tir =
    thinFilmAmplitudes(new ThinFilmLayer[]{}, Materials.BK7, Materials.AIR, cos(radians(45)), 633.0, (0, 0));
assert(fa_tir.evanescent &&abs(fa_tir.Rs - 1.0) < 1e-6 && abs(fa_tir.Rp - 1.0) < 1e-6, "Glass->Air 45deg must hit TIR");

write("All P0 checks passed.");