import optics;

// 1. signedArea on 2-node closed curved lens (was 0.0)
path two_node_lens = (0, -10) .. controls(-3, -5) and
                     (-3, 5) .. (0, 10) .. controls(3, 5) and (3, -5) .. cycle;
assert(abs(signedArea(two_node_lens)) > 30.0, "2-node curved lens area must be nonzero");

// 2. Grating align() aims +1 order at target
Element gr = Grating((0, 0), 20.0, 0.0, 500.0, reflective = false, orders = new int[]{1});
gr = align(gr, (-10, 0), (10, 10 * tan(radians(30))), lambda = 633.0);
Trace tg = system(gr).trace(Ray((-10, 0), (1, 0), lambda = 633.0));
pair dir_out = tg.last().d;
pair dir_expected = unit((10, 10 * tan(radians(30))));
assert(length(dir_out - dir_expected) < 1e-4, "Aligned grating must hit target direction");

// 3. Rotated SphericalInterface (angle = 90 deg)
Element si = SphericalInterface((0, 0), 90.0, 20.0, 15.0, Materials.BK7, Materials.AIR);
Trace ts = system(si, Materials.AIR).trace(Ray((2.0, -10.0), (0, 1)));
assert(ts.last().d.x < -0.01, "Rotated SphericalInterface must bend toward its local axis");

// 4. InterferenceInset collimated in non-vacuum medium
Ray[] col_rays = {Ray((0, -1), dir(10), medium = Materials.WATER), Ray((0, 1), dir(-10), medium = Materials.WATER)};
Inset ins = InterferenceInset(col_rays, screen_center = (10, 0), collimated = true);
assert(!ins.pic.empty(), "InterferenceInset collimated should render");

write("All P1 checks passed.");