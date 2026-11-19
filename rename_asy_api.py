#!/usr/bin/env python3
"""Migrate Asymptote API identifiers to the project's naming convention.

Usage:
    python3 rename_asy_api.py optics.asy optics_test.asy
    python3 rename_asy_api.py --dry-run examples/*.asy
    python3 rename_asy_api.py --recursive --dry-run test
    python3 rename_asy_api.py --recursive --normalize-struct-fields --normalize-local-values --backup-dir bak optics.asy optics_test.asy test

The project naming convention is:

* PascalCase for structs and other types;
* lowerCamelCase for functions and methods;
* snake_case for variables, parameters, and fields;
* UPPER_SNAKE_CASE for constants.

The predefined color helpers (``Red``, ``DarkBlue``, etc.) are intentionally
PascalCase, matching Asymptote's color and drawing helper style.

Renames apply to Asymptote code identifiers only. Comments, strings, filenames,
and standard optical notation such as ``S0``, ``rs``, and ``n1`` are preserved.
The explicit public-name map runs by default. Local variables and struct fields
are normalized only when their options are enabled, so callers can migrate the
API without a broad style rewrite.
"""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path


RENAMES = {
    # Public API parameters and fields: camelCase -> snake_case.
    "lambdaNm": "lambda_nm",
    "angleDeg": "angle_deg",
    "thetaDeg": "theta_deg",
    "psiDeg": "psi_deg",
    "chiDeg": "chi_deg",
    "phaseDeg": "phase_deg",
    "axisDeg": "axis_deg",
    "fastAxisDeg": "fast_axis_deg",
    "rotationDegPerUnit": "rotation_deg_per_unit",
    "linesPerMm": "lines_per_mm",
    "opticAxis": "optic_axis",
    "firstAxis": "first_axis",
    "screenCenter": "screen_center",
    "screenAxis": "screen_axis",
    "laserPen": "laser_pen",
    "badgeText": "badge_text",
    "alphaS": "alpha_s",
    "alphaP": "alpha_p",
    "wavelengthUm": "wavelength_um",
    "visLambda": "vis_lambda",
    "currentDs": "current_ds",
    "minDs": "min_ds",
    "maxDs": "max_ds",
    "visualLambda": "vis_lambda",
    "coherenceId": "coherence_id",
    "arcSubdivisions": "arc_subdivisions",
    "arrowHead": "arrow_head",
    "arrowPolHead": "arrow_pol_head",
    "arrowPolSize": "arrow_pol_size",
    "arrowSize": "arrow_size",
    "axisDir": "axis_dir",
    "axisNorm": "axis_norm",
    "axisOrigin": "axis_origin",
    "badgePen": "badge_pen",
    "badgeText": "badge_text",
    "boxMargin": "box_margin",
    "branchBudget": "branch_budget",
    "childrenIds": "children_ids",
    "clipVelocityEps": "clip_velocity_eps",
    "cutAtBoundary": "cut_at_boundary",
    "dashedBehind": "dashed_behind",
    "defaultChildIndex": "default_child_index",
    "defaultIndex": "default_index",
    "drawAbsorptionFading": "draw_absorption_fading",
    "drawBadges": "draw_badges",
    "drawPolarization": "draw_polarization",
    "drawReflections": "draw_reflections",
    "drawWavefronts": "draw_wavefronts",
    "epsIntegration": "eps_integration",
    "epsIntersection": "eps_intersection",
    "fadeSubpathSteps": "fade_subpath_steps",
    "flatIndex": "flat_index",
    "flatSurfaces": "flat_surfaces",
    "gratingTickAngle": "grating_tick_angle",
    "gratingTickCount": "grating_tick_count",
    "gratingTickLength": "grating_tick_length",
    "hemisphereMargin": "hemisphere_margin",
    "hitSurface": "hit_surface",
    "ignoredFlatIndices": "ignored_flat_indices",
    "imagePath": "image_path",
    "invT": "inv_t",
    "isAdaptive": "is_adaptive",
    "isBirefringent": "is_birefringent",
    "isReflection": "is_reflection",
    "isThick": "is_thick",
    "maxBisectionSteps": "max_bisection_steps",
    "maxBranches": "max_branches",
    "maxDepth": "max_depth",
    "maxOrders": "max_orders",
    "minAmplitude": "min_amplitude",
    "minDrawableLength": "min_drawable_length",
    "minIntensity": "min_intensity",
    "minLensHead": "min_lens_head",
    "nearIntersectionFactor": "near_intersection_factor",
    "parentId": "parent_id",
    "penAbsorber": "pen_absorber",
    "penLens": "pen_lens",
    "penMirror": "pen_mirror",
    "phaseModifier": "phase_modifier",
    "rayToPathExtension": "ray_to_path_extension",
    "rayTree": "ray_tree",
    "reflectedDirection": "reflected_direction",
    "reflectedWaveNormal": "reflected_wave_normal",
    "rootId": "root_id",
    "rootIds": "root_ids",
    "selfIntersectionFactor": "self_intersection_factor",
    "showDots": "show_dots",
    "showInfo": "show_info",
    "showLabels": "show_labels",
    "showP": "show_p",
    "showS": "show_s",
    "showSource": "show_source",
    "solverDenominatorEps": "solver_denominator_eps",
    "sourceMedium": "source_medium",
    "spatialIndex": "spatial_index",
    "targetMedium": "target_medium",
    "titleContentGap": "title_content_gap",
    "traceReflections": "trace_reflections",
    "useIntensityOpacity": "use_intensity_opacity",
    "visualLambda": "vis_lambda",
    "wavefrontMaxDist": "wavefront_max_dist",
    "wavefrontMaxOrders": "wavefront_max_orders",
    "wavevectorAlignmentTolerance": "wavevector_alignment_tolerance",
    "segmentIndex": "segment_index",
    "maxDist": "max_dist",
    "pLine": "p_line",
    "pDots": "p_dots",
    "laser_pen": "pen_laser",
    "arr_pen": "pen_arr",
    "arrow_pen": "pen_arrow",
    "badge_pen": "pen_badge",
    "base_pen": "pen_base",
    "dot_pen": "pen_dot",
    "draw_pen": "pen_draw",
    "eff_pen": "pen_eff",
    "line_pen": "pen_line",
    "local_pen": "pen_local",
    "slice_pen": "pen_slice",
    "spoke_pen": "pen_spoke",
    "sub_pen": "pen_sub",
    "p_line": "pen_line",
    "p_dots": "pen_dots",
    "lambda_nm_eff": "lambda_eff_nm",
    "b_pen": "pen_badge",
    "lambda_nm_test": "lambda_test_nm",
}

FUNCTION_RENAMES = {
    "ray_index": "rayIndex",
    "interface_frame": "interfaceFrame",
    "draw_polarization": "drawPolarization",
    "LowerCase": "lowerCase",
    "UpperCase": "upperCase",
    "RefractiveIndex": "refractiveIndex",
    "GRIN": "grin",
    "Absorption": "absorption",
    "Gain": "gain",
    "Dichroic": "dichroic",
    "Sellmeier": "sellmeier",
    "UniaxialOrdinary": "uniaxialOrdinary",
    "UniaxialExtraordinary": "uniaxialExtraordinary",
    "OpticalActivity": "opticalActivity",
    "FaradayRotator": "faradayRotator",
    "ComplexRefractiveIndex": "complexRefractiveIndex",
    "Refract": "refract",
    "Coating": "coating",
    "PolarizingCoating": "polarizingCoating",
    "ThinFilmStack": "thinFilmStack",
    "SingleLayerCoating": "singleLayerCoating",
    "AntiReflectionCoating": "antiReflectionCoating",
    "DistributedBraggReflector": "distributedBraggReflector",
    "DoubleRefraction": "doubleRefraction",
    "Reflect": "reflect",
    "IdealReflect": "idealReflect",
    "IdealRefract": "idealRefract",
    "Focus": "focus",
    "PhaseFocus": "phaseFocus",
    "IdealFocus": "idealFocus",
    "Diffract": "diffract",
    "PhaseShift": "phaseShift",
    "Absorb": "absorb",
    "Retard": "retard",
    "Interface": "interface",
    "SphericalInterface": "sphericalInterface",
    "Absorber": "absorber",
    "Aperture": "aperture",
    "PhaseGradient": "phaseGradient",
    "Grating": "grating",
    "Crystal": "crystal",
    "UniaxialCrystal": "uniaxialCrystal",
    "PolarizationInset": "polarizationInset",
    "InterferenceInset": "interferenceInset",
    "Transmitted": "transmittedRays",
    "Reflected": "reflectedRays",
}

# Fields whose declarations may be split across lines or whose old spelling
# overlaps a method name.  These are changed only after a member-access dot.
MEMBER_RENAMES = {
    "focalLength": "focal_length",
    "grinSolver": "grin_solver",
    "groupDir": "group_dir",
    "phaseIndex": "phase_index",
    "waveNormal": "wave_normal",
    "waveNormal_o": "wave_normal_o",
    "waveNormal_e": "wave_normal_e",
    "isBirefringent": "is_birefringent",
    "coherenceId": "coherence_id",
}


def code_mask(text: str) -> str:
    """Return same-length text with comments and quoted strings blanked out."""
    chars = list(text)
    i = 0
    state = "code"
    while i < len(text):
        char = text[i]
        nxt = text[i + 1] if i + 1 < len(text) else ""
        if state == "line_comment":
            if char == "\n":
                state = "code"
            elif char != "\r":
                chars[i] = " "
        elif state == "block_comment":
            if char == "*" and nxt == "/":
                chars[i] = chars[i + 1] = " "
                i += 1
                state = "code"
            elif char not in "\r\n":
                chars[i] = " "
        elif state == "string":
            if char == "\\":
                chars[i] = " "
                if i + 1 < len(text) and text[i + 1] not in "\r\n":
                    chars[i + 1] = " "
                    i += 1
            elif char == '"':
                chars[i] = " "
                state = "code"
            elif char not in "\r\n":
                chars[i] = " "
        elif char == "/" and nxt == "/":
            chars[i] = chars[i + 1] = " "
            i += 1
            state = "line_comment"
        elif char == "/" and nxt == "*":
            chars[i] = chars[i + 1] = " "
            i += 1
            state = "block_comment"
        elif char == '"':
            chars[i] = " "
            state = "string"
        i += 1
    return "".join(chars)


def replace_identifiers(
    text: str,
    renames: dict[str, str],
    *,
    calls_only: bool = False,
    members_only: bool = False,
) -> tuple[str, int]:
    """Replace identifier tokens in code while leaving comments and strings intact."""
    mask = code_mask(text)
    edits: list[tuple[int, int, str]] = []
    for match in re.finditer(r"\b[A-Za-z_]\w*\b", mask):
        old = match.group()
        new = renames.get(old)
        if new is None:
            continue
        if calls_only and not re.match(r"\s*\(", mask[match.end():]):
            continue
        if members_only:
            prefix = mask[:match.start()].rstrip()
            if not prefix.endswith("."):
                continue
        edits.append((match.start(), match.end(), new))
    for start, end, new in reversed(edits):
        text = text[:start] + new + text[end:]
    return text, len(edits)


def migrate(text: str, preserve_prefixes: tuple[str, ...] = ()) -> tuple[str, int]:
    renames = {
        old: new for old, new in RENAMES.items()
        if not any(old.startswith(prefix) for prefix in preserve_prefixes)
    }
    return replace_identifiers(text, renames)


def _snake_case(name: str) -> str:
    return re.sub(r"(?<!^)([A-Z])", r"_\1", name).lower()


def struct_field_renames(text: str) -> dict[str, str]:
    """Find lower-camel fields declared directly inside structs.

    Methods are excluded by tracking brace depth.  Uppercase names are left
    alone because they are normally constants or standard optical symbols.
    """
    renames: dict[str, str] = {}
    in_struct = False
    struct_depth = 0
    depth = 0
    for line in code_mask(text).splitlines():
        if re.search(r"\bstruct\s+[A-Za-z][A-Za-z0-9_]*\s*\{", line):
            in_struct = True
            struct_depth = depth + line.count("{") - line.count("}")
        elif in_struct and depth == struct_depth:
            match = re.match(
                r"\s*(?:static\s+)?[A-Za-z][A-Za-z0-9_]*(?:\[\])?\s+([a-z][A-Za-z0-9]*)\s*(?:=|;)",
                line,
            )
            if match and re.search(r"[a-z][A-Z]", match.group(1)):
                old = match.group(1)
                renames[old] = _snake_case(old)
        depth += line.count("{") - line.count("}")
        if in_struct and depth < struct_depth:
            in_struct = False
    return renames


def migrate_struct_fields(
    text: str, preserve_prefixes: tuple[str, ...] = ()
) -> tuple[str, int]:
    renames = {
        old: new for old, new in struct_field_renames(text).items()
        if not any(old.startswith(prefix) for prefix in preserve_prefixes)
    }
    return replace_identifiers(text, renames)


def migrate_member_fields(text: str) -> tuple[str, int]:
    return replace_identifiers(text, MEMBER_RENAMES, members_only=True)


def migrate_functions(text: str) -> tuple[str, int]:
    return replace_identifiers(text, FUNCTION_RENAMES, calls_only=True)


def local_value_renames(text: str) -> dict[str, str]:
    """Collect lowerCamel variable/argument names from assignments/declarations."""
    names: set[str] = set()
    # Named arguments and ordinary assignments.
    code = code_mask(text)
    names.update(re.findall(r"\b([a-z][A-Za-z0-9]*[A-Z][A-Za-z0-9]*)\b(?=\s*=)", code))
    # Typed declarations and parameters of built-in or user-defined types,
    # excluding function declarations by requiring the name not be followed by `(`.
    names.update(
        re.findall(
            r"\b[A-Za-z_]\w*(?:\s*\[\])?\s+"
            r"([a-z][A-Za-z0-9]*[A-Z][A-Za-z0-9]*)\b(?!\s*\()",
            code,
        )
    )
    return {name: _snake_case(name) for name in names if name not in {"assertTrue", "assertFalse"}}


def migrate_local_values(text: str, preserve_prefixes: tuple[str, ...] = ()) -> tuple[str, int]:
    renames = {
        old: new for old, new in local_value_renames(text).items()
        if not any(old.startswith(prefix) for prefix in preserve_prefixes)
    }
    return replace_identifiers(text, renames)


def expand_targets(paths: list[Path], recursive: bool) -> list[Path]:
    """Expand directories to Asymptote files, including git-ignored examples."""
    targets: list[Path] = []
    for path in paths:
        if path.is_dir():
            candidates = path.rglob("*.asy") if recursive else path.glob("*.asy")
            targets.extend(candidate for candidate in candidates if candidate.is_file())
        elif path.is_file():
            targets.append(path)
        else:
            raise FileNotFoundError(f"No such file or directory: {path}")
    return sorted(set(targets))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dry-run", action="store_true", help="report changes without writing files"
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="when a directory is given, scan all nested directories for .asy files",
    )
    parser.add_argument(
        "--preserve-prefix",
        action="append",
        default=[],
        help="do not rename identifiers beginning with this prefix; may be repeated",
    )
    parser.add_argument(
        "--normalize-struct-fields",
        action="store_true",
        help="rename lower-camel struct fields and their references to snake_case",
    )
    parser.add_argument(
        "--normalize-local-values",
        action="store_true",
        help="rename lower-camel local variables and named arguments to snake_case",
    )
    parser.add_argument(
        "--backup-dir",
        type=Path,
        help="copy each file here before modifying it; ignored in dry-run mode",
    )
    parser.add_argument(
        "files", nargs="+", type=Path, help="Asymptote files to migrate"
    )
    args = parser.parse_args()

    try:
        files = expand_targets(args.files, args.recursive)
    except FileNotFoundError as error:
        parser.error(str(error))

    for path in files:
        original = path.read_text()
        migrated, count = migrate(original, tuple(args.preserve_prefix))
        migrated, function_count = migrate_functions(migrated)
        count += function_count
        if args.normalize_struct_fields:
            migrated, field_count = migrate_struct_fields(
                migrated, tuple(args.preserve_prefix)
            )
            count += field_count
            migrated, member_count = migrate_member_fields(migrated)
            count += member_count
        if args.normalize_local_values:
            migrated, local_count = migrate_local_values(migrated, tuple(args.preserve_prefix))
            count += local_count
        if count and not args.dry_run:
            if args.backup_dir:
                try:
                    relative = path.resolve().relative_to(Path.cwd().resolve())
                except ValueError:
                    relative = Path(path.name)
                backup_path = args.backup_dir / relative
                backup_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, backup_path)
            path.write_text(migrated)
        action = "would update" if args.dry_run else "updated"
        print(f"{path}: {action}, {count} replacement(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
