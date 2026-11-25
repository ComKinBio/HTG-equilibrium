"""
HTGEquilibrium: the science core of the HTG equilibrium app in Julia (Phase 2,
milestone M2). Python module names are kept (htg/constants.py -> constants.jl, ...)
so the two versions are easy to compare.
"""
module HTGEquilibrium

using Clapeyron
using CSV
using EzXML
using ForwardDiff
using JSON
using Libdl
using LinearAlgebra: BLAS, ColumnNorm, diag, pinv, qr
using OrderedCollections
using Printf
using TOML
import HiGHS
import Ipopt
import NLopt
import Random

export FormulaError, parse_formula, split_charge, molar_mass, formula_summary, hill_formula,
       parse_element_list, is_plain_formula
export WaterModel, SUPCRT92Water, IAPWS95Water, WaterProperties, water_properties, water_dict
export HKFRecord, hkf_ion_g0, hkf_apparent_G
export SUPCRT92Band, supcrt92_failed, supcrt92_failure_bands, supcrt92_band_check
export ValueError, py_repr, parse_toml, read_toml, Units
export NumberListError, parse_numbers, format_numbers, parse_points
export Settings, SettingsKeyError, SettingsValueError, SETTINGS_DEFAULTS, EXPORT_LAYOUTS, EXPORT_FORMATS, load_settings,
       settings_path, export_defaults, set_export_defaults
export Salt, Feed, Case, CaseError, LegacyConversionError, load_case, case_from_dict, validate,
       temperatures, pressures, grid, chains, at, is_single_point, resolved, total_water_g_per_g_dry,
       dry_element_wt_percent, carbon_bound_messages, effective_carbon_bounds, check_carbon_bounds,
       legacy_carbon_bounds, legacy_to_v2, is_v2, fluid_phase_count, starts_from_multistart,
       BOUND_KINDS, PHASE_LABELS, WATER_MODELS, water_model_of, dump_case, to_toml
export FeedError, feed_v2_to_model, feed_summary, feed_fractions, molar_mass_function,
       species_elements_function, registry_elements, elemental_masses, ash_percent, scaled_elements, legacy_feed_to_v2
export Databases, DataTable, XMLNode, NASARecord, load_databases, normalize_name, nasa7_properties,
       nasa7_coefficients, read_nasa_lines, read_nasa_database, parse_nasa_header_elements
export water_sigma, first_complete_epcsaft_set, first_complete_pcsaft_set, first_complete_ion_epcsaft_set,
       first_complete_hkf_set, normalize_ion_name, read_hkf_candidate
export Registry, CatalogueRow, NO_SAFT_REASON, discover!, all_molecular, nasa_name, saft_parameters,
       resolve_parameters, available, formula_elements, molecular_weight
export UserDatabase, load_user_database, resolve_user_dir, default_user_dir, apply_user_layer!, hide_species!
export nasa_g0, supcrt_mineral_properties, solid_g0, hkf_records
export Model, load_model, default_data_dir, species, mineral_inventory, product_candidates
export EOSFailure, PhaseError, Phase, SolidPhase, FluidPhase, ElectrolytePhase, association_sites,
       build_pcsaft_model, gibbs, gibbs_and_mu, mu_at, lnphi, hessian, gibbs_full
export Alternative, Result, n_results, notes, for_alternative, species_table, phase_summary, SweepResult, results,
       alternatives_table, read_choices, write_choices, apply_choices, CHOICES_FILE, ALTERNATIVE_COLUMNS,
       RunTables, with_choices, flagged_points, load_saved
export read_csv, to_csv, csv_text, concat, setcolumn!, insert_column!, take_rows, nrows
export run_sweep, sweep, preflight, preflight_messages, check, point_text
export WorkerPool, Run, solve_chain, default_workers, submit, cancel, set_workers, restart, status, shutdown,
       sweep_result, live_tables, load_model_from_args
export EquilibriumError, RuntimeError, GibbsProblem, OptimizeResult, fun_and_grad, hessian_blocks,
       check_against_model, kkt, linprog, starting_points, candidate_starts, result_label, group_solutions,
       tangent_plane_test, optimize, BACKENDS, SLSQP_EVALUATIONS_PER_ITERATION, POLISH_THRESHOLD_RT
export SolveError, Stage, Setup, setup, solve, species_to_wt_percent, apply_phase_rules, carbon_bound_rows,
       carbon_bound_violation

include("constants.jl")
include("pycompat.jl")
include("toml_ordered.jl")
include("formulas.jl")
include("water.jl")
include("hkf.jl")
include("supcrt92_bands.jl")
include("units.jl")
include("number_lists.jl")
include("settings.jl")
include("case.jl")
include("feed_v2.jl")
include("databases.jl")
include("tables.jl")
include("saft.jl")
include("registry.jl")
include("user_species.jl")
include("standard_states.jl")
include("model.jl")
include("phases.jl")
include("results.jl")
include("optimizer.jl")
include("equilibrium.jl")
include("sweep.jl")
include("queue.jl")
include("cli.jl")

end # module
