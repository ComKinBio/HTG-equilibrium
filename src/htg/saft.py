"""Selection of PC-SAFT / ePC-SAFT parameter sets and SUPCRT HKF records.

Molecules: the first complete ePC-SAFT set (PROJECT_EPCSAFT.csv) is used,
otherwise the first complete PC-SAFT set (PROJECT_PCSAFT.csv, model priority
'PC-SAFT' then 'GC-PPC-SAFT'). Water's ePC-SAFT sigma is temperature dependent:
sigma(T) = 2.7927 + 10.11 exp(-0.01775 T) - 1.417 exp(-0.01146 T) Angstrom, T in K.

Ions: the first complete charged ePC-SAFT set, plus the first complete aqueous
HKF record in supcrt92.xml.
"""
import numpy as np
import pandas as pd

from .databases import normalize_name

EPC_PARAMETER_ALIASES = {
    'm': ['mseg', 'm', 'segment'],
    'sigma': ['sigma', 'sigma_a', 'sigma (a)', 'sigma (å)', 'σ'],
    'epsilon': ['u/kb', 'epsilon/k', 'epsilon/kb', 'epsilon_k', 'u_over_kb_k'],
    'epsilonAB': ['epsilonab/kb', 'association_energy_epsilonab/k', 'epsilonab_over_kb_k',
                  'epsilonab_k'],
    'kappaAB': ['kappaab', 'association_volume_kappaab'],
    'association_scheme': ['association_scheme'],
}
ION_EPC_PARAMETER_ALIASES = {
    'm': ['mseg', 'm', 'segment'],
    'sigma': ['sigma', 'sigma_a', 'sigma (a)', 'sigma (å)'],
    'epsilon': ['u/kb', 'epsilon/k', 'epsilon/kb', 'epsilon_k'],
    'dielectric_constant': ['dielectric_constant'],
}
PCSAFT_MODEL_PRIORITY = ['pc-saft', 'gc-ppc-saft']
HKF_REQUIRED_FIELDS = ['Gf', 'Hf', 'Sr', 'a1', 'a2', 'a3', 'a4', 'c1', 'c2', 'wref']


def water_sigma(T):
    """Temperature-dependent ePC-SAFT segment diameter of water [Angstrom]."""
    return float(2.7927 + 10.11 * np.exp(-0.01775 * T) - 1.417 * np.exp(-0.01146 * T))


# ---------------------------------------------------------------- molecules
def _evaluate_epc_value(species_name, quantity, raw, T):
    if pd.isna(raw):
        return np.nan
    if quantity == 'association_scheme':
        return str(raw).strip()
    try:
        return float(raw)
    except (TypeError, ValueError):
        pass
    text = str(raw).strip()
    if species_name == 'water' and quantity == 'sigma' and text.casefold().startswith('t-dependent'):
        return water_sigma(T)
    raise ValueError(f"Cannot evaluate ePC-SAFT parameter '{quantity}' for '{species_name}'.\n"
                     f"Raw value: {raw!r}")


def _extract_epc_quantity(rows, species_name, quantity, T, required=True):
    aliases = {normalize_name(a) for a in EPC_PARAMETER_ALIASES[quantity]}
    matches = rows[rows['Parameter'].map(normalize_name).isin(aliases)]
    if matches.empty:
        return None if required else np.nan
    try:
        return _evaluate_epc_value(species_name, quantity, matches.iloc[0]['Value'], T)
    except Exception:
        return None if required else np.nan


def species_aliases(species_name, aliases, formula):
    """Normalized names that identify a species in the SAFT CSV Species column."""
    names = {normalize_name(v) for v in (aliases if aliases is not None else [species_name])}
    if formula is not None:
        names.add(normalize_name(formula))
    return names


def first_complete_epcsaft_set(epcsaft, species_name, formula, aliases, T):
    """First (Model, RefID) group with finite m, sigma, epsilon; None if there is none."""
    rows = epcsaft[epcsaft['Species'].map(normalize_name).isin(
        species_aliases(species_name, aliases, formula))].copy()
    if rows.empty:
        return None
    if 'record_type' in rows.columns:
        parameter_rows = rows[rows['record_type'].astype(str).str.strip().str.casefold()
                              .eq('species_parameter')].copy()
        if not parameter_rows.empty:
            rows = parameter_rows
    rows = rows.sort_index()
    grouping = [col for col in ('Model', 'RefID') if col in rows.columns]
    if grouping:
        groups = []
        for _, key_row in rows[grouping].drop_duplicates().iterrows():
            mask = np.ones(len(rows), dtype=bool)
            for col in grouping:
                value = key_row[col]
                if pd.isna(value):
                    mask &= rows[col].isna()
                else:
                    mask &= rows[col].astype(str).eq(str(value))
            groups.append(rows.loc[mask].copy())
    else:
        groups = [rows]
    for group in groups:
        m_i = _extract_epc_quantity(group, species_name, 'm', T)
        sigma_i = _extract_epc_quantity(group, species_name, 'sigma', T)
        epsilon_i = _extract_epc_quantity(group, species_name, 'epsilon', T)
        if m_i is None or sigma_i is None or epsilon_i is None:
            continue
        if not np.all(np.isfinite([m_i, sigma_i, epsilon_i])):
            continue
        epsilonAB_i = _extract_epc_quantity(group, species_name, 'epsilonAB', T, required=False)
        kappaAB_i = _extract_epc_quantity(group, species_name, 'kappaAB', T, required=False)
        association_i = _extract_epc_quantity(group, species_name, 'association_scheme', T,
                                              required=False)
        if formula is None:
            raise KeyError(species_name)
        first_row = group.iloc[0]
        return {
            'Species': species_name, 'Formula': formula, 'Source': 'PROJECT_EPCSAFT.csv',
            'Database row': int(group.index.min()),
            'Model': str(first_row['Model']) if 'Model' in group.columns else 'ePC-SAFT',
            'RefID': (str(first_row['RefID'])
                      if 'RefID' in group.columns and pd.notna(first_row['RefID']) else ''),
            'm': float(m_i), 'sigma': float(sigma_i), 'epsilon': float(epsilon_i),
            'epsilonAB': float(epsilonAB_i) if pd.notna(epsilonAB_i) else np.nan,
            'kappaAB': float(kappaAB_i) if pd.notna(kappaAB_i) else np.nan,
            'association_scheme': str(association_i).strip() if pd.notna(association_i) else '-',
        }
    return None


def _optional_numeric(row, column):
    if column is None:
        return np.nan
    value = row[column]
    if pd.isna(value):
        return np.nan
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def first_complete_pcsaft_set(databases, species_name, formula):
    """First positive, finite PC-SAFT row matching ``formula``; None if there is none."""
    if formula is None:
        raise KeyError(f"No PC-SAFT formula mapping for '{species_name}'.")
    pcsaft, cols = databases.pcsaft, databases.pc_columns
    rows = pcsaft[pcsaft['Formula'].map(normalize_name).eq(normalize_name(formula))].copy()
    if rows.empty:
        return None
    rows = rows.sort_index()
    for preferred in PCSAFT_MODEL_PRIORITY:
        model_rows = rows[rows['Model'].astype(str).str.strip().str.casefold().eq(preferred)].copy()
        for idx, row in model_rows.iterrows():
            try:
                m_i = float(row[cols['m']])
                sigma_i = float(row[cols['sigma']])
                epsilon_i = float(row[cols['epsilon']])
            except (TypeError, ValueError, KeyError):
                continue
            if not np.all(np.isfinite([m_i, sigma_i, epsilon_i])):
                continue
            if m_i <= 0.0 or sigma_i <= 0.0 or epsilon_i <= 0.0:
                continue
            association_i = '-'
            if cols['association'] is not None and pd.notna(row[cols['association']]):
                association_i = str(row[cols['association']]).strip()
            return {
                'Species': species_name, 'Formula': formula, 'Source': 'PROJECT_PCSAFT.csv',
                'Database row': int(idx), 'Model': str(row['Model']),
                'RefID': str(row['RefID']) if 'RefID' in row.index and pd.notna(row['RefID']) else '',
                'm': m_i, 'sigma': sigma_i, 'epsilon': epsilon_i,
                'epsilonAB': _optional_numeric(row, cols['epsilonAB']),
                'kappaAB': _optional_numeric(row, cols['kappaAB']),
                'association_scheme': association_i,
            }
    return None


# ---------------------------------------------------------------- ions
def normalize_ion_name(value):
    if pd.isna(value):
        return ''
    return (str(value).strip().casefold().replace(' ', '').replace('₂', '2').replace('₃', '3')
            .replace('⁺', '+').replace('⁻', '-'))


def _numeric_epc_value(value):
    if pd.isna(value):
        return np.nan
    try:
        return float(value)
    except (TypeError, ValueError):
        raise ValueError(f'Could not convert ePC-SAFT value {value!r} to a number.') from None


def _first_epc_parameter(rows, quantity, required=True):
    for alias in (str(x).strip().casefold() for x in ION_EPC_PARAMETER_ALIASES[quantity]):
        match = rows[rows['Parameter'].astype(str).str.strip().str.casefold().eq(alias)]
        if not match.empty:
            return _numeric_epc_value(match.iloc[0]['Value'])
    if required:
        available = rows['Parameter'].dropna().astype(str).unique().tolist()
        raise ValueError(f"Required ePC-SAFT parameter '{quantity}' not found.\n"
                         f"Available parameters: {available}")
    return np.nan


def first_complete_ion_epcsaft_set(epcsaft, ion_name, saft_formula, charge):
    aliases = {normalize_ion_name(ion_name), normalize_ion_name(saft_formula)}
    rows = epcsaft[epcsaft['Species'].map(normalize_ion_name).isin(aliases)].copy()
    if rows.empty:
        return None
    rows = rows.sort_index()
    if 'record_type' in rows.columns:
        species_rows = rows[rows['record_type'].astype(str).str.strip().str.casefold()
                            .eq('species_parameter')].copy()
        if not species_rows.empty:
            rows = species_rows
    grouping = [col for col in ['RefID', 'Model'] if col in rows.columns]
    if grouping:
        parameter_sets = [g for _, g in rows.groupby(grouping, sort=False, dropna=False)]
    else:
        parameter_sets = [rows]
    for parameter_set in parameter_sets:
        try:
            m_i = _first_epc_parameter(parameter_set, 'm')
            sigma_i = _first_epc_parameter(parameter_set, 'sigma')
            epsilon_i = _first_epc_parameter(parameter_set, 'epsilon')
        except ValueError:
            continue
        if not np.all(np.isfinite([m_i, sigma_i, epsilon_i])):
            continue
        if m_i <= 0.0 or sigma_i <= 0.0 or epsilon_i <= 0.0:
            continue
        dielectric_i = _first_epc_parameter(parameter_set, 'dielectric_constant', required=False)
        first_row = parameter_set.iloc[0]
        charge_database = np.nan
        if 'Charge' in parameter_set.columns:
            values = parameter_set['Charge'].dropna()
            if not values.empty:
                charge_database = float(values.iloc[0])
        return {
            'Species': ion_name, 'Source': 'PROJECT_EPCSAFT.csv',
            'Database row': int(parameter_set.index.min()),
            'Model': (str(first_row['Model'])
                      if 'Model' in first_row.index and pd.notna(first_row['Model']) else ''),
            'RefID': (str(first_row['RefID'])
                      if 'RefID' in first_row.index and pd.notna(first_row['RefID']) else ''),
            'charge': charge, 'charge_database': charge_database,
            'm': m_i, 'sigma': sigma_i, 'epsilon': epsilon_i, 'dielectric_constant': dielectric_i,
        }
    return None


def supcrt_species_name(species_element):
    name = species_element.find('Name')
    if name is None or name.text is None:
        return ''
    return name.text.strip()


def read_hkf_candidate(species_element):
    """Charge and the aqueous HKF fields of one SUPCRT record; None if incomplete."""
    charge_element = species_element.find('Charge')
    if charge_element is None or charge_element.text is None:
        return None
    try:
        charge = float(charge_element.text.strip())
    except ValueError:
        return None
    hkf = species_element.find('./Thermo/HKF')
    if hkf is None:
        return None
    result = {'Charge': charge}
    for name in HKF_REQUIRED_FIELDS:
        child = hkf.find(name)
        if child is None or child.text is None or not child.text.strip():
            return None
        try:
            result[name] = float(child.text.strip())
        except ValueError:
            return None
    return result


def first_complete_hkf_set(supcrt_root, ion_name, supcrt_aliases):
    aliases = {normalize_ion_name(alias) for alias in supcrt_aliases}
    candidates = [node for node in supcrt_root.iter('Species')
                  if normalize_ion_name(supcrt_species_name(node)) in aliases]
    for number, candidate in enumerate(candidates):
        values = read_hkf_candidate(candidate)
        if values is None:
            continue
        result = {'Species': ion_name, 'SUPCRT name': supcrt_species_name(candidate),
                  'Candidate': number}
        result.update(values)
        for tag in ['Formula', 'MolarMass', 'Type', 'Date', 'References']:
            child = candidate.find(f'.//{tag}')
            result[tag] = child.text.strip() if child is not None and child.text is not None else ''
        return result
    return None
