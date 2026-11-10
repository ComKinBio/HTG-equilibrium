"""FROZEN copy of the model code from HTG_Model.ipynb at commit 74ac11d.

The code cells are concatenated verbatim, one "# %% cell N" block per notebook cell.
tests/reference/generate.py runs this file to produce the reference results that
the htg package must reproduce. DO NOT EDIT: the point is that it never changes.
"""

# %% cell 2
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy.linalg import qr
from scipy.optimize import linprog, minimize
R = 8.31446261815324
P0 = 100000.0
MW_WATER_KG = 0.01801528
T = 623.15
P = 25000000.0
MODEL_INITIALIZED = False

# %% cell 4
MOLECULAR_SPECIES = {'water': {'formula': 'H2O', 'nasa': 'H2O', 'clapeyron': 'water', 'elements': {'H': 2, 'O': 1}}, 'hydrogen': {'formula': 'H2', 'nasa': 'H2', 'clapeyron': 'hydrogen', 'elements': {'H': 2}}, 'carbon dioxide': {'formula': 'CO2', 'nasa': 'CO2', 'clapeyron': 'carbon dioxide', 'elements': {'C': 1, 'O': 2}}, 'carbon monoxide': {'formula': 'CO', 'nasa': 'CO', 'clapeyron': 'carbon monoxide', 'elements': {'C': 1, 'O': 1}}, 'methane': {'formula': 'CH4', 'nasa': 'CH4', 'clapeyron': 'methane', 'elements': {'C': 1, 'H': 4}}, 'methanol': {'formula': 'CH3OH', 'nasa': 'CH3OH', 'clapeyron': 'methanol', 'elements': {'C': 1, 'H': 4, 'O': 1}}, 'ammonia': {'formula': 'NH3', 'nasa': 'NH3', 'clapeyron': 'ammonia', 'elements': {'N': 1, 'H': 3}}, 'hydrogen sulfide': {'formula': 'H2S', 'nasa': 'H2S', 'clapeyron': 'hydrogen sulfide', 'elements': {'H': 2, 'S': 1}}}
GAS_SPECIES = ['water', 'hydrogen', 'carbon dioxide', 'carbon monoxide', 'methane', 'methanol', 'ammonia', 'hydrogen sulfide']
ORG_SPECIES = ['benzene', 'toluene']
ORGANIC_SPECIES = {'benzene': {'formula': 'C6H6', 'nasa': 'C6H6', 'clapeyron': 'benzene', 'elements': {'C': 6, 'H': 6}}, 'toluene': {'formula': 'C7H8', 'nasa': 'C7H8', 'clapeyron': 'toluene', 'elements': {'C': 7, 'H': 8}}}
SOLID_SPECIES = ['carbon', 'cellulose', 'Lignin1']
SOLID_DATA = {'carbon': {'formula': 'C', 'nasa': 'CHAR', 'elements': {'C': 1}}, 'cellulose': {'formula': 'C6H10O5', 'nasa': 'CELL', 'elements': {'C': 6, 'H': 10, 'O': 5}}, 'Lignin1': {'formula': 'C22H28O9', 'nasa': 'LIGH', 'elements': {'C': 22, 'H': 28, 'O': 9}}}
AQ_NEUTRAL_SPECIES = GAS_SPECIES.copy()
AQ_ION_SPECIES = ['Na+', 'HCO3-', 'CO3^2-']
ION_DATA = {'Na+': {'formula': 'Na+', 'charge': +1, 'supcrt': 'Na+', 'clapeyron': 'Na+', 'elements': {'Na': 1}}, 'HCO3-': {'formula': 'HCO3-', 'charge': -1, 'supcrt': 'HCO3-', 'clapeyron': 'HCO3-', 'elements': {'H': 1, 'C': 1, 'O': 3}}, 'CO3^2-': {'formula': 'CO3^2-', 'charge': -2, 'supcrt': 'CO3--', 'clapeyron': 'CO32-', 'elements': {'C': 1, 'O': 3}}}
ALL_MOLECULAR_DATA = {**MOLECULAR_SPECIES, **ORGANIC_SPECIES}
MOLECULAR_SPECIES['nitrogen'] = dict(formula='N2', nasa='N2', clapeyron='nitrogen', elements={'N': 2})
ALL_MOLECULAR_DATA = {**MOLECULAR_SPECIES, **ORGANIC_SPECIES}
AW = {'C': 12.011, 'H': 1.008, 'O': 15.999, 'N': 14.007, 'S': 32.06, 'Na': 22.989769}
ATOMIC_WEIGHTS = AW

# %% cell 6
def nasa_float(text):
    text = str(text).strip().replace('D', 'E').replace('d', 'e')
    if not text:
        return np.nan
    return float(text)

def nasa7_coefficients(nasa_db, species_name, T):
    if species_name not in nasa_db:
        raise KeyError(f'NASA species not found: {species_name}')
    data = nasa_db[species_name]
    T_low = data['T_low']
    T_mid = data['T_mid']
    T_high = data['T_high']
    if T < T_low or T > T_high:
        raise ValueError(f'Temperature {T:.2f} K is outside the NASA range for {species_name}: {T_low:.2f}–{T_high:.2f} K.')
    if T <= T_mid:
        coeff = data['low']
        region = 'low'
    else:
        coeff = data['high']
        region = 'high'
    return (coeff, region)

def nasa7_properties(nasa_db, species_name, T):
    T = float(T)
    if T <= 0.0:
        raise ValueError('NASA temperature must be positive.')
    a, region = nasa7_coefficients(nasa_db, species_name, T)
    a1, a2, a3, a4, a5, a6, a7 = a
    cp_R = a1 + a2 * T + a3 * T ** 2 + a4 * T ** 3 + a5 * T ** 4
    h_RT = a1 + a2 * T / 2.0 + a3 * T ** 2 / 3.0 + a4 * T ** 3 / 4.0 + a5 * T ** 4 / 5.0 + a6 / T
    s_R = a1 * np.log(T) + a2 * T + a3 * T ** 2 / 2.0 + a4 * T ** 3 / 3.0 + a5 * T ** 4 / 4.0 + a7
    cp = R * cp_R
    h = R * T * h_RT
    entropy = R * s_R
    g = h - T * entropy
    return (float(cp), float(h), float(entropy), float(g))

def read_nasa_database(path):
    nasa_db = {}
    i = 0
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        nasa_lines = f.readlines()
    while i < len(nasa_lines):
        line1 = nasa_lines[i].rstrip('\n')
        if len(line1) < 80 or line1[79:80] != '1':
            i += 1
            continue
        if i + 3 >= len(nasa_lines):
            break
        line2 = nasa_lines[i + 1].rstrip('\n')
        line3 = nasa_lines[i + 2].rstrip('\n')
        line4 = nasa_lines[i + 3].rstrip('\n')
        if not (len(line2) >= 80 and len(line3) >= 80 and (len(line4) >= 80) and (line2[79:80] == '2') and (line3[79:80] == '3') and (line4[79:80] == '4')):
            i += 1
            continue
        species_name = line1[0:18].strip()
        try:
            T_low = nasa_float(line1[45:55])
            T_high = nasa_float(line1[55:65])
            T_mid = nasa_float(line1[65:73])
        except Exception as err:
            raise ValueError(f"Could not parse NASA temperature limits for species '{species_name}'.\nHeader:\n{line1}") from err
        try:
            c2 = [nasa_float(line2[j:j + 15]) for j in range(0, 75, 15)]
            c3 = [nasa_float(line3[j:j + 15]) for j in range(0, 75, 15)]
            c4 = [nasa_float(line4[j:j + 15]) for j in range(0, 60, 15)]
            high = np.array([c2[0], c2[1], c2[2], c2[3], c2[4], c3[0], c3[1]], dtype=float)
            low = np.array([c3[2], c3[3], c3[4], c4[0], c4[1], c4[2], c4[3]], dtype=float)
        except Exception as err:
            raise ValueError(f"Could not parse NASA coefficients for species '{species_name}'.") from err
        nasa_db[species_name] = {'T_low': float(T_low), 'T_mid': float(T_mid), 'T_high': float(T_high), 'low': low, 'high': high, 'elements': parse_nasa_header_elements(line1), 'phase': line1[44:45].strip().upper()}
        i += 4
    if not nasa_db:
        raise ValueError('No NASA-7 entries found')
    return nasa_db

# %% cell 8
PCSAFT_FORMULA = {'water': 'H2O', 'hydrogen': 'H2', 'carbon dioxide': 'CO2', 'carbon monoxide': 'CO', 'methane': 'CH4', 'methanol': 'CH3OH', 'ammonia': 'NH3', 'hydrogen sulfide': 'H2S', 'benzene': 'C6H6', 'toluene': 'C7H8'}
SPECIES_ALIASES_6 = {'water': ['water', 'h2o'], 'hydrogen': ['hydrogen', 'h2', 'h₂'], 'carbon dioxide': ['carbon dioxide', 'co2', 'co₂'], 'carbon monoxide': ['carbon monoxide', 'co'], 'methane': ['methane', 'ch4', 'ch₄'], 'methanol': ['methanol', 'ch3oh', 'ch₃oh'], 'ammonia': ['ammonia', 'nh3', 'nh₃'], 'hydrogen sulfide': ['hydrogen sulfide', 'hydrogen sulphide', 'h2s', 'h₂s'], 'benzene': ['benzene', 'c6h6', 'c₆h₆'], 'toluene': ['toluene', 'c7h8', 'c₇h₈']}

def normalize_species_name_6(value):
    if pd.isna(value):
        return ''
    text = str(value).strip().casefold()
    replacements = {'₂': '2', '₃': '3', '₄': '4', '₅': '5', '₆': '6', '₇': '7', '₈': '8', '₉': '9'}
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text

def species_aliases_6(species_name):
    aliases = set()
    for value in SPECIES_ALIASES_6.get(species_name, [species_name]):
        aliases.add(normalize_species_name_6(value))
    if species_name in PCSAFT_FORMULA:
        aliases.add(normalize_species_name_6(PCSAFT_FORMULA[species_name]))
    return aliases
EPC_PARAMETER_ALIASES_6 = {'m': ['mseg', 'm', 'segment'], 'sigma': ['sigma', 'sigma_a', 'sigma (a)', 'sigma (å)', 'σ'], 'epsilon': ['u/kb', 'epsilon/k', 'epsilon/kb', 'epsilon_k', 'u_over_kb_k'], 'epsilonAB': ['epsilonab/kb', 'association_energy_epsilonab/k', 'epsilonab_over_kb_k', 'epsilonab_k'], 'kappaAB': ['kappaab', 'association_volume_kappaab'], 'association_scheme': ['association_scheme']}

def evaluate_epc_value_6(species_name, quantity, raw, T):
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
        return float(2.7927 + 10.11 * np.exp(-0.01775 * T) - 1.417 * np.exp(-0.01146 * T))
    raise ValueError(f"Cannot evaluate ePC-SAFT parameter '{quantity}' for '{species_name}'.\nRaw value: {raw!r}")

def extract_epc_quantity_6(rows, species_name, quantity, T, required=True):
    aliases = {normalize_species_name_6(a) for a in EPC_PARAMETER_ALIASES_6[quantity]}
    parameter_names = rows['Parameter'].map(normalize_species_name_6)
    matches = rows[parameter_names.isin(aliases)]
    if matches.empty:
        if required:
            return None
        return np.nan
    raw = matches.iloc[0]['Value']
    try:
        return evaluate_epc_value_6(species_name, quantity, raw, T)
    except Exception:
        if required:
            return None
        return np.nan

def first_complete_epcsaft_set(species_name, T):
    aliases = species_aliases_6(species_name)
    rows = epcsaft_raw[epcsaft_raw['Species'].map(normalize_species_name_6).isin(aliases)].copy()
    if rows.empty:
        return None
    if 'record_type' in rows.columns:
        species_parameter_rows = rows[rows['record_type'].astype(str).str.strip().str.casefold().eq('species_parameter')].copy()
        if not species_parameter_rows.empty:
            rows = species_parameter_rows
    rows = rows.sort_index()
    grouping_columns = []
    if 'Model' in rows.columns:
        grouping_columns.append('Model')
    if 'RefID' in rows.columns:
        grouping_columns.append('RefID')
    if grouping_columns:
        group_keys = rows[grouping_columns].drop_duplicates()
        candidate_groups = []
        for _, key_row in group_keys.iterrows():
            mask = np.ones(len(rows), dtype=bool)
            for col in grouping_columns:
                value = key_row[col]
                if pd.isna(value):
                    mask &= rows[col].isna()
                else:
                    mask &= rows[col].astype(str).eq(str(value))
            candidate_groups.append(rows.loc[mask].copy())
    else:
        candidate_groups = [rows]
    for group in candidate_groups:
        m_i = extract_epc_quantity_6(group, species_name, 'm', T, required=True)
        sigma_i = extract_epc_quantity_6(group, species_name, 'sigma', T, required=True)
        epsilon_i = extract_epc_quantity_6(group, species_name, 'epsilon', T, required=True)
        if m_i is None or sigma_i is None or epsilon_i is None:
            continue
        if not np.all(np.isfinite([m_i, sigma_i, epsilon_i])):
            continue
        epsilonAB_i = extract_epc_quantity_6(group, species_name, 'epsilonAB', T, required=False)
        kappaAB_i = extract_epc_quantity_6(group, species_name, 'kappaAB', T, required=False)
        association_i = extract_epc_quantity_6(group, species_name, 'association_scheme', T, required=False)
        first_row = group.iloc[0]
        return {'Species': species_name, 'Formula': PCSAFT_FORMULA[species_name], 'Source': 'PROJECT_EPCSAFT.csv', 'Database row': int(group.index.min()), 'Model': str(first_row['Model']) if 'Model' in group.columns else 'ePC-SAFT', 'RefID': str(first_row['RefID']) if 'RefID' in group.columns and pd.notna(first_row['RefID']) else '', 'm': float(m_i), 'sigma': float(sigma_i), 'epsilon': float(epsilon_i), 'epsilonAB': float(epsilonAB_i) if pd.notna(epsilonAB_i) else np.nan, 'kappaAB': float(kappaAB_i) if pd.notna(kappaAB_i) else np.nan, 'association_scheme': str(association_i).strip() if pd.notna(association_i) else '-'}
    return None

def first_existing_column_6(dataframe, candidates):
    normalized = {normalize_species_name_6(col): col for col in dataframe.columns}
    for candidate in candidates:
        key = normalize_species_name_6(candidate)
        if key in normalized:
            return normalized[key]
    return None

def optional_numeric_6(row, column):
    if column is None:
        return np.nan
    value = row[column]
    if pd.isna(value):
        return np.nan
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan

def first_complete_pcsaft_set(species_name, T):
    if species_name not in PCSAFT_FORMULA:
        raise KeyError(f"No PC-SAFT formula mapping for '{species_name}'.")
    formula = PCSAFT_FORMULA[species_name]
    formula_mask = pcsaft_raw['Formula'].map(normalize_species_name_6).eq(normalize_species_name_6(formula))
    formula_rows = pcsaft_raw[formula_mask].copy()
    if formula_rows.empty:
        return None
    formula_rows = formula_rows.sort_index()
    MODEL_PRIORITY_6 = ['pc-saft', 'gc-ppc-saft']
    for preferred_model in MODEL_PRIORITY_6:
        model_rows = formula_rows[formula_rows['Model'].astype(str).str.strip().str.casefold().eq(preferred_model)].copy()
        if model_rows.empty:
            continue
        for idx, row in model_rows.iterrows():
            try:
                m_i = float(row[PC_COL_M_6])
                sigma_i = float(row[PC_COL_SIGMA_6])
                epsilon_i = float(row[PC_COL_EPSILON_6])
            except (TypeError, ValueError, KeyError):
                continue
            if not np.all(np.isfinite([m_i, sigma_i, epsilon_i])):
                continue
            if m_i <= 0.0 or sigma_i <= 0.0 or epsilon_i <= 0.0:
                continue
            epsilonAB_i = optional_numeric_6(row, PC_COL_EPSILONAB_6)
            kappaAB_i = optional_numeric_6(row, PC_COL_KAPPAAB_6)
            association_i = '-'
            if PC_COL_ASSOCIATION_6 is not None:
                raw_assoc = row[PC_COL_ASSOCIATION_6]
                if pd.notna(raw_assoc):
                    association_i = str(raw_assoc).strip()
            return {'Species': species_name, 'Formula': formula, 'Source': 'PROJECT_PCSAFT.csv', 'Database row': int(idx), 'Model': str(row['Model']), 'RefID': str(row['RefID']) if 'RefID' in row.index and pd.notna(row['RefID']) else '', 'm': m_i, 'sigma': sigma_i, 'epsilon': epsilon_i, 'epsilonAB': epsilonAB_i, 'kappaAB': kappaAB_i, 'association_scheme': association_i}
    return None
PCSAFT_FORMULA['nitrogen'] = 'N2'
SPECIES_ALIASES_6['nitrogen'] = ['nitrogen', 'N2']

# %% cell 10
ION_SPECIES = ['Na+', 'HCO3-', 'CO3^2-']
ION_CHARGE = {'Na+': +1, 'HCO3-': -1, 'CO3^2-': -2}
ION_FORMULA = {'Na+': 'Na+', 'HCO3-': 'HCO3-', 'CO3^2-': 'CO3'}
SUPCRT_ION_ALIASES = {'Na+': ['Na+'], 'HCO3-': ['HCO3-'], 'CO3^2-': ['CO3--', 'CO3-2', 'CO32-', 'CO3^2-']}

def normalize_ion_name_7(value):
    if pd.isna(value):
        return ''
    return str(value).strip().casefold().replace(' ', '').replace('₂', '2').replace('₃', '3').replace('⁺', '+').replace('⁻', '-')

def numeric_epc_value_7(value):
    if pd.isna(value):
        return np.nan
    try:
        return float(value)
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    raise ValueError(f'Could not convert ePC-SAFT value {value!r} to a number.')
ION_EPC_PARAMETER_ALIASES = {'m': ['mseg', 'm', 'segment'], 'sigma': ['sigma', 'sigma_a', 'sigma (a)', 'sigma (å)'], 'epsilon': ['u/kb', 'epsilon/k', 'epsilon/kb', 'epsilon_k'], 'dielectric_constant': ['dielectric_constant']}

def first_epc_parameter_7(rows, quantity, required=True):
    aliases = [str(x).strip().casefold() for x in ION_EPC_PARAMETER_ALIASES[quantity]]
    for alias in aliases:
        match = rows[rows['Parameter'].astype(str).str.strip().str.casefold().eq(alias)]
        if not match.empty:
            return numeric_epc_value_7(match.iloc[0]['Value'])
    if required:
        available = rows['Parameter'].dropna().astype(str).unique().tolist()
        raise ValueError(f"Required ePC-SAFT parameter '{quantity}' not found.\nAvailable parameters: {available}")
    return np.nan

def first_complete_ion_epcsaft_set_7(ion_name):
    aliases = {normalize_ion_name_7(ion_name), normalize_ion_name_7(ION_FORMULA[ion_name])}
    rows = epcsaft_raw[epcsaft_raw['Species'].map(normalize_ion_name_7).isin(aliases)].copy()
    if rows.empty:
        return None
    rows = rows.sort_index()
    if 'record_type' in rows.columns:
        species_rows = rows[rows['record_type'].astype(str).str.strip().str.casefold().eq('species_parameter')].copy()
        if not species_rows.empty:
            rows = species_rows
    grouping_columns = [col for col in ['RefID', 'Model'] if col in rows.columns]
    if grouping_columns:
        grouped = rows.groupby(grouping_columns, sort=False, dropna=False)
        parameter_sets = [group for _, group in grouped]
    else:
        parameter_sets = [rows]
    for parameter_set in parameter_sets:
        try:
            m_i = first_epc_parameter_7(parameter_set, 'm')
            sigma_i = first_epc_parameter_7(parameter_set, 'sigma')
            epsilon_i = first_epc_parameter_7(parameter_set, 'epsilon')
        except ValueError:
            continue
        if not np.all(np.isfinite([m_i, sigma_i, epsilon_i])):
            continue
        if m_i <= 0.0 or sigma_i <= 0.0 or epsilon_i <= 0.0:
            continue
        dielectric_i = first_epc_parameter_7(parameter_set, 'dielectric_constant', required=False)
        first_row = parameter_set.iloc[0]
        charge_database = np.nan
        if 'Charge' in parameter_set.columns:
            charge_values = parameter_set['Charge'].dropna()
            if not charge_values.empty:
                charge_database = float(charge_values.iloc[0])
        return {'Species': ion_name, 'Source': 'PROJECT_EPCSAFT.csv', 'Database row': int(parameter_set.index.min()), 'Model': str(first_row['Model']) if 'Model' in first_row.index and pd.notna(first_row['Model']) else '', 'RefID': str(first_row['RefID']) if 'RefID' in first_row.index and pd.notna(first_row['RefID']) else '', 'charge': ION_CHARGE[ion_name], 'charge_database': charge_database, 'm': m_i, 'sigma': sigma_i, 'epsilon': epsilon_i, 'dielectric_constant': dielectric_i}
    return None

def xml_float_7(element, tag, required=True):
    child = element.find(f'.//{tag}')
    if child is None or child.text is None or (not child.text.strip()):
        if required:
            raise ValueError(f"Missing SUPCRT field '{tag}'.")
        return np.nan
    text = child.text.strip()
    try:
        return float(text)
    except ValueError:
        if required:
            raise ValueError(f"SUPCRT field '{tag}' is not numeric: {text!r}")
        return np.nan

def supcrt_species_name_7(species_element):
    name_element = species_element.find('Name')
    if name_element is None or name_element.text is None:
        return ''
    return name_element.text.strip()

def supcrt_candidates_7(ion_name):
    aliases = {normalize_ion_name_7(alias) for alias in SUPCRT_ION_ALIASES[ion_name]}
    matches = []
    for species_element in supcrt_root.iter('Species'):
        name = supcrt_species_name_7(species_element)
        if normalize_ion_name_7(name) in aliases:
            matches.append(species_element)
    return matches
HKF_REQUIRED_FIELDS_7 = ['Gf', 'Hf', 'Sr', 'a1', 'a2', 'a3', 'a4', 'c1', 'c2', 'wref']

def read_hkf_candidate_7(species_element):
    charge_element = species_element.find('Charge')
    if charge_element is None or charge_element.text is None:
        return None
    try:
        charge_i = float(charge_element.text.strip())
    except ValueError:
        return None
    hkf = species_element.find('./Thermo/HKF')
    if hkf is None:
        return None
    result = {'Charge': charge_i}
    for field in HKF_REQUIRED_FIELDS_7:
        child = hkf.find(field)
        if child is None or child.text is None or (not child.text.strip()):
            return None
        try:
            result[field] = float(child.text.strip())
        except ValueError:
            return None
    return result

def first_complete_hkf_set_7(ion_name):
    candidates = supcrt_candidates_7(ion_name)
    for candidate_number, candidate in enumerate(candidates):
        values = read_hkf_candidate_7(candidate)
        if values is None:
            continue
        result = {'Species': ion_name, 'SUPCRT name': supcrt_species_name_7(candidate), 'Candidate': candidate_number}
        result.update(values)
        for tag in ['Formula', 'MolarMass', 'Type', 'Date', 'References']:
            child = candidate.find(f'.//{tag}')
            result[tag] = child.text.strip() if child is not None and child.text is not None else ''
        return result
    return None

# %% cell 12
def molecular_weight_from_formula(formula):
    formula = str(formula).strip()
    tokens = re.findall('([A-Z][a-z]?)(\\d*)', formula)
    if not tokens:
        raise ValueError(f'Could not parse molecular formula: {formula}')
    molecular_weight = 0.0
    for element, number in tokens:
        if element not in ATOMIC_WEIGHTS:
            raise KeyError(f"Atomic weight missing for '{element}' in formula '{formula}'.")
        count = int(number) if number else 1
        molecular_weight += ATOMIC_WEIGHTS[element] * count
    return float(molecular_weight)

def association_sites(association_scheme, epsilonAB, kappaAB):
    scheme = str(association_scheme).strip().upper()
    if scheme == '2B':
        return (1, 1)
    if scheme == '3B':
        return (1, 2)
    if np.isfinite(epsilonAB) and np.isfinite(kappaAB) and (epsilonAB > 0.0) and (kappaAB > 0.0):
        return (1, 1)
    return (0, 0)

def build_pcsaft_model(phase_species):
    components = []
    MW = []
    segment = []
    sigma = []
    epsilon = []
    n_H = []
    n_e = []
    epsilon_assoc = {}
    bondvol = {}
    parameter_rows = []
    for species_name in phase_species:
        if species_name not in resolved_pcsaft_lookup:
            raise KeyError(f"No resolved PC-SAFT parameter set for '{species_name}'.")
        parameters = resolved_pcsaft_lookup[species_name]
        component_name = str(parameters['Formula']).strip()
        components.append(component_name)
        MW_i = molecular_weight_from_formula(component_name)
        MW.append(MW_i)
        m_i = float(parameters['m'])
        sigma_i = float(parameters['sigma'])
        epsilon_i = float(parameters['epsilon'])
        segment.append(m_i)
        sigma.append(sigma_i)
        epsilon.append(epsilon_i)
        epsilonAB_i = parameters.get('epsilonAB', np.nan)
        kappaAB_i = parameters.get('kappaAB', np.nan)
        association_scheme_i = parameters.get('association_scheme', '-')
        try:
            epsilonAB_i = float(epsilonAB_i)
        except (TypeError, ValueError):
            epsilonAB_i = np.nan
        try:
            kappaAB_i = float(kappaAB_i)
        except (TypeError, ValueError):
            kappaAB_i = np.nan
        nH_i, ne_i = association_sites(association_scheme_i, epsilonAB_i, kappaAB_i)
        n_H.append(nH_i)
        n_e.append(ne_i)
        if nH_i > 0 and ne_i > 0 and np.isfinite(epsilonAB_i) and np.isfinite(kappaAB_i):
            association_key = ((component_name, 'e'), (component_name, 'H'))
            epsilon_assoc[association_key] = float(epsilonAB_i)
            bondvol[association_key] = float(kappaAB_i)
        parameter_rows.append({'Species': species_name, 'Formula': component_name, 'Source': parameters.get('Source', '-'), 'Database row': parameters.get('Database row', np.nan), 'Model': parameters.get('Model', '-'), 'RefID': parameters.get('RefID', '-'), 'MW [g/mol]': MW_i, 'm': m_i, 'sigma [A]': sigma_i, 'epsilon/kB [K]': epsilon_i, 'epsilonAB/kB [K]': epsilonAB_i, 'kappaAB': kappaAB_i, 'association_scheme': association_scheme_i, 'n_H': nH_i, 'n_e': ne_i})
    MW = np.asarray(MW, dtype=float)
    segment = np.asarray(segment, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    epsilon = np.asarray(epsilon, dtype=float)
    n_H = np.asarray(n_H, dtype=int)
    n_e = np.asarray(n_e, dtype=int)
    for name, values in {'Mw': MW, 'segment': segment, 'sigma': sigma, 'epsilon': epsilon}.items():
        if len(values) != len(phase_species):
            raise ValueError(f'{name} length mismatch.')
        if np.any(~np.isfinite(values)):
            raise ValueError(f'{name} contains non-finite values.')
    parameter_table = pd.DataFrame(parameter_rows)
    userlocations = {'Mw': MW, 'segment': segment, 'sigma': sigma, 'epsilon': epsilon}
    has_association = bool(np.any(n_H > 0) or np.any(n_e > 0))
    if has_association:
        userlocations['n_H'] = n_H
        userlocations['n_e'] = n_e
        if epsilon_assoc:
            userlocations['epsilon_assoc'] = epsilon_assoc
        if bondvol:
            userlocations['bondvol'] = bondvol
    if 'CO2' in components and 'CH4' in components:
        k_matrix = np.zeros((len(components), len(components)))
        i = components.index('CO2')
        j = components.index('CH4')
        k_matrix[i, j] = k_matrix[j, i] = 0.065
        userlocations['k'] = k_matrix
    model = cl.PCSAFT(components, userlocations=userlocations)
    return (model, parameter_table)

def nasa_standard_chemical_potentials(phase_species):
    mu0 = np.zeros(len(phase_species), dtype=float)
    rows = []
    for i, species_name in enumerate(phase_species):
        nasa_name = get_nasa_name(species_name)
        if nasa_name not in nasa_db:
            raise KeyError(f"NASA species '{nasa_name}' for '{species_name}' is not present in nasa_db.")
        cp, h, entropy, g = nasa7_properties(nasa_db, nasa_name, T)
        mu0[i] = float(g)
        rows.append({'Species': species_name, 'NASA species': nasa_name, 'g0 NASA [kJ/mol]': g / 1000.0, 'T low [K]': nasa_db[nasa_name]['T_low'], 'T high [K]': nasa_db[nasa_name]['T_high']})
    return (mu0, pd.DataFrame(rows))

def get_nasa_name(species_name):
    registry = {**MOLECULAR_SPECIES, **ORGANIC_SPECIES, **SOLID_DATA}
    if species_name not in registry:
        raise KeyError(f'No NASA mapping for {species_name}')
    return registry[species_name]['nasa']

def solid_standard_chemical_potentials(solid_species):
    mu0_solid = np.zeros(len(solid_species), dtype=float)
    rows = []
    for i, species_name in enumerate(solid_species):
        nasa_name = get_nasa_name(species_name)
        if nasa_name not in nasa_db:
            raise KeyError(f"NASA species '{nasa_name}' for solid '{species_name}' is not present in nasa_db.")
        cp, h, entropy, g = nasa7_properties(nasa_db, nasa_name, T)
        mu0_solid[i] = float(g)
        rows.append({'Species': species_name, 'NASA species': nasa_name, 'Phase': 'pure solid', 'cp [J/mol/K]': cp, 'h [kJ/mol]': h / 1000.0, 's [J/mol/K]': entropy, 'g0 NASA [kJ/mol]': g / 1000.0, 'T low [K]': nasa_db[nasa_name]['T_low'], 'T high [K]': nasa_db[nasa_name]['T_high']})
    return (mu0_solid, pd.DataFrame(rows))

# %% cell 14
def supcrt_water(t, pressure_pa):
    props = ['rho', 'epsilon', 'QBorn', 'XBorn', 'YBorn', 'alpha', 'daldT', 'beta']
    out = pch.water(props, T=float(t), P=float(pressure_pa) / 100000.0)
    return {key: float(np.asarray(out[key]).reshape(-1)[0]) for key in props}

def hkf_standard_ions(t, pressure_pa):
    records = ion_hkf_table.set_index('Species')
    rows = []
    for species in AQ_ION_SPECIES:
        d = records.loc[species]
        rows.append({'name': str(d['SUPCRT name']), 'state': 'aq', 'G': float(d['Gf']), 'H': float(d['Hf']), 'S': float(d['Sr']), 'Cp': np.nan, 'V': np.nan, 'a1.a': float(d['a1']) * 10, 'a2.b': float(d['a2']) / 100, 'a3.c': float(d['a3']), 'a4.d': float(d['a4']) / 10000.0, 'c1.e': float(d['c1']), 'c2.f': float(d['c2']) / 10000.0, 'omega.lambda': float(d['wref']) / 100000.0, 'z.T': float(d['Charge'])})
    params = pd.DataFrame(rows)
    result, returned_water = pch.hkf(property=['G'], parameters=params, T=float(t), P=float(pressure_pa) / 100000.0, contrib=['n', 's', 'o'], H2O_props=['rho'], water_model='SUPCRT92')
    values = np.array([float(np.asarray(result[i]['G']).reshape(-1)[0]) * 4.184 for i in range(len(rows))])
    if not np.all(np.isfinite(values)):
        raise ValueError('Nonfinite HKF ion standard Gibbs energies')
    # CHNOSZ HKF returns apparent G anchored to Gf(Tr). NASA uses
    # Hf(Tr)-Tr*Sr at Tr. Shift the constant reference for EACH ion.
    # Never subtract T*Sr here: HKF already includes the entropy evolution.
    offsets = np.array([(float(records.loc[s, 'Gf']) - float(records.loc[s, 'Hf'])
                         + 298.15 * float(records.loc[s, 'Sr'])) * 4.184
                        for s in AQ_ION_SPECIES])
    return values - offsets

def first_existing_value(mapping, names, default=np.nan):
    for name in names:
        if name in mapping:
            value = mapping[name]
            if value is not None and (not pd.isna(value)):
                return value
    return default

def required_float(value, parameter, species):
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{species}: invalid {parameter} = {value!r}') from exc
    if not np.isfinite(result):
        raise ValueError(f'{species}: missing {parameter}')
    return result

def julia_string_vector(values):
    escaped = ['"' + str(value).replace('\\', '\\\\').replace('"', '\\"') + '"' for value in values]
    return jl.seval('String[' + ', '.join(escaped) + ']')

def julia_int_vector(values):
    return jl.seval('Int64[' + ', '.join((str(int(value)) for value in values)) + ']')

def build_aqueous_short_range():
    global AQUEOUS_SPECIES, aqueous_components, aqueous_charges
    global aqueous_pcsaft_userlocations, aqueous_pcsaft_model, aqueous_pcsaft_parameter_table
    global aqueous_components_julia, aqueous_charge_julia, aqueous_references_julia
    AQUEOUS_SPECIES = list(AQ_NEUTRAL_SPECIES) + list(ION_SPECIES)
    aqueous_components = [AQUEOUS_FORMULA.get(species, species) for species in AQUEOUS_SPECIES]
    aqueous_charges = [0 if species in AQ_NEUTRAL_SPECIES else int(ION_CHARGE[species]) for species in AQUEOUS_SPECIES]
    ion_lookup = {str(row['Species']): row.to_dict() for _, row in ion_epcsaft_table.iterrows()}
    ION_MW = {name: sum((AW[e] * count for e, count in data['elements'].items())) for name, data in ION_DATA.items()}
    Mw = []
    segment = []
    sigma = []
    epsilon = []
    n_H = []
    n_e = []
    parameter_rows = []
    epsilon_assoc = {}
    bondvol = {}
    for index, species in enumerate(AQUEOUS_SPECIES):
        if species in AQ_NEUTRAL_SPECIES:
            if species not in resolved_pcsaft_lookup:
                raise KeyError(f'No resolved molecular PC-SAFT parameters for {species}')
            parameters = resolved_pcsaft_lookup[species]
            elements = ALL_MOLECULAR_DATA[species]['elements']
            registry_Mw = sum((AW[element] * count for element, count in elements.items()))
            Mw_i = first_existing_value(parameters, ['Mw', 'MW', 'MW_g_mol', 'MW [g/mol]'], default=registry_Mw)
        else:
            if species not in ion_lookup:
                raise KeyError(f'No resolved ion ePC-SAFT parameters for {species}')
            parameters = ion_lookup[species]
            Mw_i = first_existing_value(parameters, ['Mw', 'MW', 'MW_g_mol', 'MW [g/mol]'], default=ION_MW[species])
        m_i = first_existing_value(parameters, ['m', 'segment'])
        sigma_i = first_existing_value(parameters, ['sigma', 'sigma_A', 'sigma [A]'])
        epsilon_i = first_existing_value(parameters, ['epsilon', 'u_over_kB_K', 'epsilon/kB [K]'])
        epsilonAB_i = first_existing_value(parameters, ['epsilonAB', 'epsilon_assoc', 'epsilonAB/kB [K]'])
        kappaAB_i = first_existing_value(parameters, ['kappaAB', 'bondvol'])
        association_scheme_i = first_existing_value(parameters, ['association_scheme', 'Association'], default='-')
        source_i = first_existing_value(parameters, ['Source', 'SourceDatabase'], default='earlier cells')
        Mw_i = required_float(Mw_i, 'Mw', species)
        m_i = required_float(m_i, 'm', species)
        sigma_i = required_float(sigma_i, 'sigma', species)
        epsilon_i = required_float(epsilon_i, 'epsilon', species)
        Mw.append(Mw_i)
        segment.append(m_i)
        sigma.append(sigma_i)
        epsilon.append(epsilon_i)
        has_association_parameters = pd.notna(epsilonAB_i) and pd.notna(kappaAB_i) and np.isfinite(float(epsilonAB_i)) and np.isfinite(float(kappaAB_i))
        if has_association_parameters:
            epsilonAB_i = float(epsilonAB_i)
            kappaAB_i = float(kappaAB_i)
            nH_i, ne_i = association_sites(association_scheme_i, epsilonAB_i, kappaAB_i)
        else:
            nH_i, ne_i = (0, 0)
        n_H.append(int(nH_i))
        n_e.append(int(ne_i))
        component = aqueous_components[index]
        if nH_i > 0 and ne_i > 0:
            key = ((component, 'H'), (component, 'e'))
            epsilon_assoc[key] = epsilonAB_i
            bondvol[key] = float(kappaAB_i)
        parameter_rows.append({'Species': species, 'Component': component, 'Charge': aqueous_charges[index], 'Source': source_i, 'Mw [g/mol]': Mw_i, 'm': m_i, 'sigma [Å]': sigma_i, 'epsilon/kB [K]': epsilon_i, 'Association': association_scheme_i, 'n_H': int(nH_i), 'n_e': int(ne_i)})
    Mw = np.asarray(Mw, dtype=float)
    segment = np.asarray(segment, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    epsilon = np.asarray(epsilon, dtype=float)
    n_H = np.asarray(n_H, dtype=int)
    n_e = np.asarray(n_e, dtype=int)
    aqueous_pcsaft_parameter_table = pd.DataFrame(parameter_rows)
    for name, values in [('Mw', Mw), ('segment', segment), ('sigma', sigma), ('epsilon', epsilon), ('n_H', n_H), ('n_e', n_e)]:
        if len(values) != len(AQUEOUS_SPECIES):
            raise ValueError(f'{name} has length {len(values)}; expected {len(AQUEOUS_SPECIES)}.')
    aqueous_pcsaft_userlocations = {'Mw': Mw, 'segment': segment, 'sigma': sigma, 'epsilon': epsilon, 'n_H': n_H, 'n_e': n_e, 'epsilon_assoc': epsilon_assoc, 'bondvol': bondvol}
    aqueous_pcsaft_model = cl.PCSAFT(aqueous_components, userlocations=aqueous_pcsaft_userlocations, verbose=False)
    aqueous_components_julia = julia_string_vector(aqueous_components)
    aqueous_charge_julia = julia_int_vector(aqueous_charges)
    aqueous_references_julia = jl.seval('String[]')

# %% cell 15
def molecular_gibbs(n, species, mu0, phase):
    if np.sum(n) <= 1e-13:
        return 0.0
    positive = np.maximum(n, 1e-14)
    model, _ = build_pcsaft_model(species)
    phis = np.asarray(cl.fugacity_coefficient(model, P, T, positive, phase=phase), dtype=float).reshape(-1)
    if np.any(~np.isfinite(phis)) or np.any(phis <= 0):
        raise ValueError(f'Invalid PC-SAFT fugacities: phase={phase}, T={T - 273.15:g} C, P={P / 100000.0:g} bar, species={species}, x={(positive / positive.sum()).tolist()}')
    mu = mu0 + R * T * np.log(positive / positive.sum() * P / P0 * phis)
    return float(np.dot(n, mu))

def aqueous_gibbs(n):
    if n[0] <= 1e-09:
        if np.sum(n[1:]) > 1e-10:
            raise ValueError('Aqueous solutes without water')
        return 0.0
    phi = phi_electrolyte(n)
    q = np.maximum(n, 1e-14)
    x = q / q.sum()
    xw = x[0]
    mu_neutral = mu0_aq_neutral + R * T * np.log(x[:len(AQ_NEUTRAL_SPECIES)] * P / P0 * phi[:len(AQ_NEUTRAL_SPECIES)])
    molal = np.maximum(n[-len(AQ_ION_SPECIES):], 1e-14) / (MW_WATER_KG * n[0])
    mu_ions = mu0_ion + R * T * (np.log(molal) + np.log(xw) + np.log(phi[-len(AQ_ION_SPECIES):]) - lnphi_infinity[-len(AQ_ION_SPECIES):])
    return float(np.dot(n[:len(AQ_NEUTRAL_SPECIES)], mu_neutral) + np.dot(n[-len(AQ_ION_SPECIES):], mu_ions))

# %% cell 17
def feed_from_species(species_moles):
    """Convert a dry-feed molar composition to elemental wt%; added water is separate."""
    registry = {**ALL_MOLECULAR_DATA, **SOLID_DATA, **ION_DATA}
    atoms = {}
    for name, amount in species_moles.items():
        if name not in registry or not np.isfinite(amount) or amount < 0:
            raise ValueError(f'Invalid dry-feed species or amount: {name}={amount}')
        for element, count in registry[name]['elements'].items():
            atoms[element] = atoms.get(element, 0.0) + amount * count
    masses = {e: count * AW[e] for e, count in atoms.items()}
    total = sum(masses.values())
    if total <= 0 or masses.get('C', 0) <= 0:
        raise ValueError('Dry feed must contain positive mass and carbon')
    return {e: mass / total * 100 for e, mass in masses.items()}

def _refresh_state(temperature_C, pressure_bar, electrolyte=False, phase_species=None):
    global T, P, resolved_pcsaft_lookup, water_operating, dielectric_constant
    global aqueous_rsp_model, aqueous_ion_model, aqueous_epcsaft_model
    global mu0_ion, mu0_aq_neutral, lnphi_infinity
    T = float(temperature_C) + 273.15
    P = float(pressure_bar) * 100000.0
    selected = list(dict.fromkeys([s for phase, names in (phase_species or {}).items() if phase != 'solid' for s in names if s not in ION_DATA] + (AQ_NEUTRAL_SPECIES if electrolyte else [])))
    resolved_pcsaft_lookup = resolve_molecular_parameters(T, selected)
    if not electrolyte:
        return
    water_operating = supcrt_water(T, P)
    dielectric_constant = water_operating['epsilon']
    if not (np.isfinite(water_operating['rho']) and water_operating['rho'] >= 350 and np.isfinite(dielectric_constant) and (dielectric_constant > 1)):
        raise ValueError(f"State outside SUPCRT aqueous validity checks: rho={water_operating['rho']}, epsilon={dielectric_constant}")
    build_aqueous_short_range()
    aqueous_rsp_model = cl.ConstRSP(float(dielectric_constant))
    aqueous_ion_model = cl.hsdDH(aqueous_components_julia, aqueous_rsp_model, aqueous_references_julia)
    aqueous_epcsaft_model = cl.ePCSAFT(aqueous_components_julia, aqueous_charge_julia, aqueous_pcsaft_model.idealmodel, aqueous_pcsaft_model, aqueous_ion_model, aqueous_references_julia)
    mu0_ion = hkf_standard_ions(T, P)
    mu0_aq_neutral, _ = nasa_standard_chemical_potentials(AQ_NEUTRAL_SPECIES)
    lnphi_infinity = infinite_dilution_reference()

# %% cell 19
def solve_equilibrium(settings):
    """Solve one carbon-normalized feed with explicit phase/species and state settings."""
    settings, phase_allocation = prepare_phase_allocation(settings)
    if not MODEL_INITIALIZED:
        raise RuntimeError('Call initialize_model(database_paths) first')
    state_keys = ('T', 'P', 'resolved_pcsaft_lookup', 'water_operating', 'dielectric_constant', 'aqueous_rsp_model', 'aqueous_ion_model', 'aqueous_epcsaft_model', 'aqueous_pcsaft_model', 'aqueous_pcsaft_userlocations', 'aqueous_pcsaft_parameter_table', 'aqueous_components', 'aqueous_charges', 'aqueous_components_julia', 'aqueous_charge_julia', 'aqueous_references_julia', 'mu0_ion', 'mu0_aq_neutral', 'lnphi_infinity', 'AQUEOUS_SPECIES', 'AQ_NEUTRAL_SPECIES', 'AQ_ION_SPECIES', 'ION_SPECIES', 'AQUEOUS_FORMULA')
    saved = {name: globals()[name] for name in state_keys if name in globals()}
    absent = [name for name in state_keys if name not in globals()]
    stage = 'feed and initial composition'
    try:
        TEMPERATURE_C = float(settings['temperature_C'])
        PRESSURE_BAR = float(settings['pressure_bar'])
        FEED_ELEMENT_WT_PERCENT = dict(settings.get('feed_element_wt_percent', {}))
        if settings.get('feed_species_moles') is not None:
            FEED_ELEMENT_WT_PERCENT = feed_from_species(settings['feed_species_moles'])
        WATER_G_PER_G_DRY_FEED = float(settings.get('water_g_per_g_dry_feed', 0.0))
        ADDED_SALTS = settings.get('added_salts', {})
        CASE_PHASE_SPECIES = settings['phase_species']
        configure_aqueous_species(CASE_PHASE_SPECIES)
        GAS_CARBON_FRACTION = settings.get('gas_carbon_fraction')
        GAS_CARBON_MIN_FRACTION = settings.get('gas_carbon_min_fraction')
        if GAS_CARBON_MIN_FRACTION is not None:
            GAS_CARBON_MIN_FRACTION = float(GAS_CARBON_MIN_FRACTION)
            if not np.isfinite(GAS_CARBON_MIN_FRACTION) or not 0 <= GAS_CARBON_MIN_FRACTION <= 1:
                raise ValueError('Minimum gas-carbon fraction must be in [0, 1]')
            if GAS_CARBON_FRACTION is not None:
                raise ValueError('Use either gas_carbon_fraction (equality) or gas_carbon_min_fraction (lower bound), not both')
        SOLID_CARBON_FRACTION = settings.get('solid_carbon_fraction')
        SOLID_CARBON_TARGET = settings.get('solid_carbon_target', 'CHAR')
        MINIMUM_AMOUNTS = settings.get('minimum_amounts', {})
        AQUEOUS_WATER_FLOOR = float(settings.get('aqueous_water_floor', 1e-08))
        MAX_ITERATIONS = int(settings.get('max_iterations', 5000))
        AQUEOUS_MODEL = settings.get('aqueous_model', 'electrolyte')
        if AQUEOUS_MODEL not in ('molecular', 'electrolyte'):
            raise ValueError('aqueous_model must be molecular or electrolyte')
        if AQUEOUS_MODEL == 'molecular' and any((s in ION_DATA for s in CASE_PHASE_SPECIES.get('aqueous', []))):
            raise ValueError('Ions require the electrolyte aqueous model')
        LIQUID_INITIALIZATION = bool(settings.get('liquid_initialization', False))
        INITIAL_GAS_WATER_FRACTION = float(settings.get('initial_gas_water_fraction', 1e-06))
        if not np.isfinite(INITIAL_GAS_WATER_FRACTION) or INITIAL_GAS_WATER_FRACTION < 0:
            raise ValueError('Invalid initial gas-water fraction')
        INITIAL_PHASE_AMOUNTS = settings.get('initial_phase_amounts')
        if not np.isfinite(TEMPERATURE_C) or TEMPERATURE_C <= -273.15:
            raise ValueError('Invalid temperature.')
        if not np.isfinite(PRESSURE_BAR) or PRESSURE_BAR <= 0:
            raise ValueError('Pressure must be finite and positive.')
        if not np.isfinite(WATER_G_PER_G_DRY_FEED) or WATER_G_PER_G_DRY_FEED < 0:
            raise ValueError('Added water must be finite and nonnegative.')
        if any((not np.isfinite(value) or value < 0 for value in FEED_ELEMENT_WT_PERCENT.values())):
            raise ValueError('Invalid feed elemental composition.')
        if FEED_ELEMENT_WT_PERCENT.get('C', 0.0) <= 0:
            raise ValueError('This carbon-normalized calculation requires positive dry-feed carbon.')
        if sum(FEED_ELEMENT_WT_PERCENT.values()) > 100.0 + 1e-08:
            raise ValueError('Feed elemental percentages exceed 100%.')
        if CASE_PHASE_SPECIES.get('aqueous') and AQUEOUS_MODEL == 'electrolyte' and (not np.isfinite(AQUEOUS_WATER_FLOOR) or AQUEOUS_WATER_FLOOR <= 1e-09):
            raise ValueError("AQUEOUS_WATER_FLOOR must exceed the existing aqueous model's 1e-9 mol solvent cutoff.")
        valid_phases = ('gas', 'aqueous', 'organic liquid', 'solid')
        if set(CASE_PHASE_SPECIES) - set(valid_phases):
            raise ValueError('Unknown phase in CASE_PHASE_SPECIES.')
        manual_phase_species = {phase: list(CASE_PHASE_SPECIES.get(phase, [])) for phase in valid_phases}
        for phase, names in manual_phase_species.items():
            if len(names) != len(set(names)):
                raise ValueError(f'Duplicate species in {phase}.')

        def manual_species_data(phase, species):
            if phase == 'solid':
                registry = SOLID_DATA
            elif phase == 'aqueous' and species in ION_DATA:
                registry = ION_DATA
            else:
                registry = {**MOLECULAR_SPECIES, **ORGANIC_SPECIES}
            if species not in registry:
                raise ValueError(f"{phase}: '{species}' is not registered. Available species: {', '.join(registry)}")
            if phase == 'aqueous' and species not in AQUEOUS_SPECIES:
                raise ValueError(f'{species}: not supported by the introduced aqueous model. Available: {AQUEOUS_SPECIES}')
            return registry[species]
        manual_pairs = [(phase, species) for phase, names in manual_phase_species.items() for species in names]
        if not manual_pairs:
            raise ValueError('Select at least one species.')
        for phase, species in manual_pairs:
            manual_species_data(phase, species)
        manual_indices = {phase: np.array([j for j, (candidate_phase, _) in enumerate(manual_pairs) if candidate_phase == phase], dtype=int) for phase in valid_phases}
        manual_elements = sorted(set(FEED_ELEMENT_WT_PERCENT) | {'H', 'O'} | {element for salt in ADDED_SALTS.values() for element in salt['elements']} | {element for phase, species in manual_pairs for element in manual_species_data(phase, species)['elements']})
        missing_atomic_weights = [element for element in manual_elements if element not in AW]
        if missing_atomic_weights:
            raise ValueError('Missing atomic weights: ' + ', '.join(missing_atomic_weights))
        manual_A = np.array([[manual_species_data(phase, species)['elements'].get(element, 0.0) for phase, species in manual_pairs] for element in manual_elements], dtype=float)
        manual_C_mol_per_g = FEED_ELEMENT_WT_PERCENT['C'] / 100.0 / AW['C']
        manual_b = np.array([FEED_ELEMENT_WT_PERCENT.get(element, 0.0) / 100.0 / AW[element] / manual_C_mol_per_g for element in manual_elements], dtype=float)
        water_MW = 2.0 * AW['H'] + AW['O']
        water_mol_per_C = WATER_G_PER_G_DRY_FEED / water_MW / manual_C_mol_per_g
        manual_b[manual_elements.index('H')] += 2.0 * water_mol_per_C
        manual_b[manual_elements.index('O')] += water_mol_per_C
        for salt_name, salt in ADDED_SALTS.items():
            amount = float(salt['mol_per_g_dry_feed'])
            if not np.isfinite(amount) or amount < 0:
                raise ValueError(f'Invalid amount for {salt_name}.')
            for element, count in salt['elements'].items():
                if not np.isfinite(count) or count < 0:
                    raise ValueError(f'Invalid elemental count for {salt_name}.')
                manual_b[manual_elements.index(element)] += amount * count / manual_C_mol_per_g
        manual_carbon = manual_A[manual_elements.index('C')]
        feed_carbon = manual_b[manual_elements.index('C')]
        manual_gas_carbon = manual_carbon * np.array([phase == 'gas' for phase, _ in manual_pairs], dtype=float)
        manual_solid_carbon = np.zeros(len(manual_pairs))
        USE_SOLID_CARBON_CONSTRAINT = SOLID_CARBON_FRACTION is not None
        if USE_SOLID_CARBON_CONSTRAINT:
            if SOLID_CARBON_TARGET == 'CHAR':
                if ('solid', 'carbon') not in manual_pairs:
                    raise ValueError("Select solid 'carbon' to constrain CHAR.")
                j = manual_pairs.index(('solid', 'carbon'))
                manual_solid_carbon[j] = manual_carbon[j]
            elif SOLID_CARBON_TARGET == 'all solids':
                manual_solid_carbon = manual_carbon * np.array([phase == 'solid' for phase, _ in manual_pairs], dtype=float)
            else:
                raise ValueError("SOLID_CARBON_TARGET must be 'CHAR' or 'all solids'.")
        manual_charge = np.array([ION_DATA[species]['charge'] if phase == 'aqueous' and species in ION_DATA else 0.0 for phase, species in manual_pairs], dtype=float)
        # Detect missing carriers before EOS calls; the full LP then checks
        # stoichiometric feasibility including counterions and carbon targets.
        missing = [e for e, b, row in zip(manual_elements, manual_b, manual_A)
                   if b > 0 and not np.any(row > 0)]
        if missing:
            raise ValueError('No allowed product carrier for feed elements: ' + ', '.join(missing)
                             + '. Add database-backed neutral gas species or solids using ion_product_species.')
        constraint_rows = list(manual_A) + [manual_charge]
        constraint_targets = list(manual_b) + [0.0]
        carbon_constraints = []
        if GAS_CARBON_FRACTION is not None:
            fraction = float(GAS_CARBON_FRACTION)
            if not np.isfinite(fraction) or not 0 <= fraction <= 1:
                raise ValueError('Gas-carbon fraction must be in [0, 1].')
            constraint_rows.append(manual_gas_carbon)
            constraint_targets.append(fraction * feed_carbon)
            carbon_constraints.append((manual_gas_carbon, fraction))
        if USE_SOLID_CARBON_CONSTRAINT:
            fraction = float(SOLID_CARBON_FRACTION)
            if not np.isfinite(fraction) or not 0 <= fraction <= 1:
                raise ValueError('Solid-carbon fraction must be in [0, 1].')
            if GAS_CARBON_FRACTION is not None and float(GAS_CARBON_FRACTION) + fraction > 1.0 + 1e-10:
                raise ValueError('Gas-carbon plus solid-carbon fractions exceed one.')
            constraint_rows.append(manual_solid_carbon)
            constraint_targets.append(fraction * feed_carbon)
            carbon_constraints.append((manual_solid_carbon, fraction))
        # Optional prescribed carbon allocation to individual phases.
        phase_fractions = settings.get('phase_carbon_fractions', {})
        if set(phase_fractions) - set(valid_phases):
            raise ValueError('Unknown phase in phase_carbon_fractions')
        fractions = {phase: float(value) for phase, value in phase_fractions.items()}
        if any(not np.isfinite(value) or not 0 <= value <= 1 for value in fractions.values()):
            raise ValueError('Phase-carbon fractions must be finite and in [0, 1]')
        if sum(fractions.values()) > 1 + 1e-10:
            raise ValueError('Phase-carbon fractions exceed one')
        if 'gas' in fractions and GAS_CARBON_FRACTION is not None and not np.isclose(fractions['gas'], GAS_CARBON_FRACTION, atol=1e-10, rtol=0):
            raise ValueError('Conflicting gas-carbon constraints')
        for phase, fraction in fractions.items():
            allocation_row = manual_carbon * np.array([p == phase for p, _ in manual_pairs], dtype=float)
            constraint_rows.append(allocation_row)
            constraint_targets.append(fraction * feed_carbon)
            carbon_constraints.append((allocation_row, fraction))
        manual_A_full = np.asarray(constraint_rows, dtype=float)
        manual_b_full = np.asarray(constraint_targets, dtype=float)
        lower = np.zeros(len(manual_pairs))
        for pair, value in MINIMUM_AMOUNTS.items():
            if pair not in manual_pairs:
                raise ValueError(f'Lower-bound species not selected: {pair}')
            if not np.isfinite(value) or value < 0:
                raise ValueError(f'Invalid lower bound for {pair}.')
            lower[manual_pairs.index(pair)] = value
        if manual_phase_species['aqueous']:
            if 'water' not in manual_phase_species['aqueous']:
                raise ValueError('The introduced aqueous model requires solvent water.')
            water_index = manual_pairs.index(('aqueous', 'water'))
            lower[water_index] = max(lower[water_index], AQUEOUS_WATER_FLOOR)
        fixed_zero = np.zeros(len(manual_pairs), dtype=bool)
        # Species containing an element absent from the feed must be zero.
        for row, inventory in zip(manual_A, manual_b):
            if inventory == 0.0:
                fixed_zero |= row > 0
        for allocation_row, fraction in carbon_constraints:
            if fraction == 0.0:
                fixed_zero |= allocation_row > 0
            elif fraction == 1.0:
                fixed_zero |= (manual_carbon > 0) & (allocation_row == 0)
        if np.any(lower[fixed_zero] > 0):
            raise ValueError('A minimum amount conflicts with a carbon constraint.')
        free = np.flatnonzero(~fixed_zero)
        if not len(free):
            raise ValueError('No free species remain.')

        def expand_manual_amounts(x):
            amounts = np.zeros(len(manual_pairs))
            amounts[free] = x
            return amounts
        reduced_A_full = manual_A_full[:, free]
        _, R_qr, pivots = qr(reduced_A_full.T, pivoting=True, mode='economic')
        diagonal = np.abs(np.diag(R_qr))
        rank = int(np.sum(diagonal > 1e-10 * max(1.0, float(diagonal.max()) if diagonal.size else 0.0)))
        selected_rows = np.sort(pivots[:rank])
        manual_A_opt = reduced_A_full[selected_rows]
        manual_b_opt = manual_b_full[selected_rows]
        manual_bounds = [(float(lower[j]), None) for j in free]
        start_objective = np.zeros(len(free))
        if settings.get('gas_water_initialization', False):
            for k, j in enumerate(free):
                if manual_pairs[j] == ('gas', 'water'):
                    start_objective[k] = -1.0
        start_bounds = list(manual_bounds)
        if LIQUID_INITIALIZATION and ('aqueous', 'water') in manual_pairs:
            for k, j in enumerate(free):
                if manual_pairs[j] == ('aqueous', 'water'):
                    start_objective[k] = -1.0
                elif manual_pairs[j] == ('gas', 'water'):
                    lower_bound, _ = start_bounds[k]
                    cap = max(lower_bound, INITIAL_GAS_WATER_FRACTION * water_mol_per_C)
                    start_bounds[k] = (lower_bound, cap)
        # Lower-bound carbon allocation is an INEQUALITY throughout the LP,
        # optimizer and final validation. It is never QR-reduced as an equality.
        gas_min_row = manual_gas_carbon[free]
        gas_min_target = (GAS_CARBON_MIN_FRACTION or 0.0) * feed_carbon
        gas_min_lp = ({'A_ub': -gas_min_row[None, :],
                       'b_ub': np.array([-gas_min_target])}
                      if GAS_CARBON_MIN_FRACTION is not None else {})
        feasible = linprog(start_objective, A_eq=manual_A_opt, b_eq=manual_b_opt,
                           bounds=start_bounds, method='highs', **gas_min_lp)
        if not feasible.success:
            raise ValueError('Selected feed, species and constraints are infeasible: ' + feasible.message)
        initial_amounts = expand_manual_amounts(feasible.x)
        if INITIAL_PHASE_AMOUNTS is not None:
            initial_amounts = np.array([INITIAL_PHASE_AMOUNTS.get(pair, 0.0) for pair in manual_pairs], dtype=float)
            if set(INITIAL_PHASE_AMOUNTS) - set(manual_pairs):
                raise ValueError('Initial composition contains an unselected phase/species')
            if np.any(~np.isfinite(initial_amounts)) or np.any(initial_amounts < lower) or np.any(initial_amounts[fixed_zero] != 0):
                raise ValueError('Initial composition violates species bounds')
            feasible.x = initial_amounts[free].copy()
        if np.max(np.abs(manual_A_full @ initial_amounts - manual_b_full)) > 1e-08:
            raise ValueError('The selected species cannot satisfy all original balances and carbon constraints.')
        if GAS_CARBON_MIN_FRACTION is not None and manual_gas_carbon @ initial_amounts < gas_min_target - 1e-8:
            raise ValueError('Initial composition violates minimum gas-carbon fraction')
        stage = 'thermochemistry / aqueous validity'
        _refresh_state(TEMPERATURE_C, PRESSURE_BAR, electrolyte=bool(manual_phase_species['aqueous']) and AQUEOUS_MODEL == 'electrolyte', phase_species=CASE_PHASE_SPECIES)
        if manual_phase_species['aqueous'] and AQUEOUS_MODEL == 'electrolyte':
            _aq_positions = np.array([AQUEOUS_SPECIES.index(s) for s in manual_phase_species['aqueous']], dtype=int)
        manual_mu0 = {}
        for phase, names in manual_phase_species.items():
            if not names:
                manual_mu0[phase] = np.array([])
            elif phase == 'aqueous' and AQUEOUS_MODEL == 'electrolyte':
                manual_mu0[phase] = np.array([])
            elif phase == 'solid':
                manual_mu0[phase], _ = solid_standard_chemical_potentials(names)
            else:
                manual_mu0[phase], _ = nasa_standard_chemical_potentials(names)

        def manual_phase_G(amounts):
            terms = {}
            for phase, names in manual_phase_species.items():
                indices = manual_indices[phase]
                phase_amounts = amounts[indices]
                if not names or phase_amounts.sum() <= 0.0:
                    terms[phase] = 0.0
                elif phase == 'aqueous':
                    if AQUEOUS_MODEL == 'molecular':
                        terms[phase] = float(molecular_gibbs(phase_amounts, names, manual_mu0[phase], 'liquid'))
                    else:
                        full_aqueous = np.zeros(len(AQUEOUS_SPECIES))
                        full_aqueous[_aq_positions] = phase_amounts
                        terms[phase] = float(aqueous_gibbs(full_aqueous))
                elif phase == 'solid':
                    terms[phase] = float(phase_amounts @ manual_mu0[phase])
                else:
                    phase_type = 'vapour' if phase == 'gas' else 'liquid'
                    terms[phase] = float(molecular_gibbs(phase_amounts, names, manual_mu0[phase], phase_type))
            if not np.all(np.isfinite(list(terms.values()))):
                raise ValueError('Nonfinite phase Gibbs energy.')
            return terms

        def manual_total_G(amounts):
            if np.min(amounts) < -1e-10:
                raise ValueError('Negative species amount.')
            terms = manual_phase_G(np.maximum(amounts, 0.0))
            return float(sum(terms.values()))
        stage = 'initial Gibbs evaluation'
        initial_G = manual_total_G(initial_amounts)
        
        
        energy_scale = float(
            settings.get(
                'optimizer_energy_scale',
                max(10000.0, abs(initial_G))
            )
        )
        if not np.isfinite(energy_scale) or energy_scale <= 0:
            raise ValueError(
                'optimizer_energy_scale must be finite and positive'
            )
        
        eos_errors = []

        def manual_objective(x):
            amounts = expand_manual_amounts(x)
            try:
                return (manual_total_G(amounts) - initial_G) / energy_scale
            except ValueError as exc:
                message = str(exc)
                if not any((text in message for text in ('Invalid PC-SAFT fugacities', 'Invalid ePC-SAFT fugacity coefficients', 'Nonfinite phase Gibbs energy'))):
                    raise
                if len(eos_errors) < 10:
                    eos_errors.append(message)
                return float(1000000.0 + np.sum((x - feasible.x) ** 2))
        stage = 'optimization / final validation'
        optimizer_constraints = [{'type':'eq', 'fun':lambda x: manual_A_opt @ x - manual_b_opt,
                                  'jac':lambda x: manual_A_opt}]
        if GAS_CARBON_MIN_FRACTION is not None:
            optimizer_constraints.append({'type':'ineq',
                'fun':lambda x: gas_min_row @ x - gas_min_target,
                'jac':lambda x: gas_min_row})
        # settings['optimizer'] selects an htg_optimizer backend (exact gradients,
        # IPOPT, multistart); None keeps the original finite-difference SLSQP.
        OPTIMIZER_SETTINGS = settings.get('optimizer')
        optimizer_diagnostics = {'method': 'slsqp-fd'}
        if OPTIMIZER_SETTINGS is None:
            optimized = minimize(manual_objective, feasible.x, method='SLSQP', bounds=manual_bounds, constraints=optimizer_constraints, options={'maxiter': MAX_ITERATIONS, 'ftol': float(settings.get('optimizer_ftol', 1e-09)), 'eps': float(settings.get('optimizer_eps', 1e-07)), 'disp': False})
        else:
            import htg_optimizer
            use_min = GAS_CARBON_MIN_FRACTION is not None
            optimizer_problem = htg_optimizer.GibbsProblem(
                ns=globals(), pairs=manual_pairs, free=free, phase_species=manual_phase_species,
                phase_indices=manual_indices, mu0=manual_mu0, aqueous_model=AQUEOUS_MODEL,
                aq_positions=_aq_positions if manual_phase_species['aqueous'] and AQUEOUS_MODEL == 'electrolyte' else None,
                A_eq=manual_A_opt, b_eq=manual_b_opt,
                ineq_row=gas_min_row if use_min else None, ineq_target=gas_min_target if use_min else None,
                lower=[low for low, _ in manual_bounds], x_start=feasible.x,
                G0=initial_G, scale=energy_scale, model_total_G=manual_total_G)
            optimized, optimizer_diagnostics = htg_optimizer.optimize(
                optimizer_problem, {'max_iterations': MAX_ITERATIONS, **OPTIMIZER_SETTINGS})
            eos_errors.extend(optimizer_problem.eos_errors[:max(0, 10 - len(eos_errors))])
        manual_amounts = expand_manual_amounts(optimized.x)
        manual_residual = float(np.max(np.abs(manual_A_full @ manual_amounts - manual_b_full)))
        gas_min_violation = (max(0.0, float(gas_min_target - manual_gas_carbon @ manual_amounts))
                             if GAS_CARBON_MIN_FRACTION is not None else 0.0)
        if gas_min_violation > 1e-7:
            raise RuntimeError(f'Minimum gas-carbon allocation violated by {gas_min_violation:.3e} mol C/mol dry-feed C')
        if not optimized.success or np.any(~np.isfinite(manual_amounts)) or np.min(manual_amounts) < -1e-09 or (manual_residual > 1e-07):
            raise RuntimeError(f'{optimized.message}; constraint residual={manual_residual:.3e}; recorded EOS errors={eos_errors}')
        manual_amounts = np.maximum(manual_amounts, 0.0)
        manual_G_terms = manual_phase_G(manual_amounts)
        manual_final_G = float(sum(manual_G_terms.values()))
        if manual_final_G > initial_G + 0.0001:
            raise RuntimeError('Final Gibbs energy exceeds the feasible start.')
        return {'gas_carbon_fraction_actual': float(manual_gas_carbon @ manual_amounts / feed_carbon), 'gas_carbon_min_fraction': GAS_CARBON_MIN_FRACTION, 'gas_carbon_min_violation': gas_min_violation, 'phase_allocation': phase_allocation, 'charge_residual': float(manual_charge @ manual_amounts), 'temperature_C': TEMPERATURE_C, 'pressure_bar': PRESSURE_BAR, 'amounts': dict(zip(manual_pairs, map(float, manual_amounts))), 'initial_amounts': dict(zip(manual_pairs, map(float, initial_amounts))), 'gibbs_J_per_mol_feed_C': manual_final_G, 'phase_gibbs_J_per_mol_feed_C': manual_G_terms, 'balance_residual': manual_residual, 'iterations': int(optimized.nit), 'eos_trial_errors': list(eos_errors), 'feed_C_mol_per_g_dry': manual_C_mol_per_g, 'elements': manual_elements, 'feed_atoms': manual_b.copy(), 'product_atoms': manual_A @ manual_amounts, 'optimizer': optimizer_diagnostics}
    except Exception as exc:
        raise RuntimeError(f'{stage}: {type(exc).__name__}: {exc}') from exc
    finally:
        globals().update(saved)
        for name in absent:
            globals().pop(name, None)

# %% cell 21
AW.update({'He':4.002602,'Li':6.94,'Be':9.0121831,'B':10.81,'F':18.998403163,
 'Ne':20.1797,'Mg':24.305,'Al':26.9815385,'Si':28.085,'P':30.973761998,
 'Cl':35.45,'Ar':39.948,'K':39.0983,'Ca':40.078,'Sc':44.955908,'Ti':47.867,
 'V':50.9415,'Cr':51.9961,'Mn':54.938044,'Fe':55.845,'Co':58.933194,
 'Ni':58.6934,'Cu':63.546,'Zn':65.38,'Ga':69.723,'Ge':72.630,'As':74.921595,
 'Se':78.971,'Br':79.904,'Kr':83.798,'Rb':85.4678,'Sr':87.62,'Y':88.90584,
 'Zr':91.224,'Nb':92.90637,'Mo':95.95,'Ru':101.07,'Rh':102.90550,
 'Pd':106.42,'Ag':107.8682,'Cd':112.414,'In':114.818,'Sn':118.710,
 'Sb':121.760,'Te':127.60,'I':126.90447,'Xe':131.293,'Cs':132.90545196,
 'Ba':137.327,'La':138.90547,'Ce':140.116,'Pr':140.90766,'Nd':144.242,
 'Sm':150.36,'Eu':151.964,'Gd':157.25,'Tb':158.92535,'Dy':162.500,
 'Ho':164.93033,'Er':167.259,'Tm':168.93422,'Yb':173.045,'Lu':174.9668,
 'Hf':178.49,'Ta':180.94788,'W':183.84,'Re':186.207,'Os':190.23,
 'Ir':192.217,'Pt':195.084,'Au':196.966569,'Hg':200.592,'Tl':204.38,
 'Pb':207.2,'Bi':208.98040,'Th':232.0377,'U':238.02891})
_DEFAULT_AQ_NEUTRALS = list(AQ_NEUTRAL_SPECIES)
_DEFAULT_AQ_IONS = list(AQ_ION_SPECIES)
_BASE_MOLECULAR_SPECIES = {s: dict(d) for s, d in MOLECULAR_SPECIES.items()}
_BASE_ORGANIC_SPECIES = {s: dict(d) for s, d in ORGANIC_SPECIES.items()}
_BASE_SOLID_DATA = {s: dict(d) for s, d in SOLID_DATA.items()}
_BASE_ION_DATA = {s: dict(d) for s, d in ION_DATA.items()}


def parse_nasa_header_elements(header):
    elements = {}
    for offset in range(24, 44, 5):
        field = header[offset:offset + 5]
        symbol = field[:2].strip().capitalize()
        if not symbol:
            continue
        try:
            count = float(field[2:].strip())
        except ValueError:
            continue
        if count:
            elements[symbol] = elements.get(symbol, 0.0) + count
    return elements


def formula_elements(formula):
    text = str(formula).strip()
    if not re.fullmatch(r'(?:[A-Z][a-z]?\d*)+', text):
        return None
    result = {}
    for symbol, number in re.findall(r'([A-Z][a-z]?)(\d*)', text):
        result[symbol] = result.get(symbol, 0.0) + int(number or 1)
    return result


def _ion_identity(name, charge):
    text = re.sub(r'\s+', '', str(name)).replace('^', '')
    # The explicit charge avoids treating the 2 in VO2+ as a divalent charge.
    suffixes = [f'{abs(int(charge))}{"+" if charge > 0 else "-"}',
                f'{"+" if charge > 0 else "-"}{abs(int(charge))}',
                ('+' if charge > 0 else '-') * abs(int(charge))]
    if text.upper() in ('NA+','K+','LI+','RB+','CS+','CL-','F-','BR-','I-'):
        return text.rstrip('+-').casefold(), int(charge)
    for suffix in suffixes:
        if text.endswith(suffix):
            return text[:-len(suffix)].casefold(), int(charge)
    return text.rstrip('+-').casefold(), int(charge)


def discover_database_species(species_overrides=None):
    """Build registries and a complete coverage catalogue from the supplied records.

    Overrides map a CSV species name to {'nasa': 'exact NASA name'} where
    an isomer or a nonstandard identifier cannot be resolved unambiguously.
    No NASA/EOS parameter values are invented for incomplete species.
    """
    global ALL_MOLECULAR_DATA, database_species_catalog, ion_epcsaft_table, ion_hkf_table
    overrides = {str(k).strip().casefold(): dict(v) for k, v in (species_overrides or {}).items()}
    original_molecules = {**MOLECULAR_SPECIES, **ORGANIC_SPECIES}
    nasa_by_formula = {}
    for nasa_name, data in nasa_db.items():
        key = tuple(sorted(data.get('elements', {}).items()))
        if key:
            nasa_by_formula.setdefault(key, []).append(nasa_name)
    pc = pcsaft_raw[pcsaft_raw['record_type'].astype(str).eq('pure_component')]
    epc = epcsaft_raw[epcsaft_raw['record_type'].astype(str).eq('species_parameter')
                      & epcsaft_raw['Charge'].fillna(0).eq(0)]
    records = []
    for frame in (pc, epc):
        for _, row in frame.iterrows():
            if pd.isna(row.get('Species')) or pd.isna(row.get('Formula')):
                continue
            name = str(row.Species).strip().casefold()
            if name not in {r[0] for r in records}:
                records.append((name, str(row.Formula).strip()))
    catalogue = []
    for name, formula in records:
        elements = formula_elements(formula)
        override = overrides.get(name, {})
        mapped = override.get('nasa') or original_molecules.get(name, {}).get('nasa')
        if not mapped:
            exact = next((key for key in nasa_db if key.casefold() in (name, formula.casefold())), None)
            if exact:
                mapped = exact
            elif elements:
                candidates = nasa_by_formula.get(tuple(sorted(elements.items())), [])
                same_formula_names = [other for other, f in records if formula_elements(f) == elements]
                if len(candidates) == 1 and len(same_formula_names) == 1:
                    mapped = candidates[0]
        if override.get('nasa') and mapped not in nasa_db:
            raise ValueError(f'NASA override for {name} not found: {mapped}')
        if mapped in nasa_db:
            elements = nasa_db[mapped].get('elements') or elements
        if not elements:
            catalogue.append(dict(Species=name, Formula=formula, Phase='molecular', Ready=False,
                                  Reason='Elemental composition unresolved', NASA=mapped))
            continue
        PCSAFT_FORMULA[name] = formula
        SPECIES_ALIASES_6[name] = [name, formula]
        data = first_complete_epcsaft_set(name, T) or first_complete_pcsaft_set(name, T)
        ready = bool(mapped in nasa_db and data and all(e in AW for e in elements))
        reason = ('ready' if ready else 'NASA mapping missing/ambiguous' if mapped not in nasa_db
                  else 'Complete SAFT parameters missing' if not data else 'Atomic weight missing')
        catalogue.append(dict(Species=name, Formula=formula, Phase='molecular', Ready=ready,
                              Reason=reason, NASA=mapped))
        if ready:
            entry = dict(formula=formula, nasa=mapped, elements=elements, clapeyron=name)
            MOLECULAR_SPECIES[name] = entry
            ORGANIC_SPECIES[name] = entry
    # Every NASA entry remains selectable as a dry-feed species. Molecular
    # EOS selection requires a matching SAFT record; condensed NASA entries
    # and explicitly marked overrides can be selected as pure solids.
    covered = {entry.get('nasa') for entry in MOLECULAR_SPECIES.values()}
    for name, data in nasa_db.items():
        elements = data.get('elements')
        if not elements or not all(e in AW for e in elements):
            continue
        key = name if name not in original_molecules else name.casefold()
        entry = dict(formula=name, nasa=name, elements=elements)
        if name not in covered:
            MOLECULAR_SPECIES.setdefault(key, entry)
            catalogue.append(dict(Species=key, Formula=name, Phase='molecular', Ready=False,
                Reason='NASA available; no resolved SAFT mapping. Supply species_overrides or EOS data.', NASA=name))
        if data.get('phase') == 'S' or overrides.get(key.casefold(), {}).get('phase') == 'solid':
            SOLID_DATA.setdefault(key, entry)
            catalogue.append(dict(Species=key, Formula=name, Phase='solid', Ready=True, Reason='NASA solid', NASA=name))
    for name, entry in SOLID_DATA.items():
        catalogue.append(dict(Species=name, Formula=entry['formula'], Phase='solid', Ready=entry['nasa'] in nasa_db, Reason='NASA solid', NASA=entry['nasa']))
    # Join ionic EOS rows to exact SUPCRT charge/formula identities.
    xml = list(supcrt_root.iter('Species'))
    ready_ion_epc, ready_ion_hkf = [], []
    charged = epcsaft_raw[epcsaft_raw['record_type'].astype(str).eq('species_parameter')
                           & epcsaft_raw['Charge'].fillna(0).ne(0)]
    for name, rows in charged.groupby('Species', sort=False):
        name = str(name).strip()
        charge = int(rows.Charge.dropna().iloc[0])
        formula = str(rows.Formula.dropna().iloc[0])
        ION_CHARGE[name] = charge
        ION_FORMULA[name] = formula
        matches = [node for node in xml if node.findtext('Type') == 'Aqueous'
                   and float(node.findtext('Charge') or 0) == charge
                   and any(_ion_identity(candidate, charge) == _ion_identity(name, charge)
                           for candidate in (node.findtext('Name'), node.findtext('Formula')))]
        SUPCRT_ION_ALIASES[name] = [node.findtext('Name') for node in matches]
        epc_data = first_complete_ion_epcsaft_set_7(name)
        hkf_data = first_complete_hkf_set_7(name) if matches else None
        node = next((node for node in matches if read_hkf_candidate_7(node)), None)
        elements = ({e: float(v) for e, v in re.findall(r'([A-Z][a-z]?)\(([-+\d.eE]+)\)', node.findtext('Elements') or '')}
                    if node is not None else None)
        ready = bool(epc_data and hkf_data and elements and all(e in AW for e in elements))
        catalogue.append(dict(Species=name, Formula=formula, Phase='aqueous ion', Ready=ready,
                              Reason='ready' if ready else 'Complete matching ionic EOS/HKF/element data missing', NASA=None))
        if ready:
            ION_DATA[name] = dict(formula=formula, charge=charge, supcrt=hkf_data['SUPCRT name'], elements=elements)
            ready_ion_epc.append(epc_data); ready_ion_hkf.append(hkf_data)
    ion_epcsaft_table = pd.DataFrame(ready_ion_epc)
    ion_hkf_table = pd.DataFrame(ready_ion_hkf)
    # Include XML-only species in the catalogue, even when this EOS has no
    # matching parameter set or condensed-phase adapter.
    for node in xml:
        name, kind = node.findtext('Name'), node.findtext('Type')
        if any(row['Species'] == name and row['Phase'] == 'aqueous ion' for row in catalogue):
            continue
        catalogue.append(dict(Species=name, Formula=node.findtext('Formula'), Phase=f'SUPCRT {kind}',
            Ready=False, Reason='SUPCRT reference present; matching EOS/reference adapter required', NASA=None))
    ALL_MOLECULAR_DATA = {**MOLECULAR_SPECIES, **ORGANIC_SPECIES}
    database_species_catalog = pd.DataFrame(catalogue).drop_duplicates(['Species', 'Phase'])
    register_supcrt_minerals()


def resolve_molecular_parameters(temperature_K, species=None):
    resolved = {}
    names = list(species) if species is not None else database_species_catalog.loc[
        database_species_catalog.Phase.eq('molecular') & database_species_catalog.Ready, 'Species'].tolist()
    for name in dict.fromkeys(names):
        if name not in ALL_MOLECULAR_DATA:
            raise ValueError(f'{name}: not found in the molecular database registry')
        data = first_complete_epcsaft_set(name, temperature_K) or first_complete_pcsaft_set(name, temperature_K)
        if data is None:
            raise ValueError(f'{name}: NASA data alone is insufficient for PC-SAFT; complete SAFT parameters/mapping required')
        data['Mw'] = sum(AW[e] * v for e, v in ALL_MOLECULAR_DATA[name]['elements'].items())
        if name == 'water' and data.get('association_scheme') in ('-', '', None):
            data['association_scheme'] = '2B'
        resolved[name] = data
    return resolved


def initialize_model(database_paths, temperature_C=350.0, pressure_bar=250.0, species_overrides=None):
    """Read all four databases, discover complete species, and apply explicit mappings."""
    global cl, jl, pch, nasa_db, pcsaft_raw, epcsaft_raw, supcrt_root
    global PC_COL_M_6, PC_COL_SIGMA_6, PC_COL_EPSILON_6
    global PC_COL_EPSILONAB_6, PC_COL_KAPPAAB_6, PC_COL_ASSOCIATION_6
    global MOLECULAR_SPECIES, ORGANIC_SPECIES, SOLID_DATA, ION_DATA, ALL_MOLECULAR_DATA
    global resolved_pcsaft_lookup, T, P, MODEL_INITIALIZED
    import pyclapeyron as cl
    from juliacall import Main as jl
    import pychnosz as pch
    paths = {key: Path(database_paths[key]) for key in ('nasa','pcsaft','epcsaft','supcrt')}
    for key, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f'{key} database not found: {path}')
    nasa_db = read_nasa_database(paths['nasa'])
    pcsaft_raw = pd.read_csv(paths['pcsaft'])
    epcsaft_raw = pd.read_csv(paths['epcsaft'])
    supcrt_root = ET.parse(paths['supcrt']).getroot()
    PC_COL_M_6, PC_COL_SIGMA_6, PC_COL_EPSILON_6 = 'm','sigma_A','u_over_kB_K'
    PC_COL_EPSILONAB_6 = first_existing_column_6(pcsaft_raw, ['epsilonAB_over_kB_K','epsilonAB_K','epsilonAB/kB','epsilonAB'])
    PC_COL_KAPPAAB_6 = first_existing_column_6(pcsaft_raw, ['kappaAB','association_volume_kappaAB'])
    PC_COL_ASSOCIATION_6 = first_existing_column_6(pcsaft_raw, ['Association_scheme','association_scheme','Association'])
    T, P = float(temperature_C) + 273.15, float(pressure_bar) * 1e5
    MOLECULAR_SPECIES = {s:dict(d) for s,d in _BASE_MOLECULAR_SPECIES.items()}
    ORGANIC_SPECIES = {s:dict(d) for s,d in _BASE_ORGANIC_SPECIES.items()}
    SOLID_DATA = {s:dict(d) for s,d in _BASE_SOLID_DATA.items()}
    ION_DATA = {s:dict(d) for s,d in _BASE_ION_DATA.items()}
    ALL_MOLECULAR_DATA = {**MOLECULAR_SPECIES, **ORGANIC_SPECIES}
    discover_database_species(species_overrides)
    resolved_pcsaft_lookup = resolve_molecular_parameters(T)
    configure_aqueous_species({'aqueous': _DEFAULT_AQ_NEUTRALS + _DEFAULT_AQ_IONS})
    MODEL_INITIALIZED = True
    return available_species()


def available_species(include_unavailable=False):
    if not globals().get('MODEL_INITIALIZED') or 'database_species_catalog' not in globals():
        return pd.DataFrame()
    catalogue = database_species_catalog.copy()
    if include_unavailable:
        return catalogue
    return catalogue[catalogue.Ready].reset_index(drop=True)


def configure_aqueous_species(phase_species):
    """Retain baseline probe components and add every selected supported solute."""
    global AQ_NEUTRAL_SPECIES, AQ_ION_SPECIES, ION_SPECIES, AQUEOUS_SPECIES, AQUEOUS_FORMULA
    selected = phase_species.get('aqueous', [])
    AQ_NEUTRAL_SPECIES = list(dict.fromkeys(_DEFAULT_AQ_NEUTRALS + [s for s in selected if s not in ION_DATA]))
    AQ_ION_SPECIES = list(dict.fromkeys(_DEFAULT_AQ_IONS + [s for s in selected if s in ION_DATA]))
    ION_SPECIES = list(AQ_ION_SPECIES)
    AQUEOUS_SPECIES = AQ_NEUTRAL_SPECIES + AQ_ION_SPECIES
    AQUEOUS_FORMULA = {s: ALL_MOLECULAR_DATA[s]['formula'] for s in AQ_NEUTRAL_SPECIES}
    AQUEOUS_FORMULA.update({s: ION_DATA[s]['formula'] for s in AQ_ION_SPECIES})


def _balanced_ion_probe_weights():
    z = np.array([ION_DATA[s]['charge'] for s in AQ_ION_SPECIES], dtype=float)
    weights = np.array([{'Na+':3., 'HCO3-':1., 'CO3^2-':1.}.get(s, 1.) for s in AQ_ION_SPECIES])
    imbalance = float(z @ weights)
    candidates = np.flatnonzero(z * imbalance < 0)
    if abs(imbalance) > 1e-14:
        if not candidates.size:
            raise ValueError('Both positive and negative ionic carriers are needed for EOS probes')
        j = candidates[0]; weights[j] += -imbalance / z[j]
    return weights


def phi_electrolyte(n):
    n = np.asarray(n, dtype=float).copy()
    if n.shape != (len(AQUEOUS_SPECIES),) or np.any(n < 0):
        raise ValueError('Invalid aqueous composition vector')
    split = len(AQ_NEUTRAL_SPECIES)
    z = np.array([ION_DATA[s]['charge'] for s in AQ_ION_SPECIES], dtype=float)
    imbalance = float(n[split:] @ z)
    smoothing = np.hypot(imbalance, 2e-14)
    pos, neg = np.flatnonzero(z > 0), np.flatnonzero(z < 0)
    if not pos.size or not neg.size:
        raise ValueError('Aqueous EOS probe requires positive and negative ionic carriers')
    n[split + pos[0]] += (smoothing - imbalance) / (2 * z[pos[0]])
    n[split + neg[0]] += (smoothing + imbalance) / (-2 * z[neg[0]])
    n[:split] += 1e-14
    n[split:] += _balanced_ion_probe_weights() * 1e-14
    phi = np.asarray(cl.fugacity_coefficient(aqueous_epcsaft_model, P, T, n, phase='liquid'), dtype=float).reshape(-1)
    if len(phi) != len(n) or np.any(~np.isfinite(phi)) or np.any(phi <= 0):
        raise ValueError('Invalid ePC-SAFT fugacity coefficients')
    return phi


def infinite_dilution_reference():
    logs = []
    split = len(AQ_NEUTRAL_SPECIES)
    for scale in (1e-5, 1e-6, 1e-7):
        n = np.zeros(len(AQUEOUS_SPECIES)); n[AQUEOUS_SPECIES.index('water')] = 10.
        n[split:] = scale * _balanced_ion_probe_weights()
        logs.append(np.log(phi_electrolyte(n)))
    x = np.sqrt(np.array([1e-5, 1e-6, 1e-7]))
    intercept = np.array([np.polyfit(x, np.array(logs)[:,j], 1)[1] for j in range(len(AQUEOUS_SPECIES))])
    estimate = logs[-1][split:] - (logs[-2][split:] - logs[-1][split:]) * x[-1] / (x[-2] - x[-1])
    if np.max(np.abs(intercept[split:] - estimate)) > .02:
        raise ValueError('Ion dilution extrapolation did not converge')
    return intercept

# %% cell 23
def prepare_phase_allocation(settings):
    """Return a copied case with the aqueous cut-off applied before EOS setup.

    No ion is renamed or discarded: elemental conservation determines the
    neutral products. Candidate names must refer to existing database records.
    The cut-off is a modelling rule, not a calculated mixture critical point.
    """
    from copy import deepcopy
    case = deepcopy(settings)
    phases = case['phase_species']
    disabled = case.get('aqueous_enabled', True) is False
    cutoff = case.get('aqueous_cutoff_C', case.get('aqueous_cutoff_temperature_C'))
    if 'aqueous_cutoff_C' in case and 'aqueous_cutoff_temperature_C' in case:
        if case['aqueous_cutoff_C'] != case['aqueous_cutoff_temperature_C']:
            raise ValueError('Conflicting aqueous cut-off settings; use aqueous_cutoff_C')
    if cutoff is not None:
        cutoff = float(cutoff)
        if not np.isfinite(cutoff) or cutoff <= -273.15:
            raise ValueError('aqueous_cutoff_temperature_C must be finite and above absolute zero')
        disabled = disabled or float(case['temperature_C']) >= cutoff
    # These are eligible equilibrium products at BOTH sides of the cut-off.
    # Do not force precipitation only after removal of the aqueous phase.
    additions = case.get('ion_product_species', {})
    if set(additions) - {'gas', 'solid'}:
        raise ValueError('ion_product_species supports only gas and solid')
    for phase, names in additions.items():
        registry = SOLID_DATA if phase == 'solid' else ALL_MOLECULAR_DATA
        for name in names:
            if name in ION_DATA or name not in registry:
                raise ValueError(f'Unknown or charged {phase} ion product: {name}')
            if phase == 'solid':
                record = nasa_db.get(registry[name]['nasa'], {})
                if registry[name].get('adapter') != 'supcrt_mineral' and record.get('phase') != 'S':
                    raise ValueError(f'{name}: ion-product solid requires a NASA solid or supported SUPCRT mineral record')
        phases[phase] = list(dict.fromkeys(phases.get(phase, []) + list(names)))
    removed = list(phases.get('aqueous', [])) if disabled else []
    if disabled:
        phases['aqueous'] = []
        # Ignore only aqueous bounds/initial guesses. Retain all feed atoms.
        case['minimum_amounts'] = {p: v for p, v in case.get('minimum_amounts', {}).items() if p[0] != 'aqueous'}
        # A supplied composition cannot be projected across phases without
        # rebalancing; use the existing feasible LP initialization instead.
        case['initial_phase_amounts'] = None
        case['liquid_initialization'] = False
        allocation = case.get('phase_carbon_fractions', {})
        if float(allocation.get('aqueous', 0.0)) != 0.0:
            raise ValueError('Aqueous cut-off conflicts with nonzero aqueous phase_carbon_fractions; revise carbon targets explicitly')
    for phase, names in phases.items():
        if phase != 'aqueous':
            forbidden = [s for s in names if s in ION_DATA]
            if forbidden:
                raise ValueError(f'Ions are only allowed in aqueous: {phase}: {forbidden}')
    report = {'aqueous_disabled': disabled,
              'aqueous_cutoff_temperature_C': cutoff,
              'removed_aqueous_species': removed,
              'effective_phase_species': deepcopy(phases),
              'initial_composition_reset': disabled and settings.get('initial_phase_amounts') is not None}
    return case, report


def ion_product_candidates(elements=None, temperature_C=None):
    """Neutral molecular gas records and pure solids with existing thermo data."""
    target = set(elements or ['Na', 'K', 'Cs', 'Li', 'Cl', 'Mg'])
    rows = []
    for phase, registry in [('gas', ALL_MOLECULAR_DATA), ('solid', SOLID_DATA)]:
        for name, d in registry.items():
            if name in ION_DATA or not target.intersection(d['elements']):
                continue
            reason = ''
            if d.get('adapter') == 'supcrt_mineral':
                high = d['supcrt_mineral']['Tmax']
                if temperature_C is not None and not 298.15 <= temperature_C+273.15 <= high:
                    reason = 'Outside SUPCRT mineral temperature range'
                source = 'user SUPCRT XML'
            else:
                record = nasa_db.get(d.get('nasa'), {})
                if record.get('phase') != ('S' if phase == 'solid' else 'G'):
                    continue
                source = 'user NASA database'
                if temperature_C is not None and not record['T_low'] <= temperature_C+273.15 <= record['T_high']:
                    reason = 'Outside NASA temperature range'
            if phase == 'gas':
                try:
                    resolve_molecular_parameters(temperature_C+273.15 if temperature_C is not None else T, [name])
                except (ValueError, KeyError) as exc:
                    reason = str(exc)
            rows.append({'Species':name, 'Phase':phase, 'Elements':dict(d['elements']),
                         'Source':source, 'Ready':not bool(reason), 'Reason':reason})
    return pd.DataFrame(rows, columns=['Species','Phase','Elements','Source','Ready','Reason'])

# %% cell 26
def register_supcrt_minerals():
    """Discover pure minerals from the CURRENT XML; no external data are added.

    This adapter implements the single-range Maier-Kelley model only.
    Polymorphic records are reported as unavailable, rather than extrapolated.
    XML <Thermo/HKF> mineral records are not aqueous HKF records.
    """
    global database_species_catalog
    rows = []
    for node in supcrt_root.iter('Species'):
        if node.findtext('Type') != 'Mineral':
            continue
        name = node.findtext('Name')
        hkf = node.find('./Thermo/HKF')
        reason, params = '', {}
        try:
            if hkf is None:
                raise ValueError('Mineral HKF data missing')
            if int(hkf.findtext('NumPhaseTrans') or 0):
                raise ValueError('Polymorphic mineral: transition adapter not implemented')
            for tag in ('Gf', 'Hf', 'Sr', 'Vr', 'a', 'b', 'c', 'Tmax'):
                params[tag] = float(hkf.findtext(tag))
            if not np.all(np.isfinite(list(params.values()))) or params['Vr'] <= 0 or params['Tmax'] < 298.15:
                raise ValueError('Invalid/missing mineral parameters')
            for tag in ('Gf', 'Hf', 'Sr', 'a', 'b', 'c'):
                if 'cal' not in hkf.find(tag).get('units', ''):
                    raise ValueError(f'Unexpected units for mineral {tag}')
            if hkf.find('Vr').get('units') != 'cm3/mol':
                raise ValueError('Expected Vr in cm3/mol (do not use mislabeled MolarVolume)')
            elements = {e: float(v) for e,v in re.findall(r'([A-Z][a-z]?)\(([-+\d.eE]+)\)', node.findtext('Elements') or '')}
            if not elements or any(e not in AW for e in elements):
                raise ValueError('Element composition or atomic weight missing')
            # Retain exact XML spelling and reference identifiers.
            SOLID_DATA[name] = dict(formula=node.findtext('Formula'), elements=elements,
                adapter='supcrt_mineral', nasa=None, supcrt_mineral=params,
                source='supcrt XML supplied to initialize_model',
                references=hkf.findtext('References'), date=hkf.findtext('Date'),
                hydrate=('*' in (node.findtext('Formula') or '') and 'H2O' in (node.findtext('Formula') or '')))
        except (TypeError, ValueError) as exc:
            reason = str(exc)
        rows.append(dict(Species=name, Formula=node.findtext('Formula'), Phase='solid',
                         Ready=not bool(reason), Reason=reason or 'SUPCRT mineral: single-range Maier-Kelley', NASA=None))
    # Replace misleading XML-only rows; keep NASA entries and aqueous inventory.
    names = {r['Species'] for r in rows}
    keep = ~((database_species_catalog.Species.isin(names)) &
             database_species_catalog.Phase.isin(['SUPCRT Mineral', 'solid']))
    database_species_catalog = pd.concat([database_species_catalog[keep], pd.DataFrame(rows)], ignore_index=True)
    return pd.DataFrame(rows)


def supcrt_mineral_properties(name, temperature_K, pressure_Pa):
    """Hf/S reference compatible with NASA; Cp=a+b*T+c/T**2.

    H(T)=Hf(Tr)+integral(Cp dT), S(T)=Sr+integral(Cp/T dT).
    G=H-T*S+Vr*(P-Pr), using a constant mineral volume.
    Cal -> J: 4.184; cm3 -> m3: 1e-6; Tr=298.15 K, Pr=1e5 Pa.
    Tmax is an equation limit, NOT a hydrate stability temperature.
    """
    d = SOLID_DATA[name]['supcrt_mineral']
    t, tr, p = float(temperature_K), 298.15, float(pressure_Pa)
    if not np.isfinite(t) or not tr <= t <= d['Tmax']:
        raise ValueError(f"{name}: T={t:g} K outside supported [{tr:g}, {d['Tmax']:g}] K")
    if not np.isfinite(p) or p <= 0:
        raise ValueError('Invalid mineral pressure')
    a,b,c = [d[k] * 4.184 for k in ('a','b','c')]
    dh = a*(t-tr) + b*(t*t-tr*tr)/2 - c*(1/t-1/tr)
    ds = a*np.log(t/tr) + b*(t-tr) - c*(t**-2-tr**-2)/2
    h, entropy = d['Hf']*4.184 + dh, d['Sr']*4.184 + ds
    pv = d['Vr'] * 1e-6 * (p-1e5)
    return {'G': h-t*entropy+pv, 'H':h+pv, 'S':entropy,
            'Cp':a+b*t+c/t**2, 'V_m3_per_mol':d['Vr']*1e-6}


def mineral_inventory(temperature_C=None, elements=None):
    """Inspect XML-derived mineral records, including hydrate composition."""
    rows = []
    for name, d in SOLID_DATA.items():
        if d.get('adapter') != 'supcrt_mineral':
            continue
        if elements and not set(elements).intersection(d['elements']):
            continue
        high = d['supcrt_mineral']['Tmax']
        rows.append({'Species':name, 'Formula':d['formula'], 'Elements':dict(d['elements']),
                     'Hydrate':d['hydrate'], 'Tmax_K':high,
                     'Within_T_range':temperature_C is None or 298.15 <= temperature_C+273.15 <= high,
                     'References_in_XML':d['references'], 'Date_in_XML':d['date']})
    return pd.DataFrame(rows)


def solid_standard_chemical_potentials(solid_species):
    values, rows = [], []
    for name in solid_species:
        d = SOLID_DATA[name]
        if d.get('adapter') == 'supcrt_mineral':
            props = supcrt_mineral_properties(name, T, P)
            value = props['G']
            rows.append({'Species':name, 'Source':'user SUPCRT XML', 'Phase':'pure solid',
                         'g0 [J/mol]':value, 'cp [J/mol/K]':props['Cp'],
                         'T high [K]':d['supcrt_mineral']['Tmax']})
        else:
            cp,h,entropy,value = nasa7_properties(nasa_db, get_nasa_name(name), T)
            rows.append({'Species':name, 'Source':'user NASA database', 'Phase':'pure solid',
                         'g0 [J/mol]':value, 'cp [J/mol/K]':cp,
                         'T high [K]':nasa_db[d['nasa']]['T_high']})
        values.append(value)
    return np.asarray(values), pd.DataFrame(rows)

# %% cell 28
SUPPORTS_GAS_CARBON_MIN_FRACTION = True

