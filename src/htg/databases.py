"""Readers for the four thermodynamic databases.

    nasa      NASA-7 polynomial records (Thermdat_WasteMech.txt)
    pcsaft    PC-SAFT pure-component parameters (PROJECT_PCSAFT.csv)
    epcsaft   ePC-SAFT species parameters, neutral and ionic (PROJECT_EPCSAFT.csv)
    supcrt    SUPCRT92 aqueous HKF and mineral records (supcrt92.xml)

The readers only parse files; interpreting the records is done in
``standard_states``, ``saft`` and ``species``.
"""
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .constants import R

DATABASE_FILES = {
    'nasa': 'Thermdat_WasteMech.txt',
    'pcsaft': 'PROJECT_PCSAFT.csv',
    'epcsaft': 'PROJECT_EPCSAFT.csv',
    'supcrt': 'supcrt92.xml',
}


@dataclass
class Databases:
    nasa: dict
    pcsaft: pd.DataFrame
    epcsaft: pd.DataFrame
    supcrt: ET.Element
    paths: dict = field(default_factory=dict)
    # PROJECT_PCSAFT.csv column names, resolved once.
    pc_columns: dict = field(default_factory=dict)

    @classmethod
    def load(cls, data_dir=None, paths=None):
        """Read all databases from ``data_dir`` (default file names) or explicit ``paths``."""
        if paths is None:
            data_dir = Path(data_dir)
            paths = {key: data_dir / name for key, name in DATABASE_FILES.items()}
        paths = {key: Path(paths[key]) for key in DATABASE_FILES}
        for key, path in paths.items():
            if not path.is_file():
                raise FileNotFoundError(f'{key} database not found: {path}')
        pcsaft = pd.read_csv(paths['pcsaft'])
        return cls(
            nasa=read_nasa_database(paths['nasa']),
            pcsaft=pcsaft,
            epcsaft=pd.read_csv(paths['epcsaft']),
            supcrt=ET.parse(paths['supcrt']).getroot(),
            paths=paths,
            pc_columns={
                'm': 'm', 'sigma': 'sigma_A', 'epsilon': 'u_over_kB_K',
                'epsilonAB': first_existing_column(pcsaft, ['epsilonAB_over_kB_K', 'epsilonAB_K',
                                                            'epsilonAB/kB', 'epsilonAB']),
                'kappaAB': first_existing_column(pcsaft, ['kappaAB', 'association_volume_kappaAB']),
                'association': first_existing_column(pcsaft, ['Association_scheme',
                                                              'association_scheme', 'Association']),
            },
        )


def normalize_name(value):
    """Case-folded name with subscript digits replaced, for matching CSV entries."""
    if pd.isna(value):
        return ''
    text = str(value).strip().casefold()
    replacements = {'₂': '2', '₃': '3', '₄': '4', '₅': '5', '₆': '6', '₇': '7', '₈': '8', '₉': '9'}
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


def first_existing_column(dataframe, candidates):
    normalized = {normalize_name(col): col for col in dataframe.columns}
    for candidate in candidates:
        key = normalize_name(candidate)
        if key in normalized:
            return normalized[key]
    return None


# ---------------------------------------------------------------- NASA-7
def nasa_float(text):
    text = str(text).strip().replace('D', 'E').replace('d', 'e')
    if not text:
        return np.nan
    return float(text)


def parse_nasa_header_elements(header):
    elements = {}
    for offset in range(24, 44, 5):
        field_ = header[offset:offset + 5]
        symbol = field_[:2].strip().capitalize()
        if not symbol:
            continue
        try:
            count = float(field_[2:].strip())
        except ValueError:
            continue
        if count:
            elements[symbol] = elements.get(symbol, 0.0) + count
    return elements


def read_nasa_database(path):
    """Parse 4-line NASA-7 records into {name: {T_low, T_mid, T_high, low, high, elements, phase}}."""
    nasa_db = {}
    i = 0
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    while i < len(lines):
        line1 = lines[i].rstrip('\n')
        if len(line1) < 80 or line1[79:80] != '1':
            i += 1
            continue
        if i + 3 >= len(lines):
            break
        line2 = lines[i + 1].rstrip('\n')
        line3 = lines[i + 2].rstrip('\n')
        line4 = lines[i + 3].rstrip('\n')
        if not (len(line2) >= 80 and len(line3) >= 80 and len(line4) >= 80
                and line2[79:80] == '2' and line3[79:80] == '3' and line4[79:80] == '4'):
            i += 1
            continue
        species_name = line1[0:18].strip()
        try:
            T_low = nasa_float(line1[45:55])
            T_high = nasa_float(line1[55:65])
            T_mid = nasa_float(line1[65:73])
        except Exception as err:
            raise ValueError(f"Could not parse NASA temperature limits for species '{species_name}'."
                             f"\nHeader:\n{line1}") from err
        try:
            c2 = [nasa_float(line2[j:j + 15]) for j in range(0, 75, 15)]
            c3 = [nasa_float(line3[j:j + 15]) for j in range(0, 75, 15)]
            c4 = [nasa_float(line4[j:j + 15]) for j in range(0, 60, 15)]
            high = np.array([c2[0], c2[1], c2[2], c2[3], c2[4], c3[0], c3[1]], dtype=float)
            low = np.array([c3[2], c3[3], c3[4], c4[0], c4[1], c4[2], c4[3]], dtype=float)
        except Exception as err:
            raise ValueError(f"Could not parse NASA coefficients for species '{species_name}'.") from err
        nasa_db[species_name] = {'T_low': float(T_low), 'T_mid': float(T_mid),
                                 'T_high': float(T_high), 'low': low, 'high': high,
                                 'elements': parse_nasa_header_elements(line1),
                                 'phase': line1[44:45].strip().upper()}
        i += 4
    if not nasa_db:
        raise ValueError('No NASA-7 entries found')
    return nasa_db


def nasa7_coefficients(nasa_db, species_name, T):
    if species_name not in nasa_db:
        raise KeyError(f'NASA species not found: {species_name}')
    data = nasa_db[species_name]
    if T < data['T_low'] or T > data['T_high']:
        raise ValueError(f"Temperature {T:.2f} K is outside the NASA range for {species_name}: "
                         f"{data['T_low']:.2f}–{data['T_high']:.2f} K.")
    if T <= data['T_mid']:
        return data['low'], 'low'
    return data['high'], 'high'


def nasa7_properties(nasa_db, species_name, T):
    """(cp [J/mol/K], h [J/mol], s [J/mol/K], g = h - T s [J/mol]) from the NASA-7 polynomials."""
    T = float(T)
    if T <= 0.0:
        raise ValueError('NASA temperature must be positive.')
    a, _ = nasa7_coefficients(nasa_db, species_name, T)
    a1, a2, a3, a4, a5, a6, a7 = a
    cp_R = a1 + a2 * T + a3 * T ** 2 + a4 * T ** 3 + a5 * T ** 4
    h_RT = a1 + a2 * T / 2.0 + a3 * T ** 2 / 3.0 + a4 * T ** 3 / 4.0 + a5 * T ** 4 / 5.0 + a6 / T
    s_R = a1 * np.log(T) + a2 * T + a3 * T ** 2 / 2.0 + a4 * T ** 3 / 3.0 + a5 * T ** 4 / 4.0 + a7
    cp = R * cp_R
    h = R * T * h_RT
    entropy = R * s_R
    g = h - T * entropy
    return float(cp), float(h), float(entropy), float(g)
