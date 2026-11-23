"""Chemical formulas: parsing, charge suffixes and the short summary shown in the app.

    parse_formula('CH3OH')            -> ({'C': 1, 'H': 4, 'O': 1}, None)
    parse_formula('CO3^2-')           -> ({'C': 1, 'O': 3}, -2)
    parse_formula('CaSO4*2H2O')       -> ({'Ca': 1, 'S': 1, 'O': 6, 'H': 4}, None)
    formula_summary('CH3OH')          -> 'C 1 · H 4 · O 1 · 32.04 g/mol'
"""
import re

from .constants import ATOMIC_WEIGHTS

_STATE_SUFFIX = re.compile(r'\((aq|g|s|l|cr)\)$', re.IGNORECASE)
_HYDRATE_SPLIT = re.compile(r'[*·•]')
_TOKEN = re.compile(r'([A-Z][a-z]?|\(|\)|\[|\]|\d+(?:\.\d+)?)')


def split_charge(text, charge=None):
    """(neutral part, charge or None) of a formula such as 'Na+', 'CO3--', 'Mg+2' or 'SO4^2-'.

    With ``charge`` given (for example from the Charge field of an ion), the
    suffix matching that charge is removed, which keeps the 2 of 'VO2+' in the
    formula. Without it, a single trailing sign after a digit counts as charge 1.
    """
    text = re.sub(r'\s+', '', str(text)).replace('⁺', '+').replace('⁻', '-')
    if charge is not None and int(charge) != 0:
        z = int(charge)
        sign = '+' if z > 0 else '-'
        n = abs(z)
        for suffix in (f'^{n}{sign}', f'^{sign}{n}', f'^{sign}', f'{sign}{n}', f'{n}{sign}', sign * n):
            if suffix == f'{n}{sign}' and n == 1:
                continue
            if text.endswith(suffix) and len(text) > len(suffix):
                return text[:-len(suffix)], z
        return text, None
    match = re.search(r'\^(\d*)([+-])$', text) or re.search(r'\^([+-])(\d*)$', text)
    if match:
        digits = next((g for g in match.groups() if g and g.isdigit()), '1')
        sign = next(g for g in match.groups() if g in '+-')
        return text[:match.start()], int(digits) * (1 if sign == '+' else -1)
    match = re.search(r'([+-])(\d+)$', text)
    if match:
        return text[:match.start()], int(match.group(2)) * (1 if match.group(1) == '+' else -1)
    match = re.search(r'(\++|-+)$', text)
    if match and len(text) > len(match.group(1)):
        signs = match.group(1)
        return text[:match.start()], len(signs) * (1 if signs[0] == '+' else -1)
    return text, None


def _parse_part(text):
    """Elements of a formula without charge or hydrate dots; raises ValueError."""
    tokens = _TOKEN.findall(text)
    if ''.join(tokens) != text:
        bad = re.sub(r'([A-Z][a-z]?|\(|\)|\[|\]|\d+(?:\.\d+)?)', '', text)
        raise ValueError(f'Formula not understood: unexpected "{bad[:5]}".')
    stack = [{}]
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if token in '([':
            stack.append({})
        elif token in ')]':
            if len(stack) == 1:
                raise ValueError('Formula not understood: unmatched bracket.')
            group = stack.pop()
            count = 1.0
            if i + 1 < len(tokens) and tokens[i + 1][0].isdigit():
                count = float(tokens[i + 1])
                i += 1
            for symbol, n in group.items():
                stack[-1][symbol] = stack[-1].get(symbol, 0.0) + n * count
        elif token[0].isdigit():
            raise ValueError('Formula not understood: a number must follow an element or bracket.')
        else:
            count = 1.0
            if i + 1 < len(tokens) and tokens[i + 1][0].isdigit():
                count = float(tokens[i + 1])
                i += 1
            stack[-1][token] = stack[-1].get(token, 0.0) + count
        i += 1
    if len(stack) != 1:
        raise ValueError('Formula not understood: unmatched bracket.')
    return stack[0]


def parse_formula(text, charge=None):
    """(elements {symbol: count} in order of appearance, charge or None).

    Accepts brackets, hydrate dots ('CaSO4*2H2O'), state suffixes ('(aq)') and
    charge suffixes. Raises ValueError with a plain message.
    """
    raw = str(text or '').strip()
    if not raw:
        raise ValueError('Enter a formula.')
    raw = _STATE_SUFFIX.sub('', raw)
    neutral, found = split_charge(raw, charge)
    if not neutral:
        raise ValueError('Formula not understood.')
    elements = {}
    for part in _HYDRATE_SPLIT.split(neutral):
        if not part:
            continue
        match = re.match(r'^(\d+(?:\.\d+)?)(?=[A-Z(\[])', part)
        factor = float(match.group(1)) if match else 1.0
        body = part[match.end():] if match else part
        for symbol, n in _parse_part(body).items():
            elements[symbol] = elements.get(symbol, 0.0) + n * factor
    if not elements:
        raise ValueError('Formula not understood.')
    unknown = [symbol for symbol in elements if symbol not in ATOMIC_WEIGHTS]
    if unknown:
        raise ValueError(f'Unknown element: {", ".join(unknown)}.')
    elements = {k: (int(v) if float(v).is_integer() else v) for k, v in elements.items()}
    return elements, found


def is_plain_formula(text):
    """True for formulas of element symbols and integer counts only ('C2H6O', not 'Ca(OH)2')."""
    return bool(re.fullmatch(r'(?:[A-Z][a-z]?\d*)+', str(text).strip()))


def molar_mass(elements):
    return float(sum(ATOMIC_WEIGHTS[e] * n for e, n in elements.items()))


def elements_text(elements):
    return ' · '.join(f'{e} {n:g}' for e, n in elements.items())


def formula_summary(text, charge=None):
    """'C 1 · H 4 · O 1 · 32.04 g/mol'; for ions the charge is added ('· charge −2')."""
    elements, found = parse_formula(text, charge)
    parts = [elements_text(elements), f'{molar_mass(elements):.2f} g/mol']
    z = charge if charge else found
    if z:
        parts.append(f'charge {"+" if z > 0 else "−"}{abs(int(z))}')
    return ' · '.join(parts)


def hill_formula(elements):
    """Formula in Hill order (C, H, then alphabetical): {'O': 1, 'C': 1, 'H': 4} -> 'CH4O'."""
    def part(symbol):
        n = elements[symbol]
        return symbol + ('' if n == 1 else f'{n:g}')
    keys = list(elements)
    if 'C' in elements:
        order = ['C'] + (['H'] if 'H' in elements else []) + sorted(k for k in keys if k not in ('C', 'H'))
    else:
        order = sorted(keys)
    return ''.join(part(k) for k in order)


def parse_element_list(text):
    """Element symbols from 'C, H, O' (also spaces or semicolons); raises ValueError for unknown ones."""
    symbols = [s.strip() for s in re.split(r'[,;\s]+', str(text or '')) if s.strip()]
    symbols = [s[:1].upper() + s[1:].lower() for s in symbols]
    unknown = [s for s in symbols if s not in ATOMIC_WEIGHTS]
    if unknown:
        raise ValueError(f'Unknown element: {", ".join(unknown)}.')
    return list(dict.fromkeys(symbols))
