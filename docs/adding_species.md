# Adding species

Most species are added by adding rows to the databases in `data/`; no code
changes are needed. When the model loads, it scans the databases and makes a
species available when its data are complete.

Check the result with:

```sh
uv run htg species --search <name>          # ready species
uv run htg species --all --search <name>    # also incomplete ones, with the reason
```

The `Reason` column tells you what is missing.

## Fluid molecules (gas, organic liquid, aqueous neutral)

A molecule needs two things:

1. **Equation-of-state parameters**: a row in `PROJECT_PCSAFT.csv`
   (`record_type = pure_component`) or rows in `PROJECT_EPCSAFT.csv`
   (`record_type = species_parameter`, `Charge = 0`). Required: `Species`
   (the name you will use in case files), `Formula`, segment number `m`,
   segment diameter `sigma` [Å] and energy `epsilon/k` [K]. Associating
   molecules also need the association scheme (e.g. `2B`), `epsilonAB/k` and
   `kappaAB`. Keep the literature reference in `RefID`/`Authors`/`DOI`. If both
   files contain the molecule, the ePC-SAFT set is used.

2. **Standard-state thermochemistry**: a NASA-7 record in
   `Thermdat_WasteMech.txt` with the same elemental composition. The record is
   four lines in the standard NASA format (name, elements, phase letter in
   column 45, temperature range; then 14 coefficients).

The molecule's name is matched to a NASA record automatically when exactly one
record has its elemental formula and no other SAFT molecule shares that
formula. **Isomers** (e.g. acetone and propanal, both C3H6O) cannot be matched
this way. Name the record explicitly in `data/species.toml`:

```toml
[nasa_names]
acetone = "CH3COCH3"
```

The left side is the `Species` name from the SAFT file (case-insensitive), the
right side the NASA record name.

### Binary interaction parameters

PC-SAFT k_ij values for the gas and organic-liquid phases are listed at the end
of `data/species.toml`, by equation-of-state formula:

```toml
[[binary_interaction]]
components = ["CO2", "CH4"]
kij = 0.065
```

## Aqueous ions

An ion needs:

1. Rows in `PROJECT_EPCSAFT.csv` with `record_type = species_parameter` and a
   nonzero `Charge`: segment number, diameter, energy (and the ion's
   dielectric constant, if the parameter set has one).
2. An aqueous HKF record in `supcrt92.xml` with the same charge. The ion is
   matched by name or formula (e.g. `Mg2+`, `Mg++`).

Ions can only be used in the aqueous phase with `model = "electrolyte"`.

## Pure solids

Either of these works:

- a NASA-7 record with phase letter `S` in `Thermdat_WasteMech.txt`; the record
  name becomes the species name, or
- a mineral in `supcrt92.xml` (Maier-Kelley heat capacity, single temperature
  range). Minerals with phase transitions are listed but not yet supported.

To give a NASA solid a friendlier name (like `carbon` for record `CHAR`), add
it to `data/species.toml`:

```toml
[solid.carbon]
formula = "C"
nasa = "CHAR"
elements = { C = 1 }
```

## Feed species

A dry feed can be given as moles of species (`species_moles` in a case file),
for example `{ cellulose = 1.0 }`. Any species with a NASA record can be
used this way, including ones without SAFT parameters.

## After adding species

1. `uv run htg species --search <name>` should list the species as ready.
2. Use it in a case file and run `uv run htg check` on the case.
3. Run `uv run pytest`. The reference tests use a frozen copy of the original
   databases (`tests/reference/databases/`), so new species do not affect them.
   `tests/test_shipped_data.py` checks that `data/` and the example cases still
   load.

Do not change existing rows that the example cases use (for example the water
parameters) unless you intend to change the model; that changes all results.
