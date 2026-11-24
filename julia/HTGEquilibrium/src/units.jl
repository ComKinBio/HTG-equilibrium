# Display units for temperature, pressure and carbon bounds (Python: htg/units.py).
#
# Case files always store temperature in °C, pressure in bar (absolute) and
# carbon bounds as fractions 0-1 (the *file units*). The interface shows and
# accepts other units; every conversion goes through this module.
#
#     temperature: "C" (default), "K"        ("°C" is another spelling of "C")
#     pressure:    "bar" (default), "MPa", "Pa", "atm"
#
#     Units.to_file(573.15, "K")        -> 300.0   (display unit -> °C or bar)
#     Units.from_file(240, "MPa")       -> 24.0    (°C or bar -> display unit)
#     Units.format_number(2.4e7, "Pa")  -> "2.4e7"
#     Units.display(240, "Pa")          -> "2.4e7" (from_file, then format_number)
#     Units.field_label("Temperatures", "K") -> "Temperatures [K]"
#     Units.percent_to_fraction(90)     -> 0.9
#
# A submodule, as in Python (`units.display`), because several names (display,
# label, quantity) are too generic for the package namespace. The conversions
# work on numbers and, element by element, on arrays. Errors are `UnitError`s
# (a `ValueError`) with the Python messages.
module Units

using Printf
import ..ValueError, ..py_format_g

"""An unknown unit, or a unit of the wrong quantity (Python: ValueError)."""
struct UnitError <: ValueError
    msg::String
end

const TEMPERATURE_UNITS = ("C", "K")
const PRESSURE_UNITS = ("bar", "MPa", "Pa", "atm")
const UNITS = (TEMPERATURE_UNITS..., PRESSURE_UNITS...)
const FILE_UNITS = Dict("temperature" => "C", "pressure" => "bar")
const LABELS = Dict("C" => "°C", "K" => "K", "bar" => "bar", "MPa" => "MPa", "Pa" => "Pa", "atm" => "atm")
const _LABEL_ORDER = ("C", "K", "bar", "MPa", "Pa", "atm")   # Python dict order of LABELS

const ZERO_CELSIUS_K = 273.15
const BAR_PER_ATM = 1.01325
const BAR_PER_MPA = 10.0
const PA_PER_BAR = 1e5
const _SCIENTIFIC = ("Pa",)
const SIGNIFICANT_FIGURES = 6

"""The canonical key of `unit` ("°C" -> "C"); UnitError for an unknown unit."""
function unit_key(unit)
    key = unit == "°C" ? "C" : unit
    if !(key in UNITS)
        throw(UnitError("Unknown unit \"$unit\". Use one of: " *
                            join((LABELS[u] for u in _LABEL_ORDER), ", ") * "."))
    end
    return String(key)
end

"""\"temperature\" or \"pressure\"."""
quantity(unit) = unit_key(unit) in TEMPERATURE_UNITS ? "temperature" : "pressure"

"""The unit as shown: "C" -> "°C", "MPa" -> "MPa"."""
label(unit) = LABELS[unit_key(unit)]

"""\"Temperatures\", \"K\" -> \"Temperatures [K]\"."""
field_label(name, unit) = "$name [$(label(unit))]"

# ---------------------------------------------------------------- conversions
function temperature_to_C(value, unit)
    unit = unit_key(unit)
    unit == "C" && return value
    unit == "K" && return value .- ZERO_CELSIUS_K
    throw(UnitError("\"$unit\" is not a temperature unit."))
end

function temperature_from_C(value, unit)
    unit = unit_key(unit)
    unit == "C" && return value
    unit == "K" && return value .+ ZERO_CELSIUS_K
    throw(UnitError("\"$unit\" is not a temperature unit."))
end

function pressure_to_bar(value, unit)
    unit = unit_key(unit)
    unit == "bar" && return value
    unit == "MPa" && return value .* BAR_PER_MPA
    unit == "Pa" && return value ./ PA_PER_BAR
    unit == "atm" && return value .* BAR_PER_ATM
    throw(UnitError("\"$unit\" is not a pressure unit."))
end

function pressure_from_bar(value, unit)
    unit = unit_key(unit)
    unit == "bar" && return value
    unit == "MPa" && return value ./ BAR_PER_MPA
    unit == "Pa" && return value .* PA_PER_BAR
    unit == "atm" && return value ./ BAR_PER_ATM
    throw(UnitError("\"$unit\" is not a pressure unit."))
end

"""A value in `unit` -> °C (temperature units) or bar (pressure units)."""
to_file(value, unit) = quantity(unit) == "temperature" ? temperature_to_C(value, unit) : pressure_to_bar(value, unit)

"""A value in °C or bar -> `unit` (°C for temperature units, bar for pressure units)."""
from_file(value, unit) = quantity(unit) == "temperature" ? temperature_from_C(value, unit) : pressure_from_bar(value, unit)

"""Carbon bounds: 90 (%) -> 0.9."""
percent_to_fraction(percent) = percent ./ 100.0

"""Carbon bounds: 0.9 -> 90 (%)."""
fraction_to_percent(fraction) = fraction .* 100.0

# ---------------------------------------------------------------- display
_strip(mantissa) = occursin('.', mantissa) ? rstrip(rstrip(mantissa, '0'), '.') : mantissa

_fmt(spec, value) = Printf.format(Printf.Format(spec), value)


"""
    format_number(value, unit=nothing, digits=SIGNIFICANT_FIGURES) -> String

A number already in `unit` as text: up to `digits` significant figures, no
trailing zeros. Pa (and very large or small numbers) in scientific notation
without '+' or leading zeros: 2.4e7, 1.5e-5.

300.0 -> "300", 0.1 + 0.2 -> "0.3", 1234567 -> "1.23457e6".
"""
function format_number(value, unit=nothing, digits=SIGNIFICANT_FIGURES)
    value = value isa AbstractString ? parse(Float64, value) : Float64(value)
    if !isfinite(value)
        return isnan(value) ? "nan" : value > 0 ? "inf" : "-inf"   # Python str(float)
    end
    value == 0 && return "0"
    scientific = unit !== nothing && unit_key(unit) in _SCIENTIFIC
    text = scientific ? _fmt("%.$(digits - 1)e", value) : py_format_g(value, digits)
    occursin('e', text) || return String(_strip(text))
    mantissa, exponent = split(text, 'e')
    return "$(_strip(mantissa))e$(parse(Int, exponent))"
end

"""[300, 350] -> "300, 350" (numbers already in `unit`)."""
format_list(values, unit=nothing, digits=SIGNIFICANT_FIGURES) =
    join((format_number(v, unit, digits) for v in values), ", ")

"""A value in file units (°C or bar) as text in `unit`: display(240, "MPa") -> "24"."""
display(value, unit, digits=SIGNIFICANT_FIGURES) = format_number(from_file(value, unit), unit, digits)

end # module Units
