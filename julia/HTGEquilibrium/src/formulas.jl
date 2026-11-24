# Chemical formulas: parsing, charge suffixes and the short summary shown in the
# app (Python: htg/formulas.py). Errors are `FormulaError`s with the same plain
# messages as the Python ValueErrors.
#
#     parse_formula("CH3OH")       -> (["C"=>1, "H"=>4, "O"=>1], nothing)
#     parse_formula("CO3^2-")      -> (["C"=>1, "O"=>3], -2)
#     parse_formula("CaSO4*2H2O")  -> (["Ca"=>1, "S"=>1, "O"=>6, "H"=>4], nothing)
#     formula_summary("CH3OH")     -> "C 1 · H 4 · O 1 · 32.04 g/mol"
#
# Element counts are Float64 in an OrderedDict (order of appearance, as Python's dict).

struct FormulaError <: Exception
    msg::String
end
Base.showerror(io::IO, e::FormulaError) = print(io, e.msg)

const Elements = OrderedDict{String,Float64}

const _STATE_SUFFIX = r"\((aq|g|s|l|cr)\)$"i
const _HYDRATE_SPLIT = r"[*·•]"
const _TOKEN = r"([A-Z][a-z]?|\(|\)|\[|\]|\d+(?:\.\d+)?)"

_g(x::Real) = py_format_g(x, 6)   # Python f"{x:g}" (Julia's @sprintf("%g") differs, see py_format_g)

"""
    split_charge(text, charge=nothing) -> (neutral part, charge or nothing)

For formulas such as "Na+", "CO3--", "Mg+2" or "SO4^2-". With `charge` given
(for example from the Charge field of an ion), the suffix matching that charge is
removed, which keeps the 2 of "VO2+" in the formula. Without it, a single trailing
sign after a digit counts as charge 1.
"""
function split_charge(text, charge=nothing)
    text = replace(replace(string(text), r"\s+" => ""), '⁺' => '+', '⁻' => '-')
    if charge !== nothing && Int(charge) != 0
        z = Int(charge)
        sign = z > 0 ? "+" : "-"
        n = abs(z)
        for suffix in ("^$n$sign", "^$sign$n", "^$sign", "$sign$n", "$n$sign", sign^n)
            suffix == "$n$sign" && n == 1 && continue
            if endswith(text, suffix) && length(text) > length(suffix)
                return first(text, length(text) - length(suffix)), z
            end
        end
        return text, nothing
    end
    m = match(r"\^(\d*)([+-])$", text)
    m === nothing && (m = match(r"\^([+-])(\d*)$", text))
    if m !== nothing
        groups = [g for g in m.captures if g !== nothing && !isempty(g)]
        digits = something(findfirst(g -> all(isdigit, g), groups), 0)
        count = digits == 0 ? 1 : parse(Int, groups[digits])
        sign = groups[findfirst(g -> g in ("+", "-"), groups)]
        return text[1:prevind(text, m.offset)], count * (sign == "+" ? 1 : -1)
    end
    m = match(r"([+-])(\d+)$", text)
    if m !== nothing
        return text[1:prevind(text, m.offset)], parse(Int, m[2]) * (m[1] == "+" ? 1 : -1)
    end
    m = match(r"(\++|-+)$", text)
    if m !== nothing && length(text) > length(m[1])
        signs = m[1]
        return text[1:prevind(text, m.offset)], length(signs) * (signs[1] == '+' ? 1 : -1)
    end
    return text, nothing
end

# Elements of a formula without charge or hydrate dots.
function _parse_part(text::AbstractString)
    tokens = [m.match for m in eachmatch(_TOKEN, text)]
    if join(tokens) != text
        bad = replace(text, _TOKEN => "")
        throw(FormulaError("Formula not understood: unexpected \"$(first(bad, 5))\"."))
    end
    stack = [Elements()]
    i = 1
    while i <= length(tokens)
        token = tokens[i]
        if token in ("(", "[")
            push!(stack, Elements())
        elseif token in (")", "]")
            length(stack) == 1 && throw(FormulaError("Formula not understood: unmatched bracket."))
            group = pop!(stack)
            count = 1.0
            if i < length(tokens) && isdigit(tokens[i+1][1])
                count = parse(Float64, tokens[i+1])
                i += 1
            end
            for (symbol, n) in group
                stack[end][symbol] = get(stack[end], symbol, 0.0) + n * count
            end
        elseif isdigit(token[1])
            throw(FormulaError("Formula not understood: a number must follow an element or bracket."))
        else
            count = 1.0
            if i < length(tokens) && isdigit(tokens[i+1][1])
                count = parse(Float64, tokens[i+1])
                i += 1
            end
            stack[end][token] = get(stack[end], token, 0.0) + count
        end
        i += 1
    end
    length(stack) != 1 && throw(FormulaError("Formula not understood: unmatched bracket."))
    return stack[1]
end

"""
    parse_formula(text, charge=nothing) -> (elements, charge or nothing)

Elements in order of appearance. Accepts brackets, hydrate dots ("CaSO4*2H2O"),
state suffixes ("(aq)") and charge suffixes. Throws `FormulaError`.
"""
function parse_formula(text, charge=nothing)
    raw = strip(text === nothing ? "" : string(text))
    isempty(raw) && throw(FormulaError("Enter a formula."))
    raw = replace(raw, _STATE_SUFFIX => "")
    neutral, found = split_charge(raw, charge)
    isempty(neutral) && throw(FormulaError("Formula not understood."))
    elements = Elements()
    for part in split(neutral, _HYDRATE_SPLIT)
        isempty(part) && continue
        m = match(r"^(\d+(?:\.\d+)?)(?=[A-Z(\[])", part)
        factor = m === nothing ? 1.0 : parse(Float64, m[1])
        body = m === nothing ? part : part[nextind(part, 0, length(m.match) + 1):end]
        for (symbol, n) in _parse_part(body)
            elements[symbol] = get(elements, symbol, 0.0) + n * factor
        end
    end
    isempty(elements) && throw(FormulaError("Formula not understood."))
    unknown = [s for s in keys(elements) if !haskey(ATOMIC_WEIGHTS, s)]
    isempty(unknown) || throw(FormulaError("Unknown element: $(join(unknown, ", "))."))
    return elements, found
end

"True for formulas of element symbols and integer counts only (\"C2H6O\", not \"Ca(OH)2\")."
is_plain_formula(text) = occursin(r"^(?:[A-Z][a-z]?\d*)+$", strip(string(text)))

molar_mass(elements) = pysum(ATOMIC_WEIGHTS[e] * n for (e, n) in elements)

elements_text(elements) = join(("$e $(_g(n))" for (e, n) in elements), " · ")

"\"C 1 · H 4 · O 1 · 32.04 g/mol\"; for ions the charge is added (\"· charge −2\")."
function formula_summary(text, charge=nothing)
    elements, found = parse_formula(text, charge)
    parts = [elements_text(elements), @sprintf("%.2f g/mol", molar_mass(elements))]
    z = (charge !== nothing && charge != 0) ? charge : found
    if z !== nothing && z != 0
        push!(parts, "charge $(z > 0 ? "+" : "−")$(abs(Int(z)))")
    end
    return join(parts, " · ")
end

"Formula in Hill order (C, H, then alphabetical): O=>1, C=>1, H=>4 -> \"CH4O\"."
function hill_formula(elements)
    part(symbol) = (n = elements[symbol]; symbol * (n == 1 ? "" : _g(n)))
    ks = collect(keys(elements))
    order = if haskey(elements, "C")
        vcat(["C"], haskey(elements, "H") ? ["H"] : String[], sort([k for k in ks if !(k in ("C", "H"))]))
    else
        sort(ks)
    end
    return join(part(k) for k in order)
end

"Element symbols from \"C, H, O\" (also spaces or semicolons); throws `FormulaError` for unknown ones."
function parse_element_list(text)
    symbols = [strip(s) for s in split(text === nothing ? "" : string(text), r"[,;\s]+") if !isempty(strip(s))]
    symbols = [uppercase(first(s, 1)) * lowercase(s[nextind(s, 1):end]) for s in symbols]
    unknown = [s for s in symbols if !haskey(ATOMIC_WEIGHTS, s)]
    isempty(unknown) || throw(FormulaError("Unknown element: $(join(unknown, ", "))."))
    return unique(symbols)
end
