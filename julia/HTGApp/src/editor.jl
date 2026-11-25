# The Case tab's form, on the server and without the page (Python: gui/editor.py
# CaseEditor). The page (case_tab.jl) shows `form_view(editor)` and the field texts,
# and sends every change here; everything the form computes is done by the
# functions of case_form.jl and the core, on every change, so the whole form can
# be tested without a browser (test/test_case_tab.jl).
#
# The form writes case-file format 2 only. Text fields are parsed on every
# change. A field that cannot be parsed keeps its text, shows the message under
# it and leaves the case data unchanged; any problem disables Save, Save as and
# Add to queue and is named in the status line. Old files that cannot be
# converted (legacy_blocker) are not shown in the form: they keep the old format
# and can only be checked, queued or saved as they are.

"""Order of the form's sections, for "the first problem" in the status line."""
const SECTIONS = ("name", "conditions", "feed", "products", "carbon", "aqueous", "solver", "advanced")

"""A field of the form: its label (for the status line), the label shown, how its text is read and stored."""
struct Field
    label::String
    shown::String
    parse::Union{Nothing,Function}   # text -> value, or a ValueError (the message goes under the field)
    store::Union{Nothing,Function}   # value -> writes the case data
    order::Int                       # creation order (Python: the element id)
    placeholder::String
end

"""
    CaseEditor(state; async_bands=true, on_change=nothing)

The form of the Case tab over `state.document`, one per page (client). Every
operation (`edit!`, `set!`, `open!`, ...) updates the case data, the field
messages and the derived numbers; `form_view(editor)` is what the page shows and
`editor.texts` the text of each field.

- `notifications`: (text, type) messages for the page to show and empty.
- `texts_changed`: the editor changed field texts itself (rebuilt part of the
  form, converted a value); the page sends `texts` again and clears it.
- `on_change()`: called after a background step changed the view (the SUPCRT92
  band check, the Check button).
"""
mutable struct CaseEditor
    state::State
    lock::ReentrantLock
    fields::Dict{String,Field}
    texts::Dict{String,String}
    parse_errors::Dict{String,String}
    feed_errors::Dict{String,String}
    feed_message::Union{Nothing,String}
    feed_summary::Any
    feed_stash::Dict{String,Any}
    feed_mode::String
    built::Bool
    order::Int
    share_units::Dict{String,String}
    water_unit::String
    species_items::Vector{Vector{Any}}   # [name, value, unit]
    minimum_items::Vector{Vector{Any}}   # [phase, species, value]
    weight_items::Vector{Vector{Any}}    # [element, value]
    carbon::Dict{String,Dict{String,Any}}
    element_list::Vector{String}
    new_element_error::Union{Nothing,String}
    new_weight_error::Union{Nothing,String}
    case_files::Vector{Tuple{String,String}}
    texts_changed::Bool
    notifications::Vector{Tuple{String,String}}
    check::Dict{String,Any}
    saveas::Dict{String,Any}
    band_points::Vector{Tuple{Float64,Float64}}
    band_hits::Vector{Tuple{Float64,Float64}}
    band_generation::Int
    band_check::Function
    band_warned::Bool
    async_bands::Bool
    on_change::Any
    elements_cache::Any
end

function CaseEditor(state::State; async_bands=true, on_change=nothing, band_check=supcrt92_band_check)
    ed = CaseEditor(state, ReentrantLock(), Dict(), Dict(), Dict(), Dict(), nothing, nothing, Dict(), "elements",
                    false, 0, Dict("solids" => "%", "ash" => "%"), "kg", [], [], [], Dict(), String[], nothing,
                    nothing, [], false, [], Dict{String,Any}("open" => false), Dict{String,Any}("open" => false),
                    [], [], 0, band_check, false, async_bands, on_change, nothing)
    lock(ed.lock) do
        ed.case_files = case_files(state)
        build!(ed)
        refresh!(ed)
        ed.texts_changed = false
    end
    return ed
end

current_doc(ed::CaseEditor) = ed.state.document
doc_data(ed::CaseEditor) = ed.state.document.data
t_unit(ed::CaseEditor) = ed.state.settings["temperature_unit"]
p_unit(ed::CaseEditor) = ed.state.settings["pressure_unit"]
notify!(ed::CaseEditor, text, type="info") = push!(ed.notifications, (text, type))

"""Case files below the cases folder: (absolute path, path relative to the folder), sorted by path."""
function case_files(state::State)
    root = state.cases_dir
    isdir(root) || return Tuple{String,String}[]
    files = String[]
    for (dir, _, names) in walkdir(root), name in names
        endswith(name, ".toml") && push!(files, joinpath(dir, name))
    end
    return [(f, relpath(f, root)) for f in sort!(files)]
end

"""`molar_mass(name)` of the species data, or `nothing` while they are loading."""
function molar_mass_fn(ed::CaseEditor)
    elements = species_elements(ed)
    return elements === nothing ? nothing : molar_mass_function(elements)
end

"""`species_elements(name)` of the loaded model (cached), or `nothing` while it is loading."""
function species_elements(ed::CaseEditor)
    model = ed.state.model
    model === nothing && return nothing
    cached = ed.elements_cache
    cached !== nothing && cached[1] === model && return cached[2]
    elements = registry_elements(model.registry)
    ed.elements_cache = (model, elements)
    return elements
end

# ---------------------------------------------------------------- files
"""Open the case file at `path` (a message if it cannot be read)."""
function open!(ed::CaseEditor, path)
    doc = try
        open_document(path; molar_mass=molar_mass_fn(ed))
    catch exc
        exc isa InterruptException && rethrow()
        notify!(ed, "Could not open $path: $(one_line(exc))", "negative")
        return false
    end
    ed.state.document = doc
    reload!(ed)
    return true
end

"""A new case from the template."""
function new!(ed::CaseEditor)
    ed.state.document = new_document()
    reload!(ed)
    return true
end

function reload!(ed::CaseEditor)
    ed.feed_stash = Dict()
    build!(ed)
    refresh!(ed)
    mark_saved!(current_doc(ed))   # what building the form tidied up is not an unsaved change
    ed.texts_changed = true
    return ed
end

"""Save to the file it came from; false when it has none yet (then the page asks for a path)."""
function save!(ed::CaseEditor)
    current_doc(ed).path === nothing && return false
    write!(ed, current_doc(ed).path)
    return true
end

"""The path Save as suggests: the case name in the cases folder."""
function default_save_path(ed::CaseEditor)
    name = string(something(get(doc_data(ed), "name", nothing), ""))
    return joinpath(ed.state.cases_dir, (isempty(name) ? "case" : name) * ".toml")
end

"""Write the case to `path`; true when it was saved."""
function write!(ed::CaseEditor, path)
    problem = first_problem(ed)
    if problem !== nothing
        notify!(ed, "Not saved: $problem", "negative")
        return false
    end
    saved = try
        save_document(current_doc(ed), path)
    catch exc
        exc isa InterruptException && rethrow()
        notify!(ed, "Not saved: $(one_line(exc))", "negative")
        return false
    end
    notify!(ed, "Saved $saved", "positive")
    ed.case_files = case_files(ed.state)
    refresh!(ed)
    return true
end

# ---------------------------------------------------------------- status and actions
"""Set `data[path...] = value` (`nothing` removes it) and update the form."""
function set!(ed::CaseEditor, path, value)
    putpath!(doc_data(ed), path, value)
    refresh!(ed)
end

"""Problems with the form, first one first (field problems in form order, then the feed)."""
function problems(ed::CaseEditor)
    out = String[]
    ed.built || return out
    rank(key) = (something(findfirst(==(first(split(key, '.'))), SECTIONS), length(SECTIONS) + 1), ed.fields[key].order)
    for key in sort!(collect(keys(ed.fields)); by=rank)
        message = get(ed.parse_errors, key, get(ed.feed_errors, key, nothing))
        message === nothing || push!(out, problem_text(ed.fields[key].label, message))
    end
    if ed.feed_message !== nothing && !(ed.feed_message in values(ed.feed_errors))
        push!(out, ed.feed_message)
    end
    return out
end

"""The first problem, or `nothing` when the case can be saved and queued."""
function first_problem(ed::CaseEditor)
    list = problems(ed)
    isempty(list) || return list[1]
    try
        document_case(current_doc(ed))
    catch exc
        exc isa InterruptException && rethrow()
        return one_line(exc)
    end
    return nothing
end

"""(status line text, ok)."""
function status_line(ed::CaseEditor)
    list = problems(ed)
    if isempty(list)
        try
            c = document_case(current_doc(ed))
            return ready_text(length(grid(c)), length(chains(c))), true
        catch exc
            exc isa InterruptException && rethrow()
            list = [one_line(exc)]
        end
    end
    return "✗ $(list[1])", false
end

"""
Queue the case (Python: CaseEditor._queue); returns the run, or `nothing` with a
message (also while the species data are loading).
"""
function queue!(ed::CaseEditor)
    problem = first_problem(ed)
    if problem !== nothing
        notify!(ed, "Not queued: $problem", "negative")
        return nothing
    end
    c = document_case(current_doc(ed))
    run = try
        submit(ed.state, c)
    catch exc
        exc isa InterruptException && rethrow()
        notify!(ed, "Not queued: $(one_line(exc))", "negative")
        return nothing
    end
    notify!(ed, "Queued $(c.name): $(point_count_text(run.total))", "positive")
    return run
end

"""
The case for the Check button, or `nothing` with a message (a problem in the
form, or the species data still loading).
"""
function check_case(ed::CaseEditor)
    problem = first_problem(ed)
    if problem !== nothing
        notify!(ed, problem, "negative")
        return nothing
    end
    if ed.state.model === nothing
        notify!(ed, LOADING_MESSAGE, "warning")
        return nothing
    end
    return document_case(current_doc(ed))
end

"""Show the result of the checks (`messages` from `preflight_messages`)."""
function check_done!(ed::CaseEditor, c, messages)
    ed.check = Dict{String,Any}("open" => true, "ok" => isempty(messages),
                                "title" => isempty(messages) ? "All $(point_count_text(length(grid(c)))) passed the checks." :
                                           "Problems found",
                                "messages" => messages[1:min(end, 50)],
                                "more" => length(messages) > 50 ? "… and $(length(messages) - 50) more." : "")
    return ed
end

# ---------------------------------------------------------------- fields
"""A text field parsed on every change: `parse(text)` -> value or a ValueError (shown under the field),
then `store(value)` writes the case data."""
function field!(ed::CaseEditor, key, label, text, parse, store; shown=label, placeholder="")
    ed.order += 1
    ed.fields[key] = Field(label, shown, parse, store, ed.order, placeholder)
    ed.texts[key] = string(text)
    return key
end

"""A select (or other choice) whose problems are shown like a field's."""
function choice_field!(ed::CaseEditor, key, label)
    ed.order += 1
    ed.fields[key] = Field(label, label, nothing, nothing, ed.order, "")
    return key
end

"""The text of a field changed (typed by the user)."""
function edit!(ed::CaseEditor, key, text)
    text = text === nothing ? "" : string(text)
    ed.texts[key] = text
    if haskey(ed.fields, key) && ed.fields[key].parse !== nothing
        apply!(ed, key)
    elseif key == "notes"
        set!(ed, ["description"], text)
    elseif key == "measured.label"
        commit_measured!(ed)
    end
    return ed
end

"""The editor sets a field's text itself (a converted value): stored as if typed."""
function set_text!(ed::CaseEditor, key, text)
    ed.texts_changed = true
    edit!(ed, key, text)
end

"""Parse the field's text and store it (or record the problem), then update the form."""
function apply!(ed::CaseEditor, key)
    f = ed.fields[key]
    value = try
        f.parse(get(ed.texts, key, ""))
    catch exc
        exc isa HTGEquilibrium.ValueError || rethrow()
        ed.parse_errors[key] = exc.msg
        refresh!(ed)
        return ed
    end
    delete!(ed.parse_errors, key)
    f.store(value)
    refresh!(ed)
    return ed
end

"""Drop the fields (and their problems) of a part of the form that is rebuilt."""
function forget!(ed::CaseEditor, prefix)
    for key in [k for k in keys(ed.fields) if startswith(k, prefix)]
        delete!(ed.fields, key)
        delete!(ed.parse_errors, key)
        delete!(ed.texts, key)
    end
end

"""A number field that may be left empty (= the default, shown as placeholder)."""
function optional!(ed::CaseEditor, key, label, path, default; limits...)
    field!(ed, key, label, number_text(getpath(doc_data(ed), path)), text -> parse_number(text; limits...),
           value -> putpath!(doc_data(ed), path, value); placeholder=default)
end

# ---------------------------------------------------------------- the form
function build!(ed::CaseEditor)
    ed.fields = Dict()
    ed.texts = Dict()
    ed.parse_errors = Dict()
    ed.feed_errors = Dict()
    ed.feed_message = nothing
    ed.feed_summary = nothing
    ed.built = false
    ed.new_element_error = ed.new_weight_error = nothing
    ed.texts["new_element"] = ed.texts["new_weight"] = ""
    current_doc(ed).legacy_blocker === nothing || return ed
    d = doc_data(ed)
    field!(ed, "name.name", "Name (also the results folder name)", string(get(d, "name", "")), strip,
           value -> putpath!(doc_data(ed), ["name"], String(value)))
    ed.texts["notes"] = string(get(d, "description", ""))
    build_conditions!(ed)
    build_feed!(ed)
    build_carbon!(ed)
    build_aqueous!(ed)
    build_solver!(ed)
    build_advanced!(ed)
    ed.built = true
    return ed
end

# ---------------------------------------------------------------- temperature and pressure
function build_conditions!(ed::CaseEditor)
    forget!(ed, "conditions.")
    T, P = t_unit(ed), p_unit(ed)
    conditions = get!(doc_data(ed), "conditions") do
        OrderedDict{String,Any}()
    end
    if get(conditions, "points", nothing) !== nothing
        text = join(("$(Units.display(t, T)), $(Units.display(p, P))" for (t, p) in conditions["points"]), "\n")
        field!(ed, "conditions.points", "(T, P) points, one per line [$(Units.label(T)), $(Units.label(P))]", text,
               text -> parse_points_field(ed, text), value -> putpath!(doc_data(ed), ["conditions", "points"], value))
    else
        field!(ed, "conditions.temperatures", Units.field_label("Temperatures", T),
               format_numbers(get(conditions, "temperature_C", nothing), T), text -> parse_list(text, t_unit(ed)),
               value -> putpath!(doc_data(ed), ["conditions", "temperature_C"], scalar(value)))
        field!(ed, "conditions.pressures", Units.field_label("Pressures", P),
               format_numbers(get(conditions, "pressure_bar", nothing), P), text -> parse_list(text, p_unit(ed)),
               value -> putpath!(doc_data(ed), ["conditions", "pressure_bar"], scalar(value)))
    end
end

function parse_list(text, unit)
    values = parse_numbers(text; unit)
    isempty(values) && throw(FieldError("Give at least one value."))
    return values
end

function parse_points_field(ed::CaseEditor, text)
    points = parse_points(text, t_unit(ed), p_unit(ed))
    isempty(points) && throw(FieldError("Give at least one point."))
    return points
end

"""Temperature or pressure unit (a preference; `key` "temperature_unit" or "pressure_unit")."""
function set_unit!(ed::CaseEditor, key, unit)
    (unit === nothing || unit == ed.state.settings[key]) && return ed
    ed.state.settings[key] = unit
    build_conditions!(ed)
    if key == "temperature_unit" && haskey(ed.fields, "aqueous.cutoff")
        f = ed.fields["aqueous.cutoff"]
        ed.fields["aqueous.cutoff"] = Field(f.label, Units.field_label("Remove at and above", unit), f.parse, f.store,
                                            f.order, f.placeholder)
        haskey(ed.parse_errors, "aqueous.cutoff") ||
            (ed.texts["aqueous.cutoff"] = format_numbers(getpath(doc_data(ed), ("aqueous", "cutoff_C")), unit))
    end
    ed.texts_changed = true
    refresh!(ed)
end

"""The Input formats table open or closed (a preference)."""
function set_input_formats!(ed::CaseEditor, open)
    ed.state.settings["input_formats_open"] = Bool(open)
    return ed
end

"""\"grid\" (every T at every P) or \"points\" (a list of (T, P) points)."""
function switch_conditions!(ed::CaseEditor, mode)
    conditions = get!(doc_data(ed), "conditions") do
        OrderedDict{String,Any}()
    end
    if mode == "points" && !haskey(conditions, "points")
        temperatures = _values(get(conditions, "temperature_C", nothing))
        pressures = _values(get(conditions, "pressure_bar", nothing))
        isempty(temperatures) && (temperatures = Any[nothing])
        isempty(pressures) && (pressures = Any[nothing])
        points = Any[Any[t, p] for p in pressures for t in temperatures if t !== nothing && p !== nothing]
        empty!(conditions)
        conditions["points"] = isempty(points) ? Any[Any[300.0, 240.0]] : points
    elseif mode == "grid" && haskey(conditions, "points")
        points = pop!(conditions, "points")
        points = HTGEquilibrium.py_truth(points) ? points : Any[Any[300.0, 240.0]]
        conditions["temperature_C"] = sort!(unique(Any[p[1] for p in points]))
        conditions["pressure_bar"] = scalar(sort!(unique(Any[p[2] for p in points])))
    end
    build_conditions!(ed)
    ed.texts_changed = true
    refresh!(ed)
end

# ---------------------------------------------------------------- feed
feed_data(ed::CaseEditor) = get!(doc_data(ed), "feed") do
    OrderedDict{String,Any}()
end

function build_feed!(ed::CaseEditor)
    forget!(ed, "feed.")
    ed.feed_mode = feed_mode(feed_data(ed))
    ed.feed_mode == "species" ? species_feed!(ed) : elemental_feed!(ed)
end

"""\"elements\" (elemental analysis) or \"species\" (species amounts); the other mode's feed is kept to switch back."""
function switch_feed!(ed::CaseEditor, mode)
    feed = feed_data(ed)
    current = feed_mode(feed)
    (mode === nothing || mode == current) && return ed
    ed.feed_stash[current] = OrderedDict{String,Any}(k => deepcopy(v) for (k, v) in feed if k != "salts")
    new = deepcopy(get(ed.feed_stash, mode, nothing) === nothing ? default_feed(mode, ed.state.feed_species) :
                   ed.feed_stash[mode])
    HTGEquilibrium.py_truth(get(feed, "salts", nothing)) && (new["salts"] = feed["salts"])
    doc_data(ed)["feed"] = new
    build_feed!(ed)
    ed.texts_changed = true
    refresh!(ed)
end

function masses(ed::CaseEditor)
    try
        return elemental_masses(feed_data(ed))
    catch exc
        exc isa HTGEquilibrium.ValueError || rethrow()
        return nothing
    end
end

function elemental_feed!(ed::CaseEditor)
    feed = feed_data(ed)
    field!(ed, "feed.total", "1 · Total mass [kg]", number_text(get(feed, "total_kg", nothing)),
           text -> parse_number(text; required=true, above=0), value -> putpath!(feed_data(ed), ["total_kg"], value))
    share!(ed, "solids", "2 · Solids", true)
    share!(ed, "ash", "3 · Ash", false)
    build_elements!(ed)
end

"""Solids or ash: a value with a kg | % switch."""
function share!(ed::CaseEditor, name, label, required)
    entry = something(get(feed_data(ed), name, nothing), Dict())
    ed.share_units[name] = entry isa AbstractDict ? string(get(entry, "unit", "%")) : "%"
    limits = required ? (; above=0) : (; minimum=0)
    store(value) = putpath!(feed_data(ed), [name], value === nothing ? nothing :
                            OrderedDict{String,Any}("value" => value, "unit" => ed.share_units[name]))
    field!(ed, "feed.$name", label, number_text(entry isa AbstractDict ? get(entry, "value", nothing) : nothing),
           text -> parse_number(text; required, limits...), store)
end

"""The hint under the solids or ash field."""
function share_hint(ed::CaseEditor, name)
    unit = ed.share_units[name]
    name == "solids" && return unit == "%" ? "% of total mass" : "kg"
    return (unit == "%" ? "% of solids" : "kg") * ". Empty = no ash."
end

"""Switch solids or ash between "%" and "kg"; the value is converted when it can be."""
function share_unit!(ed::CaseEditor, name, unit)
    old = ed.share_units[name]
    unit == old && return ed
    ed.share_units[name] = unit
    m = masses(ed)
    whole = m === nothing ? nothing : name == "solids" ? m["total"] : m["solids"]
    entry = get(feed_data(ed), name, nothing)
    current = entry isa AbstractDict ? get(entry, "value", nothing) : nothing
    key = "feed.$name"
    (haskey(ed.parse_errors, key) || current === nothing) && return apply!(ed, key)
    converted = convert_share(current, old, unit, whole)
    converted == current && return apply!(ed, key)
    set_text!(ed, key, number_text(converted))   # stores in the new unit
end

function build_elements!(ed::CaseEditor)
    forget!(ed, "feed.element:")
    elements = something(get(feed_data(ed), "elements_wt_percent", nothing), OrderedDict{String,Any}())
    ed.element_list = String[string(e) for e in keys(elements)]
    for element in ed.element_list
        field!(ed, "feed.element:$element", "Element $element", number_text(elements[element]),
               text -> parse_number(text; required=true, minimum=0),
               value -> (get!(feed_data(ed), "elements_wt_percent") do
                             OrderedDict{String,Any}()
                         end[element] = value); shown="")
    end
end

"""Add an element to the elemental analysis (`text`: the symbol typed)."""
function add_element!(ed::CaseEditor, text)
    elements = get!(feed_data(ed), "elements_wt_percent") do
        OrderedDict{String,Any}()
    end
    symbol = element_symbol(text)
    ed.texts["new_element"] = string(something(text, ""))
    ed.new_element_error = element_error(symbol, keys(elements))
    ed.new_element_error === nothing || return refresh!(ed)
    elements[symbol] = 0.0
    build_elements!(ed)
    ed.texts_changed = true
    refresh!(ed)
end

function remove_element!(ed::CaseEditor, element)
    elements = get(feed_data(ed), "elements_wt_percent", nothing)
    elements isa AbstractDict && delete!(elements, element)
    build_elements!(ed)
    ed.texts_changed = true
    refresh!(ed)
end

# species amounts
function species_molar_mass(ed::CaseEditor, name)
    mm = molar_mass_fn(ed)
    mm === nothing && return nothing
    try
        return mm(name)
    catch exc
        exc isa Union{HTGEquilibrium.ValueError,ArgumentError,KeyError} || rethrow()
        return nothing
    end
end

function species_feed!(ed::CaseEditor)
    feed = feed_data(ed)
    species = get!(feed, "species") do
        OrderedDict{String,Any}()
    end
    if !haskey(species, "water")
        species = feed["species"] = merge(OrderedDict{String,Any}("water" => OrderedDict{String,Any}("value" => 0.0, "unit" => "kg")),
                                          species)
    end
    ed.species_items = Vector{Any}[Any[name, entry isa AbstractDict ? get(entry, "value", nothing) : nothing,
                                       entry isa AbstractDict ? string(get(entry, "unit", "kg")) : "kg"]
                                   for (name, entry) in species if name != "water"]
    water = species["water"]
    ed.water_unit = water isa AbstractDict ? string(get(water, "unit", "kg")) : "kg"
    store_water(value) = (feed_data(ed)["species"]["water"] = OrderedDict{String,Any}("value" => value, "unit" => ed.water_unit))
    field!(ed, "feed.water", "Water", number_text(water isa AbstractDict ? get(water, "value", nothing) : nothing),
           text -> parse_number(text; required=true, minimum=0), store_water; shown="Amount")
    build_species_rows!(ed)
end

"""Water in "kg" or "mol"; the amount is converted."""
function water_unit!(ed::CaseEditor, unit)
    old = ed.water_unit
    unit == old && return ed
    ed.water_unit = unit
    current = get(feed_data(ed)["species"]["water"], "value", nothing)
    converted = convert_species_amount(current, old, unit, water_molar_mass(get(doc_data(ed), "atomic_weights", nothing)))
    (haskey(ed.parse_errors, "feed.water") || converted == current) && return apply!(ed, "feed.water")
    set_text!(ed, "feed.water", number_text(converted))
end

function build_species_rows!(ed::CaseEditor)
    forget!(ed, "feed.species")
    for (i, (name, value, unit)) in enumerate(ed.species_items)
        k = i - 1   # 0-based in the keys, as Python
        choice_field!(ed, "feed.species-name:$k", "Species")
        field!(ed, "feed.species:$k", "Species $(isempty(name) ? i : name)", number_text(value),
               text -> parse_number(text; required=true, minimum=0), v -> species_value!(ed, i, v); shown="Amount")
    end
    commit_species!(ed)
end

"""Row `k` (0-based) of the species table: its species."""
function species_name!(ed::CaseEditor, k, name)
    ed.species_items[k+1][1] = something(name, "")
    commit_species!(ed)
    refresh!(ed)
end

species_value!(ed::CaseEditor, i, value) = (ed.species_items[i][2] = value; commit_species!(ed))

"""Row `k` (0-based): "kg" or "mol"; the amount is converted with the species' molar mass."""
function species_unit!(ed::CaseEditor, k, unit)
    name, value, old = ed.species_items[k+1]
    unit == old && return ed
    ed.species_items[k+1][3] = unit
    converted = convert_species_amount(value, old, unit, isempty(name) ? nothing : species_molar_mass(ed, name))
    key = "feed.species:$k"
    (haskey(ed.parse_errors, key) || converted == value) && return apply!(ed, key)
    set_text!(ed, key, number_text(converted))
end

function add_species!(ed::CaseEditor)
    push!(ed.species_items, Any["", 1.0, "kg"])
    build_species_rows!(ed)
    ed.texts_changed = true
    refresh!(ed)
end

function remove_species!(ed::CaseEditor, k)
    deleteat!(ed.species_items, k + 1)
    build_species_rows!(ed)
    ed.texts_changed = true
    refresh!(ed)
end

"""Write the species rows to the feed; unnamed or repeated species are problems of their row."""
function commit_species!(ed::CaseEditor)
    species = get!(feed_data(ed), "species") do
        OrderedDict{String,Any}()
    end
    out = OrderedDict{String,Any}("water" => get(species, "water", OrderedDict{String,Any}("value" => 0.0, "unit" => "kg")))
    for (i, (name, value, unit)) in enumerate(ed.species_items)
        key = "feed.species-name:$(i - 1)"
        if isempty(name)
            ed.parse_errors[key] = "Choose a species or remove the row."
        elseif haskey(out, name)
            ed.parse_errors[key] = "$name is already in the feed."
        else
            delete!(ed.parse_errors, key)
            value === nothing || (out[name] = OrderedDict{String,Any}("value" => value, "unit" => unit))
        end
    end
    feed_data(ed)["species"] = out
end

"""Check the feed: the summary, or the first problem (under its field when it has one)."""
function update_feed!(ed::CaseEditor)
    feed = feed_data(ed)
    mode = feed_mode(feed)
    elements = mode == "species" ? species_elements(ed) : nothing
    summary, field, message = if mode == "species" && elements === nothing
        nothing, "feed", LOADING_MESSAGE
    else
        try
            check_feed(feed, get(doc_data(ed), "atomic_weights", nothing), elements)
        catch exc
            exc isa InterruptException && rethrow()
            nothing, "feed", one_line(exc)
        end
    end
    ed.feed_summary = summary
    ed.feed_errors = Dict()
    ed.feed_message = message
    message === nothing && return ed
    key = get(Dict("total" => "feed.total", "solids" => "feed.solids", "ash" => "feed.ash",
                   "species:water" => "feed.water"), field, nothing)
    if startswith(field, "element:")
        key = "feed." * field
    elseif startswith(field, "species:") && key === nothing
        name = field[length("species:")+1:end]
        index = findfirst(item -> item[1] == name, ed.species_items)
        key = index === nothing ? nothing : "feed.species:$(index - 1)"
    end
    key !== nothing && haskey(ed.fields, key) && (ed.feed_errors[key] = message)
    return ed
end

# ---------------------------------------------------------------- products
"""The products allowed in a phase."""
set_products!(ed::CaseEditor, phase, names) = set!(ed, ["phases", phase], Any[string(n) for n in something(names, [])])

# ---------------------------------------------------------------- carbon split
function build_carbon!(ed::CaseEditor)
    forget!(ed, "carbon.")
    for phase in CARBON_PHASES
        choice, first_, second = carbon_choice(getpath(doc_data(ed), ("constraints", "carbon", phase)))
        row = ed.carbon[phase] = Dict{String,Any}("choice" => choice, "first" => first_, "second" => second,
                                                  "range" => choice == CARBON_RANGE)
        name = PHASE_NAMES[phase]
        percent_text(fraction) = fraction isa Real ? number_text(Units.fraction_to_percent(fraction)) : ""
        parser(slot) = function (text)
            (row["choice"] == "free" || (slot == "second" && row["choice"] != CARBON_RANGE)) && return nothing
            return parse_number(text; required=true, minimum=0, maximum=100)
        end
        storer(slot) = function (value)
            row[slot] = value === nothing ? nothing : tidy(Units.percent_to_fraction(value))
            commit_carbon!(ed, phase)
        end
        field!(ed, "carbon.$phase", "$name carbon [%]", percent_text(first_), parser("first"), storer("first");
               shown=choice == CARBON_RANGE ? "min %" : "%")
        field!(ed, "carbon.$phase.max", "$name carbon maximum [%]", percent_text(second), parser("second"), storer("second");
               shown="max %")
    end
end

function commit_carbon!(ed::CaseEditor, phase)
    row = ed.carbon[phase]
    kind = row["choice"]
    needed = kind == CARBON_RANGE ? [row["first"], row["second"]] : [row["first"]]
    bound = kind == "free" || any(isnothing, needed) ? nothing : carbon_bound(kind, row["first"], row["second"])
    putpath!(doc_data(ed), ["constraints", "carbon", phase], bound)
end

"""free, exactly, at least, at most (or between, when the file has a range) for a phase."""
function carbon_choice!(ed::CaseEditor, phase, choice)
    ed.carbon[phase]["choice"] = choice
    apply!(ed, "carbon.$phase.max")
    apply!(ed, "carbon.$phase")
end

# ---------------------------------------------------------------- aqueous phase
function build_aqueous!(ed::CaseEditor)
    store(value) = putpath!(doc_data(ed), ["aqueous", "cutoff_C"],
                            value === nothing ? nothing : tidy(Units.to_file(value, t_unit(ed))))
    field!(ed, "aqueous.cutoff", Units.field_label("Remove at and above", t_unit(ed)),
           format_numbers(getpath(doc_data(ed), ("aqueous", "cutoff_C")), t_unit(ed)), text -> parse_number(text), store)
end

"""Why the water model cannot be chosen now (`nothing` when it can)."""
function water_model_reason(d)
    getpath(d, ("aqueous", "enabled"), true) === false && return WATER_MODEL_NO_AQUEOUS
    getpath(d, ("aqueous", "model"), "electrolyte") == "molecular" && return WATER_MODEL_MOLECULAR
    return nothing
end

"""SUPCRT92 (the default: no key in the file) or IAPWS-95."""
set_water_model!(ed::CaseEditor, name) = set!(ed, ["aqueous", "water_model"], name == "SUPCRT92" ? nothing : name)

# ---------------------------------------------------------------- solver
function build_solver!(ed::CaseEditor)
    starts = getpath(doc_data(ed), ("solver", "starts"), 1)
    field!(ed, "solver.starts", "Initial guesses", number_text(starts),
           text -> parse_number(text; required=true, minimum=1, integer=true),
           value -> putpath!(doc_data(ed), ["solver", "starts"], value))
end

# ---------------------------------------------------------------- advanced settings
function build_advanced!(ed::CaseEditor)
    optional!(ed, "advanced.gas_water", "Water in the gas at the start [fraction]",
              ["initial_guess", "initial_gas_water_fraction"], "1e-6"; minimum=0, maximum=1)
    optional!(ed, "advanced.water_floor", "Smallest amount of aqueous water [mol / mol C]", ["aqueous", "water_floor"],
              "1e-8"; minimum=0)
    optional!(ed, "advanced.max_iterations", "Maximum iterations", ["solver", "max_iterations"], "5000";
              integer=true, minimum=1)
    optional!(ed, "advanced.seed", "Random seed for starting points", ["solver", "seed"], "0"; integer=true, minimum=0)
    optional!(ed, "advanced.energy_scale", "Energy scale [J / mol C]", ["solver", "energy_scale"], "automatic"; above=0)
    optional!(ed, "advanced.composition_tol", "Composition [mol / mol C]", ["solver", "distinct_composition_tol"],
              "1e-4"; minimum=0)
    optional!(ed, "advanced.gibbs_tol", "Gibbs energy [J / mol C]", ["solver", "distinct_gibbs_tol"], "1"; minimum=0)
    minimum = something(getpath(doc_data(ed), ("constraints", "minimum_amounts")), OrderedDict())
    ed.minimum_items = Vector{Any}[Any[string(p), string(s), v] for (p, names) in minimum if names isa AbstractDict
                                   for (s, v) in names]
    build_minimum_rows!(ed)
    weights = get(doc_data(ed), "atomic_weights", nothing)
    ed.weight_items = Vector{Any}[Any[string(e), v] for (e, v) in (weights isa AbstractDict ? weights : ())]
    build_weight_rows!(ed)
    ed.texts["measured.label"] = string(something(getpath(doc_data(ed), ("comparison", "label")), ""))
end

# minimum amounts
function build_minimum_rows!(ed::CaseEditor)
    forget!(ed, "advanced.min")
    for (i, (_, _, value)) in enumerate(ed.minimum_items)
        choice_field!(ed, "advanced.min-species:$(i - 1)", "Minimum amounts")
        field!(ed, "advanced.min:$(i - 1)", "Minimum amounts", number_text(value),
               text -> parse_number(text; required=true, minimum=0), v -> limit_part!(ed, i - 1, 3, v; refresh=false);
               shown="")
    end
    commit_limits!(ed)
end

"""Row `k` (0-based) of Minimum amounts: part 1 phase (clears the species), 2 species, 3 value
(`limit_phase!`, `limit_species!`)."""
function limit_part!(ed::CaseEditor, k, part, value; refresh=true)
    ed.minimum_items[k+1][part] = part == 3 ? value : something(value, "")
    if part == 1
        ed.minimum_items[k+1][2] = ""
        build_minimum_rows!(ed)
        ed.texts_changed = true
    end
    commit_limits!(ed)
    refresh && refresh!(ed)
    return ed
end

limit_phase!(ed::CaseEditor, k, phase) = limit_part!(ed, k, 1, phase)
limit_species!(ed::CaseEditor, k, species) = limit_part!(ed, k, 2, species)

function add_limit!(ed::CaseEditor)
    push!(ed.minimum_items, Any["gas", "", 1e-8])
    build_minimum_rows!(ed)
    ed.texts_changed = true
    refresh!(ed)
end

function remove_limit!(ed::CaseEditor, k)
    deleteat!(ed.minimum_items, k + 1)
    build_minimum_rows!(ed)
    ed.texts_changed = true
    refresh!(ed)
end

function commit_limits!(ed::CaseEditor)
    out = OrderedDict{String,Any}()
    for (i, (phase, species, value)) in enumerate(ed.minimum_items)
        key = "advanced.min-species:$(i - 1)"
        if isempty(phase) || isempty(species)
            ed.parse_errors[key] = "Choose a phase and a species, or remove the row."
        elseif haskey(get(out, phase, Dict()), species)
            ed.parse_errors[key] = "$species ($(lowercase(PHASE_NAMES[phase]))) is listed twice."
        else
            delete!(ed.parse_errors, key)
            value === nothing || (get!(out, phase) do
                OrderedDict{String,Any}()
            end[species] = value)
        end
    end
    putpath!(doc_data(ed), ["constraints", "minimum_amounts"], isempty(out) ? nothing : out)
end

# atomic weights
function build_weight_rows!(ed::CaseEditor)
    forget!(ed, "advanced.aw")
    for (i, (element, value)) in enumerate(ed.weight_items)
        field!(ed, "advanced.aw:$(i - 1)", "Atomic weight of $element", value === nothing ? "" : number_text(value),
               text -> parse_number(text; required=true, above=0), v -> weight_value!(ed, i, v); shown="")
    end
end

function commit_weights!(ed::CaseEditor)
    weights = OrderedDict{String,Any}(e => v for (e, v) in ed.weight_items if v !== nothing)
    putpath!(doc_data(ed), ["atomic_weights"], isempty(weights) ? nothing : weights)
end

weight_value!(ed::CaseEditor, i, value) = (ed.weight_items[i][2] = value; commit_weights!(ed))

"""Add an atomic-weight override (`text`: the symbol typed); it starts at the standard value."""
function add_weight!(ed::CaseEditor, text)
    symbol = element_symbol(text)
    ed.texts["new_weight"] = string(something(text, ""))
    ed.new_weight_error = element_error(symbol, [e for (e, _) in ed.weight_items])
    ed.new_weight_error === nothing || return refresh!(ed)
    push!(ed.weight_items, Any[symbol, get(HTGEquilibrium.ATOMIC_WEIGHTS, symbol, nothing)])
    commit_weights!(ed)
    build_weight_rows!(ed)
    ed.texts_changed = true
    refresh!(ed)
end

function remove_weight!(ed::CaseEditor, k)
    deleteat!(ed.weight_items, k + 1)
    commit_weights!(ed)
    build_weight_rows!(ed)
    ed.texts_changed = true
    refresh!(ed)
end

# measured data
"""CSV files to choose from, as (value, label): none, the files below the cases folder, the current one."""
function measured_options(ed::CaseEditor)
    files = comparison_files(ed.state, base_dir(current_doc(ed)))
    current = getpath(doc_data(ed), ("comparison", "file"))
    options = Tuple{String,String}[("", "(none)")]
    append!(options, files)
    current isa AbstractString && !isempty(current) && !any(o -> o[1] == current, files) && push!(options, (current, current))
    return options
end

"""The measured-data file or its unit changed (`part` "file" or "quantity")."""
function set_measured!(ed::CaseEditor, part, value)
    comparison = getpath(doc_data(ed), ("comparison",))
    current = comparison isa AbstractDict ? comparison : Dict()
    file = part == "file" ? something(value, "") : string(get(current, "file", ""))
    quantity = part == "quantity" ? value : get(current, "quantity", "mol_per_kg_dry_feed")
    commit_measured!(ed; file, quantity)
end

function commit_measured!(ed::CaseEditor; file=nothing, quantity=nothing)
    comparison = getpath(doc_data(ed), ("comparison",))
    current = comparison isa AbstractDict ? comparison : Dict()
    file = something(file, string(get(current, "file", "")))
    quantity = something(quantity, get(current, "quantity", "mol_per_kg_dry_feed"))
    label = get(ed.texts, "measured.label", "")
    if isempty(file)
        set!(ed, ["comparison"], nothing)
    else
        value = OrderedDict{String,Any}("file" => file, "quantity" => quantity)
        isempty(strip(label)) || (value["label"] = label)
        set!(ed, ["comparison"], value)
    end
end

# ---------------------------------------------------------------- selects and checkboxes
"""
The selects and checkboxes that set one value of the case: name => (path, value
-> what to store). A default is stored as `nothing` (no key in the file), as Python.
"""
const CHOICES = Dict{String,Tuple{Vector{String},Function}}(
    "aqueous_enabled" => (["aqueous", "enabled"], v -> v === true ? nothing : false),
    "aqueous_model" => (["aqueous", "model"], identity),
    "solid_target" => (["constraints", "solid_carbon_target"], v -> v == "char" ? nothing : v),
    "method" => (["solver", "method"], identity),
    "continuation" => (["sweep", "continuation"], v -> v === true ? nothing : false),
    "initialization" => (["initial_guess", "initialization"], v -> v == "auto" ? nothing : v),
    "stability" => (["solver", "stability_test"], v -> get(Dict("on" => true, "off" => false), v, nothing)),
)

"""A select or checkbox of `CHOICES` changed."""
function set_choice!(ed::CaseEditor, name, value)
    path, store = CHOICES[name]
    set!(ed, path, store(value))
end

# ---------------------------------------------------------------- derived numbers
"""Update what depends on the data: the feed check and the SUPCRT92 band check."""
function refresh!(ed::CaseEditor)
    ed.built || return ed
    update_feed!(ed)
    update_bands!(ed)
    return ed
end

"""
Check the case's temperatures against the SUPCRT92 failure bands when they
changed. A new pressure costs about 0.3 s, so the check runs in a background
task (`async_bands`); its result is shown when it is still current.
"""
function update_bands!(ed::CaseEditor)
    points = band_points(doc_data(ed))
    points == ed.band_points && return ed
    ed.band_points = points
    ed.band_generation += 1
    generation = ed.band_generation
    if isempty(points)
        ed.band_hits = []
        return ed
    end
    function run()
        hits = try
            band_hits(points; check=ed.band_check)
        catch exc
            exc isa InterruptException && rethrow()
            ed.band_warned || @warn "The SUPCRT92 band check failed; no band notice is shown" exception = exc
            ed.band_warned = true
            Tuple{Float64,Float64}[]
        end
        lock(ed.lock) do
            generation == ed.band_generation || return
            ed.band_hits = hits
            ed.async_bands && ed.on_change !== nothing && ed.on_change()
        end
    end
    ed.async_bands ? errormonitor(Threads.@spawn run()) : run()
    return ed
end

# ---------------------------------------------------------------- what the page shows
_options(dict) = [Dict("value" => k, "label" => v) for (k, v) in dict]

"""
    form_view(ed) -> Dict

Everything the Case tab shows except the field texts (`ed.texts`): file bar,
status line, notices, field labels and messages, the derived numbers, and the
state of the selects, toggles and dialogs.
"""
function form_view(ed::CaseEditor)
    doc, d = current_doc(ed), doc_data(ed)
    status, ok = status_line(ed)
    errors = Dict{String,String}()
    for key in keys(ed.fields)
        message = get(ed.parse_errors, key, get(ed.feed_errors, key, nothing))
        message === nothing || (errors[key] = message)
    end
    T, P = t_unit(ed), p_unit(ed)
    out = Dict{String,Any}(
        "ready" => true,
        "file" => Dict("status" => file_status(doc.path, dirty(doc)),
                       "selected" => doc.path !== nothing && any(f -> f[1] == doc.path, ed.case_files) ? doc.path : nothing,
                       "cases" => [Dict("value" => p, "label" => l) for (p, l) in ed.case_files]),
        "status" => status, "ok" => ok, "loading" => ed.state.model === nothing,
        "notices" => [n for n in notices(doc) if n != SALTS_NOTICE],
        "blocked" => doc.legacy_blocker !== nothing, "salts" => has_salts(doc),
        "labels" => Dict(k => f.shown for (k, f) in ed.fields),
        "placeholders" => Dict(k => f.placeholder for (k, f) in ed.fields if !isempty(f.placeholder)),
        "errors" => errors,
        "new_element_error" => ed.new_element_error, "new_weight_error" => ed.new_weight_error,
        "check" => ed.check, "saveas" => ed.saveas)
    ed.built || return out
    conditions = get(d, "conditions", Dict())
    out["conditions"] = Dict(
        "mode" => get(conditions, "points", nothing) !== nothing ? "points" : "grid",
        "t_unit" => T, "p_unit" => P,
        "count" => point_count_text(point_count(conditions, getpath(d, ("aqueous", "cutoff_C")))),
        "formats_open" => ed.state.settings["input_formats_open"],
        "band" => band_notice(ed.band_hits, T, P))
    feed = Dict{String,Any}("mode" => ed.feed_mode, "caption" => ed.feed_mode == "species" ? FEED_CAPTION_SPECIES : FEED_CAPTION,
                            "general" => ed.feed_message !== nothing && isempty(ed.feed_errors) ? ed.feed_message : nothing)
    summary = ed.feed_summary
    feed["summary"] = summary === nothing ? nothing :
                      [[row["name"], format_kg(row["kg"]), format_mol(row["mol"])] for row in summary["rows"]]
    feed["summary_message"] = summary === nothing ? something(ed.feed_message, "") : ""
    if ed.feed_mode == "elements"
        used, entered = try
            scaled_elements(something(get(feed_data(ed), "elements_wt_percent", nothing), OrderedDict()))
        catch exc
            exc isa FeedError || rethrow()
            Dict(), nothing
        end
        feed["units"] = Dict("solids" => ed.share_units["solids"], "ash" => ed.share_units["ash"])
        feed["hints"] = Dict("solids" => share_hint(ed, "solids"), "ash" => share_hint(ed, "ash"))
        feed["elements"] = ed.element_list
        feed["used"] = Dict(e => haskey(used, e) ? format_percent(used[e]) : "—" for e in ed.element_list)
        feed["entered_sum"] = entered === nothing ? "—" : format_percent(entered)
        feed["used_sum"] = isempty(used) ? "—" : format_percent(100.0)
        feed["scale_notice"] = entered === nothing ? nothing : entered_sum_notice(entered)
    else
        feed["water_unit"] = ed.water_unit
        feed["species"] = [Dict("name" => name, "unit" => unit) for (name, _, unit) in ed.species_items]
    end
    out["feed"] = feed
    out["products"] = Dict(phase => Any[string(s) for s in _values(getpath(d, ("phases", phase)))] for phase in keys(PHASE_NAMES))
    out["carbon"] = Dict(phase => Dict("choice" => row["choice"],
                                       "options" => _options(row["range"] ? merge(CARBON_CHOICES, Dict(CARBON_RANGE => CARBON_RANGE_LABEL)) :
                                                             CARBON_CHOICES))
                         for (phase, row) in ed.carbon)
    out["carbon_target"] = getpath(d, ("constraints", "solid_carbon_target"), "char")
    out["aqueous"] = Dict("enabled" => getpath(d, ("aqueous", "enabled"), true) !== false,
                          "model" => getpath(d, ("aqueous", "model"), "electrolyte"),
                          "water_model" => getpath(d, ("aqueous", "water_model"), "SUPCRT92"),
                          "water_model_reason" => water_model_reason(d),
                          "conflict" => conflict_notice(d, T))
    starts = getpath(d, ("solver", "starts"), 1)
    out["solver"] = Dict("method" => getpath(d, ("solver", "method"), "ipopt"),
                         "starts_tooltip" => initial_guesses_tooltip(starts),
                         "continuation" => getpath(d, ("sweep", "continuation"), true) !== false)
    stability = getpath(d, ("solver", "stability_test"))
    out["advanced"] = Dict(
        "initialization" => getpath(d, ("initial_guess", "initialization"), "auto"),
        "stability" => stability === true ? "on" : stability === false ? "off" : "default",
        "minimum" => [Dict("phase" => p, "species" => s) for (p, s, _) in ed.minimum_items],
        "weights" => [e for (e, _) in ed.weight_items],
        "measured" => Dict("file" => something(getpath(d, ("comparison", "file")), ""),
                           "quantity" => getpath(d, ("comparison", "quantity"), "mol_per_kg_dry_feed"),
                           "options" => [Dict("value" => v, "label" => l) for (v, l) in measured_options(ed)]))
    return out
end

"""Species names a select offers: `options` and the current ones (Python: options + current)."""
with_current(options, current) = unique(vcat(options, [string(c) for c in current if !isempty(string(c))]))

"""The message of `exc` on one line (as Python's str(exc): no "ArgumentError: " in front)."""
function one_line(exc)
    text = hasproperty(exc, :msg) && getproperty(exc, :msg) isa AbstractString ? exc.msg : sprint(showerror, exc)
    return first(join(split(strip(text), '\n'; keepempty=false), " "), 500)
end
