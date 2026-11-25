# The Case tab on the page (Python: gui/editor.py CaseEditor, the NiceGUI part).
#
# The form itself is the CaseEditor (editor.jl), one per page (client), over the
# shared `STATE[].document`. The page shows three fields of the model:
#
#   case_view     form_view(editor): labels, messages, derived numbers, choices
#   case_texts    the text of every field (sent only when the editor changed a
#                 text itself, so typing is never overwritten by an older text)
#   case_options  the long species lists (sent when the species data change)
#
# Every change in the page is one event `case_op` {op, args}; the server applies
# it to the editor and sends the new view. Genie handles the messages of one page
# in order, so the view after the last keystroke is the last one sent. Typing
# sends no field update of its own (no 300 ms debounce): the text goes with the
# event.

"""
A connected page's Case tab: its editor (once built) and the queue of its
events. The events are handled in order by a task of their own, not in the
page's message handler: a handler that waits (the first events compile code,
slowly while the solver's warm-up compiles) would hold up the page's other
messages and its keep-alive, and the page would lose its connection.
"""
mutable struct CasePage
    editor::Union{Nothing,CaseEditor}
    queue::Channel{Any}
end

"""
The Case tab of each connected page (a reloaded page gets a new one over the
same document). Keyed by identity: Stipple hashes a model by its field
values, so a `Dict` or `WeakKeyDict` would lose the page once a field changes.
The entry goes when the page's event task ends.
"""
const CASE_PAGES = IdDict{Any,CasePage}()
const CASE_PAGES_LOCK = ReentrantLock()

case_page(model) = lock(() -> get(CASE_PAGES, model, nothing), CASE_PAGES_LOCK)
case_editor(model) = (page = case_page(model); page === nothing ? nothing : page.editor)

"""A page connected: build its Case tab and handle its events, in a task of its own."""
function case_start!(model)
    page = CasePage(nothing, Channel{Any}(Inf))
    old = lock(CASE_PAGES_LOCK) do
        old = get(CASE_PAGES, model, nothing)
        CASE_PAGES[model] = page
        old
    end
    old === nothing || close(old.queue)
    errormonitor(Threads.@spawn begin
        try
            case_connect!(model, page)
        catch exc
            exc isa InterruptException && rethrow()
            @error "The Case tab could not start" exception = (exc, catch_backtrace())
        end
        for event in page.queue
            event === :species ? case_species_changed!(model) : case_event!(model, event)
        end
    end)
    return page
end

"""The page went away: its event task ends and the page is forgotten."""
function case_stop!(model)
    page = lock(() -> pop!(CASE_PAGES, model, nothing), CASE_PAGES_LOCK)
    page === nothing || close(page.queue)
    return nothing
end

"""Queue an event of the page (`:species`: the species data changed)."""
function case_enqueue!(model, event)
    page = case_page(model)
    if page === nothing || !isopen(page.queue)
        # An event of a page that is still (or again) connected starts its tab again.
        Stipple.isconnected(model) === false && return nothing
        page = case_start!(model)
    end
    put!(page.queue, event)
    return nothing
end

"""Make the page's editor and send the form."""
function case_connect!(model, page::CasePage)
    STATE[] === nothing && return nothing
    ed = CaseEditor(STATE[])
    ed.on_change = () -> push_case!(model, ed)
    page.editor = ed
    push_options!(model)
    lock(ed.lock) do
        ed.texts_changed = true
        push_case!(model, ed)
    end
    return ed
end

"""The species lists for the selects with search (per phase, and the feed species)."""
function push_options!(model)
    state = STATE[]
    model.case_options[] = Dict{String,Any}("phases" => Dict{String,Any}(k => v for (k, v) in state.phase_options),
                                            "feed" => state.feed_species, "loaded" => state.model !== nothing)
end

"""Send the view (and the texts when the editor changed them) and show the editor's messages."""
function push_case!(model, ed::CaseEditor)
    lock(ed.lock) do
        if ed.texts_changed
            model.case_texts[] = Dict{String,Any}(ed.texts)
            ed.texts_changed = false
        end
        v = form_view(ed)
        model.case_view[] = v
        for (text, type) in ed.notifications
            notify(model, text, Symbol(type); multiLine=type == "negative")
        end
        empty!(ed.notifications)
    end
    return nothing
end

"""The species data were loaded (or reloaded): new lists, and a species feed can be checked now."""
function case_species_changed!(model)
    ed = case_editor(model)
    ed === nothing && return
    push_options!(model)
    lock(ed.lock) do
        refresh!(ed)
        push_case!(model, ed)
    end
end

# ---------------------------------------------------------------- events
const CASE_OPS = Dict{String,Function}(
    "edit" => (ed, a) -> edit!(ed, a[1], a[2]),
    "new" => (ed, a) -> new!(ed),
    "open" => (ed, a) -> a[1] === nothing || open!(ed, a[1]),
    "conditions_mode" => (ed, a) -> switch_conditions!(ed, a[1]),
    "temperature_unit" => (ed, a) -> set_unit!(ed, "temperature_unit", a[1]),
    "pressure_unit" => (ed, a) -> set_unit!(ed, "pressure_unit", a[1]),
    "formats" => (ed, a) -> set_input_formats!(ed, a[1]),
    "feed_mode" => (ed, a) -> switch_feed!(ed, a[1]),
    "share_unit" => (ed, a) -> share_unit!(ed, a[1], a[2]),
    "water_unit" => (ed, a) -> water_unit!(ed, a[1]),
    "species_name" => (ed, a) -> species_name!(ed, Int(a[1]), a[2]),
    "species_unit" => (ed, a) -> species_unit!(ed, Int(a[1]), a[2]),
    "add_species" => (ed, a) -> add_species!(ed),
    "remove_species" => (ed, a) -> remove_species!(ed, Int(a[1])),
    "add_element" => (ed, a) -> add_element!(ed, get(ed.texts, "new_element", "")),
    "remove_element" => (ed, a) -> remove_element!(ed, a[1]),
    "products" => (ed, a) -> set_products!(ed, a[1], a[2]),
    "carbon_choice" => (ed, a) -> carbon_choice!(ed, a[1], a[2]),
    "choice" => (ed, a) -> set_choice!(ed, a[1], a[2]),
    "water_model" => (ed, a) -> set_water_model!(ed, a[1]),
    "add_limit" => (ed, a) -> add_limit!(ed),
    "limit_phase" => (ed, a) -> limit_phase!(ed, Int(a[1]), a[2]),
    "limit_species" => (ed, a) -> limit_species!(ed, Int(a[1]), a[2]),
    "remove_limit" => (ed, a) -> remove_limit!(ed, Int(a[1])),
    "add_weight" => (ed, a) -> add_weight!(ed, get(ed.texts, "new_weight", "")),
    "remove_weight" => (ed, a) -> remove_weight!(ed, Int(a[1])),
    "measured" => (ed, a) -> set_measured!(ed, a[1], a[2]),
    "close_check" => (ed, a) -> (ed.check = Dict{String,Any}("open" => false)),
    "saveas_cancel" => (ed, a) -> (ed.saveas = Dict{String,Any}("open" => false)),
)

"""One event of the Case tab: {op, args}."""
function case_event!(model, event)
    ed = case_editor(model)
    ed === nothing && return nothing
    op = string(get(event, "op", ""))
    args = Any[get(event, "args", Any[])...]
    try
        if haskey(CASE_OPS, op)
            lock(() -> CASE_OPS[op](ed, args), ed.lock)
        elseif op == "save"
            saved = lock(() -> save!(ed), ed.lock)
            saved || save_as!(model, ed)
        elseif op == "save_as"
            save_as!(model, ed)
        elseif op == "saveas_ok"
            lock(ed.lock) do
                ed.saveas = Dict{String,Any}("open" => false)
                write!(ed, with_extension(strip(get(ed.texts, "saveas.path", ""))))
            end
        elseif op == "queue"
            run = lock(() -> queue!(ed), ed.lock)
            run === nothing || (model.tab[] = "queue")
        elseif op == "check"
            check!(model, ed)
        else
            @warn "Unknown Case tab event" op
        end
    catch exc
        exc isa InterruptException && rethrow()
        @error "The Case tab could not do this ($op)" exception = (exc, catch_backtrace())
        lock(() -> notify!(ed, "Something went wrong: $(one_line(exc))", "negative"), ed.lock)
    end
    push_case!(model, ed)
    return nothing
end

"""A path typed or chosen without an extension gets .toml (so the case appears in the Case file list)."""
with_extension(path) = isempty(path) || !isempty(splitext(path)[2]) ? String(path) : path * ".toml"

"""
Save as: the window's native save dialog (owner decision PR #3); in the browser
(no native dialog) a dialog with the path as text, as in Python.
"""
function save_as!(model, ed::CaseEditor)
    default = lock(() -> default_save_path(ed), ed.lock)
    if NATIVE_SAVE[] === nothing
        lock(ed.lock) do
            ed.saveas = Dict{String,Any}("open" => true)
            ed.texts["saveas.path"] = default
            ed.texts_changed = true
        end
        return nothing
    end
    # The dialog waits for the user: not in the page's message handler.
    errormonitor(@async begin
        path = save_dialog(Dict("title" => "Save case file as", "defaultPath" => default,
                                "filters" => [Dict("name" => "Case files", "extensions" => ["toml"])]))
        if path !== nothing
            lock(() -> write!(ed, with_extension(path)), ed.lock)
            push_case!(model, ed)
        end
    end)
    return nothing
end

"""Check: the preflight checks of every point, in the background, with a spinner (Python: _check)."""
function check!(model, ed::CaseEditor)
    c = lock(() -> check_case(ed), ed.lock)
    c === nothing && return nothing
    model_ = STATE[].model
    Base.run(model, "this.caseCheckStart()")
    errormonitor(Threads.@spawn begin
        try
            messages = preflight_messages(model_, c)
            lock(() -> check_done!(ed, c, messages), ed.lock)
        catch exc
            exc isa InterruptException && rethrow()
            lock(() -> notify!(ed, "The check failed: $(one_line(exc))", "negative"), ed.lock)
        finally
            Base.run(model, "this.caseCheckEnd()")
            push_case!(model, ed)
        end
    end)
    return nothing
end

# ---------------------------------------------------------------- the template
# Raw Vue/Quasar markup (in-DOM template: kebab-case props, explicit closing tags).

"""Escape `s` for an HTML attribute in double quotes."""
_attr(s) = replace(string(s), "&" => "&amp;", "\"" => "&quot;", "<" => "&lt;")
"""`s` as a JavaScript string literal (single quotes)."""
_js(s) = "'" * replace(string(s), "\\" => "\\\\", "'" => "\\'") * "'"
"""Options for a q-select or q-btn-toggle, as a JavaScript array in an attribute."""
_opts(dict) = _attr(JSON.json([Dict("value" => k, "label" => v) for (k, v) in dict]))

"""
A text field of the form; `key` is a JavaScript expression giving the field's key.
Its label, message and placeholder come from the view.
"""
function _field(key; class="", textarea=false, cell=false, hint=nothing, mark=nothing, extra="")
    tag = "q-input"
    props = """:model-value="case_texts[$key]" @update:model-value="v => caseEdit($key, v)" :label="case_view.labels[$key] || undefined" """ *
            """:error="!!case_view.errors[$key]" :error-message="case_view.errors[$key]" """ *
            """:placeholder="case_view.placeholders[$key]" :stack-label="!!case_view.placeholders[$key]" """ *
            """:data-key="$key" """
    textarea && (props *= """type="textarea" autogrow """)
    cell && (props *= "dense hide-bottom-space ")
    hint === nothing || (props *= """:hint="$hint" """)
    return """<$tag $props class="$class" $extra></$tag>"""
end

"""A select whose value is in the view; `op` is sent with `args...` and the new value."""
function _select(value, op; options, label=nothing, class="", args="", mark=nothing, extra="")
    lbl = label === nothing ? "" : """label="$(_attr(label))" """
    m = mark === nothing ? "" : """data-mark="$mark" """
    pre = isempty(args) ? "" : args * ", "
    return """<q-select :model-value="$value" @update:model-value="v => caseOp('$op', $(pre)v)" :options="$options" """ *
           """emit-value map-options $lbl class="$class" $m $extra></q-select>"""
end

"""A select with search over a list of names (the filter runs in the page)."""
function _search_select(value, op; list, id, label=nothing, class="", args="", multiple=false, mark=nothing, footer=nothing,
                        error=nothing)
    lbl = label === nothing ? "" : """label="$(_attr(label))" """
    m = mark === nothing ? "" : """data-mark="$mark" """
    pre = isempty(args) ? "" : args * ", "
    mode = multiple ? "multiple use-chips use-input" : "use-input fill-input hide-selected"
    err = error === nothing ? "" : """:error="!!case_view.errors[$error]" :error-message="case_view.errors[$error]" """
    slot = footer === nothing ? "" :
           """<template v-slot:after-options><div class="text-caption text-grey-7 q-pa-sm">$(_attr(footer))</div></template>"""
    return """<q-select :model-value="$value" @update:model-value="v => caseOp('$op', $(pre)v)" """ *
           """:options="case_filtered[$id] || $list" @filter="(val, update) => caseFilter(val, update, $id, $list)" """ *
           """$mode input-debounce="0" $lbl $err class="$class" $m>$slot</q-select>"""
end

_toggle(value, op; options, args="", mark=nothing, class="") =
    """<q-btn-toggle :model-value="$value" @update:model-value="v => caseOp('$op', $(isempty(args) ? "" : args * ", ")v)" """ *
    """:options="$options" dense no-caps toggle-color="primary" class="$class" """ *
    (mark === nothing ? "" : """data-mark="$mark" """) * "></q-btn-toggle>"

_card(title, body; caption=nothing, class="htg-full") =
    """<q-card class="htg-card $class"><div class="text-subtitle1 text-weight-medium">$(_attr(title))</div>""" *
    (caption === nothing ? "" : """<div class="text-caption text-grey-7">$(_attr(caption))</div>""") * body * "</q-card>"

_notice(text; mark="case-notice", extra="") =
    """<div class="row items-center no-wrap bg-amber-1 q-pa-sm rounded-borders htg-full q-gutter-x-sm" data-mark="$mark" $extra>""" *
    """<q-icon name="info" color="warning" size="sm"></q-icon><div class="text-body2">$text</div></div>"""

_caption(text, class="") = """<div class="text-caption text-grey-7 $class">$(_attr(text))</div>"""

_remove(op, arg, mark="") =
    """<q-btn color="primary" flat dense round size="sm" icon="close" @click="caseOp('$op', $arg)" class="q-mt-sm" $mark></q-btn>"""

function _file_bar()
    """
<div class="row items-center no-wrap q-gutter-x-sm htg-full">
  <q-select :model-value="case_view.file.selected" @update:model-value="v => caseOp('open', v)" :options="case_filtered.cases || case_view.file.cases"
     @filter="(val, update) => caseFilterLabels(val, update, 'cases', case_view.file.cases)" use-input fill-input hide-selected input-debounce="0"
     emit-value map-options label="Case file" class="htg-w72" data-mark="case-file"></q-select>
  <q-btn color="primary" flat no-caps icon="note_add" label="New" @click="caseOp('new')" data-mark="case-new"></q-btn>
  <q-btn color="primary" flat no-caps icon="save" label="Save" @click="caseOp('save')" :disable="!case_view.ok" data-mark="case-save"></q-btn>
  <q-btn color="primary" flat no-caps icon="save_as" label="Save as…" @click="caseOp('save_as')" :disable="!case_view.ok" data-mark="case-save-as"></q-btn>
  <div class="text-caption text-grey-7" data-mark="case-file-status">{{ case_view.file.status }}</div>
  <q-space></q-space>
  <q-btn color="primary" outline no-caps icon="fact_check" label="Check" @click="caseOp('check')" data-mark="case-check"></q-btn>
  <q-btn color="primary" no-caps icon="play_arrow" label="Add to queue" @click="caseOp('queue')" :disable="!case_view.ok" data-mark="case-queue"></q-btn>
</div>
<div :class="'text-body2 ' + (case_view.ok ? 'text-positive' : 'text-negative')" data-mark="case-status">{{ case_view.status }}</div>
<template v-for="text in case_view.notices" :key="text">$(_notice("{{ text }}"))</template>
"""
end

function _name_card()
    _card("", """
<div class="row no-wrap items-start htg-full q-gutter-x-md">
  $(_field("'name.name'"; class="htg-w80"))
  <q-input :model-value="case_texts.notes" @update:model-value="v => caseEdit('notes', v)" label="Notes" type="textarea" autogrow rows="1" class="col-grow" data-key="notes"></q-input>
</div>""") |> s -> replace(s, """<div class="text-subtitle1 text-weight-medium"></div>""" => "")
end

function _conditions_card()
    formats = join(("<tr><td>$(_attr(f))</td><td><code>$(_attr(t))</code></td><td>$(_attr(g))</td></tr>" for (f, t, g) in INPUT_FORMATS), "")
    _card("Temperature and pressure", """
<div class="row items-center htg-gap16">
  $(_toggle("case_view.conditions.mode", "conditions_mode"; options=_opts(MODE_TOGGLE), mark="case-conditions-mode"))
  $(_select("case_view.conditions.t_unit", "temperature_unit"; options=_opts(OrderedDict(u => Units.label(u) for u in Units.TEMPERATURE_UNITS)),
            label="Temperature unit", class="htg-w36", mark="case-temperature-unit"))
  $(_select("case_view.conditions.p_unit", "pressure_unit"; options=_opts(OrderedDict(u => Units.label(u) for u in Units.PRESSURE_UNITS)),
            label="Pressure unit", class="htg-w36", mark="case-pressure-unit"))
</div>
<div class="row items-start no-wrap htg-gap16 htg-full">
  <template v-if="case_view.conditions.mode == 'points'">$(_field("'conditions.points'"; textarea=true, class="col-grow"))</template>
  <template v-else>
    $(_field("'conditions.temperatures'"; class="col-grow"))
    $(_field("'conditions.pressures'"; class="htg-w80"))
  </template>
  <div class="text-body2 q-mt-md text-no-wrap" data-mark="case-point-count">{{ case_view.conditions.count }}</div>
</div>
<div v-if="case_view.conditions.band" class="htg-full">$(_notice("{{ case_view.conditions.band }}"; mark="case-band"))</div>
<q-expansion-item dense label="Input formats" :model-value="case_view.conditions.formats_open"
   @update:model-value="v => { case_view.conditions.formats_open = v; caseOp('formats', v) }" class="htg-full" data-mark="case-input-formats">
  <q-markup-table dense flat><thead><tr><th class="text-left">Format</th><th class="text-left">Type</th><th class="text-left">Gives</th></tr></thead>
  <tbody>$formats</tbody></q-markup-table>
</q-expansion-item>""")
end

function _elemental()
    """
$(_field("'feed.total'"; class="htg-w56", hint=_js("Solids + water")))
<template v-for="name in ['solids', 'ash']" :key="name">
  <div class="row items-start no-wrap q-gutter-x-sm">
    $(_field("'feed.' + name"; class="htg-w56", hint="case_view.feed.hints[name]"))
    $(_toggle("case_view.feed.units[name]", "share_unit"; args="name", options=_opts(OrderedDict("kg" => "kg", "%" => "%")),
              class="q-mt-sm"))
  </div>
</template>
<div class="row items-baseline no-wrap htg-full q-mt-sm htg-gap24">
  <div>4 · Elements [mass %, ash-free solids]</div><q-space></q-space>$(_caption(ELEMENTS_CAPTION))
</div>
<div class="htg-elements">
  <div class="text-caption text-grey-7">Element</div><div class="text-caption text-grey-7">Entered</div>
  <div class="text-caption text-grey-7">Used</div><div></div>
  <template v-for="e in case_view.feed.elements" :key="e">
    <div>{{ e }}</div>
    $(_field("'feed.element:' + e"; cell=true, class="htg-w28"))
    <div :data-used="e">{{ case_view.feed.used[e] }}</div>
    <q-btn color="primary" flat dense round size="sm" icon="close" @click="caseOp('remove_element', e)" :data-remove="e"></q-btn>
  </template>
  <div class="text-weight-medium">Sum</div><div class="text-weight-medium" data-mark="case-entered-sum">{{ case_view.feed.entered_sum }}</div>
  <div class="text-weight-medium">{{ case_view.feed.used_sum }}</div><div></div>
</div>
<div class="row items-start q-gutter-x-sm">
  <q-input :model-value="case_texts.new_element" @update:model-value="v => caseEdit('new_element', v)" label="Element symbol"
     :error="!!case_view.new_element_error" :error-message="case_view.new_element_error" class="htg-w32" data-key="new_element"></q-input>
  <q-btn color="primary" flat dense no-caps label="+ Add element" @click="caseOp('add_element')" class="q-mt-sm" data-mark="case-add-element"></q-btn>
</div>
<div v-if="case_view.feed.scale_notice" class="text-body2 text-warning" data-mark="case-scale-notice">{{ case_view.feed.scale_notice }}</div>
"""
end

function _species_amounts()
    """
<div class="row items-start no-wrap q-gutter-x-sm">
  <div class="htg-w64 q-mt-md">Water (always present)</div>
  $(_field("'feed.water'"; class="htg-w28"))
  $(_toggle("case_view.feed.water_unit", "water_unit"; options=_opts(OrderedDict("mol" => "mol", "kg" => "kg")), class="q-mt-sm",
            mark="case-feed.water-unit"))
</div>
<template v-for="(row, i) in case_view.feed.species" :key="i + ':' + row.name">
  <div class="row items-start no-wrap q-gutter-x-sm">
    $(_search_select("row.name || null", "species_name"; list="case_options.feed", id="'species'", label="Species", class="htg-w64",
                     args="i", footer=SPECIES_FOOTER, error="'feed.species-name:' + i"))
    $(_field("'feed.species:' + i"; class="htg-w28"))
    $(_toggle("row.unit", "species_unit"; args="i", options=_opts(OrderedDict("mol" => "mol", "kg" => "kg")), class="q-mt-sm"))
    <q-btn color="primary" flat dense round size="sm" icon="close" @click="caseOp('remove_species', i)" class="q-mt-sm"></q-btn>
  </div>
</template>
<q-btn color="primary" flat dense no-caps label="+ Add species" @click="caseOp('add_species')" data-mark="case-add-species"></q-btn>
"""
end

function _feed_card()
    _card("Feed", """
<div class="text-caption text-grey-7">{{ case_view.feed.caption }}</div>
$(_toggle("case_view.feed.mode", "feed_mode"; options=_opts(FEED_TOGGLE), mark="case-feed-mode"))
<div class="row items-start no-wrap htg-full htg-gap40">
  <div class="column items-start htg-gap8">
    <template v-if="case_view.feed.mode == 'species'">$(_species_amounts())</template>
    <template v-else>$(_elemental())</template>
    <div v-if="case_view.feed.general" class="text-caption text-negative">{{ case_view.feed.general }}</div>
  </div>
  <div class="column q-ml-auto">
    <div class="text-weight-medium">Feed summary</div>
    <div v-if="case_view.feed.summary" class="htg-summary" data-mark="case-feed-summary">
      <div></div><div class="text-caption text-grey-7 text-right">kg</div><div class="text-caption text-grey-7 text-right">mol</div>
      <template v-for="row in case_view.feed.summary" :key="row[0]">
        <div>{{ row[0] }}</div><div class="text-right">{{ row[1] }}</div><div class="text-right">{{ row[2] }}</div>
      </template>
    </div>
    <div v-else class="text-caption text-grey-7 htg-w64">{{ case_view.feed.summary_message }}</div>
  </div>
</div>""")
end

function _products_card()
    selects = join((_search_select("case_view.products.$phase", "products"; list="(case_options.phases.$phase || [])",
                                   id=_js(phase), label=label, class="htg-full", args=_js(phase), multiple=true,
                                   mark="case-phase-$phase")
                    for (phase, label) in PHASE_NAMES), "\n")
    _card("Possible products", selects; caption=PRODUCTS_CAPTION)
end

function _carbon_card()
    rows = """
<template v-for="phase in $(_attr(JSON.json(collect(CARBON_PHASES))))" :key="phase">
  <div class="row items-start no-wrap q-gutter-x-sm">
    <div class="q-mt-md htg-w28">{{ $(_attr(JSON.json(PHASE_NAMES)))[phase] }}</div>
    <q-select :model-value="case_view.carbon[phase].choice" @update:model-value="v => caseOp('carbon_choice', phase, v)"
       :options="case_view.carbon[phase].options" emit-value map-options class="htg-w28" :data-mark="'case-carbon-' + phase"></q-select>
    <div v-show="case_view.carbon[phase].choice != 'free'">$(_field("'carbon.' + phase"; class="htg-w24"))</div>
    <div v-show="case_view.carbon[phase].choice == 'range'">$(_field("'carbon.' + phase + '.max'"; class="htg-w24"))</div>
  </div>
</template>
$(_select("case_view.carbon_target", "choice"; args="'solid_target'", options=_opts(SOLID_CARBON), label="Solid carbon means", class="htg-w64"))
"""
    _card("Carbon split", rows; caption=CARBON_CAPTION, class="col htg-basis0")
end

function _aqueous_card()
    _card("Aqueous phase", """
<q-checkbox :model-value="case_view.aqueous.enabled" @update:model-value="v => caseOp('choice', 'aqueous_enabled', v)"
   label="Include an aqueous phase" data-mark="case-aqueous-enabled"></q-checkbox>
$(_select("case_view.aqueous.model", "choice"; args="'aqueous_model'", options=_opts(AQUEOUS_MODELS), label="Model", class="htg-full",
          mark="case-aqueous-model"))
<div class="column htg-gap4">
  <div class="text-caption text-grey-7">$(_attr(WATER_MODEL_LABEL))</div>
  <div>
    $(_toggle("case_view.aqueous.water_model", "water_model"; options=_opts(WATER_MODEL_CHOICES), mark="case-water-model")
      |> s -> replace(s, "dense no-caps" => """dense no-caps :disable="!!case_view.aqueous.water_model_reason" """))
    <q-tooltip class="text-body2 htg-tooltip">{{ case_view.aqueous.water_model_reason ? case_view.aqueous.water_model_reason + ' ' : '' }}$(_attr(WATER_MODEL_TOOLTIP))</q-tooltip>
  </div>
  <div v-if="case_view.aqueous.water_model_reason" class="text-caption text-grey-7" data-mark="case-water-model-reason">{{ case_view.aqueous.water_model_reason }}</div>
</div>
$(_field("'aqueous.cutoff'"; class="htg-full"))
$(_caption(CUTOFF_CAPTION))
<div v-if="case_view.aqueous.conflict" class="text-body2 text-warning" data-mark="case-conflict">{{ case_view.aqueous.conflict }}</div>
"""; class="col htg-basis0")
end

function _solver_card()
    _card("Solver", """
<div class="htg-full">
  $(_select("case_view.solver.method", "choice"; args="'method'", options=_opts(METHODS), label="Method", class="htg-full", mark="case-method"))
  <q-tooltip class="text-body2 htg-tooltip">$(_attr(METHOD_TOOLTIP))</q-tooltip>
</div>
<div class="htg-full">
  $(_field("'solver.starts'"; class="htg-full"))
  <q-tooltip class="text-body2 htg-tooltip" data-mark="case-starts-tooltip">{{ case_view.solver.starts_tooltip }}</q-tooltip>
</div>
<div>
  <q-checkbox :model-value="case_view.solver.continuation" @update:model-value="v => caseOp('choice', 'continuation', v)" label="Continuation"></q-checkbox>
  <q-tooltip class="text-body2 htg-tooltip">$(_attr(CONTINUATION_TOOLTIP))</q-tooltip>
</div>
"""; class="col htg-basis0")
end

_optional_field(key, caption=nothing) = _field(_js(key); class="htg-w72") * (caption === nothing ? "" : _caption(caption, "htg-w72"))

function _advanced()
    """
<q-expansion-item icon="tune" label="Advanced settings" caption="$(_attr(ADVANCED_CAPTION))" class="htg-full" data-mark="case-advanced">
<div class="q-pa-md">
<div class="row no-wrap items-start htg-gap32">
  <div class="column htg-gap4">
    <div class="text-weight-medium">Starting point</div>
    $(_select("case_view.advanced.initialization", "choice"; args="'initialization'", options=_opts(INITIAL_COMPOSITIONS),
              label="Initial composition", class="htg-w72"))
    $(_caption(INITIAL_COMPOSITION_CAPTION, "htg-w72"))
    $(_optional_field("advanced.gas_water"))
    $(_optional_field("advanced.water_floor"))
  </div>
  <div class="column htg-gap4">
    <div class="text-weight-medium">Solver details</div>
    $(_optional_field("advanced.max_iterations"))
    $(_optional_field("advanced.seed", SEED_CAPTION))
    $(_select("case_view.advanced.stability", "choice"; args="'stability'", options=_opts(STABILITY_TEST), label="Phase stability test",
              class="htg-w72"))
    $(_optional_field("advanced.energy_scale", ENERGY_SCALE_CAPTION))
  </div>
  <div class="column htg-gap4">
    <div class="text-weight-medium">When are two results different?</div>
    $(_caption(DISTINCT_CAPTION, "htg-w72"))
    $(_optional_field("advanced.composition_tol", COMPOSITION_TOL_CAPTION))
    $(_optional_field("advanced.gibbs_tol", GIBBS_TOL_CAPTION))
  </div>
</div>
<div class="row no-wrap items-start htg-gap32 q-mt-md">
  <div class="column htg-gap4">
    <div class="text-weight-medium">Minimum amounts [mol / mol C]</div>
    <div v-if="case_view.advanced.minimum.length" class="htg-minimum">
      <div class="text-caption text-grey-7">Phase</div><div class="text-caption text-grey-7">Species</div>
      <div class="text-caption text-grey-7">Minimum</div><div></div>
      <template v-for="(row, i) in case_view.advanced.minimum" :key="i + ':' + row.phase + ':' + row.species">
        $(_select("row.phase || null", "limit_phase"; args="i", options=_opts(PHASE_NAMES), class="htg-w36"))
        $(_search_select("row.species || null", "limit_species"; list="(case_options.phases[row.phase] || [])", id="'min-' + row.phase",
                         class="htg-w56", args="i", error="'advanced.min-species:' + i"))
        $(_field("'advanced.min:' + i"; cell=true, class="htg-w28"))
        <q-btn color="primary" flat dense round size="sm" icon="close" @click="caseOp('remove_limit', i)" class="q-mt-sm"></q-btn>
      </template>
    </div>
    <q-btn color="primary" flat dense no-caps label="+ Add limit" @click="caseOp('add_limit')" data-mark="case-add-limit"></q-btn>
  </div>
  <div class="column htg-gap4">
    <div class="text-weight-medium">Atomic weights</div>
    $(_caption(ATOMIC_WEIGHTS_CAPTION, "htg-w72"))
    <div v-if="case_view.advanced.weights.length" class="htg-weights">
      <div class="text-caption text-grey-7">Element</div><div class="text-caption text-grey-7">g/mol</div><div></div>
      <template v-for="(e, i) in case_view.advanced.weights" :key="i + ':' + e">
        <div class="q-mt-sm">{{ e }}</div>
        $(_field("'advanced.aw:' + i"; cell=true, class="htg-w28"))
        <q-btn color="primary" flat dense round size="sm" icon="close" @click="caseOp('remove_weight', i)" class="q-mt-sm"></q-btn>
      </template>
    </div>
    <div class="row items-start q-gutter-x-sm">
      <q-input :model-value="case_texts.new_weight" @update:model-value="v => caseEdit('new_weight', v)" label="Element symbol"
         :error="!!case_view.new_weight_error" :error-message="case_view.new_weight_error" class="htg-w32" data-key="new_weight"></q-input>
      <q-btn color="primary" flat dense no-caps label="+ Add element" @click="caseOp('add_weight')" class="q-mt-sm" data-mark="case-add-weight"></q-btn>
    </div>
  </div>
  <div class="column htg-gap4">
    <div class="text-weight-medium">Measured data to show in plots</div>
    $(_select("case_view.advanced.measured.file", "measured"; args="'file'", options="case_view.advanced.measured.options",
              label="CSV file", class="htg-w80"))
    $(_select("case_view.advanced.measured.quantity", "measured"; args="'quantity'", options=_opts(QUANTITY_CHOICES),
              label="Its unit", class="htg-w80"))
    <q-input :model-value="case_texts['measured.label']" @update:model-value="v => caseEdit('measured.label', v)"
       label="Label in the legend" class="htg-w80"></q-input>
    $(_caption(MEASURED_CAPTION, "htg-w80"))
  </div>
</div>
</div>
</q-expansion-item>
"""
end

function _dialogs()
    """
<q-dialog :model-value="!!case_view.check.open" @update:model-value="v => v || caseOp('close_check')">
  <q-card class="htg-card" style="width: 48rem; max-width: 100%" data-mark="case-check-dialog">
    <div :class="'text-h6 ' + (case_view.check.ok ? 'text-positive' : 'text-negative')">{{ case_view.check.title }}</div>
    <div class="column htg-gap4">
      <div v-for="m in case_view.check.messages" class="text-body2">{{ m }}</div>
      <div v-if="case_view.check.more" class="text-body2">{{ case_view.check.more }}</div>
    </div>
    <div class="text-caption">$(_attr(CHECK_CAPTION))</div>
    <q-btn color="primary" flat label="Close" @click="caseOp('close_check')"></q-btn>
  </q-card>
</q-dialog>
<q-dialog :model-value="!!case_view.saveas.open" @update:model-value="v => v || caseOp('saveas_cancel')">
  <q-card class="htg-card" style="width: 36rem; max-width: 100%" data-mark="case-saveas-dialog">
    <div class="text-h6">Save case file as</div>
    <q-input :model-value="case_texts['saveas.path']" @update:model-value="v => caseEdit('saveas.path', v)" label="File" class="htg-full"
       data-key="saveas.path"></q-input>
    <div class="row q-gutter-x-sm">
      <q-btn label="Save" color="primary" @click="caseOp('saveas_ok')"></q-btn>
      <q-btn color="primary" flat label="Cancel" @click="caseOp('saveas_cancel')"></q-btn>
    </div>
  </q-card>
</q-dialog>
"""
end

"""The Case tab (the panel's content)."""
function case_tab_html()
    form = """
<div v-if="!case_view.blocked" class="column htg-gap12 htg-full">
  $(_name_card())
  $(_conditions_card())
  <div v-if="case_view.salts" class="htg-full">$(_notice(_attr(SALTS_NOTICE)))</div>
  $(_feed_card())
  $(_products_card())
  <div class="row no-wrap items-stretch htg-gap12 htg-full">
    $(_carbon_card())
    $(_aqueous_card())
    $(_solver_card())
  </div>
  $(_advanced())
</div>"""
    return """
<div v-if="case_view.ready" class="column htg-gap8 htg-case" data-mark="case-tab">
$(_file_bar())
$form
$(_dialogs())
</div>
<div v-else class="text-body1 text-grey-7 htg-case-starting">Starting the Case tab.</div>
"""
end

# Layout as NiceGUI's (cards with 16 px padding and gap, rows, Tailwind widths).
const CASE_CSS = """
<style>
.htg-card { padding: 16px; display: flex; flex-direction: column; align-items: flex-start; gap: 16px; }
.htg-full { width: 100%; }
.htg-basis0 { flex-basis: 0; }
.htg-gap4 { gap: 4px; } .htg-gap8 { gap: 8px; } .htg-gap12 { gap: 12px; } .htg-gap16 { gap: 16px; }
.htg-gap24 { gap: 24px; } .htg-gap32 { gap: 32px; } .htg-gap40 { gap: 40px; }
.htg-w24 { width: 6rem; } .htg-w28 { width: 7rem; } .htg-w32 { width: 8rem; } .htg-w36 { width: 9rem; }
.htg-w56 { width: 14rem; } .htg-w64 { width: 16rem; } .htg-w72 { width: 18rem; } .htg-w80 { width: 20rem; }
.htg-elements { display: grid; grid-template-columns: 3rem 7rem 5rem 2rem; column-gap: 16px; align-items: center; }
.htg-summary { display: grid; grid-template-columns: auto auto auto; column-gap: 24px; }
.htg-minimum { display: grid; grid-template-columns: 9rem 14rem 7rem 2rem; column-gap: 8px; align-items: start; }
.htg-weights { display: grid; grid-template-columns: 3rem 8rem 2rem; column-gap: 8px; align-items: start; }
.htg-tooltip { max-width: 24rem; }
.htg-case .q-btn-toggle { box-shadow: 0 1px 3px rgba(0,0,0,.2); }
</style>"""

@methods """
caseEdit(key, value) {
  this.case_texts[key] = value;
  this.handle_event({op: 'edit', args: [key, value]}, 'case_op');
},
caseOp(op, ...args) {
  this.handle_event({op: op, args: args}, 'case_op');
},
caseFilter(val, update, id, list) {
  update(() => {
    const n = (val || '').toLowerCase();
    this.case_filtered[id] = n === '' ? null : list.filter(v => v.toLowerCase().indexOf(n) > -1);
  });
},
caseFilterLabels(val, update, id, list) {
  update(() => {
    const n = (val || '').toLowerCase();
    this.case_filtered[id] = n === '' ? null : list.filter(o => o.label.toLowerCase().indexOf(n) > -1);
  });
},
caseCheckStart() {
  this.caseCheckDismiss = this.\$q.notify({spinner: true, message: 'Checking every point…', timeout: 0, group: false});
},
caseCheckEnd() {
  if (this.caseCheckDismiss) { this.caseCheckDismiss(); this.caseCheckDismiss = null; }
}
"""

@client_data begin
    case_filtered = Dict{String,Any}()
    caseCheckDismiss = nothing
end

@event :case_op begin
    case_enqueue!(__model__, event)
end
