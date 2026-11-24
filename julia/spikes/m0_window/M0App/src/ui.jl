
# The page. Hooks set by shell.jl so the handlers can call Electron's native dialogs
# (Electron is only used in shell.jl). Each takes a Dict of dialog options and returns
# a path (String) or nothing.
const NATIVE_OPEN = Ref{Any}(nothing)
const NATIVE_SAVE = Ref{Any}(nothing)
# Called once per client when the page is connected (shell.jl uses it to time start-up).
const ON_READY = Ref{Any}(nothing)

const SPECIES = vcat(
    ["water", "hydrogen", "methane", "ethane", "propane", "n-butane", "carbon dioxide",
     "carbon monoxide", "nitrogen", "ammonia", "hydrogen sulfide", "phenol", "cellulose",
     "glucose", "acetic acid", "formic acid", "methanol", "ethanol", "graphite", "sodium chloride",
     "Na+", "Cl-", "HCO3-", "CO3-2", "OH-", "H+", "K+", "SO4-2", "HS-", "NH4+"],
    ["n-alkane C$(n)H$(2n + 2)" for n in 5:60],
    ["aromatic C$(n)H$(n)" for n in 6:40])

const RUNS = ["run $(i): fig$(70 + i) $(300 + 10i) °C" for i in 1:25]

to_kelvin(t, unit) = unit == "K" ? t : t + 273.15
from_kelvin(t, unit) = unit == "K" ? t : t - 273.15

function curves(unit)
    tk = collect(573.15:10.0:873.15)
    x = from_kelvin.(tk, unit)
    h2 = @. 0.05 + 0.4 / (1 + exp(-(tk - 720) / 25))
    ch4 = @. 0.3 / (1 + exp(-(tk - 650) / 30))
    co2 = @. 0.25 + 0.0003 * (tk - 573.15)
    [PlotData(x = x, y = y, name = n, plot = StipplePlotly.Charts.PLOT_TYPE_SCATTER, mode = "lines+markers")
     for (n, y) in (("H2", h2), ("CH4", ch4), ("CO2", co2))]
end

make_layout(unit) = PlotLayout(
    title = PlotLayoutTitle(text = "Gas composition"),
    xaxis = [PlotLayoutAxis(xy = "x", title_text = "Temperature ($(unit))")],
    yaxis = [PlotLayoutAxis(xy = "y", title_text = "Mole fraction")],
    height = 360)

function species_table()
    DataTable(DataFrame(
        Species = SPECIES[1:12],
        Phase = repeat(["gas", "aqueous", "organic liquid"], 4),
        var"Amount (mol)" = round.(collect(range(0.1, 12.0; length = 12)) .^ 1.3; digits = 3),
        var"Mass (kg)" = round.(collect(range(0.001, 0.4; length = 12)); digits = 4)))
end

# A tiny SVG as a data URL: stands in for the PNG preview of the save-plot dialog.
preview_svg(unit) = "data:image/svg+xml;base64," * base64encode("""
<svg xmlns="http://www.w3.org/2000/svg" width="320" height="120"><rect width="320" height="120" fill="#eef"/>
<text x="20" y="65" font-size="20" font-family="sans-serif">Preview, axis in $(unit)</text></svg>""")

@app begin
    # header, tabs, status line updated by a server-side timer
    @in tab = "case"
    @out status = "Workers: 0 running"
    @private ticks = 0

    # select with search (single) and multi-select with chips
    @in species = "methane"
    @out species_all = SPECIES
    @in runs = String[]
    @out runs_all = RUNS

    # number field with unit switch
    @in temperature = 350.0
    @in temp_unit = "°C"
    @private prev_unit = "°C"
    @out unit_options = [Dict("label" => "°C", "value" => "°C"), Dict("label" => "K", "value" => "K")]

    # other form controls
    @in method = "ipopt"
    @out method_options = [Dict("label" => "IPOPT", "value" => "ipopt"),
                           Dict("label" => "SLSQP", "value" => "slsqp")]
    @in continuation = true
    @in starts = 4
    @out starts_tooltip = "4 initial guesses per point"
    @in notes = "Case notes, several lines.\nSecond line."
    @in chip_gas = true
    @in chip_aqueous = false
    @in chip_organic = true
    @out progress = 0.35
    @out log_text = "12:00:01 run 1 started\n12:00:04 point 1/31 done"

    # table, dialog, plot, image
    @out results = species_table()
    @in show_dialog = false
    @in dialog_ok = false
    @out plot_data::Vector{PlotData} = curves("°C")
    @out plot_layout::PlotLayout = make_layout("°C")
    @out preview_src = preview_svg("°C")

    # buttons handled on the server
    @in notify_btn = false
    @in check_btn = false
    @in native_open_btn = false
    @in native_save_btn = false
    @out chosen_path = "(none)"
    @out upload_info = "(no upload yet)"

    @onchange isready begin
        isready || return
        ON_READY[] === nothing || ON_READY[](__model__)
        # status line: a server-side timer pushes a new value every second
        @async while isready && Stipple.isconnected(__model__) !== false
            ticks += 1
            status = "Workers: $(ticks % 4) running · $(ticks) s"
            sleep(1)
        end
    end

    @onchange temp_unit begin
        temp_unit == prev_unit && return
        temperature = round(from_kelvin(to_kelvin(temperature, prev_unit), temp_unit); digits = 2)
        prev_unit = temp_unit
        plot_data = curves(temp_unit)
        plot_layout = make_layout(temp_unit)
        preview_src = preview_svg(temp_unit)
    end

    @onchange starts begin
        starts_tooltip = "$(starts) initial guess$(starts == 1 ? "" : "es") per point"
    end

    @onbutton notify_btn begin
        notify(__model__, "Saved case file.", :positive)
    end

    @onbutton check_btn begin
        # persistent notification with a spinner, dismissed from the server
        Base.run(__model__, "this.startCheck()")
        sleep(1.5)
        Base.run(__model__, "this.endCheck()")
        notify(__model__, "All points checked.", :info)
    end

    @onbutton dialog_ok begin
        show_dialog = false
        notify(__model__, "Dialog confirmed on the server.")
    end

    @onbutton native_open_btn begin
        f = NATIVE_OPEN[]
        if f === nothing
            notify(__model__, "No native dialogs (not running in Electron).", :warning)
        else
            path = f(Dict("title" => "Import species", "properties" => ["openFile"],
                          "filters" => [Dict("name" => "Species files", "extensions" => ["csv", "toml", "dat"]),
                                        Dict("name" => "All files", "extensions" => ["*"])]))
            chosen_path = path === nothing ? "(cancelled)" : "open: " * path
        end
    end

    @onbutton native_save_btn begin
        f = NATIVE_SAVE[]
        if f === nothing
            notify(__model__, "No native dialogs (not running in Electron).", :warning)
        else
            path = f(Dict("title" => "Save plot", "defaultPath" => "plot.csv",
                          "filters" => [Dict("name" => "CSV table", "extensions" => ["csv"])]))
            if path === nothing
                chosen_path = "(cancelled)"
            else
                write(path, "T,H2\n300,0.1\n")
                chosen_path = "saved: " * path
            end
        end
    end

    @onchange fileuploads begin
        isempty(fileuploads) && return
        name = get(fileuploads, "name", "?")
        tmp = get(fileuploads, "path", "")
        n = isfile(tmp) ? filesize(tmp) : -1
        upload_info = "uploaded $(name): $(n) bytes (server copy $(tmp))"
        notify(__model__, "Uploaded $(name)")
        fileuploads = Dict{AbstractString,AbstractString}()
    end
end

# Client-side helpers: filtering for the searchable selects (Quasar's @filter),
# the spinner notification, and a browser download (Blob + <a download>).
@methods """
filterSpecies(val, update) {
  update(() => {
    const n = val.toLowerCase();
    this.species_filtered = n === '' ? this.species_all : this.species_all.filter(v => v.toLowerCase().indexOf(n) > -1);
  })
},
filterRuns(val, update) {
  update(() => {
    const n = val.toLowerCase();
    this.runs_filtered = n === '' ? this.runs_all : this.runs_all.filter(v => v.toLowerCase().indexOf(n) > -1);
  })
},
startCheck() {
  this.dismissCheck = this.\$q.notify({spinner: true, message: 'Checking every point…', timeout: 0, group: false});
},
endCheck() {
  if (this.dismissCheck) { this.dismissCheck(); this.dismissCheck = null; }
},
browserDownload() {
  const blob = new Blob(['T,H2\\n300,0.1\\n'], {type: 'text/csv'});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = 'browser-download.csv';
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}
"""

@client_data begin
    species_filtered = SPECIES
    runs_filtered = RUNS
    dismissCheck = nothing
end

const FORM_ROW = "items-center q-gutter-md q-mb-md"

function case_tab()
    [
        row(class = FORM_ROW, [
            select(:species, options = :species_filtered, label = "Species (type to search)",
                   useinput = true, var"fill-input" = true, var"hide-selected" = true, inputdebounce = 0,
                   clearable = true, dense = true, style = "width: 18rem",
                   @on(:filter, "filterSpecies")),
            numberfield("Temperature", :temperature, dense = true, style = "width: 10rem",
                        var":rules" = "[v => v !== null && v !== '' || 'Enter a number']"),
            btntoggle(:temp_unit, options = :unit_options, dense = true, nocaps = true,
                      var"toggle-color" = "primary"),
            numberfield("Initial guesses", :starts, dense = true, style = "width: 9rem",
                        [tooltip("{{ starts_tooltip }}", class = "text-body2")]),
        ]),
        row(class = FORM_ROW, [
            select(:runs, options = :runs_filtered, label = "Runs to compare", multiple = true,
                   usechips = true, useinput = true, inputdebounce = 0, dense = true,
                   style = "min-width: 24rem", @on(:filter, "filterRuns")),
        ]),
        row(class = FORM_ROW, [
            # NiceGUI's ui.radio is a q-option-group
            quasar(:option__group, fieldname = :method, var":options" = "method_options",
                   inline = true, dense = true),
            radio("IPOPT (q-radio)", :method, val = "ipopt", dense = true),
            radio("SLSQP (q-radio)", :method, val = "slsqp", dense = true),
            checkbox("Continuation", :continuation, dense = true),
            icon("info", color = "warning", [tooltip("A tooltip on an icon")]),
        ]),
        row(class = FORM_ROW, [
            span("Phases: "),
            chip("Gas", selected = :chip_gas, clickable = true, color = "blue-1", textcolor = "blue-10"),
            chip("Aqueous", selected = :chip_aqueous, clickable = true, color = "blue-1", textcolor = "blue-10"),
            chip("Organic liquid", selected = :chip_organic, clickable = true, color = "blue-1", textcolor = "blue-10"),
            badge("3 running", color = "orange"),
        ]),
        expansionitem(label = "Advanced settings", icon = "tune",
                      caption = "Solver options, limits and weights", class = "q-mb-md", dense = true,
                      [card(flat = true, bordered = true, [card_section([
                          textarea("Notes", :notes, autogrow = true, dense = true, outlined = true),
                          # CSS grid as in NiceGUI's ui.grid
                          htmldiv(class = "m0-grid", style = "display:grid; grid-template-columns: 3rem 7rem 5rem; gap: 0 1rem; margin-top: 8px",
                                  ["C", "50.0", "%", "H", "6.0", "%", "O", "44.0", "%"]),
                      ])])]),
        row(class = FORM_ROW, [
            btn("Notify", @click(:notify_btn), nocaps = true, outline = true),
            btn("Check (spinner notification)", @click(:check_btn), nocaps = true, outline = true),
            btn("Open dialog", @click("show_dialog = true"), nocaps = true, color = "primary",
                [tooltip("Opens a q-dialog")]),
        ]),
        dialog(:show_dialog, persistent = true, [
            card(style = "min-width: 22rem", [
                card_section(htmldiv(class = "text-h6", "Remove species?")),
                card_section("This removes methane from your database."),
                card_actions(align = "right", [
                    btn("Cancel", flat = true, nocaps = true, @click("show_dialog = false")),
                    btn("Remove", color = "negative", nocaps = true, @click(:dialog_ok)),
                ]),
            ]),
        ]),
        htmldiv(class = "text-subtitle2 q-mt-lg", "Long content below, to check the header stays fixed:"),
        [htmldiv(class = "q-py-sm", "Filler line $(i)") for i in 1:40]...,
    ]
end

function queue_tab()
    [
        quasar(:linear__progress, var":value" = "progress", size = "10px", class = "q-mb-md"),
        quasar(:list, bordered = true, separator = true, dense = true, [
            item([itemsection("run 1: fig78 380 °C"), itemsection(side = true, [badge("running", color = "orange")])]),
            item([itemsection("run 2: fig79 400 °C"), itemsection(side = true, [badge("finished", color = "green")])]),
        ]),
        htmldiv(class = "text-subtitle2 q-mt-md", "Log (ui.log equivalent):"),
        quasar(:scroll__area, style = "height: 120px; border: 1px solid #ccc", [
            Html.pre("{{ log_text }}", class = "q-ma-sm text-caption"),
        ]),
    ]
end

function results_tab()
    [
        row(class = FORM_ROW, [
            btntoggle(:temp_unit, options = :unit_options, dense = true, nocaps = true),
            span("Unit switch also redraws the plot and the preview image."),
        ]),
        plot(:plot_data, layout = :plot_layout, class = "q-mb-md"),
        table(:results, dense = true, flat = true, bordered = true, title = "Species",
              pagination = DataTablePagination(rows_per_page = 0)),
        htmldiv(class = "text-subtitle2 q-mt-md", "Save-plot preview (q-img):"),
        imageview(src = :preview_src, style = "max-width: 320px", class = "q-mb-md"),
    ]
end

function files_tab()
    [
        htmldiv(class = "text-subtitle1", "1. Quasar uploader (browser picker, file content posted to the server)"),
        uploader(label = "Choose file", autoupload = true, maxfiles = 1, class = "q-mb-sm", ref = "uploader"),
        htmldiv("{{ upload_info }}", class = "q-mb-lg"),
        htmldiv(class = "text-subtitle1", "2. Electron native dialogs (called from Julia via Electron.jl)"),
        row(class = FORM_ROW, [
            btn("Open file… (native)", @click(:native_open_btn), nocaps = true, outline = true),
            btn("Save… (native)", @click(:native_save_btn), nocaps = true, outline = true),
        ]),
        htmldiv("{{ chosen_path }}", class = "q-mb-lg"),
        htmldiv(class = "text-subtitle1", "3. Browser download (Blob + a[download])"),
        btn("Download CSV (browser)", @click("browserDownload()"), nocaps = true, outline = true),
    ]
end

function ui()
    [
        quasar(:layout, view = "hHh lpR fFf", [
            quasar(:header, elevated = true, class = "bg-primary text-white", [
                row(class = "items-center q-px-md q-pt-sm", [
                    htmldiv(class = "text-h6", "HTG equilibrium (M0 spike)"),
                    space(),
                    htmldiv(class = "text-body2", "{{ status }}"),
                ]),
                tabgroup(:tab, align = "left", [
                    tab(name = "case", icon = "edit_note", label = "Case", nocaps = true),
                    tab(name = "queue", icon = "queue", label = "Queue", nocaps = true),
                    tab(name = "results", icon = "insights", label = "Results", nocaps = true),
                    tab(name = "files", icon = "folder_open", label = "Files", nocaps = true),
                ]),
            ]),
            quasar(:page__container, [
                quasar(:page, class = "q-pa-md", [
                    tabpanels(:tab, animated = false, [
                        tabpanel(name = "case", case_tab()),
                        tabpanel(name = "queue", queue_tab()),
                        tabpanel(name = "results", results_tab()),
                        tabpanel(name = "files", files_tab()),
                    ]),
                ]),
            ]),
        ]),
    ]
end

# Stipple's default layout loads the Material Icons font from fonts.googleapis.com,
# so icons are missing offline. Serve a local copy instead when one is available
# (M0_ICON_FONT, else the copy that NiceGUI ships in the Python app's venv).
function icon_font()
    p = get(ENV, "M0_ICON_FONT", "")
    isfile(p) && return p
    dir = @__DIR__
    while dirname(dir) != dir   # walk up to the repository that has the Python venv
        lib = joinpath(dir, ".venv", "lib")
        if isdir(lib)
            for py in readdir(lib; join = true)
                f = joinpath(py, "site-packages", "nicegui", "static", "fonts", "0c19a63c7528cc1a.woff2")
                isfile(f) && return f
            end
        end
        dir = dirname(dir)
    end
    nothing
end

const ICON_CSS = """
<style>
@font-face { font-family: 'Material Icons'; font-style: normal; font-weight: 400;
  src: url(/m0/material-icons.woff2) format('woff2'); }
.material-icons { font-family: 'Material Icons'; font-weight: normal; font-style: normal; font-size: 24px;
  line-height: 1; letter-spacing: normal; text-transform: none; display: inline-block; white-space: nowrap;
  word-wrap: normal; direction: ltr; -webkit-font-feature-settings: 'liga'; -webkit-font-smoothing: antialiased; }
</style>"""

function register()
    StippleUI.enable_uploads()
    font = icon_font()
    if font === nothing
        @warn "No local Material Icons font; icons load from fonts.googleapis.com (need internet)"
        @page("/", ui, layout = DEFAULT_LAYOUT(title = "HTG equilibrium (M0)", footer = ""))
    else
        bytes = read(font)
        route("/m0/material-icons.woff2") do
            Genie.Router.HTTP.Response(200, ["Content-Type" => "font/woff2"], bytes)
        end
        # a Vector: DEFAULT_LAYOUT joins head_content, so a String would be split into characters
        head = [ICON_CSS, Stipple.Layout.stipplecore_css()...]
        @page("/", ui, layout = DEFAULT_LAYOUT(title = "HTG equilibrium (M0)", footer = "",
                                                core_theme = false, head_content = head))
    end
end

