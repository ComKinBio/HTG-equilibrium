# The page (Python: gui/app.py build_page): a fixed header with the title, the
# status line and the five tabs; one panel per tab below it.
#
# Offline: Stipple's default layout loads the Material Icons font and its theme
# fonts from fonts.googleapis.com. The page is registered with Stipple's core
# theme and theme files turned off (`core_theme = false`, `include_themes =
# false`): no Google font, no `.theme` file written into the working
# directory. The icon font comes from the `material_icons` artifact
# (Artifacts.toml, the material-icons npm package) and is served by this server.

"""The app's shared state; set by `run_app` before the server starts."""
const STATE = Ref{Union{Nothing,State}}(nothing)
"""Called once per client when its page is connected (`shell.jl` uses it to time the start-up)."""
const ON_READY = Ref{Any}(nothing)

"""The tabs: (name, label, Material icon), as in Python's build_page."""
const TABS = (("case", "Case", "edit_note"), ("queue", "Queue", "queue"), ("results", "Results", "insights"),
              ("compare", "Compare", "compare_arrows"), ("species", "Species", "science"))

"""Text of a tab panel until its screen is ported."""
placeholder_text(label) = "The $label tab is not in the Julia version yet. For now, use the Python app (uv run htg gui)."

const STATUS_PERIOD = 1.0   # s; Python: ui.timer(1.0, ...)

@app begin
    @in tab = "case"
    @out status = "Loading the species data"
    # The Case tab (case_tab.jl): the form's view, its field texts, the species lists.
    @out case_view = Dict{String,Any}("ready" => false)
    @out case_texts = Dict{String,Any}()
    @out case_options = Dict{String,Any}("phases" => Dict{String,Any}(), "feed" => Any[])

    @onchange isready begin
        isready || return
        ON_READY[] === nothing || ON_READY[](__model__)
        case_start!(__model__)
        # The status line: the server pushes a new text every second (Python: ui.timer).
        # The same loop tells the Case tab when the species data have been loaded.
        @async begin
            loaded = STATE[] === nothing ? nothing : STATE[].model
            while isready && Stipple.isconnected(__model__) !== false
                text = STATE[] === nothing ? "" : status_text(STATE[])
                text == status || (status = text)
                if STATE[] !== nothing && STATE[].model !== loaded
                    loaded = STATE[].model
                    case_enqueue!(__model__, :species)
                end
                sleep(STATUS_PERIOD)
            end
            case_stop!(__model__)
        end
    end
end

function ui()
    [
        quasar(:layout, view = "hHh lpR fFf", [
            # The tabs are in the (fixed) header so they stay visible when a long tab is scrolled.
            quasar(:header, elevated = true, class = "bg-primary text-white column items-stretch", [
                row(class = "items-center q-px-md q-pt-sm", [
                    htmldiv(class = "text-h6", "HTG equilibrium"),
                    space(),
                    htmldiv(class = "text-body2 htg-status", "{{ status }}"),
                ]),
                tabgroup(:tab, class = "full-width", [tab(name = name, label = label, icon = icon)
                                                     for (name, label, icon) in TABS]),
            ]),
            quasar(:page__container, [
                quasar(:page, [
                    tabpanels(:tab, animated = false, class = "htg-panels",
                              [tabpanel(name = name, [name == "case" ? case_tab_html() :
                                                      htmldiv(class = "text-body1 htg-placeholder", placeholder_text(label))])
                               for (name, label, _) in TABS]),
                ]),
            ]),
        ]),
    ]
end

"""The Material Icons font (woff2) from the `material_icons` artifact."""
icon_font_path() = joinpath(artifact"material_icons", "package", "iconfont", "material-icons.woff2")

const FONT_ROUTE = "/htg/material-icons.woff2"

# Material Icons (as Google's CSS defines it), NiceGUI's colours (NiceGUI's
# default primary is #5898d4, not Quasar's) and the max-width panel container
# (Python: w-full max-w-screen-2xl mx-auto, 1536 px).
const PAGE_CSS = """
<style>
@font-face { font-family: 'Material Icons'; font-style: normal; font-weight: 400; font-display: block;
  src: url($FONT_ROUTE) format('woff2'); }
.material-icons { font-family: 'Material Icons'; font-weight: normal; font-style: normal; font-size: 24px;
  line-height: 1; letter-spacing: normal; text-transform: none; display: inline-block; white-space: nowrap;
  word-wrap: normal; direction: ltr; -webkit-font-feature-settings: 'liga'; font-feature-settings: 'liga';
  -webkit-font-smoothing: antialiased; }
html:root { --q-primary: #5898d4; --q-secondary: #26a69a; --q-accent: #9c27b0; --q-dark: #1d1d1d;
  --q-dark-page: #121212; --q-positive: #21ba45; --q-negative: #c10015; --q-info: #31ccec; --q-warning: #f2c037; }
.htg-panels { width: 100%; max-width: 1536px; margin-left: auto; margin-right: auto; }
</style>"""

"""Register the font route and the page at "/"."""
function register_page()
    Genie.Secrets.secret_token!()   # a random one per session (Genie warns without one)
    font = read(icon_font_path())
    route(FONT_ROUTE) do
        Genie.Router.HTTP.Response(200, ["Content-Type" => "font/woff2", "Cache-Control" => "max-age=86400"], font)
    end
    # The stylesheets Stipple's theme would add (Quasar's CSS from StippleUI), except its core
    # theme and its default theme file (theme-default-light.css @imports a Google font and
    # restyles Quasar's colours). A Vector: DEFAULT_LAYOUT joins head_content.
    styles = [f() for f in Stipple.Layout.THEMES[] if f !== Stipple.Layout.coretheme]
    links = String[s for x in styles for s in (x isa AbstractString ? [x] : x)]
    head = vcat(filter(s -> !occursin("/css/themes/", s) && !occursin("://", s), links), [PAGE_CSS, CASE_CSS])
    @page("/", ui, core_theme = false,
          layout = DEFAULT_LAYOUT(title = "HTG equilibrium", footer = "", core_theme = false,
                                  include_themes = false, head_content = head))
    return nothing
end
