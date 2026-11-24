# M0 spike (throw-away): a Genie/Stipple page with the Quasar components the Phase 1
# app uses, shown in an Electron.jl window. A package, so a sysimage can contain it
# the way it would contain the real app.
module M0App

using Stipple, Stipple.ReactiveTools
using StippleUI
# Genie.Renderer.Html exports HTML elements with the same names; take the Quasar ones.
using StippleUI.Selects: select
using StippleUI.Tables: table
using StippleUI: textarea, radio, checkbox, icon, badge, chip, item, itemsection, space, tooltip, dialog,
                 card, card_section, card_actions, uploader, imageview
using StipplePlotly
using DataFrames: DataFrame
using Base64: base64encode
using Sockets
using Electron
using Genie
using JSON

# launch time (Unix seconds) that all [M0] time stamps are relative to
const T0 = Ref(time())
stamp(msg) = (println("[M0] $(lpad(round(time() - T0[]; digits = 2), 6)) s  ", msg); flush(stdout))

include("ui.jl")        # the page: model, handlers, view
include("shell.jl")     # server + Electron window + clean shut-down, memory report
include("selftest.jl")  # drives the page and checks each component

end # module
