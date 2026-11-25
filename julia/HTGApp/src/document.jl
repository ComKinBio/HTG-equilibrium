# Case files being edited (Python: gui/document.py CaseDocument).
#
# The form edits `data`, plain nested OrderedDicts with the structure of a case
# file. Owner decision (PR #8): saving writes a fresh file with the core's
# `dump_case`; comments and the layout of the original file are not kept (Python
# patched the original text with tomlkit, `update_document`, which is not ported).

const LEGACY_NOTICE = "This case uses the old format. Saving will update it to the new format."
const SALTS_NOTICE = "This case has added salts. The form can't edit them; they are kept when you save."

const _Data = OrderedDict{String,Any}

"""
    CaseDocument(text; path=nothing, molar_mass=nothing)
    CaseDocument.open(path; molar_mass=nothing)     # open_document(path)
    new_document(name="new_case")

A case file being edited: `data` (nested OrderedDicts) and the file's `path`.
Old (format-1) files are converted to format 2 when opened (`data` holds the
new keys). Flags for the form:

- `converted_from_legacy`: true until the converted file is saved;
- `legacy_blocker`: why an old file was left in the old format (or `nothing`);
- `legacy_notes`: what was changed when it was opened (a dropped solver and its
  settings, owner decision PR #10; the core's `case.legacy_notes` texts);
- `has_salts(doc)`: the feed has [feed.salts] (kept, not editable).

`snapshot` is the cleaned data as last opened or saved (`dirty(doc)`); it lives
in the document so that a reloaded page still knows about unsaved changes.
`molar_mass(name) -> g/mol` is needed to convert an old species-moles feed.
"""
mutable struct CaseDocument
    data::_Data
    path::Union{Nothing,String}
    converted_from_legacy::Bool
    legacy_blocker::Union{Nothing,String}
    legacy_notes::Vector{String}
    snapshot::Any
end

function CaseDocument(text::AbstractString; path=nothing, molar_mass=nothing)
    data = _Data(parse_toml(text))
    converted, blocker = false, nothing
    if !is_v2(data)
        try
            data = _Data(legacy_to_v2(data, molar_mass))
            converted = true
        catch exc
            exc isa Union{HTGEquilibrium.ValueError,ArgumentError} || rethrow()   # LegacyConversionError, or a broken old file
            blocker = sprint(showerror, exc)
        end
    end
    notes = drop_removed_solver!(data)
    doc = CaseDocument(data, path === nothing ? nothing : abspath(path), converted, blocker, notes, nothing)
    mark_saved!(doc)
    return doc
end

"""Open the case file at `path`."""
open_document(path; molar_mass=nothing) = CaseDocument(read(path, String); path, molar_mass)

"""A new case from the template, named `name`."""
function new_document(name="new_case")
    doc = CaseDocument(read(HTGEquilibrium.CASE_TEMPLATE, String))
    doc.data["name"] = name
    mark_saved!(doc)
    return doc
end

"""
Owner decision (PR #10): the solvers "slsqp-fd" and "trust-constr" and the
[solver] keys fd_ftol and fd_eps are gone. Remove them from `data` (the form
then shows IPOPT, the default) and return the notes the core gives for them.
"""
function drop_removed_solver!(data)
    notes = String[]
    solver = get(data, "solver", nothing)
    solver isa AbstractDict || return notes
    method = get(solver, "method", nothing)
    if method isa AbstractString && method in HTGEquilibrium.LEGACY_SOLVERS
        delete!(solver, "method")
        push!(notes, HTGEquilibrium.legacy_solver_note(method))
    end
    dropped = [k for k in HTGEquilibrium.LEGACY_SOLVER_KEYS if haskey(solver, k)]
    if !isempty(dropped)
        foreach(k -> delete!(solver, k), dropped)
        push!(notes, HTGEquilibrium.legacy_solver_keys_note(dropped))
    end
    return notes
end

"""True when the feed has added salts ([feed.salts])."""
function has_salts(doc::CaseDocument)
    feed = get(doc.data, "feed", nothing)
    return feed isa AbstractDict && HTGEquilibrium.py_truth(get(feed, "salts", nothing))
end

"""Notices for the Case tab, in the order to show them."""
function notices(doc::CaseDocument)
    out = String[]
    doc.converted_from_legacy && push!(out, LEGACY_NOTICE)
    doc.legacy_blocker === nothing || push!(out, doc.legacy_blocker)
    append!(out, doc.legacy_notes)
    has_salts(doc) && push!(out, SALTS_NOTICE)
    return out
end

"""Folder that relative paths in the case refer to."""
base_dir(doc::CaseDocument) = doc.path === nothing ? pwd() : dirname(resolve_path(doc.path))

"""The Case described by `data` (throws with a readable message)."""
document_case(doc::CaseDocument) = case_from_data(doc.data, doc.path)

_stem(path) = splitext(basename(path))[1]

"""Build a Case from case-file data; the comparison file is made absolute, relative to `path`."""
function case_from_data(data, path=nothing)
    c = case_from_dict(clean(data); name=path === nothing ? nothing : _stem(path))
    if HTGEquilibrium.py_truth(c.comparison)
        base = path === nothing ? pwd() : dirname(abspath(path))
        file = abspath(joinpath(base, c.comparison["file"]))
        c.comparison["file"] = ispath(file) ? realpath(file) : file
    end
    return c
end

"""
The TOML text of `data`: a fresh file in the core's layout (`dump_case`), the
comparison file as written in `data` (relative). Comments are not kept (PR #8).
"""
function document_text(doc::CaseDocument, path=doc.path)
    return dump_case(case_from_dict(clean(doc.data); name=path === nothing ? nothing : _stem(path)))
end

"""
    save_document(doc, path=nothing) -> path

Write the case to `path` (default: where it came from). Refuses a case that
cannot be loaded. A comparison file stays reachable from a new folder.
"""
function save_document(doc::CaseDocument, path=nothing)
    path === nothing && doc.path === nothing && throw(FieldError("Choose a file name."))
    path = resolve_path(something(path, doc.path))
    document_case(doc)   # refuse to save something that cannot be loaded
    comparison = get(doc.data, "comparison", nothing)
    if comparison isa AbstractDict && HTGEquilibrium.py_truth(get(comparison, "file", nothing)) &&
       dirname(path) != base_dir(doc)
        target = normpath(joinpath(base_dir(doc), comparison["file"]))
        comparison["file"] = replace(relpath(target, dirname(path)), '\\' => '/')
    end
    content = document_text(doc, path)
    mkpath(dirname(path))
    write(path, content)
    doc.path = path
    doc.converted_from_legacy = false
    empty!(doc.legacy_notes)
    mark_saved!(doc)
    return path
end

"""Remember the current data as saved (for `dirty`)."""
mark_saved!(doc::CaseDocument) = (doc.snapshot = plain(clean(doc.data)); doc)

"""True when `data` differs from the file as last opened or saved."""
dirty(doc::CaseDocument) = plain(clean(doc.data)) != doc.snapshot

"""Copy without `nothing` values, NaN numbers and empty tables (lists are kept, also empty ones)."""
function clean(data::AbstractDict)
    out = _Data()
    for (key, value) in data
        if value isa AbstractDict
            value = clean(value)
            isempty(value) && continue
        elseif value === nothing || (value isa AbstractFloat && isnan(value))
            continue
        end
        out[string(key)] = deepcopy(value)
    end
    return out
end

"""Nested Dicts (key order does not matter) and Vectors, for comparing data."""
plain(x::AbstractDict) = Dict{String,Any}(string(k) => plain(v) for (k, v) in x)
plain(x::AbstractVector) = Any[plain(v) for v in x]
plain(x) = x
