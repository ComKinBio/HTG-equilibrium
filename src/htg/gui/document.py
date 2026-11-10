"""Case files as editable TOML documents.

The GUI edits a plain dict with the structure of a case file. Saving merges it
into the original document with tomlkit, so comments and layout of values the
user did not change are kept.
"""
import copy
import math
import os
from pathlib import Path

import tomlkit
from tomlkit.items import InlineTable, Table

from ..case import case_from_dict

TEMPLATE = Path(__file__).resolve().parents[1] / 'templates' / 'case_template.toml'


class CaseDocument:
    """A case file being edited: ``data`` (plain dict) plus the original document."""

    def __init__(self, text, path=None):
        self.doc = tomlkit.parse(text)
        self.data = self.doc.unwrap()
        self.path = Path(path) if path else None

    @classmethod
    def open(cls, path):
        return cls(Path(path).read_text(encoding='utf-8'), path)

    @classmethod
    def new(cls, name='new_case'):
        document = cls(TEMPLATE.read_text(encoding='utf-8'))
        document.data['name'] = name
        return document

    def base_dir(self):
        """Folder that relative paths in the case refer to."""
        return self.path.resolve().parent if self.path else Path.cwd()

    def case(self):
        """The Case described by ``data`` (raises ValueError with a readable message)."""
        return case_from_data(self.data, self.path)

    def text(self):
        """TOML text of ``data``, keeping comments of the original document."""
        update_document(self.doc, clean(self.data))
        return tomlkit.dumps(self.doc)

    def save(self, path=None):
        if path is None and self.path is None:
            raise ValueError('Choose a file name.')
        path = Path(path or self.path).resolve()
        self.case()  # refuse to save something that cannot be loaded
        comparison = self.data.get('comparison') or {}
        if comparison.get('file') and path.parent != self.base_dir():
            # keep the comparison file reachable from the new location
            target = (self.base_dir() / comparison['file']).resolve()
            comparison['file'] = os.path.relpath(target, path.parent).replace(os.sep, '/')
        text = self.text()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding='utf-8')
        self.path = path
        return path


def case_from_data(data, path=None):
    """Build a Case from case-file data; the comparison file is relative to ``path``."""
    path = Path(path) if path else None
    case = case_from_dict(clean(data), name=path.stem if path else None)
    if case.comparison:
        base = path.parent if path else Path.cwd()
        case.comparison['file'] = str((base / case.comparison['file']).resolve())
    return case


def clean(data):
    """Copy without None values, NaN numbers and empty tables (except [phases] lists)."""
    out = {}
    for key, value in data.items():
        if isinstance(value, dict):
            value = clean(value)
            if not value:
                continue
        elif value is None or (isinstance(value, float) and math.isnan(value)):
            continue
        out[key] = copy.deepcopy(value)
    return out


ROOT_INLINE = {'atomic_weights'}


def update_document(container, data, root=None):
    """Make ``container`` (a tomlkit document or table) hold ``data``.

    Unchanged values keep their formatting and comments; changed values keep
    their trailing comment; inline tables stay inline. New sections at the top
    level become [tables]; new tables of plain values below it become inline.
    """
    root = isinstance(container, tomlkit.TOMLDocument) if root is None else root
    for key in list(container.keys()):
        if key not in data:
            del container[key]
    for key, value in data.items():
        old = container.get(key)
        if isinstance(value, dict):
            if isinstance(old, Table):
                update_document(old, value, root=False)
                continue
            inline = _flat(value) and (not root or key in ROOT_INLINE)
            if isinstance(old, InlineTable) or (old is None and inline):
                if old is not None and old.unwrap() == value:
                    continue
                container[key] = _inline(value)
            else:
                table = tomlkit.table(is_super_table=not _flat(value) and all(
                    isinstance(v, dict) for v in value.values()))
                update_document(table, value, root=False)
                if not table.is_super_table():
                    table.add(tomlkit.nl())
                container[key] = table
            continue
        if old is not None and _unwrap(old) == value:
            continue
        item = tomlkit.item(value)
        if isinstance(value, list) and len(value) > 6:
            item.multiline(True)
        comment = getattr(getattr(old, 'trivia', None), 'comment', '')
        if comment:
            item.comment(comment.lstrip('#').strip())
        container[key] = item


def _flat(value):
    return all(not isinstance(v, dict) for v in value.values())


def _inline(value):
    table = tomlkit.inline_table()
    for k, v in value.items():
        table[k] = _inline(v) if isinstance(v, dict) else v
    return table


def _unwrap(item):
    return item.unwrap() if hasattr(item, 'unwrap') else item
