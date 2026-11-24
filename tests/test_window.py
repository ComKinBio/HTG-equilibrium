"""Choosing between the app's own window and the browser (no window is opened here)."""
import os

from htg.gui import window

NO_DISPLAY = '''Failed to create wl_display (No such file or directory)
qt.qpa.xcb: could not connect to display
qt.qpa.plugin: From 6.5.0, xcb-cursor0 or libxcb-cursor0 is needed to load the Qt xcb platform plugin.
'''


def test_explain_qt_failures(monkeypatch):
    monkeypatch.setenv('DISPLAY', ':0')
    assert 'no display' in window.explain_qt_failure(NO_DISPLAY)
    message = window.explain_qt_failure(
        'qt.qpa.plugin: From 6.5.0, xcb-cursor0 or libxcb-cursor0 is needed to load the Qt xcb platform plugin.')
    assert 'libxcb-cursor0 is missing' in message and 'package manager' in message
    message = window.explain_qt_failure('ImportError: libEGL.so.1: cannot open shared object file')
    assert 'libEGL.so.1 is missing' in message
    message = window.explain_qt_failure("ModuleNotFoundError: No module named 'PyQt6'")
    assert 'PyQt6' in message and 'uv sync' in message
    assert window.explain_qt_failure('Segmentation fault\n') == 'Qt could not start (Segmentation fault).'
    monkeypatch.delenv('DISPLAY')
    monkeypatch.delenv('WAYLAND_DISPLAY', raising=False)
    assert 'no display' in window.explain_qt_failure('anything')


def test_browser_flag_skips_the_check(monkeypatch, capsys):
    monkeypatch.setattr(window, 'native_window_problem', lambda: 1 / 0)
    assert window.native_or_browser(browser=True) is False
    assert capsys.readouterr().out == ''


def test_fallback_prints_one_line(monkeypatch, capsys):
    monkeypatch.setattr(window, 'native_window_problem', lambda: 'no display was found. Run htg gui from a desktop session.')
    assert window.native_or_browser() is False
    out = capsys.readouterr().out
    assert out.count('\n') == 1
    assert 'no display was found' in out and 'browser' in out
    monkeypatch.setattr(window, 'native_window_problem', lambda: None)
    assert window.native_or_browser() is True


def test_quiet_qt_hides_wayland_text_input_warnings(monkeypatch):
    monkeypatch.delenv('QT_LOGGING_RULES', raising=False)
    window.quiet_qt()
    assert os.environ['QT_LOGGING_RULES'] == 'qt.qpa.wayland.textinput=false'
    monkeypatch.setenv('QT_LOGGING_RULES', '*.debug=true')
    window.quiet_qt()
    assert os.environ['QT_LOGGING_RULES'] == '*.debug=true'
