"""Can the app open in its own window (NiceGUI native mode, pywebview)?

Checked before ``ui.run``: NiceGUI opens the window in a child process and stops
the whole app when that process dies, so a window that cannot start must be
detected first. On Linux pywebview uses Qt (installed with ``pywebview[qt]``);
the check starts a Qt application with the web engine in a separate Python
process, because a missing display or system library makes Qt abort the process.
On Windows and macOS pywebview uses the system web view and no check is made.
"""
import importlib.util
import os
import re
import subprocess
import sys

_QT_PROBE = (
    'import sys\n'
    'import webview.platforms.qt\n'
    'from qtpy.QtWidgets import QApplication\n'
    'QApplication(sys.argv[:1])\n'
)
_LINUX_PACKAGES = 'libegl1, libxcb-cursor0, libxkbcommon-x11-0 and libnss3'


def use_qt():
    """Make pywebview use Qt on Linux (unless PYWEBVIEW_GUI says otherwise)."""
    return sys.platform.startswith('linux') and 'PYWEBVIEW_GUI' not in os.environ


def native_window_problem(timeout=30):
    """None if a native window can open; otherwise one line saying what is missing and how to fix it."""
    if importlib.util.find_spec('webview') is None:
        return 'pywebview is not installed. Run "uv sync" to install it.'
    if not use_qt():
        return None
    try:
        result = subprocess.run([sys.executable, '-c', _QT_PROBE], capture_output=True, text=True,
                                timeout=timeout)
    except subprocess.TimeoutExpired:
        return f'Qt did not start within {timeout} s.'
    except OSError as exc:
        return f'Qt could not be checked ({exc}).'
    if result.returncode == 0:
        return None
    return explain_qt_failure(result.stderr)


def explain_qt_failure(stderr):
    """One line for the output of a failed Qt start."""
    text = stderr or ''
    module = re.search(r"No module named '([\w.]+)'", text)
    if module:
        return f'the Python package {module.group(1)} is missing. Run "uv sync" to install it.'
    if ('could not connect to display' in text or 'Failed to create wl_display' in text
            or not (os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY'))):
        return 'no display was found. Run htg gui from a desktop session.'
    library = re.search(r'(lib[\w+-]+?(?:\.so[\d.]*)?)(?::| is needed| cannot)', text)
    if library:
        return (f'the system library {library.group(1)} is missing. Install it with your package manager '
                f'(on Debian or Ubuntu the usual Qt packages are {_LINUX_PACKAGES}).')
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    reason = lines[-1][:150] if lines else 'Qt stopped without a message'
    return f'Qt could not start ({reason}).'


def native_or_browser(browser=False):
    """True to open the native window. Prints one line and returns False when it cannot open."""
    if browser:
        return False
    problem = native_window_problem()
    if problem is None:
        return True
    print(f'Cannot open a separate window: {problem} Opening the app in the browser instead '
          '(htg gui --browser does this directly).')
    return False
