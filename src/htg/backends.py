"""Lazy access to the heavy external libraries.

pyclapeyron (Clapeyron.jl through juliacall) starts Julia on first import, which
takes from a few seconds to minutes. Importing htg stays fast; Julia starts the
first time an equation of state is needed.
"""
import functools


@functools.cache
def clapeyron():
    import pyclapeyron
    return pyclapeyron


@functools.cache
def julia():
    from juliacall import Main
    return Main


@functools.cache
def chnosz():
    import pychnosz
    return pychnosz


def julia_string_vector(values):
    escaped = ['"' + str(v).replace('\\', '\\\\').replace('"', '\\"') + '"' for v in values]
    return julia().seval('String[' + ', '.join(escaped) + ']')


def julia_int_vector(values):
    return julia().seval('Int64[' + ', '.join(str(int(v)) for v in values) + ']')
