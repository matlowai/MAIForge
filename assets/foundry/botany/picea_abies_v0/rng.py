"""Keyed deterministic draws.

A draw is a pure function of the fixture seed and a key built from stable node
ids. It does not depend on traversal order, on how many draws came before it, or
on which detail tier asked, so near, mid and far compilers reading the same graph
see identical values.
"""
from hashlib import blake2b
import math

_MAX = (1 << 64) - 1


def _bits(seed, key):
    return int.from_bytes(blake2b(f'{seed}|{key}'.encode(), digest_size=8).digest(), 'big')


def unit(seed, key):
    """Uniform in [0,1) for this exact key."""
    return _bits(seed, key) / (_MAX + 1)


def between(seed, key, low, high):
    return low + (high - low) * unit(seed, key)


def bernoulli(seed, key, probability):
    return unit(seed, key) < probability


def pick(seed, key, options, masses):
    """Pick an option by its mass. Masses need not sum to exactly one."""
    total = math.fsum(masses)
    if total <= 0:
        raise ValueError('pick needs positive mass')
    target = unit(seed, key) * total
    cumulative = 0.0
    for option, mass in zip(options, masses):
        cumulative += mass
        if target < cumulative:
            return option
    return options[-1]
