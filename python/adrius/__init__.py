# Copyright (c) 2025 InGenifold Research LLC. MIT License.
"""Adrius: Diophantine approximation library (Python bindings, EigenBackend)."""

from ._adrius import (
    # Functions
    illl,
    lll_reduce,
    cf_expansion,
    cf_convergents,
    best_rational,
    simultaneous_approx,
    # Result types
    ILLLResult,
    LLLResult,
    SimultApproxResult,
    Rational,
    # Exceptions
    Error,
    DomainError,
    ConvergenceError,
    PrecisionError,
)

__all__ = [
    "illl",
    "lll_reduce",
    "cf_expansion",
    "cf_convergents",
    "best_rational",
    "simultaneous_approx",
    "ILLLResult",
    "LLLResult",
    "SimultApproxResult",
    "Rational",
    "Error",
    "DomainError",
    "ConvergenceError",
    "PrecisionError",
]
