// Copyright (c) 2025 InGenifold Research LLC. MIT License.
//
// pybind11 bindings for Adrius — EigenBackend (double precision).
//
// Include order matters: pybind11/pybind11.h first, then Eigen (via adrius),
// then pybind11/eigen.h so pybind11 sees the Eigen types.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <adrius/adrius.hpp>

#include <pybind11/eigen.h>

#include <span>
#include <sstream>
#include <string>
#include <vector>

namespace py = pybind11;
using B = adrius::EigenBackend;

// Convert a 1-D py::array_t<double> (or any array-like) to std::span<const double>.
// The forcecast flag coerces Python lists and non-contiguous arrays to a
// temporary contiguous buffer; span_alpha borrows from that buffer, which stays
// alive for the duration of the function call.
static std::span<const double>
to_span(const py::array_t<double, py::array::c_style | py::array::forcecast>& arr)
{
    const auto buf = arr.request();
    if (buf.ndim != 1)
        throw adrius::DomainError("alpha must be a 1-D array");
    return {static_cast<const double*>(buf.ptr), static_cast<std::size_t>(buf.shape[0])};
}

PYBIND11_MODULE(_adrius, m)
{
    m.doc() = R"doc(
Adrius: Diophantine approximation library (Python bindings, EigenBackend).

Functions
---------
illl                 Iterated LLL — sequence of simultaneous approximations
lll_reduce           LLL lattice basis reduction
cf_expansion         Continued-fraction partial quotients (eager)
cf_convergents       CF convergents p_k/q_k
best_rational        Best rational p/q with q ≤ max_denominator
simultaneous_approx  One-shot simultaneous Diophantine approximation

Types
-----
ILLLResult, LLLResult, SimultApproxResult, Rational

Exceptions
----------
Error (base), DomainError, ConvergenceError, PrecisionError
)doc";

    // ── Exceptions ────────────────────────────────────────────────────────
    // Register base class first so derived exceptions inherit from it in Python.
    static auto exc_base =
        py::register_exception<adrius::Error>(m, "Error", PyExc_RuntimeError);
    py::register_exception<adrius::DomainError>(m,
        "DomainError", exc_base.ptr());
    py::register_exception<adrius::ConvergenceError>(m,
        "ConvergenceError", exc_base.ptr());
    py::register_exception<adrius::PrecisionError>(m,
        "PrecisionError", exc_base.ptr());

    // ── Rational<int64_t> ─────────────────────────────────────────────────
    using Rat = adrius::Rational<std::int64_t>;
    py::class_<Rat>(m, "Rational",
        R"doc(
Exact rational number stored in lowest terms (denominator > 0).

Attributes
----------
numerator   : int
denominator : int

Methods
-------
to_float()  -> float
is_integer() -> bool
)doc")
        .def(py::init<std::int64_t, std::int64_t>(), py::arg("numerator"), py::arg("denominator"))
        .def(py::init<std::int64_t>(), py::arg("n"))
        .def_property_readonly("numerator",   &Rat::numerator)
        .def_property_readonly("denominator", &Rat::denominator)
        .def("to_float",    &Rat::to_double)
        .def("is_integer",  &Rat::is_integer)
        .def("is_zero",     &Rat::is_zero)
        .def("__float__",   &Rat::to_double)
        .def("__int__",     [](const Rat& r) {
            if (!r.is_integer())
                throw py::value_error("Rational is not an integer");
            return r.numerator();
        })
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * py::self)
        .def(py::self / py::self)
        .def(-py::self)
        .def(py::self == py::self)
        .def(py::self <  py::self)
        .def(py::self <= py::self)
        .def(py::self >  py::self)
        .def(py::self >= py::self)
        .def("__repr__", [](const Rat& r) {
            std::ostringstream oss;
            oss << "Rational(" << r.numerator() << ", " << r.denominator() << ")";
            return oss.str();
        })
        .def("__str__", [](const Rat& r) {
            std::ostringstream oss;
            oss << r;
            return oss.str();
        });

    // ── ILLLResult ────────────────────────────────────────────────────────
    py::class_<adrius::ILLLResult<B>>(m, "ILLLResult",
        R"doc(
Result of illl().

Attributes
----------
relations  : list of list of int — [[q, p1, …, pn], …], one per iteration
quality    : list of float       — max_i |q*alpha_i - p_i| per iteration
iterations : int                 — number of relations found
)doc")
        .def_readonly("relations",  &adrius::ILLLResult<B>::relations)
        .def_readonly("quality",    &adrius::ILLLResult<B>::quality)
        .def_readonly("iterations", &adrius::ILLLResult<B>::iterations)
        .def("__repr__", [](const adrius::ILLLResult<B>& r) {
            return "<ILLLResult iterations=" + std::to_string(r.iterations) + ">";
        });

    // ── LLLResult ─────────────────────────────────────────────────────────
    py::class_<adrius::LLLResult<B>>(m, "LLLResult",
        R"doc(
Result of lll_reduce().

Attributes
----------
reduced_basis : np.ndarray (float64, shape n×n) — LLL-reduced basis (columns)
transform     : np.ndarray (int64,   shape n×n) — unimodular U: reduced = basis @ U
swap_count    : int                              — number of Lovász swaps performed
)doc")
        .def_readonly("reduced_basis", &adrius::LLLResult<B>::reduced_basis)
        .def_readonly("transform",     &adrius::LLLResult<B>::transform)
        .def_readonly("swap_count",    &adrius::LLLResult<B>::swap_count)
        .def("__repr__", [](const adrius::LLLResult<B>& r) {
            return "<LLLResult swap_count=" + std::to_string(r.swap_count) + ">";
        });

    // ── SimultApproxResult ────────────────────────────────────────────────
    py::class_<adrius::SimultApproxResult<B>>(m, "SimultApproxResult",
        R"doc(
Result of simultaneous_approx().

Attributes
----------
denominator : int         — common denominator q
numerators  : list of int — integer vector [p1, …, pn]
quality     : float       — max_i |q*alpha_i - p_i|
)doc")
        .def_readonly("denominator", &adrius::SimultApproxResult<B>::denominator)
        .def_readonly("numerators",  &adrius::SimultApproxResult<B>::numerators)
        .def_readonly("quality",     &adrius::SimultApproxResult<B>::quality)
        .def("__repr__", [](const adrius::SimultApproxResult<B>& r) {
            return "<SimultApproxResult q=" + std::to_string(r.denominator)
                 + " quality=" + std::to_string(static_cast<double>(r.quality)) + ">";
        });

    // ── illl ──────────────────────────────────────────────────────────────
    m.def("illl",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> alpha,
           std::size_t   max_iterations,
           std::int64_t  max_denominator,
           double        epsilon,
           double        quality_tol)
        {
            adrius::ILLLParams params;
            params.max_iterations  = max_iterations;
            params.max_denominator = max_denominator;
            params.epsilon         = epsilon;
            params.quality_tol     = quality_tol;
            return adrius::illl<B>(to_span(alpha), params);
        },
        py::arg("alpha"),
        py::arg("max_iterations")  = 30,
        py::arg("max_denominator") = 1'000'000LL,
        py::arg("epsilon")         = 0.5,
        py::arg("quality_tol")     = 1e-12,
        R"doc(
Iterated LLL simultaneous Diophantine approximation (Bosma & Smeets 2010).

For each iteration k, finds integers q_k, p_{k,1}, …, p_{k,n} such that
max_i |q_k * alpha_i - p_{k,i}| decreases and q_k grows.  The Dirichlet
coefficient q_k^{1/n} * quality_k is bounded by a constant depending only on n.

Parameters
----------
alpha : array-like of float
    Target real numbers (alpha_1, …, alpha_n).
max_iterations : int, default 30
    Maximum outer ILLL iterations.
max_denominator : int, default 1_000_000
    Stop when q >= max_denominator; 0 = no limit.
epsilon : float, default 0.5
    Controls the initial lattice scale c_0 = epsilon^{n+1}.
quality_tol : float, default 1e-12
    Stop when max|q*alpha_i - p_i| < quality_tol.

Returns
-------
ILLLResult

Raises
------
DomainError
    If alpha is empty or epsilon is outside (0, 1).

Examples
--------
>>> import adrius, numpy as np
>>> r = adrius.illl([np.sqrt(2), np.pi], max_denominator=10_000)
>>> for k in range(r.iterations):
...     q, *ps = r.relations[k]
...     print(q, ps, r.quality[k])
)doc");

    // ── lll_reduce ────────────────────────────────────────────────────────
    m.def("lll_reduce",
        [](Eigen::MatrixXd basis, double delta, double eta, std::size_t max_iter)
        {
            adrius::LLLParams params;
            params.delta    = delta;
            params.eta      = eta;
            params.max_iter = max_iter;
            return adrius::lll_reduce<B>(basis, params);
        },
        py::arg("basis"),
        py::arg("delta")    = 0.75,
        py::arg("eta")      = 0.51,
        py::arg("max_iter") = 0,
        R"doc(
LLL lattice basis reduction.

Parameters
----------
basis : np.ndarray, shape (n, n), dtype float
    Input basis matrix — each **column** is a basis vector (Eigen convention).
    A row-major numpy array is automatically transposed to column-major.
delta : float, default 0.75
    Lovász condition parameter delta in (0.25, 1.0).
eta : float, default 0.51
    Size-reduction threshold (theoretical minimum 0.5).
max_iter : int, default 0
    Maximum column swaps; 0 = unlimited.

Returns
-------
LLLResult
    .reduced_basis : np.ndarray (float64) — LLL-reduced basis (columns)
    .transform     : np.ndarray (int64)   — unimodular U: reduced = basis @ U
    .swap_count    : int

Raises
------
DomainError
    If the basis is linearly dependent.
ConvergenceError
    If max_iter swaps are exceeded.

Examples
--------
>>> import adrius, numpy as np
>>> basis = np.array([[1, -1, 3], [1, 0, 5], [1, 2, 6]], dtype=float)
>>> r = adrius.lll_reduce(basis)
>>> np.allclose(basis @ r.transform, r.reduced_basis)
True
)doc");

    // ── cf_expansion ─────────────────────────────────────────────────────
    m.def("cf_expansion",
        [](double x, std::size_t max_depth) {
            return adrius::cf_expansion(x, adrius::CFParams{.max_depth = max_depth});
        },
        py::arg("x"),
        py::arg("max_depth") = 64,
        R"doc(
Compute the continued-fraction partial quotients of x.

Parameters
----------
x : float
    The real number to expand.
max_depth : int, default 64
    Maximum number of partial quotients.

Returns
-------
list of int
    Partial quotients [a_0; a_1, a_2, …].

Examples
--------
>>> adrius.cf_expansion(3.14159265358979, max_depth=6)
[3, 7, 15, 1, 292, 1]
)doc");

    // ── cf_convergents ───────────────────────────────────────────────────
    m.def("cf_convergents",
        [](double x, std::size_t max_depth) {
            const auto convs = adrius::cf_convergents(
                x, adrius::CFParams{.max_depth = max_depth});
            std::vector<std::pair<std::int64_t, std::int64_t>> out;
            out.reserve(convs.size());
            for (const auto& c : convs)
                out.emplace_back(c.p, c.q);
            return out;
        },
        py::arg("x"),
        py::arg("max_depth") = 64,
        R"doc(
Compute the convergents p_k/q_k of the CF expansion of x.

Parameters
----------
x : float
    The real number.
max_depth : int, default 64
    Maximum number of convergents.

Returns
-------
list of (int, int)
    Pairs (p_k, q_k) for k = 0, 1, 2, …

Examples
--------
>>> adrius.cf_convergents(3.14159265358979, max_depth=4)
[(3, 1), (22, 7), (333, 106), (355, 113)]
)doc");

    // ── best_rational ────────────────────────────────────────────────────
    m.def("best_rational",
        [](double x, std::int64_t max_denominator) {
            return adrius::best_rational(
                x, adrius::RationalApproxParams{.max_denominator = max_denominator});
        },
        py::arg("x"),
        py::arg("max_denominator") = 1'000'000LL,
        R"doc(
Find the best rational approximation p/q to x with q <= max_denominator.

Uses CF convergents and Farey semi-convergents (Hardy & Wright, Theorem 171).
Guaranteed to return the fraction that minimises |x - p/q| among all p/q
with denominator <= max_denominator.

Parameters
----------
x : float
    The real number to approximate.
max_denominator : int, default 1_000_000
    Maximum allowed denominator.

Returns
-------
Rational
    .numerator   : int
    .denominator : int
    .to_float()  : float

Raises
------
DomainError
    If max_denominator < 1.

Examples
--------
>>> r = adrius.best_rational(3.14159265, max_denominator=1000)
>>> str(r)
'355/113'
)doc");

    // ── simultaneous_approx ───────────────────────────────────────────────
    m.def("simultaneous_approx",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> alpha,
           double scale,
           double delta,
           double eta)
        {
            adrius::LLLParams lll_params;
            lll_params.delta = delta;
            lll_params.eta   = eta;
            return adrius::simultaneous_approx<B>(to_span(alpha), scale, lll_params);
        },
        py::arg("alpha"),
        py::arg("scale"),
        py::arg("delta") = 0.75,
        py::arg("eta")   = 0.51,
        R"doc(
One-shot simultaneous Diophantine approximation via a single LLL call.

Finds integers q, p_1, …, p_n such that max_i |q*alpha_i - p_i| ~ 1/scale^{1/n}.
For an iterated sequence with improving quality use illl() instead.

Parameters
----------
alpha : array-like of float
    Target real numbers (alpha_1, …, alpha_n).
scale : float
    Lattice scale N > 0.  Larger N gives smaller errors at higher cost.
    Expected quality ~ 1/N^{1/n} for n targets.
delta : float, default 0.75
    LLL Lovász parameter.
eta : float, default 0.51
    LLL size-reduction threshold.

Returns
-------
SimultApproxResult
    .denominator : int
    .numerators  : list of int  — [p_1, …, p_n]
    .quality     : float        — max_i |q*alpha_i - p_i|

Raises
------
DomainError
    If alpha is empty or scale <= 0.

Examples
--------
>>> import adrius, numpy as np
>>> r = adrius.simultaneous_approx([np.sqrt(2), np.pi], scale=1e6)
>>> print(r.denominator, r.numerators, r.quality)
)doc");
}
