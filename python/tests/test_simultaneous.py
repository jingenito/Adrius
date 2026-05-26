# Copyright (c) 2025 InGenifold Research LLC. MIT License.

import math
import numpy as np
import pytest
import adrius


class TestSimultaneousApprox:
    def test_result_has_expected_attributes(self):
        r = adrius.simultaneous_approx([math.sqrt(2)], scale=10.0)
        assert hasattr(r, 'denominator')
        assert hasattr(r, 'numerators')
        assert hasattr(r, 'quality')

    def test_numerators_length_matches_input(self):
        alpha = [math.sqrt(2), math.pi]
        r = adrius.simultaneous_approx(alpha, scale=10.0)
        assert len(r.numerators) == len(alpha)

    def test_quality_matches_actual_error_when_nonzero(self):
        # If q > 0, quality must equal max_i |q*alpha_i - p_i|.
        alpha = [math.sqrt(2), math.pi]
        r = adrius.simultaneous_approx(alpha, scale=10.0)
        if r.denominator > 0:
            q = r.denominator
            actual = max(abs(q * alpha[i] - r.numerators[i]) for i in range(len(alpha)))
            assert actual == pytest.approx(r.quality, abs=1e-12)

    def test_repr_contains_type_name(self):
        r = adrius.simultaneous_approx([math.sqrt(2)], scale=10.0)
        assert "SimultApproxResult" in repr(r)

    def test_numpy_array_input_accepted(self):
        alpha = np.array([math.sqrt(2), math.pi])
        r = adrius.simultaneous_approx(alpha, scale=10.0)
        assert len(r.numerators) == 2

    def test_empty_alpha_raises_domain_error(self):
        with pytest.raises(adrius.DomainError):
            adrius.simultaneous_approx([], scale=100.0)

    def test_zero_scale_raises_domain_error(self):
        with pytest.raises(adrius.DomainError):
            adrius.simultaneous_approx([math.sqrt(2)], scale=0.0)

    def test_negative_scale_raises_domain_error(self):
        with pytest.raises(adrius.DomainError):
            adrius.simultaneous_approx([math.sqrt(2)], scale=-1.0)


class TestExceptionHierarchy:
    def test_domain_error_is_subclass_of_error(self):
        assert issubclass(adrius.DomainError, adrius.Error)

    def test_convergence_error_is_subclass_of_error(self):
        assert issubclass(adrius.ConvergenceError, adrius.Error)

    def test_precision_error_is_subclass_of_error(self):
        assert issubclass(adrius.PrecisionError, adrius.Error)

    def test_error_is_subclass_of_runtime_error(self):
        assert issubclass(adrius.Error, RuntimeError)

    def test_domain_error_is_catchable_as_error(self):
        with pytest.raises(adrius.Error):
            adrius.illl([])

    def test_domain_error_is_catchable_as_runtime_error(self):
        with pytest.raises(RuntimeError):
            adrius.illl([])
