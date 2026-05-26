# Copyright (c) 2025 InGenifold Research LLC. MIT License.

import math
import numpy as np
import pytest
import adrius


class TestILLL:
    def test_result_has_expected_attributes(self):
        r = adrius.illl([math.sqrt(2)], max_iterations=3)
        assert hasattr(r, 'relations')
        assert hasattr(r, 'quality')
        assert hasattr(r, 'iterations')

    def test_iterations_matches_relation_and_quality_length(self):
        r = adrius.illl([math.sqrt(2), math.pi], max_iterations=6)
        assert r.iterations == len(r.relations)
        assert r.iterations == len(r.quality)

    def test_relation_width_is_n_plus_1(self):
        # Each relation row is [q, p_1, ..., p_n]
        n = 2
        r = adrius.illl([math.sqrt(2), math.pi], max_iterations=5)
        for row in r.relations:
            assert len(row) == n + 1

    def test_denominators_non_decreasing(self):
        r = adrius.illl([math.sqrt(2), math.pi], max_iterations=8)
        qs = [row[0] for row in r.relations]
        assert all(qs[i] <= qs[i + 1] for i in range(len(qs) - 1))

    def test_quality_non_increasing(self):
        r = adrius.illl([math.sqrt(2), math.pi], max_iterations=8)
        assert all(r.quality[i] >= r.quality[i + 1]
                   for i in range(r.iterations - 1))

    def test_quality_matches_actual_error(self):
        # For each relation [q, p_1, ..., p_n], quality = max_i |q*alpha_i - p_i|
        alpha = [math.sqrt(2), math.pi]
        r = adrius.illl(alpha, max_iterations=10)
        for k in range(r.iterations):
            q = r.relations[k][0]
            ps = r.relations[k][1:]
            actual_quality = max(abs(q * alpha[i] - ps[i]) for i in range(len(alpha)))
            assert actual_quality == pytest.approx(r.quality[k], abs=1e-11)

    def test_one_dimensional_quality(self):
        alpha = [math.sqrt(2)]
        r = adrius.illl(alpha, max_iterations=5, max_denominator=1000)
        for k in range(r.iterations):
            q, p = r.relations[k]
            assert abs(q * math.sqrt(2) - p) == pytest.approx(r.quality[k], abs=1e-12)

    def test_max_denominator_bound_respected(self):
        # The relation that triggers the stop is stored as the final entry, so
        # at most one q can be >= max_denominator; all preceding ones are strictly below.
        r = adrius.illl([math.sqrt(2), math.pi], max_denominator=100)
        qs = [row[0] for row in r.relations]
        assert all(q < 100 for q in qs[:-1])
        assert sum(1 for q in qs if q >= 100) <= 1

    def test_pi_approximation_22_over_7(self):
        # Bosma-Smeets should recover 22/7 as an early denominator for pi
        r = adrius.illl([math.pi], max_iterations=20, max_denominator=10_000)
        denominators = [row[0] for row in r.relations]
        assert 7 in denominators

    def test_repr_contains_iterations(self):
        r = adrius.illl([math.sqrt(2)], max_iterations=3)
        assert "ILLLResult" in repr(r)
        assert "iterations" in repr(r)

    def test_empty_alpha_raises_domain_error(self):
        with pytest.raises(adrius.DomainError):
            adrius.illl([])

    def test_epsilon_above_1_raises_domain_error(self):
        with pytest.raises(adrius.DomainError):
            adrius.illl([math.sqrt(2)], epsilon=1.5)

    def test_epsilon_zero_raises_domain_error(self):
        with pytest.raises(adrius.DomainError):
            adrius.illl([math.sqrt(2)], epsilon=0.0)

    def test_numpy_array_input_accepted(self):
        alpha = np.array([math.sqrt(2), math.pi])
        r = adrius.illl(alpha, max_iterations=5)
        assert r.iterations > 0
