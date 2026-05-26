# Copyright (c) 2025 InGenifold Research LLC. MIT License.

import math
import pytest
import adrius


class TestCFExpansion:
    def test_pi_known_partial_quotients(self):
        # pi = [3; 7, 15, 1, 292, 1, ...]
        assert adrius.cf_expansion(math.pi, max_depth=6) == [3, 7, 15, 1, 292, 1]

    def test_integer_terminates_immediately(self):
        result = adrius.cf_expansion(5.0, max_depth=20)
        assert result == [5]

    def test_one_half(self):
        assert adrius.cf_expansion(0.5, max_depth=10) == [0, 2]

    def test_sqrt2_periodic_twos(self):
        # sqrt(2) = [1; 2, 2, 2, ...]
        result = adrius.cf_expansion(math.sqrt(2), max_depth=5)
        assert result[0] == 1
        assert all(a == 2 for a in result[1:])

    def test_depth_limit_respected(self):
        result = adrius.cf_expansion(math.pi, max_depth=3)
        assert len(result) <= 3

    def test_all_quotients_positive(self):
        result = adrius.cf_expansion(math.sqrt(3), max_depth=10)
        assert result[0] >= 0
        assert all(a >= 1 for a in result[1:])


class TestCFConvergents:
    def test_pi_known_convergents(self):
        # pi = [3; 7, 15, 1, ...] → convergents 3/1, 22/7, 333/106, 355/113
        result = adrius.cf_convergents(math.pi, max_depth=4)
        assert result == [(3, 1), (22, 7), (333, 106), (355, 113)]

    def test_each_convergent_satisfies_dirichlet_bound(self):
        # |x - p/q| < 1/q² for every convergent (Hardy & Wright, Thm 171)
        x = math.pi
        for p, q in adrius.cf_convergents(x, max_depth=8):
            assert abs(x - p / q) < 1.0 / (q * q) + 1e-15

    def test_denominators_strictly_increasing(self):
        convs = adrius.cf_convergents(math.sqrt(2), max_depth=8)
        qs = [q for _, q in convs]
        assert qs == sorted(set(qs))  # strictly increasing (no repeats)

    def test_length_matches_expansion(self):
        x = math.e
        depth = 7
        exp = adrius.cf_expansion(x, max_depth=depth)
        convs = adrius.cf_convergents(x, max_depth=depth)
        assert len(convs) == len(exp)

    def test_first_convergent_is_floor(self):
        x = 2.718281828
        convs = adrius.cf_convergents(x, max_depth=5)
        p0, q0 = convs[0]
        assert q0 == 1
        assert p0 == math.floor(x)
