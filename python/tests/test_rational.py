# Copyright (c) 2025 InGenifold Research LLC. MIT License.

import math
import pytest
import adrius


class TestRational:
    def test_numerator_and_denominator(self):
        r = adrius.Rational(3, 4)
        assert r.numerator == 3
        assert r.denominator == 4

    def test_reduces_to_lowest_terms(self):
        r = adrius.Rational(6, 4)
        assert r.numerator == 3
        assert r.denominator == 2

    def test_negative_denominator_normalized(self):
        r = adrius.Rational(3, -4)
        assert r.numerator == -3
        assert r.denominator == 4

    def test_from_single_integer(self):
        r = adrius.Rational(7)
        assert r.numerator == 7
        assert r.denominator == 1

    def test_zero(self):
        r = adrius.Rational(0)
        assert r.is_zero()
        assert r.numerator == 0
        assert r.denominator == 1

    def test_to_float(self):
        assert adrius.Rational(1, 4).to_float() == pytest.approx(0.25)

    def test_float_dunder(self):
        assert float(adrius.Rational(1, 3)) == pytest.approx(1 / 3)

    def test_int_dunder_on_integer(self):
        assert int(adrius.Rational(6, 2)) == 3

    def test_int_dunder_on_non_integer_raises(self):
        with pytest.raises(ValueError):
            int(adrius.Rational(1, 3))

    def test_is_integer_true(self):
        assert adrius.Rational(4, 2).is_integer()

    def test_is_integer_false(self):
        assert not adrius.Rational(1, 3).is_integer()

    def test_addition(self):
        assert adrius.Rational(1, 2) + adrius.Rational(1, 3) == adrius.Rational(5, 6)

    def test_addition_to_zero(self):
        assert adrius.Rational(1, 2) + adrius.Rational(-1, 2) == adrius.Rational(0)

    def test_subtraction(self):
        assert adrius.Rational(3, 4) - adrius.Rational(1, 4) == adrius.Rational(1, 2)

    def test_multiplication(self):
        assert adrius.Rational(2, 3) * adrius.Rational(3, 4) == adrius.Rational(1, 2)

    def test_division(self):
        assert adrius.Rational(1, 2) / adrius.Rational(1, 4) == adrius.Rational(2)

    def test_negation(self):
        assert -adrius.Rational(3, 5) == adrius.Rational(-3, 5)
        assert -adrius.Rational(-3, 5) == adrius.Rational(3, 5)

    def test_equality(self):
        assert adrius.Rational(2, 4) == adrius.Rational(1, 2)

    def test_less_than(self):
        assert adrius.Rational(1, 3) < adrius.Rational(1, 2)

    def test_less_than_or_equal(self):
        assert adrius.Rational(1, 3) <= adrius.Rational(1, 3)
        assert adrius.Rational(1, 4) <= adrius.Rational(1, 3)

    def test_greater_than(self):
        assert adrius.Rational(2, 3) > adrius.Rational(1, 2)

    def test_greater_than_or_equal(self):
        assert adrius.Rational(1, 2) >= adrius.Rational(1, 2)
        assert adrius.Rational(2, 3) >= adrius.Rational(1, 2)

    def test_repr(self):
        assert repr(adrius.Rational(3, 4)) == "Rational(3, 4)"

    def test_str(self):
        assert str(adrius.Rational(3, 4)) == "3/4"

    def test_str_integer(self):
        # operator<< omits /1 for whole numbers
        assert str(adrius.Rational(5)) == "5"


class TestBestRational:
    def test_pi_355_over_113(self):
        r = adrius.best_rational(3.14159265, max_denominator=1000)
        assert r.numerator == 355
        assert r.denominator == 113

    def test_half_exact(self):
        r = adrius.best_rational(0.5, max_denominator=10)
        assert r.numerator == 1
        assert r.denominator == 2

    def test_exact_integer(self):
        r = adrius.best_rational(3.0, max_denominator=100)
        assert r.numerator == 3
        assert r.denominator == 1

    def test_denominator_does_not_exceed_limit(self):
        r = adrius.best_rational(math.sqrt(2), max_denominator=50)
        assert r.denominator <= 50

    def test_result_is_rational_type(self):
        r = adrius.best_rational(math.pi, max_denominator=100)
        assert isinstance(r, adrius.Rational)

    def test_invalid_denominator_raises(self):
        with pytest.raises(adrius.DomainError):
            adrius.best_rational(1.0, max_denominator=0)
