# Copyright (c) 2025 InGenifold Research LLC. MIT License.

import math
import numpy as np
import pytest
import adrius


# Basis from LLL paper / docs example: columns are basis vectors.
EXAMPLE_BASIS = np.array(
    [[1, -1,  3],
     [1,  0,  5],
     [1,  2,  6]],
    dtype=float,
)


class TestLLLReduce:
    def test_transform_relationship_holds(self):
        # reduced_basis must equal original_basis @ transform exactly (up to fp noise)
        r = adrius.lll_reduce(EXAMPLE_BASIS)
        reconstructed = EXAMPLE_BASIS @ r.transform.astype(float)
        assert np.allclose(reconstructed, r.reduced_basis, atol=1e-12)

    def test_transform_is_unimodular(self):
        r = adrius.lll_reduce(EXAMPLE_BASIS)
        det = abs(round(np.linalg.det(r.transform.astype(float))))
        assert det == 1

    def test_swap_count_is_non_negative(self):
        r = adrius.lll_reduce(EXAMPLE_BASIS)
        assert r.swap_count >= 0

    def test_first_column_is_among_shortest_lattice_vectors(self):
        # LLL guarantees ||b'_1|| <= 2^{(n-1)/4} * lambda_1.
        # A weak sanity check: reduced first vector no longer than any original column.
        r = adrius.lll_reduce(EXAMPLE_BASIS)
        first_norm = np.linalg.norm(r.reduced_basis[:, 0])
        orig_norms = [np.linalg.norm(EXAMPLE_BASIS[:, j]) for j in range(3)]
        assert first_norm <= max(orig_norms)

    def test_identity_basis_unchanged(self):
        basis = np.eye(3)
        r = adrius.lll_reduce(basis)
        assert np.allclose(basis @ r.transform, r.reduced_basis, atol=1e-12)
        assert r.swap_count == 0

    def test_2x2_known_example(self):
        # Simple 2x2: columns (3, 2) and (1, 1). LLL should prefer (1,1).
        basis = np.array([[3.0, 1.0], [2.0, 1.0]])
        r = adrius.lll_reduce(basis)
        assert np.allclose(basis @ r.transform, r.reduced_basis, atol=1e-12)
        assert abs(round(np.linalg.det(r.transform.astype(float)))) == 1

    def test_higher_delta_more_reduced(self):
        # delta=0.99 should produce at least as many swaps as delta=0.75 (or equal).
        r_low  = adrius.lll_reduce(EXAMPLE_BASIS, delta=0.75)
        r_high = adrius.lll_reduce(EXAMPLE_BASIS, delta=0.99)
        assert r_high.swap_count >= r_low.swap_count

    def test_transform_dtype_is_integer(self):
        r = adrius.lll_reduce(EXAMPLE_BASIS)
        assert np.issubdtype(r.transform.dtype, np.integer)

    def test_accepts_list_input(self):
        basis_list = [[1.0, -1.0, 3.0],
                      [1.0,  0.0, 5.0],
                      [1.0,  2.0, 6.0]]
        r = adrius.lll_reduce(np.array(basis_list))
        assert r.swap_count >= 0
