## Summary
<!-- Brief description of what this PR changes and why -->

## Type of Change
<!-- Mark the relevant option with an "x" -->
- [ ] Feature (new functionality)
- [ ] Fix (bug fix)
- [ ] Documentation (docs/guides)
- [ ] Refactoring (no behavior change)
- [ ] Tests (test suite improvements)

## Algorithm Changes
<!-- Only required if this PR implements or modifies an algorithm -->
<!-- If no algorithm changes, remove this section -->
- [ ] New algorithm documentation in `docs/algorithm-name.md`
- [ ] Research papers referenced in PR body with specific sections
- [ ] Examples added to README.md
- [ ] Test suite covers hand-computed examples and bounds

## Research Papers
<!-- List any papers referenced by this PR -->
<!-- Include specific sections/equations used (e.g., Lemma 3.4, eq. 28) -->
- **Citation**: Authors. "Title." Journal/Conference YYYY.
  - Paper reference: arXiv link or DOI
  - Sections/equations used: (e.g., Algorithm 1, Lemma 3.4, eq. 28)

Example:
- **Bosma, W. & Smeets, I.** (2010). "Finding simultaneous Diophantine approximations with prescribed quality." arXiv:1001.4455.
  - Sections: §2.1 (lattice construction), Lemma 3.4 (scale decay), Theorem 3.5 (bounds)

## Testing
<!-- Proof that tests pass locally -->
- [ ] New tests added (if applicable)
- [ ] All tests pass locally: `ctest --test-dir build --output-on-failure`
- [ ] No compiler errors or warnings

**Test results**:
```
[Paste output from: ctest --test-dir build --output-on-failure]
```

Example:
```
100% tests passed, 0 tests failed out of 46

Total Test time (real) =   0.53 sec
```

## Documentation
- [ ] CLAUDE.md updated if adding new style rules (unlikely)
- [ ] DEVELOPMENT.md updated if adding new workflows (unlikely)
- [ ] Algorithm documentation added if new algorithm (required for feature)
- [ ] README.md updated if new public API (required for feature)
- [ ] Code comments explain *why*, not *what* (CLAUDE.md §Code Style Rules)

## Checklist
- [ ] Branch created off latest `main`: `git checkout main && git pull && git checkout -b feature/...`
- [ ] Commit messages follow format: imperative mood, reference papers if applicable
- [ ] All tests pass
- [ ] No breaking changes to public API (unless intentional with explanation)
- [ ] Code follows CLAUDE.md style rules

## Related Issues
<!-- Link to any related issues: Closes #123 -->

---

**Note**: This PR description uses the template from `.github/pull_request_template.md`. For detailed guidance, see:
- **CLAUDE.md** — Code style, API design, and commit conventions
- **DEVELOPMENT.md** — Workflow and testing procedures
