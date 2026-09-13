# ReGIFT 0.1.0

## Frozen r4 release

- Package the final frozen r4 estimator (base 0.0.1.9007 plus conditional donor-jackknife state shrinkage).
- Apply state shrinkage through normal `regift_fit()` and `regift()` calls; external loaders and function-body patching are no longer needed.
- Preserve the base-model one-sweep penalty calibration in all high-level and grouped fitting entry points.
- Set the primary `regift()` preset to K=5, H=2, lambda_fraction=1/32, lambda_Delta=3, max_iter=100 and tol=1e-6. HPAP analyses specify K=10 explicitly.
- Include the frozen fitting, held-out projection, inference, grouped tuning and working-response corrections that were absent from version 0.0.1.
- Retain the original algorithmic scope and fallback behavior; unsupported state-shrinkage designs retain the base anchors and report a diagnostic.

# ReGIFT 0.0.1.9007 (unreleased paired-anchor candidate)

* In paired designs, initialize the common condition anchor conditional on donor and state baselines.
* Retain within-donor baseline removal when rare states have fewer than three complete donor-specific contrasts; jointly identifiable directions from partial pairs remain usable.
* Record within-donor versus between-donor anchor estimation explicitly. Keep the 9006 nonpaired branch and unsupported-direction pooling rule.
* Preserve the declared robust summary for states with sufficient complete pairs. The scalar response reliability and continuous-response optimizer are unchanged.
* Synthetic invariance checks pass; no universal benchmark superiority is claimed.

# ReGIFT 0.0.1.9006 (unreleased condition-anchor candidate)

* Separate baseline expression from condition slopes when initializing global and state-specific response anchors. Unequal numbers of unpaired donors must not turn a state baseline into a condition response.
* Use weighted centering with an unpenalized intercept and the existing numerical ridge on slopes. Keep biological-replicate weights and the declared condition contrast.
* With state anchors enabled, initialize the common response conditional on state intercepts. In incomplete states, estimate only supported deviation directions and explicitly pool unidentified directions to the common response; save ranks, identified projections and condition support. Globally unidentifiable effects still fail. The paired donor-specific intercept and optional refit paths are retained.
* Analytical null, nonzero-effect and baseline-shift checks pass. Real-data superiority remains unestablished.

# ReGIFT 0.0.1.9005 (unreleased development candidate)

* Match donor sizes to rowsum row names in the working response, making unequal-donor normalization invariant to cell row order.
* Propagate fixed model settings and row-level strata consistently through grouped tuning, lambda probes and final fitting; record the actual settings.
* Batch equivalent held-out latent solves by condition; preserve scalar projection results.
* Skip absent states in response extraction to avoid zero-extent matrix warnings.
* Preserve all failed, capped and unfavorable development outcomes. Formal superiority and resource limits have not been established.

# ReGIFT 0.0.1.9004 (unreleased candidate)

* Preserve frozen state anchors in held-out projection and reuse canonical response calibration.
* Require known state annotations for state-anchored held-out targets.
* No real-data benchmark superiority has been established.

# ReGIFT 0.0.1

* Initial development release of the CPU-based ReGIFT model.
* Added the high-level `regift()` workflow and response-table summaries.
* Added paired-donor inference, held-out donor projection, and grouped tuning.
* Added a compact, deterministic HPAP pancreas subset as `regift_example`.
