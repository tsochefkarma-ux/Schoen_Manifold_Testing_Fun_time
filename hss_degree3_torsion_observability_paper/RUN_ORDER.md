# Reproduction order

Activate SageMath and run from the bundle root.

1. `28e_specialize_first_worker.sage` through `run_calibration_v4_specialize_first.sh`
2. `29_exact_polar_dependency_certificate.sage`
3. `30_degree3_bootstrap_dimension_census.sage`
4. `31_shared_bijacobi_torsion_model.sage`
5. `32_degree3_A3_B3_hss_torsion_candidate_v2.sage 40`
6. `33_degree3_evaluation_kernel_and_sector_lift.sage`
7. `34_degree3_certified_15d_quotient_v3.sage`
8. `35b_degree3_torsion_modularity_audit.sage`
9. `36_degree2_full_space_baseline_and_degree3_selectivity_ledger.sage`
10. `37_degree3_cm_jet_extraction.sage`
11. `37b_degree3_linear_curvature_response.sage` (prototype)
12. `37b2_degree3_response_rank_growth_audit.sage` (authoritative scalar-rank audit)
13. `37c_degree3_information_loss_and_observable_stack.sage`
14. `37d_degree3_minimal_observable_lift.sage`
15. `37e_degree3_parity_block_mixing_certificate_v2.sage`
16. `37f_degree3_exchange_parity_decomposition.sage`
17. `38_degree2_tensorial_parity_calibration.sage`
18. `39_degree3_torsion_recombination_constraint_census.sage`
19. `40_degree3_exact_reconstruction_atlas.sage`
20. `41_degree3_minimal_external_data_bridge.sage`

Recommended long-run settings:

```bash
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

Use `tmux` for Cells 37A-37D.
