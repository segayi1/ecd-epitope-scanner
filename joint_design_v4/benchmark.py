"""Input contract for deferred retrospective validation; no fabricated results."""
from __future__ import annotations


def benchmark_contract():
    return {
        "schema_version": "1.0", "status": "INTERFACE_ONLY_NOT_EXECUTED",
        "required_fields": ["complex_id", "antigen_structure_path", "antigen_chain", "structure_sha256", "residue_mapping", "hidden_truth_path", "hidden_truth_sha256"],
        "mapping_fields": ["chain_id", "auth_seq_id", "insertion_code", "uniprot_pos"],
        "hidden_truth_fields": ["antibody_chains", "contact_definition", "contact_cutoff_a", "antigen_contact_residues"],
        "design_input_policy": "Antibody coordinates and contact labels must be absent from antigen-only design inputs and selection.",
        "tiers": ["PD1", "PDL1", "CTLA4"], "budget_curve": [4,5,6,7,8],
        "metrics": ["top1_localization", "top3_localization", "geodesic_localization_error_a", "weighted_jaccard", "information_per_construct"],
        "interpretation": "Geometric simulated responses do not validate experimental SNIPR reporter or affinity effects.",
    }


def validate_benchmark_case(case):
    missing = set(benchmark_contract()["required_fields"]) - set(case)
    if missing:
        raise ValueError("Missing benchmark fields: " + ", ".join(sorted(missing)))
    if case["antigen_structure_path"] == case["hidden_truth_path"]:
        raise ValueError("Hidden truth must be stored separately from the design structure")
    mapping = case["residue_mapping"]
    if not isinstance(mapping, list) or not mapping:
        raise ValueError("Benchmark requires an explicit residue mapping")
    for row in mapping:
        if set(benchmark_contract()["mapping_fields"]) - set(row):
            raise ValueError("Incomplete benchmark residue mapping")
    return dict(case)
