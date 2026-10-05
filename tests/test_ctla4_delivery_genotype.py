import pytest
import json
from ecd_epitope_scanner.scripts.audit_ctla4_coverage_v3 import verify_candidate_pair


def test_canonical_genotype_must_match_normalized_saved_model():
    candidate={'mutations':[{'residue':{'uniprot_pos':61,'one_letter':'A','chain_id':'A',
                           'auth_seq_id':24,'insertion_code':''},'mutant_aa':'W'}]}
    mapping={61:{'normalized_position':'25','one_letter':'A','chain_id':'A',
                 'auth_seq_id':'24','insertion_code':''}}
    verify_candidate_pair(candidate,{'mutations':['A25W'],'chain':'A'},mapping)
    # A internally valid PDB/manifest for a different substitution cannot
    # certify this canonical candidate merely because its own hashes match.
    with pytest.raises(ValueError,match='genotype differ'):
        verify_candidate_pair(candidate,{'mutations':['A25R'],'chain':'A'},mapping)
    with pytest.raises(ValueError,match='mapping identity'):
        verify_candidate_pair(candidate,{'mutations':['A25W'],'chain':'A'},
                              {61:{**mapping[61],'insertion_code':'A'}})


def test_selected_csv_json_fasta_must_describe_same_construct(tmp_path):
    from ecd_epitope_scanner.scripts.audit_ctla4_coverage_v3 import verify_selected_exports
    candidate={'mutations':[{'notation':'A1W','residue':{'uniprot_pos':1},'mutant_aa':'W'}]}
    local={'candidate_id':'c','mutation_notations':['A1W'],'model_repeats':5,
           'actual_changed_area_a2':2.,'hard_scan_center_area_a2':10.,'foldx_ddg_kcal_mol':.2}
    (tmp_path/'selected_panel_v4.json').write_text(json.dumps({'selected_candidate_ids':['c']}))
    (tmp_path/'surface_perturbation_local_constructs.json').write_text(json.dumps({'constructs':[local]}))
    (tmp_path/'surface_perturbation_local_constructs.csv').write_text(
        'candidate_id,mutations_uniprot,model_repeats,actual_changed_area_a2,hard_scan_center_area_a2,foldx_ddg_kcal_mol\n'
        'c,A1W,5,2.0,10.0,0.2\n')
    fasta=tmp_path/'full_length_protein_sequences.fasta'
    fasta.write_text('>V4_01 full-length mutations=A1W\nWCDE\n')
    verify_selected_exports(tmp_path,['c'],{'c':candidate},'ACDE')
    analysis={'designs':[{'id':'V4_01','design_role':'STAGE1_PROBE',
        'mutations':[{'wt':'A','position':1,'mutant':'W'}]},
        {'id':'PAIR_01','design_role':'STAGE2_PAIRWISE','parent_design_id':'V4_01',
         'mutations':[{'wt':'A','position':1,'mutant':'W'}]}]}
    (tmp_path/'analysis.json').write_text(json.dumps(analysis))
    fasta.write_text('>V4_01 full-length mutations=A1W\nWCDE\n>PAIR_01 full-length mutations=A1W\nWCDE\n')
    verify_selected_exports(tmp_path,['c'],{'c':candidate},'ACDE')
    fasta.write_text('>V4_01 full-length mutations=A1W\nWCDE\n>PAIR_01 full-length mutations=A1W\nRCDE\n')
    with pytest.raises(ValueError,match='FASTA sequence'):
        verify_selected_exports(tmp_path,['c'],{'c':candidate},'ACDE')
    (tmp_path/'analysis.json').unlink()
    fasta.write_text('>V4_01 full-length mutations=A1W\nRCDE\n')
    with pytest.raises(ValueError,match='FASTA sequence'):
        verify_selected_exports(tmp_path,['c'],{'c':candidate},'ACDE')
    fasta.write_text('>V4_01 full-length mutations=A1W\nWCDE\n')
    (tmp_path/'selected_panel_v4.json').write_text(json.dumps({'selected_candidate_ids':['other']}))
    with pytest.raises(ValueError,match='Selected panel JSON'):
        verify_selected_exports(tmp_path,['c'],{'c':candidate},'ACDE')
