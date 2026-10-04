"""Freeze interface geometry labels using explicit mmCIF UniProt alignment.

This module has no dependency on candidate generation, calibration or design.
Labels from receptor/ligand interfaces are never called antibody epitopes.
"""
from __future__ import annotations
from pathlib import Path
import hashlib
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1
from scipy.spatial import cKDTree
from ..util import sha256_file


def label_mapping(data,accession,chain,sequence,*,declared_conflicts=None):
    refs={ref:acc for ref,db,acc in zip(data['_struct_ref.id'],data['_struct_ref.db_name'],
        data['_struct_ref.pdbx_db_accession']) if db=='UNP'}
    mapping={}
    columns=('_struct_ref_seq.ref_id','_struct_ref_seq.pdbx_strand_id',
             '_struct_ref_seq.seq_align_beg','_struct_ref_seq.seq_align_end',
             '_struct_ref_seq.db_align_beg','_struct_ref_seq.db_align_end')
    for ref,chains,left,right,start,end in zip(*(data[k] for k in columns)):
        if refs.get(ref)!=accession or chain not in chains.split(','):continue
        left,right,start,end=map(int,(left,right,start,end))
        if right-left!=end-start or not 1<=start<=end<=len(sequence):
            raise ValueError('ambiguous or out-of-range UniProt alignment')
        for label in range(left,right+1):
            position=start+label-left
            if label in mapping and mapping[label]!=position:raise ValueError('conflicting UniProt alignments')
            mapping[label]=position
    if not mapping:raise ValueError('no explicit UniProt alignment for target author chain')
    conflict_columns=('_struct_ref_seq_dif.pdbx_pdb_strand_id','_struct_ref_seq_dif.seq_num',
        '_struct_ref_seq_dif.pdbx_seq_db_accession_code','_struct_ref_seq_dif.mon_id',
        '_struct_ref_seq_dif.db_mon_id','_struct_ref_seq_dif.pdbx_seq_db_seq_num')
    declared={}
    if all(k in data for k in conflict_columns):
        for author,label,acc,actual,expected,position in zip(*(data[k] for k in conflict_columns)):
            if author==chain and acc==accession:
                declared[int(label)]=(seq1(actual),seq1(expected),int(position))
    for author,label,residue in zip(data['_atom_site.auth_asym_id'],data['_atom_site.label_seq_id'],
                                    data['_atom_site.label_comp_id']):
        if author!=chain or label in {'.','?'} or int(label) not in mapping:continue
        position=mapping[int(label)]
        if seq1(residue,custom_map={'MSE':'M'})!=sequence[position-1]:
            observed=seq1(residue,custom_map={'MSE':'M'})
            if declared_conflicts is not None and declared.get(int(label))==(observed,sequence[position-1],position):
                declared_conflicts.append({'label_seq_id':int(label),'uniprot_pos':position,
                    'observed':observed,'canonical':sequence[position-1],
                    'status':'EXCLUDED_EXPLICIT_DATABASE_SEQUENCE_CONFLICT'})
                del mapping[int(label)]
                continue
            raise ValueError(f'contact source sequence mismatch at UniProt {position}')
    return mapping


def freeze_contacts(path,accession,chain,partner_chains,sequence,name,kind,cutoff=4.5):
    if kind not in {'antibody_interface_geometry','receptor_interface_geometry','ligand_interface_geometry'}:
        raise ValueError('explicit interface label class required')
    if chain in partner_chains or not partner_chains or cutoff<=0:raise ValueError('invalid contact partners/cutoff')
    data=MMCIF2Dict(str(path));excluded=[]
    mapping=label_mapping(data,accession,chain,sequence,declared_conflicts=excluded)
    columns=('_atom_site.auth_asym_id','_atom_site.label_seq_id','_atom_site.auth_seq_id',
             '_atom_site.pdbx_PDB_ins_code','_atom_site.label_atom_id','_atom_site.type_symbol',
             '_atom_site.Cartn_x','_atom_site.Cartn_y','_atom_site.Cartn_z')
    count=len(data[columns[0]])
    model=data.get('_atom_site.pdbx_PDB_model_num',['1']*count)
    alt=data.get('_atom_site.label_alt_id',['.']*count)
    target=[];partners=[]
    for row,model_id,alternative in zip(zip(*(data[k] for k in columns)),model,alt):
        author,label,residue,icode,atom,element,x,y,z=row
        if model_id!=model[0] or alternative not in {'.','?','A'} or element.upper() in {'H','D'} or label in {'.','?'}:continue
        xyz=np.array([float(x),float(y),float(z)])
        if not np.isfinite(xyz).all():raise ValueError('nonfinite interface coordinates')
        if author==chain and int(label) in mapping:
            target.append((mapping[int(label)],residue,'' if icode in {'.','?'} else icode,atom,xyz))
        elif author in partner_chains:partners.append(xyz)
    if not target or not partners:raise ValueError('target or partner heavy atoms missing')
    distances=cKDTree(np.array(partners)).query(np.array([row[4] for row in target]))[0]
    contacts={}
    for (position,residue,icode,atom,_),distance in zip(target,distances):
        if distance<=cutoff:
            record=contacts.setdefault(position,{'uniprot_pos':position,'author_chain':chain,
                'auth_seq_id':residue,'insertion_code':icode,'minimum_heavy_distance_a':float(distance)})
            record['minimum_heavy_distance_a']=min(record['minimum_heavy_distance_a'],float(distance))
    if not contacts:raise ValueError('specified chains have no contacts under the frozen cutoff')
    return {'accession':accession,'numbering':'UniProt precursor numbering',
        'excluded_declared_sequence_conflicts':excluded,
        'source_sha256':sha256_file(Path(path)),'sequence_sha256':hashlib.sha256(sequence.encode()).hexdigest(),
        'contacts':{name:{'uniprot_positions':sorted(contacts),'kind':kind,'source_file':Path(path).name,
            'target_author_chain':chain,'partner_author_chains':list(partner_chains),
            'heavy_atom_cutoff_a':cutoff,'mapping':'explicit mmCIF UNP alignment; canonical sequence checked',
            'mapped_contact_records':[contacts[p] for p in sorted(contacts)]}},
        'interpretation':'hidden geometric sanity labels; no functional mutation or antibody detection validation'}
