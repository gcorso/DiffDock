# Significant contribution from Ben Fry and Nick Polizzi

three_to_one = {'ALA':	'A',
                'ARG':	'R',
                'ASN':	'N',
                'ASP':	'D',
                'CYS':	'C',
                'GLN':	'Q',
                'GLU':	'E',
                'GLY':	'G',
                'HIS':	'H',
                'ILE':	'I',
                'LEU':	'L',
                'LYS':	'K',
                'MET':	'M',
                'MSE':  'M', # MSE this is almost the same AA as MET. The sulfur is just replaced by Selen
                'PHE':	'F',
                'PRO':	'P',
                'PYL':	'O',
                'SER':	'S',
                'SEC':	'U',
                'THR':	'T',
                'TRP':	'W',
                'TYR':	'Y',
                'VAL':	'V',
                'ASX':	'B',
                'GLX':	'Z',
                'XAA':	'X',
                'XLE':	'J'}


aa_name2aa_idx = {'ALA': 0, 'ARG': 1, 'ASN': 2, 'ASP': 3, 'CYS': 4, 'GLU': 5, 'GLN': 6, 'GLY': 7,
                  'HIS': 8, 'ILE': 9, 'LEU': 10, 'LYS': 11, 'MET': 12, 'PHE': 13, 'PRO': 14,
                  'SER': 15, 'THR': 16, 'TRP': 17, 'TYR': 18, 'VAL': 19, 'MSE': 12}


aa_short2long = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS', 'I': 'ILE',
                 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 'G': 'GLY', 'H': 'HIS',
                 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 'A': 'ALA', 'V': 'VAL', 'E': 'GLU',
                 'Y': 'TYR', 'M': 'MET'}


aa_short2aa_idx = {aa_short: aa_name2aa_idx[aa_long] for aa_short, aa_long in aa_short2long.items()}
aa_idx2aa_short = {aa_idx: aa_short for aa_short, aa_idx in aa_short2aa_idx.items()}
aa_long2short = {aa_long: aa_short for aa_short, aa_long in aa_short2long.items()}
aa_long2short['MSE'] = 'M'

chi = { 'C' :
        { 1: ('N'  , 'CA' , 'CB' , 'SG' )   },
        'D' :
        { 1: ('N'  , 'CA' , 'CB' , 'CG' ),
          2: ('CA' , 'CB' , 'CG' , 'OD1'), },
        'E' :
        { 1: ('N'  , 'CA' , 'CB' , 'CG' ),
          2: ('CA' , 'CB' , 'CG' , 'CD' ),
          3: ('CB' , 'CG' , 'CD' , 'OE1'), },
        'F' :
        { 1: ('N'  , 'CA' , 'CB' , 'CG' ),
          2: ('CA' , 'CB' , 'CG' , 'CD1'), },
        'H' :
        { 1: ('N'  , 'CA' , 'CB' , 'CG' ),
          2: ('CA' , 'CB' , 'CG' , 'ND1'), },
        'I' :
        { 1: ('N'  , 'CA' , 'CB' , 'CG1'),
          2: ('CA' , 'CB' , 'CG1', 'CD1'), },
        'K' :
        { 1: ('N'  , 'CA' , 'CB'  ,'CG' ),
          2: ('CA' , 'CB' , 'CG'  ,'CD' ),
          3: ('CB' , 'CG' , 'CD'  ,'CE' ),
          4: ('CG' , 'CD' , 'CE'  ,'NZ' ), },
        'L' :
        { 1: ('N'  , 'CA' , 'CB' , 'CG' ),
          2: ('CA' , 'CB' , 'CG' , 'CD1'), },
        'M' :
        { 1: ('N'  , 'CA' , 'CB'  ,'CG' ),
          2: ('CA' , 'CB' , 'CG'  ,'SD' ),
          3: ('CB' , 'CG' , 'SD'  ,'CE' ), },
        'N' :
        { 1: ('N'  , 'CA' , 'CB' , 'CG' ),
          2: ('CA' , 'CB' , 'CG' , 'OD1'), },
        'P' :
        { 1: ('N'  , 'CA' , 'CB' , 'CG' ),
          2: ('CA' , 'CB' , 'CG' , 'CD' ), },
        'Q' :
        { 1: ('N'  , 'CA' , 'CB' , 'CG' ),
          2: ('CA' , 'CB' , 'CG' , 'CD' ),
          3: ('CB' , 'CG' , 'CD' , 'OE1'), },
        'R' :
        { 1: ('N'  , 'CA' , 'CB'  ,'CG' ),
          2: ('CA' , 'CB' , 'CG'  ,'CD' ),
          3: ('CB' , 'CG' , 'CD'  ,'NE' ),
          4: ('CG' , 'CD' , 'NE'  ,'CZ' ), },
        'S' :
        { 1: ('N'  , 'CA' , 'CB' , 'OG' ), },
        'T' :
        { 1: ('N'  , 'CA' , 'CB' , 'OG1'), },
        'V' :
        { 1: ('N'  , 'CA' , 'CB' , 'CG1'), },
        'W' :
        { 1: ('N'  , 'CA' , 'CB' , 'CG' ),
          2: ('CA' , 'CB' , 'CG' , 'CD1'), },
        'Y' :
        { 1: ('N'  , 'CA' , 'CB' , 'CG' ),
          2: ('CA' , 'CB' , 'CG' , 'CD1'), },
        }


atom_order = {'G': ['N', 'CA', 'C', 'O'],
'A': ['N', 'CA', 'C', 'O', 'CB'],
'S': ['N', 'CA', 'C', 'O', 'CB', 'OG'],
'C': ['N', 'CA', 'C', 'O', 'CB', 'SG'],
'T': ['N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2'],
'P': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD'],
'V': ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2'],
'M': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE'],
'N': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2'],
'I': ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1'],
'L': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2'],
'D': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2'],
'E': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2'],
'K': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'],
'Q': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2'],
'H': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'],
'F': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
'R': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
'Y': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
'W': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE2', 'CE3', 'NE1', 'CZ2', 'CZ3', 'CH2'],
'X': ['N', 'CA', 'C', 'O']}     # unknown amino acid



amino_acid_smiles = {
    'PHE': '[NH3+]CC(=O)N[C@@H](Cc1ccccc1)C(=O)NCC(=O)O',
    'MET': 'CSCC[C@H](NC(=O)C[NH3+])C(=O)NCC(=O)O',
    'TYR': '[NH3+]CC(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)NCC(=O)O',
    'ILE': 'CC[C@H](C)[C@H](NC(=O)C[NH3+])C(=O)NCC(=O)O',
    'TRP': '[NH3+]CC(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)NCC(=O)O',
    'THR': 'C[C@@H](O)[C@H](NC(=O)C[NH3+])C(=O)NCC(=O)O',
    'CYS': '[NH3+]CC(=O)N[C@@H](CS)C(=O)NCC(=O)O',
    'ALA': 'C[C@H](NC(=O)C[NH3+])C(=O)NCC(=O)O',
    'LYS': '[NH3+]CCCC[C@H](NC(=O)C[NH3+])C(=O)NCC(=O)O',
    'PRO': '[NH3+]CC(=O)N1CCC[C@H]1C(=O)NCC(=O)O',
    'LEU': 'CC(C)C[C@H](NC(=O)C[NH3+])C(=O)NCC(=O)O',
    'GLY': '[NH3+]CC(=O)NCC(=O)NCC(=O)O',
    'ASP': '[NH3+]CC(=O)N[C@@H](CC(=O)O)C(=O)NCC(=O)O',
    'HIS': '[NH3+]CC(=O)N[C@@H](Cc1c[nH]c[nH+]1)C(=O)NCC(=O)O',
    'VAL': 'CC(C)[C@H](NC(=O)C[NH3+])C(=O)NCC(=O)O',
    'SER': '[NH3+]CC(=O)N[C@@H](CO)C(=O)NCC(=O)O',
    'ARG': 'NC(=[NH2+])NCCC[C@H](NC(=O)C[NH3+])C(=O)NCC(=O)O',
    'GLU': '[NH3+]CC(=O)N[C@@H](CCC(=O)O)C(=O)NCC(=O)O',
    'GLN': 'NC(=O)CC[C@H](NC(=O)C[NH3+])C(=O)NCC(=O)O',
    'ASN': 'NC(=O)C[C@H](NC(=O)C[NH3+])C(=O)NCC(=O)O',
 }

cg_rdkit_indices = {
    'PHE': {4: 'N', 5: 'CA', 13: 'C', 14: 'O', 6: 'CB', 7: 'CG', 8: 'CD1', 12: 'CD2', 9: 'CE1', 11: 'CE2', 10: 'CZ'},
    'MET': {5: 'N', 4: 'CA', 10: 'C', 11: 'O', 3: 'CB', 2: 'CG', 1: 'SD', 0: 'CE'},
    'TYR': {4: 'N', 5: 'CA', 14: 'C', 15: 'O', 6: 'CB', 7: 'CG', 8: 'CD1', 13: 'CD2', 9: 'CE1', 12: 'CE2', 10: 'CZ', 11: 'OH'},
    'ILE': {5: 'N', 4: 'CA', 10: 'C', 11: 'O', 2: 'CB', 1: 'CG1', 3: 'CG2', 0: 'CD1'},
    'TRP': {4: 'N', 5: 'CA', 16: 'C', 17: 'O', 6: 'CB', 7: 'CG', 8: 'CD1', 15: 'CD2', 9: 'NE1', 10: 'CE2', 14: 'CE3', 11: 'CZ2', 13: 'CZ3', 12: 'CH2'},
    'THR': {4: 'N', 3: 'CA', 9: 'C', 10: 'O', 1: 'CB', 2: 'OG1', 0: 'CG2'},
    'CYS': {4: 'N', 5: 'CA', 8: 'C', 9: 'O', 6: 'CB', 7: 'SG'},
    'ALA': {2: 'N', 1: 'CA', 7: 'C', 8: 'O', 0: 'CB'},
    'LYS': {6: 'N', 5: 'CA', 11: 'C', 12: 'O', 4: 'CB', 3: 'CG', 2: 'CD', 1: 'CE', 0: 'NZ'},
    'PRO': {4: 'N', 8: 'CA', 9: 'C', 10: 'O', 7: 'CB', 6: 'CG', 5: 'CD'},
    'LEU': {5: 'N', 4: 'CA', 10: 'C', 11: 'O', 3: 'CB', 1: 'CG', 0: 'CD1', 2: 'CD2'},
    'GLY': {4: 'N', 5: 'CA', 6: 'C', 7: 'O'},
    'ASP': {4: 'N', 5: 'CA', 10: 'C', 11: 'O', 6: 'CB', 7: 'CG', 8: 'OD1', 9: 'OD2'},
    'HIS': {4: 'N', 5: 'CA', 12: 'C', 13: 'O', 6: 'CB', 7: 'CG', 11: 'ND1', 8: 'CD2', 10: 'CE1', 9: 'NE2'},
    'VAL': {4: 'N', 3: 'CA', 9: 'C', 10: 'O', 1: 'CB', 0: 'CG1', 2: 'CG2'},
    'SER': {4: 'N', 5: 'CA', 8: 'C', 9: 'O', 6: 'CB', 7: 'OG'},
    'ARG': {8: 'N', 7: 'CA', 13: 'C', 14: 'O', 6: 'CB', 5: 'CG', 4: 'CD', 3: 'NE', 1: 'CZ', 0: 'NH1', 2: 'NH2'},
    'GLU': {4: 'N', 5: 'CA', 11: 'C', 12: 'O', 6: 'CB', 7: 'CG', 8: 'CD', 9: 'OE1', 10: 'OE2'},
    'GLN': {6: 'N', 5: 'CA', 11: 'C', 12: 'O', 4: 'CB', 3: 'CG', 1: 'CD', 2: 'OE1', 0: 'NE2'},
    'ASN': {5: 'N', 4: 'CA', 10: 'C', 11: 'O', 3: 'CB', 1: 'CG', 2: 'OD1', 0: 'ND2'}
}

aa_to_cg_indices = {aa_long2short[x]: [atom_order[aa_long2short[x]].index(aname) for aname in index_dict.values()]  for x, index_dict in cg_rdkit_indices.items()}

