PDB_id = '7nei' #@param {type:"string"}
SMILES_or_pubchem_id = 'CC(=O)C1=CC=C(C=C1)C(=O)OCCOC' #@param {type:"string"}
chain = 'A' #@param {type:"string"}
import os
import requests
import time
from random import random

def download_pdb_file(pdb_id: str) -> str:
    """Download pdb file as a string from rcsb.org"""
    PDB_DIR ="./tmp/pdb/"
    os.makedirs(PDB_DIR, exist_ok=True)

    # url or pdb_id
    if pdb_id.startswith('http'):
        url = pdb_id
        filename = url.split('/')[-1]
    else:
        url = f"http://files.rcsb.org/view/{pdb_id}.pdb"
        filename = f'{pdb_id}.pdb'

    cache_path = os.path.join(PDB_DIR, filename)
    if os.path.exists(cache_path):
        return cache_path

    pdb_req = requests.get(url)
    pdb_req.raise_for_status()
    open(cache_path, 'w').write(pdb_req.text)
    # If chain is specified, remove all other chains
    if chain is not None:
        pdb_lines = open(cache_path).readlines()
        open(cache_path, 'w').write(''.join([line for line in pdb_lines if line.startswith('ATOM') and line[21] == chain]))
    return cache_path

def download_smiles_str(pubchem_id: str, retries:int = 2) -> str:
    """Given a pubchem id, get a smiles string"""
    while True:
        req = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{pubchem_id}/property/CanonicalSMILES/CSV")
        smiles_url_csv = req.text if req.status_code == 200 else None
        if smiles_url_csv is not None:
            break
        if retries == 0:
            return None
        time.sleep(1+random())
        retries -= 1

    return smiles_url_csv.splitlines()[1].split(',')[1].strip('"').strip("'") if smiles_url_csv is not None else None




if not PDB_id or not SMILES_or_pubchem_id:
    PDB_id = "7nei"
    SMILES_or_pubchem_id = "CC(=O)C1=CC=C(C=C1)C(=O)OCCOC"
    print(f"No input supplied. Using example data: {PDB_id} and {SMILES_or_pubchem_id}")
# to run many PDB+smiles at once, fill in a list of PDB_files and smiles here...
pdb_files = [download_pdb_file(PDB_id)]
smiless = [download_smiles_str(SMILES_or_pubchem_id) if str(SMILES_or_pubchem_id).isnumeric() else SMILES_or_pubchem_id]

with open("./tmp/input_protein_ligand.csv", 'w') as out:
    out.write("protein_path,ligand\n")
    for pdb_file in pdb_files:
        for smiles in smiless:
            out.write(f"{pdb_file},{smiles}\n")




