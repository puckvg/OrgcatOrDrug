import pandas as pd
from rdkit import Chem
def get_substructures():
    smiles_df = pd.read_csv('data/SMILES_Substructure_search.csv')
    smiles_list = smiles_df['SMILES'].to_list()
    smarts_list = smiles_df['SMARTS'].to_list()

    assert len(smiles_list) == len(smarts_list)

   # smiles_list_oscar = smiles_df[smiles_df['Motif_number_OSCAR'] != 'na']['SMILES'].to_list()
    return smiles_list, smarts_list
def get_drugs_organocats():
    df = pd.read_csv("data/oscar_plus_drugbank.csv")
    drug_smiles = df[df['is_drug'] == 1]['smiles'].to_list()
    organocats_smiles = df[df['is_drug'] == 0]['smiles'].to_list()
    return drug_smiles, organocats_smiles

def isNaN(string):
    return string != string
def substruct_search(mol, fragments_list):
    # check if there are any matches in frag list
    bools = [mol.HasSubstructMatch(patt) for patt in fragments_list]
    res = any(bools)
    return res

def get_substruct_mols(smiles_list):
    mols_list = []
    for smi in smiles_list:
        if isNaN(smi):
            mol = None
        else:
            mol = Chem.MolFromSmiles(smi)
        mols_list.append(mol)
    return mols_list

def get_substruct_smarts(smiles_list, smarts_list):
    substruct_mols = get_substruct_mols(smiles_list)
    f_smarts_list = []
    for i, mol in enumerate(substruct_mols):
        if mol is None:
            smarts = smarts_list[i]
        else:
            smarts = Chem.MolToSmarts(mol)
        f_smarts_list.append(smarts)
    return f_smarts_list

if __name__ == "__main__":

    drugs_list, organocats_list = get_drugs_organocats()
    organocat_mols = [Chem.MolFromSmiles(x, sanitize=True) for x in organocats_list]

    # now substructure search
    smiles_list, smarts_list = get_substructures()
    substruct_smarts = get_substruct_smarts(smiles_list, smarts_list)

    #print(substruct_smarts)
    substruct_mols = [Chem.MolFromSmarts(s) for s in substruct_smarts]
    organocats_bools = [substruct_search(org, substruct_mols) for org in organocat_mols]

    false_organocats = []
    for i, bool in enumerate(organocats_bools):
        if not bool:
            organocat = organocats_list[i]
            false_organocats.append(organocat)

    with open('false_organocats.txt', 'w') as f:
        for organocat in false_organocats:
            f.write(organocat)
            f.write('\n')

    results = {"organocats":organocats_list, "match substructures?":organocats_bools}
    results_df = pd.DataFrame(results)
    results_df.to_csv('data/organocats_substructure_match.csv')

    drugs_mols = [Chem.MolFromSmiles(x, sanitize=True) for x in drugs_list]
    drugs_bools = [substruct_search(drug, substruct_mols) for drug in drugs_mols]

    true_drugs = []
    for i, bool in enumerate(drugs_bools):
        if bool:
            drug = drugs_list[i]
            true_drugs.append(drug)

    with open('true_drugs.txt', 'w') as f:
        for drug in true_drugs:
            f.write(drug)
            f.write('\n')

    results = {"drugs":drugs_list, "match_substructures?":drugs_bools}
    results_df = pd.DataFrame(results)
    results_df.to_csv('data/drugs_substructure_match.csv')


