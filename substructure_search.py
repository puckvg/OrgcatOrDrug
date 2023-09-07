import pandas as pd
from rdkit import Chem
def get_substructures():
    smiles_df = pd.read_csv('data/SMILES_Substructure_search.csv')
    smiles_list = smiles_df['SMILES'].to_list()

    smiles_list_oscar = smiles_df[smiles_df['Motif_number_OSCAR'] != 'na']['SMILES'].to_list()
    return smiles_list, smiles_list_oscar
def get_drugs_organocats():
    df = pd.read_csv("data/oscar_plus_drugbank.csv")
    drug_smiles = df[df['is_drug'] == 1]['smiles'].to_list()
    organocats_smiles = df[df['is_drug'] == 0]['smiles'].to_list()
    return drug_smiles, organocats_smiles

def substruct_search(mol, fragments_list):
    # check if there are any matches in frag list
    bools = [mol.HasSubstructMatch(patt) for patt in fragments_list]
    res = any(bools)
    print(res)
    return res

if __name__ == "__main__":
    smiles_list, smiles_list_oscar = get_substructures()
    drugs_list, organocats_list = get_drugs_organocats()
    organocat_mols = [Chem.MolFromSmiles(x) for x in organocats_list]

    # now substructure search
    # TODO oscar list only after
    substruct_mols = [Chem.MolFromSmiles(x) for x in smiles_list]
    substruct_smarts = [Chem.MolToSmarts(m) for m in substruct_mols]
    substruct_mols = [Chem.MolFromSmarts(s) for s in substruct_smarts]
    organocats_bools = [substruct_search(org, substruct_mols) for org in organocat_mols]


    results = {"organocats":organocats_list, "match substructures?":organocats_bools}
    results_df = pd.DataFrame(results)
    results_df.to_csv('data/organocats_substructure_match.csv')
