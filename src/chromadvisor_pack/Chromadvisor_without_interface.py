#First the convertion to the english name of  molecule into its SMILES
def get_smiles(molecule_name):
    try: # Attempt to get compounds from PubChem by name
        results = pcp.get_compounds(molecule_name, 'name')
        if results: # Extract the canonical SMILES representation from the first result
            smiles = results[0].canonical_smiles
            return smiles
        else:
            return "Molecule not found. Please try another name."
    except Exception as e:
        return "An error occurred: {}".format(str(e)) # Print an error message if an exception occurs


#Then find the functional groups of a molecule and show it on 2D

# List of the sub-structure SMARTS representing the functional groups of molecules
functional_group_smarts = {
    "anhydride": "[CX3](=[OX1])[OX2][CX3](=[OX1])",
    "carboxylic acid": "[CX3](=O)[OX2H1]",
    "carboxylic acid": "[CX3](=O)[OX1H0-,OX2H1]",
    "acyl halide":"[CX3](=[OX1])[F,Cl,Br,I]",
    "imide": "[CX3](=[OX1])[NX3H][CX3](=[OX1])",
    "imide": "[CX3](=[OX1])[NX3H0]([#6])[CX3](=[OX1])",
    "imide": "[CX3](=[OX1])[NX3H0]([NX3H0]([CX3](=[OX1]))[CX3](=[OX1]))[CX3](=[OX1])",
    "ester": "[#6][CX3](=O)[OX2H0][#6]",
    "amide": "[NX3][CX3](=[OX1])[#6]",
    "phenol": "[OX2H][cX3]:[c]",
    "enol": "[OX2H][#6X3]=[#6]",
    "sulfinic acid": "[$([#16X3](=[OX1])[OX2H,OX1H0-]),$([#16X3+]([OX1-])[OX2H,OX1H0-])]",
    "sulfinate": "[$([#16X3](=[OX1])[OX2H0]),$([#16X3+]([OX1-])[OX2H0])]",
    "sulfoxide": "[$([#16X3](=[OX1])([#6])[#6]),$([#16X3+]([OX1-])([#6])[#6])]",
    "sulfoxide": "[$([#16X3]=[OX1]),$([#16X3+][OX1-])]",

    "aldehyde": "[CX3H1](=O)[#6]",
    "ketone": "[CX3](=O)[#6]",

    "ammonia": "[NX3][CX3]=[NX3+]",
    "ether": "[OD2]([#6])[#6]",

    "sulfonamide": "[$([#16X4]([NX3])(=[OX1])(=[OX1])[#6]),$([#16X4+2]([NX3])([OX1-])([OX1-])[#6])]",
    "nitrile": "[NX1]#[CX2]",
    "enamine": "[NX3][CX3]=[CX3]",
    "enamine": "[NX3][$(C=C),$(cc)]",
    "amine": "[NX3;H2,H1;!$(NC=O)]",
    "amine": "[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6]",
    "amine": "[NX3;H2,H1;!$(NC=O)].[NX3;H2,H1;!$(NC=O)] ",
    "imine": "[CX3;$([C]([#6])[#6]),$([CH][#6])]=[NX2][#6]",
    "imine": "[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]",
    "thiol": "[#16X2H] ",
    "sulfide": "[SX2]",
    "sulfide": "[#16X2H0]",

    #Est ce qu'on garde alcool ET hydroxyl ? on peut regrouper les smarts
    "alcohol": "[OX2H][CX4]",
    "hydroxyl": "[O]", #/!\ attention juste un O !!
}

# Function to find the functional groups of a molecule
def find_functional_groups(smiles):

    functional_groups = {}
    functional_groups2 = {}

    # Initialize a dictionary to store the number of each functional group
    functional_groups_count = {key: 0 for key in functional_group_smarts.keys()}

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Error : Impossible to convert the SMILES into a molecule.")
        return None

   
    for name, smarts in functional_group_smarts.items():

        #Eliminate the functionnal group of the rdkit molecular rdkit once the fuction has detected it
        #To avoid repetition. /!\ problèmes si cycle !!

        #Anhydride : Delete the functionnal group once it has been detected to avoid the function thinking it has anhydride AND ketone AND ester etc..
        if 'anhydride' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('C(=O)OC(=O)') #Select the pattern to remove from the smiles of the molecule
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove) #Delete the pattern

        #Carboxylic acid :
        if 'carboxylic acid' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('C(=O)O')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Acyl halide : a simplifier mais là j'arrive pas. est ce que faut différencier les différents types d'acyl halide (avec Br, I, ...) quand on les compte ??
        if 'acyl halide' in functional_groups :
            pattern_to_remove = 'BrC(=O)'
            pattern_to_remove = Chem.MolFromSmarts(pattern_to_remove)
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)
            pattern_to_remove = 'ClC(=O)'
            pattern_to_remove = Chem.MolFromSmarts(pattern_to_remove)
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)
            pattern_to_remove = 'FC(=O)'
            pattern_to_remove = Chem.MolFromSmarts(pattern_to_remove)
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)
            pattern_to_remove = 'IC(=O)'
            pattern_to_remove = Chem.MolFromSmarts(pattern_to_remove)
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Imide : /!\ fonctionne pas. pas détecté ("2x amide")
        if 'imide' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('C(=O)NC(=O)')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)
       
        #Ester :  supprimé par a.carboxylique donc non détecté quand ya les 2
        if 'ester' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('C(=O)O')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Amide :
        if 'amide' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('C(=O)N')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Phenol :
        if 'phenol' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('c1ccc(O)cc1')
        #    pattern_to_remove = Chem.MolFromSmarts('[c,C]1:[c,C]:[c,C]:[c,C]:[c,C]:1-[OH0]')
        #    pattern_to_remove = Chem.MolFromSmarts('[c,C][c,C][c,C][c,C][c,C]O')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Enol :
        if 'enol' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('C=CO')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Sulfinic acid :
        if 'sulfinic acid' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('S(=O)O')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #je crois ca détecte pas amine (?!)

        #Sulfinate : jsp comment faire la diff avec sulfinic acid
        if 'sulfinate' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('S(=O)O')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Sulfoxide : utile ?
        if 'sulfoxide' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('S(=O)')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)
       
        #Aldehyde : probleme avec C(=O)C(=O) /!\
        if 'aldehyde' in functional_groups :
            # Define SMARTS to identify aldehydes and adjacent atoms
            aldehyde_pattern = Chem.MolFromSmarts('[CX3H1](=O)[#6][!#1]')

            # Find correspondences in molecule
            matches = mol.GetSubstructMatches(aldehyde_pattern)

            # Browse correspondences
            for match in matches:
                # Browse atoms in correspondence
                for atom_idx in match:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    # If atom not C, H or O, mark it for suppression
                    if atom.GetSymbol() not in ['H', 'C', 'O']:
                        atom.SetAtomicNum(0)

            # Delete atoms marked for suppression
            mol = Chem.DeleteSubstructs(mol, Chem.MolFromSmarts("[#0]"))
           
        #Ketone :
        if 'ketone' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('C(C=O)C')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Ammonia : ammonia pas détecté ??
        if 'ammonia' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('N')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Ether : utile ?
        if 'ether' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('COC')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Alcohol : utile?
        if 'alcohol' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('CO')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Hydroxyl : capte pas quand ca détecte alcohol et quand hydroxyl
        if 'hydroxyl' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('CO')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Sulfonamide : utile?
        if 'sulfonamide' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('S(=O)(=O)N')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Nitrile : utile ?
        if 'nitrile' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('CN')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Enamine : utile? pas testé
        if 'enamine' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('C=CN')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Amine : utile? meme que ammonia??? pas détecté
        if 'amine' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('N')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Thiol :
        if 'thiol' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('CS')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        pattern = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(pattern)
        functional_groups_count[name] = {'count' : len(matches), 'positions' : matches}
       
        if matches:
            #functional_groups[name] = matches
            functional_groups[name] = {
                "count": len(matches),
                "positions": matches
            }
               
    functional_groups2.update(functional_groups)

    return functional_groups2

#Polarity and log(P)
def calculate_logp_and_recommend_solvent(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None, "Invalid SMILES"
    logp = Crippen.MolLogP(molecule) # Calculate the log(P) using a RDKit module
    recommendation = ""
    if logp > 3:
        recommendation = "Use an apolar eluent such as : hexane or toluene with a bit of EtOH/acetone."
    elif 0 <= logp <= 3:
        recommendation = "Use a mix of DCM/MeOH or a mix of ethyl acetate/hexane."
    else:  # logP < 0
        recommendation = "Use a polar eluent such as : MeOH or acetone, even water."
    return logp, recommendation


#Draw of molecule 2D + 3D
def display_molecule_2d(smiles): # Function to display a representation 2D  of the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Error", "Impossible to convert the SMILES into a molecule.")
        return
    img = Draw.MolToImage(mol)
    display(img)


def generate_3d_structure(smiles):
    # Convert the SMILES in an molecular RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Error : Impossible to convert the SMILES into a molecule.")
        return None
    # Generate a 3D conformation
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)
    # Convertir la molécule en format PDB
    pdb = Chem.MolToPDBBlock(mol)
    # Visualise the molecule in 3D
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(pdb, 'pdb')
    viewer.setStyle({'stick': {}})
    viewer.zoomTo()
    # Display the visualisator in the notebook
    return viewer.show()


def on_submit(event=None):
    molecule_name = entry.get() # Retrieve the molecule in english from the entry field
    smiles = get_smiles(molecule_name) # Obtain the representation SMILES of the molecule
    if smiles:  # If the SMILES was obtained with success
        functional_groups = find_functional_groups(smiles)  # Find the functional groups of the molecule
        logp, recommendation = calculate_logp_and_recommend_solvent(smiles)  # Calculate the logP and propose an eluent

        if functional_groups:  # If the functional groups are found, they will be displayed
            functional_groups_str = ''
            for name, data in functional_groups.items():
                count = len(data["positions"])

                if name in ['ketone', 'phenol']:
                    count = int(count/2)

                positions = data['positions']

                if count != 0:
                    functional_groups_str2 = "\n".join([f"Functional group {name} found {count} times in the molecule."]) #at the positions : {positions}"]) #for name, data in functional_groups.items()])
                    functional_groups_str = '\n'.join([functional_groups_str, functional_groups_str2])
                    # Delete the hollow lines at the beggining and end of the chain
                    functional_groups_str = functional_groups_str.strip()

        else:
            print("No functional groups found in the molecule.")
    else:
        print("Error", "Molecule not found. Please try another name.")


# Using it
molecule_name = input("Enter the name of your desired molecule (in English) :")
smiles = get_smiles(molecule_name)
print(smiles)
functional_groups = find_functional_groups(smiles)
if functional_groups:
        print("Functional groups of the molecule:")
        for name, data in functional_groups.items():
            count = len(data["positions"])
            print(f"Functional group {name} found {count} times in the molecule.")

        # Calculate logP and recommend eluent
        logp, recommendation = calculate_logp_and_recommend_solvent(smiles)
        print(f"Log(P): {logp}")
        print("Recommendation:", recommendation)
else:
        print("No functional groups found in the molecule.")

# Display the molecule in 2D and in 3D in the notebook
display_molecule_2d(smiles)
generate_3d_structure(smiles)
