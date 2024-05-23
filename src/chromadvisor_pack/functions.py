"""Recommend the eluant for a chromatography based on the desired molecule."""
import pubchempy as pcp

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import FunctionalGroups
from rdkit.Chem import Crippen
from PIL import ImageTk
#Warning, does not work for.py files
import py3Dmol


def get_smiles(molecule_name):
    """
    Retrieve the SMILES representation of a molecule from its name.

    Parameters
    ----------
    molecule_name : str
        The name of the molecule in English.

    Returns
    -------
    str
        The SMILES representation of the molecule, or None if the molecule is not found.

    Examples
    --------
    >>> get_smiles("water")
    'O'
    """
    try: # Attempt to get compounds from PubChem by name
        results = pcp.get_compounds(molecule_name, 'name')
        if results: # Extract the canonical SMILES representation from the first result
            smiles = results[0].canonical_smiles
            return smiles
        else:
            return None #"Molecule not found. Please try another name."
    except Exception as e:
        return None 


# List of the sub-structures SMARTS representing the unctional groups of molecules
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
    "alcohol": "[OX2H][CX4]",
    "hydroxyl": "[O]", #/!\ attention juste un O !!
}

# Function to find the functional groups of a molecule
def find_functional_groups(smiles):
    """
    Find the functional groups of a molecule given its SMILES representation.

    Parameters
    ----------
    smiles : str
        The SMILES representation of the molecule.

    Returns
    -------
    dict
        A dictionary containing information about the detected functional groups.

    Examples
    --------
    >>> find_functional_groups("CC(=O)O")
    {'carboxylic acid': {'count': 1, 'positions': ((1, 2, 3),)}}
    """
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

        #Acyl halide : 
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

        #Imide :
        if 'imide' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('C(=O)NC(=O)')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)
       
        #Ester :
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
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Enol :
        if 'enol' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('C=CO')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Sulfinic acid :
        if 'sulfinic acid' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('S(=O)O')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Sulfoxide :
        if 'sulfoxide' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('S(=O)')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)
       
        #Aldehyde : 
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

        #Ammonia :
        if 'ammonia' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('N')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Ether :
        if 'ether' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('COC')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Alcohol :
        if 'alcohol' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('CO')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Hydroxyl :
        if 'hydroxyl' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('CO')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Sulfonamide : 
        if 'sulfonamide' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('S(=O)(=O)N')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Nitrile : 
        if 'nitrile' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('CN')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Enamine :
        if 'enamine' in functional_groups :
            pattern_to_remove = Chem.MolFromSmarts('C=CN')
            mol = Chem.DeleteSubstructs(mol, pattern_to_remove)

        #Amine : 
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


def calculate_logp_and_recommend_solvent(smiles):
    """
    Calculate the logP value of a molecule and recommend a solvent for chromatography.

    Args:
        smiles (str): The SMILES representation of the molecule.

    Returns:
        tuple: A tuple containing the logP value and the recommended solvent.

    Examples:
        >>> logp, recommendation = calculate_logp_and_recommend_solvent("CCO")
        >>> print(logp)
        0.17
        >>> print(recommendation)
        Use a mix of DCM/MeOH or a mix of ethyl acetate/hexane.
    """
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None, "Invalid SMILES"
    logp = Crippen.MolLogP(molecule)
    recommendation = ""
    if logp > 3:
        recommendation = "Use an apolar eluent such as : hexane or toluene with a bit of EtOH/acetone."
    elif 0 <= logp <= 3:
        recommendation = "Use a mix of DCM/MeOH or a mix of ethyl acetate/hexane."
    else:  # logP < 0
        recommendation = "Use a polar eluent such as : MeOH or acetone, even water."
    return logp, recommendation


# Function to display the 2D representation of a molecule
def display_molecule_2d(smiles, parent_window):
    """
    Display the 2D representation of a molecule in a Tkinter window.

    Args:
    smiles (str): The SMILES representation of the molecule.
    parent_window: The parent Tkinter window to display the molecule in.

    Returns:
    None

    Examples:
    --------
    >>> display_molecule_2d("CCO", parent_window)
    # This will display the 2D structure of ethanol in the provided Tkinter window.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        messagebox.showerror("Error", "Impossible to convert the SMILES into a molecule.")
        return
    img = Draw.MolToImage(mol)
    img_tk = ImageTk.PhotoImage(img)
    # Afficher l'image dans une étiquette
    molecule_image_label = ttk.Label(parent_window, image=img_tk)
    molecule_image_label.image = img_tk
    molecule_image_label.pack()

# Function to generate the 3D structure of a molecule
def generate_3d_structure(smiles):
    """
    Generate the 3D structure of a molecule and display it in a 3D viewer.

    Args:
    smiles (str): The SMILES representation of the molecule.

    Returns:
    None

    Examples:
    --------
    >>> generate_3d_structure("CCO")
    # This will display the 3D structure of ethanol in a 3D viewer.
    """
    # Convert the SMILES into an RDKIT molecular object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Error : Impossible to convert the SMILES into a molecule.")
        return None
    # Generate the 3D
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)
    # Convert the molecule into a format PDB
    pdb = Chem.MolToPDBBlock(mol)
    # Visualise the molecule in 3D
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(pdb, 'pdb')
    viewer.setStyle({'stick': {}})
    viewer.zoomTo()
    # Show the visualisator in the notebook (sadly not in the interface because it is not possible yet)
    return viewer.show()


# Callback function for the submit button
def on_submit(event=None):
    """
    Callback function for the submit button.

    Retrieves the molecule name from the entry field,
    obtains the SMILES representation of the molecule,
    finds its functional groups, calculates logP, and
    recommends a solvent for chromatography. Displays
    the results in a Tkinter window.

    Args:
    event: The event that triggered the callback (default is None).

    Returns:
    None
    """
    molecule_name = entry.get() # Retrieve the molecule in english from the entry field
    smiles = get_smiles(molecule_name) # Obtain the representation SMILES of the molecule
    if smiles:  # If the representation SMILES is obtained correctly
        functional_groups = find_functional_groups(smiles)  # Fin the functional groups in the molecule
        logp, recommendation = calculate_logp_and_recommend_solvent(smiles)  # Calculate the logP and propose an eluant

        if functional_groups:  # If functional groups are found, they will be showed in a tkinker interface 
            functional_groups_str = ''
            for name, data in functional_groups.items():
                count = len(data["positions"])

                if name in ['ketone', 'phenol']:
                    count = int(count/2)

                positions = data['positions']

                if count != 0:
                    functional_groups_str2 = "\n".join([f"Functional group {name} found {count} times in the molecule."]) #at the positions : {positions}"]) #for name, data in functional_groups.items()])
                    functional_groups_str = '\n'.join([functional_groups_str, functional_groups_str2])
                    # Delete the lines at the beggining and at the enf of the chain
                    functional_groups_str = functional_groups_str.strip()

            # Create a new window to display the informations and the image in 2D of the molecule
            result_window = tk.Toplevel(root)
            result_window.title("Molecule Analysis Result")
            result_window.geometry("900x400")

            # Organise the elements in the window
            # Frame creation for the informations
            info_frame = ttk.Frame(result_window)
            info_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

            # Left : informations on the functional groups, the logP and the recommendations
            num_lines = len(functional_groups_str.split('\n')) #Make the text window as large as the number of functional groups found in the mol
            functional_groups_label = ttk.Label(info_frame, text="Functional groups:")
            functional_groups_label.configure(font=("Helvetica", 14))
            functional_groups_label.pack()
            functional_groups_text = tk.Text(info_frame, height=num_lines, width=50)
            functional_groups_text.insert(tk.END, functional_groups_str)
            functional_groups_text.config(state=tk.DISABLED)
            functional_groups_text.configure(font=("Helvetica", 12))
            functional_groups_text.pack()
            # Afficher la polarité en log(P)
            formatted_logp = "{:.2f}".format(logp) # Format the logP value with two decimal places
            logp_label = ttk.Label(info_frame, text="Log(P):")
            logp_label.configure(font=("Helvetica", 14))
            logp_label.pack()
            logp_text = tk.Text(info_frame, height=1, width=50)
            logp_text.insert(tk.END, formatted_logp)
            logp_text.config(state=tk.DISABLED)
            logp_text.configure(font=("Helvetica", 12))
            logp_text.pack()

            # Show eluent recommendation
            recommendation_label = ttk.Label(info_frame, text="Recommendation :")
            recommendation_label.configure(font=("Helvetica", 14))
            recommendation_label.pack()
            recommendation_text = tk.Text(info_frame, height=1, width=50)
            recommendation_text.insert(tk.END, recommendation)
            recommendation_text.config(state=tk.DISABLED)
            recommendation_text.configure(font=("Helvetica", 12))
            recommendation_text.pack()

            # Right: image in 2D of the molecule
            # Add a “Molecule Display” title above the 2D image
            molecule_display_label = ttk.Label(result_window, text="Molecule Display:")
            molecule_display_label.configure(font=("Helvetica", 14))
            molecule_display_label.pack(side=tk.TOP, anchor="w", padx=110, pady=(20, 0))
            
            # Display the 2D image of the molecule in the interface and 3D in the notebook
            display_molecule_2d(smiles, result_window)
            generate_3d_structure(smiles)
        else:
            messagebox.showinfo("No functional groups found in the molecule.")
    else:
        messagebox.showerror("Error", "Molecule not found. Please try another name.")


# Crate a GUI window
root = tk.Tk()
root.title("Molecule analysis") # Define the title of the window

# Style for the widgets ttk
style = ttk.Style()
style.configure("TButton", font=("Helvetica", 12))
style.configure("TLabel", font=("Helvetica", 14))


# Creating a label and an entry field for entering the SMILES representation
label = ttk.Label(root, text="Enter the name of your desired molecule (in English) :")
label.configure(font=("Helvetica", 14))
label.grid(row=0, column=0, padx=10, pady=5, sticky="w")
entry = ttk.Entry(root, width=50)
entry.grid(row=0, column=1, padx=10, pady=5)

# Creating a submit button to trigger the analysis
submit_button = ttk.Button(root, text="Submit", command=on_submit)
submit_button.grid(row=1, column=0, columnspan=2, padx=10, pady=5)

# Binding the "Enter" key to the submit function
root.bind("<Return>", on_submit)

# Starting the main event loop for the GUI
root.mainloop()

