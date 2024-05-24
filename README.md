<div align="center">
  <img src="assets/Chromadvisor_logo.png" alt="Project Logo">
</div>


![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">
Chromadvisor
</h1>

<br>


(Recommend the eluant for a chromatography based on the desired molecule) A CHANGER ABSOLUMENT 

## 👩‍💻 Installation

Create a new environment, you may also give the environment a different name.
```
git clone https://github.com/Squwiddly/Chromadvisor.git
cd Chromadvisor
```

```
conda create -n chromadvisor_pack python=3.10 
```

```
conda activate chromadvisor_pack
(chromadvisor_pack) $ pip install .
```

You will also need to install tkinter, which is the module for the interface :

```
(chromadvisor_pack) $ conda install -c anaconda tk
```

If you need jupyter lab, install it 

```
(chromadvisor_pack) $ pip install jupyterlab
```


## 🔥 Usage

```python
import tkinter as tk
from tkinter import ttk
from src.chromadvisor_pack.functions import on_submit

# Wrapper function to pass entry and root to on_submit
def on_submit_wrapper(event=None):
    on_submit(entry, root)

#Below is the code to create and apply the interface to all the functions
# Crate a GUI window
root = tk.Tk()
root.title("Molecule analysis")

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
submit_button = ttk.Button(root, text="Submit", command=on_submit_wrapper)
submit_button.grid(row=1, column=0, columnspan=2, padx=10, pady=5)

# Binding the "Enter" key to the submit function
root.bind("<Return>", on_submit_wrapper)

# Starting the main event loop for the GUI
root.mainloop()
```
OR, if you don't want the interface :

```python
from src.chromadvisor_pack.functions_without_interface import get_smiles, find_functional_groups, calculate_logp_and_recommend_solvent, display_molecule_2d, generate_3d_structure

molecule_name = input("Enter the name of your desired molecule (in English) :")# Retrieve the molecule in english from the entry field
smiles = get_smiles(molecule_name)
print(smiles)
if smiles:
        functional_groups = find_functional_groups(smiles)
        logp, recommendation = calculate_logp_and_recommend_solvent(smiles)
   
        if functional_groups: # If some functional groups are found, they will be displayed
            functional_groups_str = ''
            for name, data in functional_groups.items():
                count = len(data["positions"])
                if name in ['ketone', 'phenol']:
                    count = int(count/2)
                positions = data['positions']
                if count != 0:
                    print(f"Functional group {name} found {count} times in the molecule.")
        else :
            print("No functional groups found in the molecule.")

        print(f"Log(P): {logp}")
        print("Recommendation:", recommendation)
        # Display the 2D and 3D image of the molecule in the notebook
        display_molecule_2d(smiles)
        generate_3d_structure(smiles) #erase if it does not work on a file.py
else:
        print("Error", "Molecule not found. Please try another name.")
```

## 🛠️ Development installation

Initialize Git (only for the first time). 

Note: You should have create an empty repository on `https://github.com:Squwiddly/Chromadvisor`.

```
git init
git add * 
git add .*
git commit -m "Initial commit" 
git branch -M main
git remote add origin git@github.com:Squwiddly/Chromadvisor.git 
git push -u origin main
```

Then add and commit changes as usual. 

To install the package, run

```
(chromadvisor_pack) $ pip install -e ".[test,doc]"
```

### Run tests and coverage

```
(chromadvisor_pack) $ pip install tox
(chromadvisor_pack) $ tox
(chromadvisor_pack) $ coverage html
```


