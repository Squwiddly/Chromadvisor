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

#Warning, py3Dmol does only work on .ipynb files
import py3Dmol

#For the code without the interface 
from IPython.display import display, Image
