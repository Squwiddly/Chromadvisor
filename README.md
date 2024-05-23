<div align="center">
  <img src="assets/banner.png" alt="Project Logo">
</div>


![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">
Chromadvisor
</h1>

<br>


(Recommend the eluant for a chromatography based on the desired molecule) A CHANGER ABSOLUMENT 

## 🔥 Usage

```python
from chromadvisor.functions import get_smiles, find_functional_groups, calculate_logp_and_recommend_solvent, display_molecule_2d, on_submit

#Obtain the smiles of an english name molecule
SMILES = get_smiles(?)

#Find the functional groups of the submitted molecule
functional_groups = find_functional_groups(SMILES)

#Display the 2D visualization of the molecule
2d_mol = display_molecule_2d(SMILES)

#Give the log(P) and the eluent recommended for a chromatography
logP_&_recommendation = calculate_logp_and_recommend_solvent(SMILES)

#Use it
USING THE INTERFACE ?
```

## 👩‍💻 Installation

Create a new environment, you may also give the environment a different name. 

```
conda create -n chromadvisor_pack python=3.10 
```

```
conda activate chromadvisor_pack
(conda_env) $ pip install .
```

If you need jupyter lab, install it 

```
(chromadvisor_pack) $ pip install jupyterlab
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
(conda_env) $ pip install tox
(conda_env) $ tox
```



