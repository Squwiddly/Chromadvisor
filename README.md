div style="text-align:center;">
![Project Logo](assets/banner.png)
</div>

![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">
Chromadvisor
</h1>

<br>


Recommend the eluant for a chromatography based on the desired molecule

## 🔥 Usage

```python
from mypackage import main_func

# One line to rule them all
result = main_func(data)
```

This usage example shows how to quickly leverage the package's main functionality with just one line of code (or a few lines of code). 
After importing the `main_func` (to be renamed by you), you simply pass in your `data` and get the `result` (this is just an example, your package might have other inputs and outputs). 
Short and sweet, but the real power lies in the detailed documentation.

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



