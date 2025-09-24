import tomllib
import yaml

# Get pixi config
with open("pixi.toml", "rb") as f:
    pixi = tomllib.load(f)
# get conda env config
with open("conf/env.yml", "r") as f:
    conda = yaml.safe_load(f)
_conda_deps = [i.split(' ')[0] for i in conda['dependencies']]

with open('conf/requirements.txt', 'r') as f:
    for line in f:
        if not line.startswith('git'):
            _conda_deps.append(line.strip().split(' ')[0])

_pixi_deps = list(pixi['dependencies'].keys()) + [i for i in list(pixi['pypi-dependencies'].keys()) if i != 'deeptools']
print(list(pixi['pypi-dependencies'].keys()))
for i in _pixi_deps:
    assert i in _conda_deps, f"{i} in pixi.toml not in conda env.yml"
for i in _conda_deps:
    assert i in _pixi_deps, f"{i} in conda env.yml not in pixi.toml"
