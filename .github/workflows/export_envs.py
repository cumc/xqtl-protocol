import os
import sys
import re
from pathlib import Path
from functools import reduce
from operator import concat
import yaml
import pandas as pd
import numpy as np

class IndentDumper(yaml.Dumper):
    def increase_indent(self, flow=False, indentless=False):
        return super(IndentDumper, self).increase_indent(flow, False)

container_table = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "../../container/containers.csv"))
container_table_universal = container_table[(container_table.Environment == "universal")]
container_table_envs = container_table[(container_table.Environment != "universal")]

# Hard code channels in the correct order
channels = ["dnachun", "conda-forge", "bioconda", "nodefaults"]

environment_list = [re.sub(r'\s+', "", env).split(",") for env in container_table_envs.Environment]
environments = list(np.unique(reduce(concat, environment_list, [])))

package_dict = {package: version for package, version in zip(container_table_envs.Package, container_table_envs.Version)}
environment_dict = {package: environment for package, environment in zip(container_table_envs.Package, container_table_envs.Environment)}
environment_package_dict = {environment: {} for environment in environments}
for environment in environments:
    for package in sorted(environment_dict.keys()):
        if environment in environment_dict[package]:
            environment_package_dict[environment][package] = package_dict[package]
for environment in environment_package_dict.keys():
    if environment != "dapars2":
        universal_packages = sorted(set(container_table_universal.Package) - set(environment_package_dict[environment].keys()))
        for package in universal_packages:
            environment_package_dict[environment][package] = np.nan

yaml_dict = {environment: {"name": str(environment), "channels": channels, "dependencies": []} for environment in environment_package_dict.keys()}
for environment in yaml_dict.keys():
    for package in environment_package_dict[environment].keys():
        version = environment_package_dict[environment][package]
        if pd.notna(version):
            yaml_dict[environment]["dependencies"].append(f"{package}={version}")
        else:
            yaml_dict[environment]["dependencies"].append(package)

for environment in yaml_dict.keys():
    yml_file = f"../../container/{environment}/{environment}.yml"
    yml_path = Path(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), yml_file))
    yml_path.parent.mkdir(parents=True, exist_ok=True)
    with open(yml_path, 'w') as file_handle:
        yaml.dump(yaml_dict[environment], file_handle, sort_keys = False, Dumper=IndentDumper)
