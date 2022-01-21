# Automatically generated singularity configurations

Warning: all files in this folder are generated from docker configurations under `../docker` folder. **Please do not manually edit contents in this folder**. To generate them, run

```
./release docker_to_singularity
```

under the root directory of this repo.

## Build singularity images

```
singularity build --fakeroot {filename}.sif {filename}.def
```

It will generate the `sif` file in the current directory.
