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

If fakeroot didn't work and gave the error message of:
1. `FATAL:   could not use fakeroot: no mapping entry found in /etc/subuid for {username}`

   Check if your username is in /etc/subuid

2. `ERROR  : Failed to create container namespaces`

   For debian system, do `sudo echo 1 > /proc/sys/kernel/unprivileged_userns_clone`
   
3. `FATAL:   Unable to create build: while searching for mksquashfs: exec: "mksquashfs": executable file not found in $PATH`

   Do `sudo apt-get update && sudo apt-get install squashfs-tools`
