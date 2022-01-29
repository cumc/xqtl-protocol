# Automatically generated singularity configurations

Warning: all files in this folder are generated from docker configurations under `../docker` folder. **Please do not manually edit contents in this folder**. To generate them, run

```
./release.sos docker_to_singularity
```


## Build singularity images

```
./release.sos singularity --config bioinfo.def --out-dir "./"
```

under the root directory of this repo. It will generate the `sif` file in the directory specified by `--out-dir`.

### Potential issues

If `--fakeroot` did not work with error message of one of the following:

1. `FATAL:   could not use fakeroot: no mapping entry found in /etc/subuid for {username}`
   - Solution: check if your username is in `/etc/subuid` by `cat /etc/subuid`. If not, add it in

2. `ERROR  : Failed to create container namespaces`
   - Solution: for debian system, do `sudo echo 1 > /proc/sys/kernel/unprivileged_userns_clone`
   
3. `FATAL:   Unable to create build: while searching for mksquashfs: exec: "mksquashfs": executable file not found in $PATH`
   - Solution: for debian system, do `sudo apt-get update && sudo apt-get install squashfs-tools`

If you are on a HPC please ask your admin to read these suggestions and implement them for your system.
