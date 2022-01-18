# Singularity Recipe
These files are meant to be used as templates to build a singularity sif file. Each sos module can then use the sif file as the container in place of the docker.

The command to generate the sif file is:
`singularity build {filename}.sif {filename}.def`

If the previous command returns an error due to an environment issue. The sif file can also be built remotely.
Please login (or create an account by following prompt) by 
`singularity remote login`
And then run
`singularity build --remote {filename}.sif {filename}.def`

Internet connection is required for the remote option.
