### Generate new recipes using HPC Container Maker (HPCCM)

    $ hpccm --recipe recipe_<variant>.py --format singularity --singularity-version=3.2 > Singularity.<variant>

### Build Singularity image locally (must be done from project root directory)

    $ sudo singularity build <image_name>.sif recipes/Singularity.<variant>

### Pull Singularity image from GitHub Container Registry

Latest `master` version (here OpenMP variant):

    $ singularity pull oras://ghcr.io/MRChemSoft/mrchem-nompi:latest

Tagged version (here MRChem-v1.0.2 using OpenMPI-v4.0):

    $ singularity pull oras://ghcr.io/MRChemSoft/mrchem-openmpi4.0:v1.0.2

### Run Singularity container (non MPI)

    $ singularity exec <image-name>.sif mrchem molecule

### Run Singularity container (MPI)

    $ singularity exec <image-name>.sif mrchem -D molecule
    $ mpirun singularity exec <image-name>.sif mrchem.x molecule.json
