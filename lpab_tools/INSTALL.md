# Installation guide

The commands below target a Linux workstation or an HPC cluster. Module names
such as GCC and CUDA differ between clusters; ask the local administrator when
a listed module is unavailable. Run the commands from the repository root.

This workflow normally uses two Python environments:

1. the `lpab` environment for this repository, Ladderpath, and ipSAE;
2. the AlphaFold 3 environment, which has its own tightly coupled scientific
   dependencies.

The two environments can be combined only if the resulting dependency set is
known to work. The pipeline invokes AlphaFold 3 directly as
`[alphafold_python, alphafold_script, ...]`, so a Docker-only AlphaFold 3
installation is not directly usable without adding a wrapper command.

## 1. Ladderpath (`lppack`)

The current release supports Python 3.10 or newer. The project requirement file
pins the official revision tested with this workflow:

```bash
python -m pip install -r lpab_tools/requirements.txt
python -c "import lppack, numpy; print('LPAB dependencies are ready')"
```

See the [official lppack repository](https://github.com/yuernestliu/lppack) for
newer versions and its full user guide.

The Graphviz system program is only needed for drawing laddergraphs:

```bash
mamba install -c conda-forge graphviz python-graphviz
dot -V
```

## 2. Clustal Omega

```bash
mamba install -c bioconda clustalo -y
clustalo --help
```

## 3. Rosetta 3.14

Download the licensed academic Linux package from the
[RosettaCommons download page](https://www.rosettacommons.org/downloads/academic/3.14/),
then unpack and compile it. The compiler module below is an example from the
original laboratory environment.

```bash
tar -xvf rosetta.binary.linux.release-371.tar.gz
module load gcc/13.1.0-gcc-9.4.0-hye
cd rosetta.binary.linux.release-371/main/source
mamba install scons -y
./scons.py -j10 mode=release bin
```

Add the binary directory to `PATH` (replace the example path):

```bash
export PATH="$PATH:/path/to/rosetta.binary.linux.release-371/main/source/bin"
antibody.default.linuxgccrelease -help
```

The pipeline also needs BLASTP:

```bash
mamba install -c bioconda blast -y
blastp -version
```

Rosetta is third-party licensed software and is not redistributed here.

## 4. MMseqs2

Install from Bioconda as described by the
[official MMseqs2 project](https://github.com/soedinglab/MMseqs2):

```bash
mamba install -c conda-forge -c bioconda mmseqs2 -y
mmseqs -h
```

The pipeline runs the antibody settings used in the workflow:

```bash
mmseqs createdb heavy_light.fasta sequence_db
mmseqs cluster sequence_db clusters tmp --min-seq-id 0.9 -c 0.8 --cov-mode 2
mmseqs createtsv sequence_db sequence_db clusters cluster_map.tsv
```

## 5. AlphaFold 3

Start with the [official AlphaFold 3 installation guide](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md).
The current official requirements include Linux, an NVIDIA GPU with a
supported compute capability, substantial disk space for the databases, and
model parameters obtained under AlphaFold 3's terms of use. The guide's
commands and supported versions can change, so do not treat an old lab module
recipe as a universal installation recipe.

After validating the AlphaFold 3 installation, set these four values in the
top-level `parameters.txt`:

```text
alphafold_python = /path/to/alphafold3/.venv/bin/python
alphafold_script = /path/to/alphafold3/run_alphafold.py
alphafold_model_dir = /path/to/alphafold3_models
alphafold_database_dir = /path/to/alphafold3_databases
```

The first path must point to the Python interpreter that can import the
installed AlphaFold 3 package. The last two paths must contain the approved
model parameters and the downloaded databases. Do not commit model parameters
or databases to this repository.

Verify the GPU:

```bash
python -c "import jax; print(jax.devices())"
```

Run that check in the AlphaFold 3 environment. Run `python run_lpab.py
--check` in the LPAB environment; it verifies the configured files and
executables, but an AlphaFold 3 import test still requires the selected AF3
interpreter and should be run separately when troubleshooting.
