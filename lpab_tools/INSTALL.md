# Installation guide

The commands below target a Linux workstation or an HPC cluster. Module names
such as GCC and CUDA differ between clusters; ask the local administrator when
a listed module is unavailable.

## 1. Ladderpath (`lppack`)

The current release supports Python 3.10 or newer. The project requirement file
pins the official revision tested with this workflow:

```bash
python -m pip install -r lpab_tools/requirements.txt
python -c "import lppack; print('lppack is ready')"
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

Start with the current [official AlphaFold 3 installation guide](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md).
The official guide currently targets Linux, a large sequence database, and a
supported NVIDIA GPU. Model parameters require a separate application.

The following condensed native-conda recipe records the original laboratory
setup. Package versions can become incompatible, so compare it with the
official guide before reproducing it:

```bash
conda create -n AF3_3.11 python=3.11 -y
conda activate AF3_3.11
module load gcc/12.1.0-gcc-9.4.0-jtv
module load cuda/12.0.1-gcc-12.1.0-kof
module load cmake/3.23.1-gcc-9.4.0-wx2

git clone https://github.com/google-deepmind/alphafold3.git ~/git_develop/af3
cd ~/git_develop/af3
mamba install -c bioconda hmmer -y
pip install -r dev-requirements.txt
pip install --upgrade "jax[cuda12]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
pip install . --no-deps --verbose
build_data
python run_alphafold_test.py --model_dir=./alphafold3_models
```

Download the official databases listed in the AlphaFold 3 guide and place the
approved model parameters in a private model directory. Do not commit model
parameters or databases to this repository.

Verify the GPU:

```bash
python -c "import jax; print(jax.devices())"
```

Finally, put `run_alphafold.py`, the model directory, and the database directory
paths into the top-level `parameters.txt`.
