
># NOTE:
>## This repository is no longer maintained.
>## For using PAQR, please visit the new repository at [https://github.com/zavolanlab/PAQR2](https://github.com/zavolanlab/PAQR2)

---


# Use PAQR and KAPAC

This repository provides all necessary information and scripts required to use the methods PAQR and/or KAPAC.
Further down you find the link to a self-contained pipeline that runs both tools in combination.


To get you started, clone the repository and change to the directory you specified as target-directory:
```bash
git clone https://github.com/zavolanlab/PAQR_KAPAC.git path/to/workdir
cd path/to/workdir
```

## Prepare the virtual environment
The execution of PAQR and KAPAC depends on specific software (like python, including several packages, or R). With the installation and activation of the virtual environment as it is shown below, you can ensure that both tools run properly.

### Step 1: Download Miniconda3
On Linux:
  ```bash
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh
  ```

On MacOS X:
  ```bash
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
  bash Miniconda3-latest-MacOSX-x86_64.sh
  ```

#### Step 2: Creating a new python3 environment (with all required packages/software)
**Important Update**: Please not that we experienced some difficulties with some anaconda-packages which have been removed from anaconda cloud. In order to create the python3 environment please follow these steps:
  ```bash
  conda create -n paqr_kapac -c bioconda -c ostrokach -c r -c defaults --file requirements_py3_UPDATED.txt
  ```
This will result in an environment without the R package "optparse". In order to install this package please run:
  ```bash
  conda activate paqr_kapac
  R
  ```
Now, you should find your console to be in an active R session. Here, execute:
  ```R
  install.packages("optparse")
  q()
  ```
  
If you already have a conda version with python2 installed, create the new python 3 environment like this:
  ```bash
  conda create -n paqr_kapac -c bioconda -c ostrokach --file requirements_py3_UPDATED.txt python=3
  ```
And follow the steps above to also install the required R package.

The original file `requirements_py3.txt` is deprecated.

With the installation of the environment, the following software is installed as well:
- [Python](https://www.python.org/) (v3.5.1)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) (v3.13.0)
- [NumPy](http://www.numpy.org/) (v1.11.3)
- [Graphviz](http://www.graphviz.org/) (v2.38.0)
- [PyYaml](http://pyyaml.org/) (v3.12)
- [Docutils](http://docutils.sourceforge.net/) (v0.12)
- [STAR](https://github.com/alexdobin/STAR) (v2.5.2a)
- [samtools](http://www.htslib.org/) (v1.3.1)
- [bedtools](http://bedtools.readthedocs.io/en/latest/) (v2.26.0)
- [gzip](http://www.gzip.org/) (v1.7)
- [R](https://www.r-project.org/) (v3.3.2)
- [SRA tools](https://github.com/ncbi/sra-tools) (v.2.8.2)

#### Step 3: Creating a new python2 environment to run PAQR
PAQR relies on packages that were available only for python2 at the time of implementation. To ensure that PAQR runs properly, this python2 environment is installed.
**Important Update**: Please not that due to changes in the anaconda default channels, the packages listed in the `requirements_py2.txt` are not available by default anymore. Therefore, an updated version of the file is now used to create the environment as follows:
  ```bash
  conda create -n py2_paqr python=2 -c bioconda -c anaconda --file requirements_py2_UPDATED.txt
  ```
The original file `requirements_py2.txt` is deprecated.

The following python packages are available in this python2 environment:
- [HTSeq](https://htseq.readthedocs.io/en/release_0.9.1/index.html) (v0.7.2)
- [RseQC](http://rseqc.sourceforge.net) (v2.6.4)

#### Step 4: Activate the environment
(only the python3 environment needs to be activated; the py2_paqr environment is only used internally)
  ```bash
  source activate paqr_kapac
  ```
To exit the environment (after finishing the usage of the pipeline), just execute
  ```bash
  source deactivate
  ```

