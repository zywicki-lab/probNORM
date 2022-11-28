<style>
/* requirements table */
td, th, table {
  border: none!important;
 }
 </style>

# Installing probNORM

probNORM can be installed through conda, pypi and from the repository.

!!! tip "Required"

    |              |                                                                                                                                                                                                                                                                                                      |
    |--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
    |[**Python**](https://www.python.org):       |version 3.6 or greater (Python 3 is supported). If you’re setting up Python for the first time,<br>the [Anaconda Python distribution](https://www.anaconda.com/products/distribution) is highly recommended.                                                                                                                                  |
    |**Libraries**:    |[pysam](https://pysam.readthedocs.io/en/latest/installation.html), [numpy](https://numpy.org/install/), [scipy](https://scipy.org/install/)                                                                                                                                                                                                                                                                                   |
    |[**BEDTools**](https://bedtools.readthedocs.io/en/latest/index.html):     |The version is not important, but later versions will have more features so it’s a good idea<br>to get the latest. Follow the instructions at [https://github.com/arq5x/bedtools2](https://github.com/arq5x/bedtools2) to install,<br>and make sure the programs are on your path. That is, you should be able to call bedtools<br>from any directory.|


## Conda installation

This is by far the easiest option. If you’re using the Anaconda Python distribution, then the following will install probNORM:

!!! tip ""

        conda install -c bioconda probNORM

    You can also install  probNORM together with BEDTools via conda:

        conda install -c bioconda probNORM bedtools
    <br>
    Both commands will install probNORM from the bioconda/conda channel and automatically makes sure that dependencies are installed.

Otherwise, read on for installation on other platforms and in other environments.

## PyPi installation

probNORM is on PyPI under the **probnorm** name, so you can install via pip like most Python packages.

!!! tip ""

    Depending on your Python installation, this may require admin rights:

        pip install probnorm

## GitHub development version

1. Get probNORM from repository

    !!! tip ""
        Download by clonning this repository:

            git clone https://github.com/zywicki-lab/probNORM.git

        or download a .zip package and unpack it:

            wget https://github.com/zywicki-lab/probNORM/archive/refs/heads/main.zip

            unzip main.zip

2. Install python dependencies:

    !!! tip ""
        You can install them automatically with the following command:

            pip3 install requirements.txt

        or one-by-one using standard pip installation command, i.e.

            pip3 install pysam numpy scipy

3. Install BEDtools:

    !!! tip ""
        Following command installs mentioned software system-wide.

            sudo apt-get install bedtools