# Workflow to analyze label-free data with the compomics pipeline

This workflow is based on a docker image which is available at the docker hub (https://hub.docker.com/repository/docker/veitveit/compomicsworkflow) and is downloaded automatically.
If you want to build your own docker image, be aware that it needs to be named veitveit/compomicsworkflow:dev or change the name in the configuration file. 

## Content
_Nextflow folder:_ Implementation of workflow

_Results folder:_ Results from the running the UPS data set

## Getting started

Here is a little script that can be used to test the workflow execution.
This script assumes that Docker is available on your system and is targeting Debian/Ubuntu linux distributions.

```
# Install Java
sudo apt-get install openjdk-8-jdk

# Fix docker permission issue (https://stackoverflow.com/questions/48957195/how-to-fix-docker-got-permission-denied-issue)
# Note: this allows to run docker as non-root user
sudo usermod -aG docker $USER
newgrp docker

# Fetch scripts and data
git clone https://github.com/veitveit/IS_Benchmarking.git
cd ./Compomics-Workflow/Nextflow/data

# Nextflow
cd ..
curl -s https://get.nextflow.io | bash

# Worfklow
./nextflow run main.nf -profile docker,test

```

## Run benchmarking data set

You should have successfully tested the workflow using the procedure above.

Download the raw files from PRIDE: http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD001819

Run the workflow with the following command and parameters after changing _RAWFOLDER_ to the folder where the raw files are located. You will also need to place the files _yeast_UPS.fasta_ and _pxd001819.txt_ into in the current folder. These files is given in the Results folder: https://github.com/veitveit/IS_Benchmarking/tree/master/Compomics/Results

Also adjust the parameter values _max_cpus_ and _max_memory_ to the computing power you have available.
```
nextflow run main.nf --raws 'RAWFOLDER/*.raw' --fasta yeast_UPS.fasta --miscleavages 2 --fragment_mass_tolerance 0.8 \
--precursor_mass_tolerance 5 --enzyme 'Trypsin (no P rule)' --variable_mods 'Oxidation of M,Acetylation of protein N-term' \
--experiment_design pxd001819.txt --max_cpus 8 --max_memory \
8GB -profile docker -with-report -with-trace -with-timeline


```


