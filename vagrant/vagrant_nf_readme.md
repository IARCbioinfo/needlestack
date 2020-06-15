# Needlestack test run using Vagrant and Nextflow

While running Needlestack directly on your machine with Nextflow and Docker is less complicated, the extra step of seting up a Vagrant VM is helpful for testing i.e. different Linux distros.   


 

## Prerequisites

* Vagrant
* VirtualBox

### Setup used

```
uname -a
# Linux amd-linux 5.4.43-1-MANJARO

vagrant version 
# Installed Version: 2.2.9

virtualbox --help
#Oracle VM VirtualBox VM Selector v6.1.8
#<snip>
```

While these are current stable versions for both Vagrant & VirtualBox, it should work with older versions as well.  

## Steps

```
# create a directory for Needlestack Vagrant test

mkdir needlestack_vagrant

# download the test data (size: 403M )

cd needlestack_vagrant
git clone --depth=1 https://github.com/IARCbioinfo/data_test

## Caveat:
## Linking directories/files to a Vagrant VM directory 
## (== shared to /vagrant dir inside VM) 
## will not work (broken links).

# get the Vagrant file from the dev branch of the repository

wget https://raw.githubusercontent.com/IARCbioinfo/needlestack/dev/Vagrantfile

# start the Vagrant machine & ssh into it

vagrant up 
vagrant ssh


# create the run directory
## Caveat: 
## starting nextflow inside /vagrant/data_test
## gives an error:
## "nextflow needs to be executed in a shared file system that supports file locks."
cd
mkdir nf_run_01
cd /home/vagrant/nf_run_01

# create shell file nf_testrun.sh
# with the content below

# make it executable & execute it 
chmod +x nf_testrun.sh
./nf_testrun.sh


```

### nf_testrun.sh

```
#!/bin/sh

data_dir=/vagrant/data_test

nextflow run iarcbioinfo/needlestack -with-docker \
                 --bed $data_dir/BED/TP53_all.bed \
                 --input_bams $data_dir/BAM/BAM_multiple/ \
                 --ref $data_dir/REF/17.fasta \
                 --output_vcf all_variants.vcf
```
