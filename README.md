# THIS REPO IS DEPRECATED. PLEASE USE [THIS ONE](https://github.com/sicotteh/CAVA) INSTEAD
#  
#  
#  
#  
#  

* [CAVA README](#cava-readme)
    * [1 INTRODUCTION](#1-introduction)
    * [2 PUBLICATION](#2-publication)
    * [3 DEPENDENCIES](#3-dependencies)
    * [4 INSTALLATION ON LINUX OR MAC](#4-installation-on-linux-or-mac)
    * [5 RUNNING CAVA](#5-running-cava)
    * [6 LICENCE](#6-licence)


CAVA README
==================

1 INTRODUCTION
--------------

CAVA (Clinical Annotation of VAriants) is a lightweight, fast, flexible 
and easy-to-use Next Generation Sequencing (NGS) variant annotation tool. 
It implements a clinical sequencing nomenclature (CSN), a fixed variant 
annotation consistent with the principles of the Human Genome Variation 
Society (HGVS) guidelines, optimised for automated clinical variant 
annotation of NGS data. 

Since 2017, CAVA has been maintained by a group of bioinformaticians involved 
in both research and clinical genomics, adding several key functionalities, enhancements, 
and general support along the way.

2 PUBLICATION
-------------

If you use CAVA, please cite:

Márton Münz, Elise Ruark, Anthony Renwick, Emma Ramsay, Matthew Clarke, 
Shazia Mahamdallie, Victoria Cloke, Sheila Seal, Ann Strydom, 
Gerton Lunter, Nazneen Rahman. CSN and CAVA: variant annotation tools 
for rapid, robust next-generation sequencing analysis in the clinical 
setting. Genome Medicine 7:76, doi:10.1186/s13073-015-0195-6 (2015).

Maybe some day, we'll get around to publishing what we've done to rescue this abandoned project.

3 DEPENDENCIES
--------------

To install and run CAVA you will need the following dependencies installed:
- Python 3
- GCC and GNU make
- virtualenv

It just makes sense to keep things in a virtualenv. Here's how you do it if you are
unfamiliar.

```bash 
pip install virtualenv
virtualenv cava
source cava/bin/activate
```
At this point, your terminal should change to let you know you are in a virtual environment.

4 INSTALLATION ON LINUX OR MAC
------------------------------

```bash 
pip install cava

# - or -
git clone git@github.com:Steven-N-Hart/CAVA.git
# optional to checkout release
# e.g. git checkout v.1.2.4
python setup.py install
```

If you get an error with pycurl, run the following command before running the setup.py
```bash
pip uninstall pycurl
export PYCURL_SSL_LIBRARY=nss
pip install --compile --install-option="--with-nss" --no-cache-dir pycurl  
```


5 RUNNING CAVA
--------------

Before using CAVA, you will need to create a database of transcripts for which to base your annotations from.
Details can be found in [this README](cava/ensembldb/README.md). In short, we reccomend using MANE transcripts, 
so to get started, you would simply:
```bash
# Download GTF files for either RefSeq or ENSEMBLE
wget -O data/ENST.gtf.gz ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_0.91/MANE.GRCh38.v0.91.select_ensembl_genomic.gtf.gz and
wget -O data/RefSeq.gtf.gz ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_0.91/MANE.GRCh38.v0.91.select_refseq_genomic.gtf.gz

# Separate into ENST and NM Transcripts
zcat data/ENST.gtf.gz |cut -f9|cut -f4 -d' '|grep ENST|sed 's/;//;s/\"//g'|sort -u > data/ENST.txt
zcat data/RefSeq.gtf.gz |cut -f9|cut -f4 -d' '|grep "NM_"|sed 's/;//;s/\"//g'|sort -u > data/RefSeq.txt

# Look at the options and configure appropriately
python3 MANE.py -h

```

CAVA can be run with the following simple command:

```bash
python3 CAVA.py -c config.txt -i input.vcf -o output
```

It requires three command line arguments: 
the name of the configuration file (-c), the name of the input file (-i) 
and the prefix of the output file name (-o). 

6 LICENCE
---------

CAVA is released under MIT licence (see the LICENCE file).

