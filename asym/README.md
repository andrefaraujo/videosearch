# Asymmetric Image Comparisons

This directory contains code for setting up asymmetric image comparison 
datasets (Asym-QCD and Asym-DCQ), and for using them in retrieval experiments.

Please follow the instructions below -- if any questions/issues arise, please contact Andre at afaraujo@alumni.stanford.edu

## Pre-requisites

The main portion of this github repository must be working.
Please follow the instructions outlined in https://github.com/andrefaraujo/videosearch/blob/master/README.md#quick-start

## Download of relevant data

In the following, "mypath" refers to the path you downloaded the repository to.

**Download configuration files:**

    $ cd $mypath/videosearch/asym/asym_qcd
TODO    $ wget # Download config files
    $ unzip -j asym_qcd.zip
    $ cd ../asym_dcq
TODO: download
    $ unzip -j asym_dcq.zip

**Create data folder:**

    $ cd $mypath/videosearch/asym/
    $ mkdir data

**Download SMVS dataset:**

    $ cd $mypath/videosearch/asym/data
    $ mkdir smvs
    $ cd smvs
    $ wget https://stacks.stanford.edu/file/druid:rb470rw0983/cd_covers.zip # Note: size of this file is 442MB
    $ unzip cd_covers.zip
    $ wget https://stacks.stanford.edu/file/druid:rb470rw0983/dvd_covers.zip # Note: size of this file is 439MB
    $ unzip dvd_covers.zip

**Download Holidays dataset:** 

(note: this can be slow)

    $ cd $mypath/videosearch/asym/data
    $ mkdir holidays 
    $ cd holidays
    $ wget ftp://ftp.inrialpes.fr/pub/lear/douze/data/jpg1.tar.gz # Note: size of this file is 1.1GB
    $ tar -zxvf jpg1.tar.gz
    $ wget ftp://ftp.inrialpes.fr/pub/lear/douze/data/jpg2.tar.gz # Note: size of this file is 1.6GB
    $ tar -zxvf jpg2.tar.gz

**Download MIR-FLICKR-1M:** 

(note: this can be quite slow)

(only the required parts for these experiments are downloaded)

    $ cd $mypath/videosearch/asym/data
    $ mkdir mirflickr1m
    $ cd mirflickr1m
    $ for i in 0 1 2 3 4 5; do
    $   wget http://press.liacs.nl/mirflickr/mirflickr1m/images${i}.zip # Note: size of these files is 12GB
    $   unzip images${i}.zip
    $ done

## Asym-QCD

Here, we illustrate an example using 5 clutter images per database image (ie, C=5).
This choice is reflected in the parameters of the files used below.
In the following, "mypath" refers to the path you downloaded the repository to.

**Extract query features, using 10 threads:**

    $ cd $mypath/videosearch/asym/asym_qcd
    $ ./run_sift_extraction_query.sh

**Extract database features:** 

(this can take some time)

Note: This extracts features for images from C=5 set, using 10 threads.
For other values of C, change the chosen file in the following script.

    $ ./run_sift_extraction_database.sh

**Extract global descriptors of database using simple parameters:**

    $ ./run_gd_index.sh

**Retrieve using simple parameters:**

    $ ./run_retriever_asym.sh

**Convert to scoreable format:**

    $ ./run_convert_results_format.sh

**Score:**

    $ ./run_evaluate_results.sh

