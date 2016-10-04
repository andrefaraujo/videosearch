# THIS IS ACTIVELY BEING EDITED, please wait until we're finished before using this

# Query-by-Image Video Retrieval Using Bloom Filters

This directory contain code for the Bloom Filter technique developed in our
paper ["Large-Scale Query-by-Image Video Retrieval Using Bloom Filters"](https://arxiv.org/abs/1604.07939).

Please follow the instructions below -- if any questions/issues arise, feel free to reach out.

## Pre-requisites

The main portion of this github repository must be working.
Please follow the instructions outlined [here](https://github.com/andrefaraujo/videosearch/blob/master/README.md#quick-start).
We're using the Hessian-Affine detector for these experiments, so please make sure to have the pipeline using it working (see instructions [here](https://github.com/andrefaraujo/videosearch#indexingretrievingscoring-using-hessian-affine-detector)).

## Instructions

We'll use the same images and videos from the example in the main portion of this repo.
Make sure that the Hessian-Affine+SIFT descriptors have been extracted, instructions
[here](https://github.com/andrefaraujo/videosearch#indexingretrievingscoring-using-hessian-affine-detector)

Here, we present the usage of BF-PI with LSH-B hashes (please refer to the [paper](https://arxiv.org/abs/1604.07939) for a description of
the terminology).

In the following, `mypath` refers to the path you downloaded the repository to.

**Extract Fisher point-indexed descriptors for database frames, and their binarized versions**

```bash
cd $mypath/bloom_filters/point_indexed/
./run_extract_pi_descriptors.sh
./run_binarize_pi_descriptors.sh
```

**Extract Fisher point-indexed descriptors for query images, and their binarized versions**

```bash
cd $mypath/bloom_filters/point_indexed/
./run_extract_pi_descriptors_query.sh
./run_binarize_pi_descriptors_query.sh
```

**Index and retrieve**

```bash
cd $mypath/bloom_filters/retriever/
./run_retrieve_bf.sh
```

**Scoring results:**

```bash
cd $mypath/bloom_filters/scoring/
./run_evaluate_results.sh
```

This last step gives the results: `Total Results: mAP = XXX, mP@1 = XXX`.
This illustrates the usage of this technique.

