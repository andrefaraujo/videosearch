# videosearch

**Large-Scale Video Retrieval Using Image Queries**  

By Andr&eacute; Araujo, in collaboration with Jason Chaves, David Chen and Haricharan Lakshman

Image, Video and Multimedia Systems Group, Stanford University

This project currently contains code for 
- Keyframe extraction from videos
- Shot boundary detector for videos
- SIFT descriptor extraction per image/frame
- Global descriptor (Binarized Fisher vectors) extraction per image/frame, shot or scene
- Retrieval in image or video databases using image queries
- Evaluating retrieval results based on Average Precision and Precision at 1

With these, you can reproduce the main results from the papers mentioned below, following the steps outlined in the next section.

This repository can also be useful if one is interested in searching a database of images using query images. In
that case, one can simply use the frame-based techniques described below.

Our implementation has been tested on Linux (Ubuntu) and Mac OS X.

For any questions or issues, feel free to get in touch.

## Quick start

Here we illustrate the usage of this repository's code by running through a simple example containing
4 database video clips and two image queries. This also serves as a way to make sure your code is working
properly.

**Prerequisites**: (all of these can be easily obtained for Ubuntu or OS X using 'apt-get install' 
or 'brew install' respectively)
- opencv (tested with version 2.4.0 on Ubuntu, and version 2.4.12 on OS X)
- ffmpeg (tested with version git-2012-08-24-fef9e84 on Ubuntu, and version 2.6.1 on OS X)
- pkg-config (tested with version 0.25 on Ubuntu, and version 0.28 on OS X)

**Step 1**: Clone repository (where `mypath` is the path you'll download the repository to):

```bash
cd $mypath
git clone https://github.com/andrefaraujo/videosearch.git
```

**Step 2**: Building VLFEAT library:

```bash
cd $mypath/videosearch/common/vlfeat-0.9.18/
make
```

**Step 3**: Building YAEL library:

```bash
cd $mypath/videosearch/common/yael_v260_modif/
./configure.sh
cd yael
make
```

**Step 4**: Extract keyframes from test database videos:

```bash
cd $mypath/videosearch/indexer/keyframes
./run_keyframe_extraction_test.sh
```

**Step 5**: Build shot boundary detector and extract shot boundaries for test database videos:

```bash
cd $mypath/videosearch/indexer/shot_detector
make
./run_shot_detector_test.sh
```

**Step 6**: Build SIFT extractor and extract SIFT for each keyframe in database:

```bash
cd $mypath/videosearch/indexer/local_descriptors/
make
./run_sift_extraction_test.sh
```

**Step 7**: Build global descriptor extractors and extract global descriptors per frame, shot and scene:

```bash
cd $mypath/videosearch/indexer/global_descriptors/
make
    
# Extract frame-based global descriptors (GD)
./run_frame_based_index_test.sh # extract GDs for each clip
./run_join_frame_based_index_test.sh # join all GDs in one index
    
# Extract shot-based global descriptors (GD) with mode LOC
./run_shot_based_index_mode_1_test.sh # extract GDs for each clip
./run_join_shot_based_index_mode_1_test.sh # join all GDs in one index
./run_process_shot_files_mode_1_test.sh # process auxiliary shot files for this mode

# Extract shot-based global descriptors (GD) with mode INDEP
./run_shot_based_index_mode_0_test.sh # extract GDs for each clip
./run_join_shot_based_index_mode_0_test.sh # join all GDs in one index
./run_process_shot_files_mode_0_test.sh # process auxiliary shot files for this mode
    
# Extract scene-based global descriptors (GD)
./run_scene_based_index_test.sh # extract GD for each clip
./run_join_scene_based_index_test.sh # join all GDs in one index
./run_process_scene_files_test.sh # process auxiliary scene files
./run_process_scene_rerank_files_test.sh # process auxiliary file for scene reranking
```

**Step 8**: Extract local descriptors (and optionally global descriptors) for query images (you need to do this before running retriever, which is the next step):

```bash
cd $mypath/videosearch/indexer/local_descriptors/
./run_sift_extraction_test_query.sh
# Optional: extract global descriptors
cd $mypath/videosearch/indexer/global_descriptors/
./run_query_index_test.sh
```

**Step 9**: Build retriever and run it for frame-, shot- and scene-based indexes:

```bash
cd $mypath/videosearch/retriever/
make

# Retrieve using frame-based global descriptors
./run_frame_test.sh

# Optional: Retrieve using frame-based global descriptors, using pre-computed query global descriptors
./run_frame_test_with_query_index.sh

# Retrieve using shot-based global descriptors, mode LOC
./run_shot_mode_1_test.sh

# Retrieve using shot-based global descriptors, mode INDEP
./run_shot_mode_0_test.sh

# Retrieve using scene-based global descriptors in first stage,
# then shot-based global descriptors in second stage
./run_scene_test.sh
```

**Step 10**: Evaluate retrieval results (calculate AP and p@1):

```bash
cd $mypath/videosearch/scoring/

# Evaluate frame-based results
./run_convert_frame_based_results_test.sh # converting results to scoreable format
./run_evaluate_frame_based_test.sh # calculating AP and p@1

# Optional: Evaluate frame-based results which used pre-computed query global descriptors
./run_convert_frame_based_results_test_query_index.sh # converting results to scoreable format
./run_evaluate_frame_based_test_query_index.sh # calculating AP and p@1

# Evaluate shot-based results, mode LOC
./run_convert_shot_based_mode_1_results_test.sh # converting results to scoreable format
./run_evaluate_shot_based_mode_1_test.sh # calculating AP and p@1

# Evaluate shot-based results, mode INDEP
./run_convert_shot_based_mode_0_results_test.sh # converting results to scoreable format
./run_evaluate_shot_based_mode_0_test.sh # calculating AP and p@1

# Evaluate scene-based results
./run_convert_scene_based_results_test.sh # converting results to scoreable format
./run_evaluate_scene_based_test.sh # calculating AP and p@1
```

After running the `run_evaluate_*` scripts, you should see the scores for each query and at the end the mean scores (mAP, mP@1). 
For this small example dataset, we get `mAP = 1` and `mP@1 = 1` for all of the cases illustrated above. 
You should obtain the same results if your code is working properly.
The retrieval results of frame-based experiments using pre-computed query global descriptors (the "optional" commands above) should be exactly the same as those without pre-computation.

## Indexing/Retrieving/Scoring using Hessian-Affine detector

The example above uses SIFT (DoG) detector + SIFT descriptor.
Usually, the Hessian-Affine (HA) detector provides better retrieval performance, compared to the SIFT detector.
In this section, we walk through an example using the HA detector (for the frame-based pipeline, but shot/scene-based results can also be obtained in a straightforward manner).

We'll use the HA detector from INRIA; the detector can be found [here](http://lear.inrialpes.fr/~jegou/data.php) -- download the program named "compute_descriptors_linux64" or "compute_descriptors_mac", depending on your platform.
Place the downloaded file under "indexer/local_descriptors".
Also, make sure netpbm is installed (see instructions [here](http://netpbm.sourceforge.net/getting_netpbm.php)).

**Step 1**: Extract HesAff+SIFT descriptors for the test example.
Edit the "indexer/local_descriptors/run_siftHesAff_extraction_test.sh" and "indexer/local_descriptors/run_siftHesAff_extraction_test_query.sh" by setting the netpbm path to the variable NETPBM_BIN_PATH.
Then, run:

    $  cd $mypath/videosearch/indexer/local_descriptors
    $ ./run_siftHesAff_extraction_test.sh
    $ ./run_siftHesAff_extraction_test_query.sh

**Step 2**: Build global descriptors for database frames and query images.
Edit the files "indexer/global_descriptors/run_frame_based_index_test.sh", "indexer/global_descriptors/run_join_frame_based_index_test.sh" and "indexer/global_descriptors/run_query_index_test.sh" by setting the variable LD_MODE to 1.
Then, run:

    $ cd $mypath/videosearch/indexer/global_descriptors
    $ ./run_frame_based_index_test.sh
    $ ./run_join_frame_based_index_test.sh
    $ ./run_query_index_test.sh

**Step 3**: Run retriever.
Edit the file "retriever/run_frame_test_with_query_index.sh" by setting the variable FEAT_MODE to 1.
Then, run:

    $ cd $mypath/videosearch/retriever/
    $ ./run_frame_test_with_query_index.sh

**Step 4**: Evaluate retrieval results.
Edit the files "scoring/run_convert_frame_based_results_test_query_index.sh" and "scoring/run_evaluate_frame_based_test_query_index.sh" by setting the variable FEAT_MODE to 1.
Then, run:

    $ cd $mypath/videosearch/scoring/
    $ ./run_convert_frame_based_results_test_query_index.sh
    $ ./run_evaluate_frame_based_test_query_index.sh

You should obtain "Total Results: mAP = 1.000000, mP@1 = 1.000000".

## Performing retrieval on the Stanford I2V dataset

Here we provide some scripts to use our programs and obtain results on the Stanford I2V dataset ([Dataset page](http://blackhole1.stanford.edu/vidsearch/dataset/stanfordi2v.html), [Download link](http://purl.stanford.edu/zx935qw7203)). For this to work, you need to download the dataset beforehand and follow the instructions (found [here](https://stacks.stanford.edu/file/druid:zx935qw7203/README.txt)) for setting it up.

Extracting keyframes from entire Stanford I2V dataset:

    $ cd $mypath/videosearch/stanford_i2v/indexer/
    $ python extract_database_keyframes.py # Look at script for more details and for changing parameters

After extracting frames from the dataset, you can run through steps 5 to 10 above to index, retrieve and get results on the Stanford I2V dataset.

To score results obtained with the Stanford I2V dataset, you should use a file with a specific format (as explained in the scoring/\*format\*.txt files). We provide examples of such files (scoring/example\*) and even helper conversion scripts if your system outputs results based on keyframes (scoring/convert\*). To score Scene Retrieval and Temporal Refinement results (refer to our [MMSys'15 paper](http://web.stanford.edu/~afaraujo/Araujo_et_al_MMSys_v14.pdf) for explanation of this terminology), respectively, do:

    $ cd $mypath/videosearch/scoring
    $ python evaluate_scene_retrieval.py example_scene_retrieval_results_file.txt light_dataset_public.txt 100
    $ python evaluate_temporal_refinement.py example_temporal_refinement_results_file_frames.txt light_dataset_public.txt frames

Most often, when using this dataset, one is interested in Scene Retrieval results (i.e., retrieving the correct video clips), as in the ICIP'15 paper mentioned below.

## Citation
If you use this code, please cite:

A. Araujo, J. Chaves, R. Angst and B. Girod. "Temporal Aggregation for Large-Scale Query-by-Image Video Retrieval", in Proc. ICIP, 2015 ([Paper](http://web.stanford.edu/~afaraujo/Araujo_et_al_ICIP15_v12.pdf)) ([Poster](http://web.stanford.edu/~afaraujo/2015_09_28_ICIP_poster_v3.pdf))

Bibtex:

    @inproceedings{AraujoICIP2015,
    author = {Araujo, A. and Chaves, J. and Angst, R. and Girod, B.},
    booktitle = {Proc. ICIP},
    title = {{Temporal Aggregation for Large-Scale Query-by-Image Video Retrieval}},
    year = {2015}
    }

If you use the Stanford I2V dataset, please cite:

A. Araujo, J. Chaves, D. Chen, R. Angst and B. Girod. "Stanford I2V: A News Video Dataset for Query-by-Image Experiments", in Proc. ACM Multimedia Systems, 2015 ([Paper](http://web.stanford.edu/~afaraujo/Araujo_et_al_MMSys_v14.pdf)) ([Slides](http://web.stanford.edu/~afaraujo/2015_03_19_MMSys_talk.pdf))

Bibtex:

    @inproceedings{AraujoMMSYS2015,
    author = {Araujo, A. and Chaves, J. and Chen, D. and Angst, R. and Girod, B.},
    booktitle = {Proc. ACM Multimedia Systems},
    title = {{Stanford I2V: A News Video Dataset for Query-by-Image Experiments}},
    year = {2015}
    }

## License
MIT (except for the pieces of code from other sources -- check their original licenses)

