# videosearch: Large-Scale Video Retrieval Using Images

Andr&eacute; Araujo ([stanford.edu/~afaraujo](http://stanford.edu/~afaraujo), afaraujo@stanford.edu)  
In collaboration with Jason Chaves, David Chen and Haricharan Lakshman

Image, Video and Multimedia Systems Group, Stanford University

This project currently contains code for 
- Keyframe extraction from videos
- Shot boundary detector for videos
- SIFT descriptor extraction per image/frame
- Global descriptor (Binarized Fisher vectors) extraction per frame, shot or scene
- Retrieval in image or video databases using image queries
- Scoring retrieval results based on Average Precision and Precision at 1
- Reproducing main results from the papers mentioned below (IN PROGRESS: to be concluded by Oct 18)

This repository can also be useful if one is interested in searching a database of images using query images. In
that case, one can simply use the frame-based techniques described below.

For any questions or issues, please get in touch through github or using the contact information mentioned above.

## Quick start

Here we illustrate the usage of this repository's code by running through a simple example containing
4 database video clips and two image queries. This also serves as a way to make sure your code is working
properly.

Clone repository (where "mypath" is the path you'll download the repository to):

    > cd $mypath
    > git clone https://github.com/andrefaraujo/videosearch.git

Creating VLFEAT library:

    > cd $mypath/videosearch/common/vlfeat-0.9.18/
    > make # Making vlfeat
    > sudo ln -s bin/glnxa64/libvl.so /usr/lib/libvl.so #sym-link it to your library

Test keyframe extraction:

    > cd $mypath/videosearch/indexer/keyframes
    > ./run_keyframe_extraction_test.sh

Build and test shot boundary detection:

    > cd $mypath/videosearch/indexer/shot_detector
    > make # Building shot boundary detector
    > ./run_shot_detector_test.sh

Build and test SIFT extraction:

    > cd $mypath/videosearch/indexer/local_descriptors/
    > make # Building programs for extracting, reading and writing features
    > ./run_sift_extraction_test.sh

Build global descriptors:

    > cd $mypath/videosearch/indexer/global_descriptors/
    > make # Building programs for extracting and joining indexes of global descriptors
    
    > # Extract frame-based global descriptors (GD)
    > ./run_frame_based_index_test.sh # extract GDs for each clip
    > ./run_join_frame_based_index_test.sh # join all GDs in one index
    
    > # Extract shot-based global descriptors (GD) with mode LOC
    > ./run_shot_based_index_mode_1_test.sh # extract GDs for each clip
    > ./run_join_shot_based_index_mode_1_test.sh # join all GDs in one index
    > ./run_process_shot_files_mode_1_test.sh # process auxiliary shot files for this mode

    > # Extract shot-based global descriptors (GD) with mode INDEP
    > ./run_shot_based_index_mode_0_test.sh # extract GDs for each clip
    > ./run_join_shot_based_index_mode_0_test.sh # join all GDs in one index
    > ./run_process_shot_files_mode_0_test.sh # process auxiliary shot files for this mode
    
    > # Extract scene-based global descriptors (GD)
    > ./run_scene_based_index_test.sh # extract GD for each clip
    > ./run_join_scene_based_index_test.sh # join all GDs in one index
    > ./run_process_scene_files_test.sh # process auxiliary scene files

Extract local descriptors for query images (you need to do this before running retriever, which is the next step):

    > cd $mypath/videosearch/indexer/local_descriptors/
    > ./run_sift_extraction_test_query.sh

Build and run retriever:

    > cd $mypath/videosearch/retriever/
    > make # Making library and program to do query-by-image video retrieval 

    > # Retrieve using frame-based global descriptors
    > ./run_frame_test.sh

    > # Retrieve using shot-based global descriptors, mode LOC
    > ./run_shot_mode_1_test.sh

    > # Retrieve using shot-based global descriptors, mode INDEP
    > ./run_shot_mode_0_test.sh

    TODO(andrefaraujo): retrieve in two stages using scene descriptors

Evaluate retrieval results (calculate AP and p@1):

    > cd $mypath/videosearch/scoring/

    > # Evaluate frame-based results
    > ./run_convert_frame_based_results_test.sh # converting results to scoreable format
    > ./run_evaluate_frame_based_test.sh # calculating AP and p@1

    > # Evaluate shot-based results, mode LOC
    > ./run_convert_shot_based_mode_1_results_test.sh # converting results to scoreable format
    > ./run_evaluate_shot_based_mode_1_test.sh # calculating AP and p@1

    > # Evaluate shot-based results, mode INDEP
    > ./run_convert_shot_based_mode_0_results_test.sh # converting results to scoreable format
    > ./run_evaluate_shot_based_mode_0_test.sh # calculating AP and p@1

    TODO(andrefaraujo): scene-based results

## Performing retrieval on Stanford I2V dataset

Here we provide helpful scripts to use our programs and obtain results on the Stanford I2V dataset ([Dataset page](http://blackhole1.stanford.edu/vidsearch/dataset/stanfordi2v.html), [Download link](http://purl.stanford.edu/zx935qw7203)). For this to work, you need to download the dataset beforehand and follow the instructions (found [here](https://stacks.stanford.edu/file/druid:zx935qw7203/README.txt)) for setting it up.

TODO(andrefaraujo): complete this by Oct/18

Extracting keyframes from entire Stanford I2V dataset:

    > cd $mypath/videosearch/stanford_i2v/indexer/
    > python extract_database_keyframes.py # Look at script for more details and for changing parameters

Scoring results obtained with the Stanford I2V dataset. In this case, your should use a file with a specific format (as explained in the scoring/\*format\*.txt files). We provide examples of such files (scoring/example\*) and even helper conversion scripts if your system outputs results based on keyframes (scoring/convert\*). To score Scene Retrieval and Temporal Refinement results (refer to our [MMSys'15 paper](http://web.stanford.edu/~afaraujo/Araujo_et_al_MMSys_v14.pdf) for explanation of this terminology), respectively, do:

    > cd $mypath/videosearch/scoring
    > python evaluate_scene_retrieval.py example_scene_retrieval_results_file.txt light_dataset_public.txt 100
    > python evaluate_temporal_refinement.py example_temporal_refinement_results_file_frames.txt light_dataset_public.txt frames

## Citation
If you use this code, please cite:

A. Araujo, J. Chaves, R. Angst and B. Girod. "Temporal Aggregation for Large-Scale Query-by-Image Video Retrieval", in Proc. ICIP, 2015

Bibtex:

    @inproceedings{AraujoICIP2015,
    author = {Araujo, A. and Chaves, J. and Angst, R. and Girod, B.},
    booktitle = {Proc. ICIP},
    title = {{Temporal Aggregation for Large-Scale Query-by-Image Video Retrieval}},
    year = {2015}
    }

If you use the Stanford I2V dataset, please cite:

A. Araujo, J. Chaves, D. Chen, R. Angst and B. Girod. "Stanford I2V: A News Video Dataset for Query-by-Image Experiments", in Proc. ACM Multimedia Systems, 2015

Bibtex:

    @inproceedings{AraujoMMSYS2015,
    author = {Araujo, A. and Chaves, J. and Chen, D. and Angst, R. and Girod, B.},
    booktitle = {Proc. ACM Multimedia Systems},
    title = {{Stanford I2V: A News Video Dataset for Query-by-Image Experiments}},
    year = {2015}
    }
