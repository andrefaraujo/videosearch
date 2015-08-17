# videosearch: Large-Scale Video Retrieval Using Images

Andre Araujo ([stanford.edu/~afaraujo](http://stanford.edu/~afaraujo), afaraujo@stanford.edu)

Image, Video and Multimedia Systems Group, Stanford University

This project currently contains code for 
- Keyframe extraction from videos
- Shot boundary detector for videos
- SIFT descriptor extraction per frame
- Scoring a retrieval system that uses the Stanford I2V dataset

## Quick start

Clone repository:

    > cd $mypath
    > git clone https://github.com/andrefaraujo/videosearch.git

Creating VLFEAT library:

    > cd $mypath/videosearch/common/vlfeat-0.9.18/
    > make # Making vlfeat
    > sudo ln -s bin/glnxa64/libvl.so /usr/lib/libvl.so #sym-link it to your library

Build and test SIFT extraction:

    > cd $mypath/videosearch/indexer/local_features/
    > make # Making programs for extracting, reading and writing features
    > ./test_extract

You can also test reading a ".siftb" file (you can look at test_read.cc code for an example on how to read SIFT features binary files):

    > ./test_read # under $mypath/videosearch/indexer/local_features/

A larger scale extraction of features can be done with (look at the code of the shell script for more details):

    > ./run_sift_extraction.sh # under $mypath/videosearch/indexer/local_features/ 

Test keyframe extraction:

    > cd $mypath/videosearch/indexer/keyframes
    > ./run_keyframe_extraction.sh

Build and test shot boundary detection:

    > cd $mypath/videosearch/indexer/shot_detector
    > make
    > ./run_shot_detector.sh

Extracting keyframes and features from entire Stanford I2V dataset ([Dataset page](http://blackhole1.stanford.edu/vidsearch/dataset/stanfordi2v.html), [Download link](http://purl.stanford.edu/zx935qw7203)). Note: For this to work, you need to download the dataset beforehand and follow the instructions (found [here](https://stacks.stanford.edu/file/druid:zx935qw7203/README.txt)) for setting it up.

    > cd $mypath/videosearch/stanford_i2v/indexer/
    > python extract_database_keyframes.py # Look at script for more details and for changing parameters

Scoring results obtained with the Stanford I2V dataset. In this case, your should use a file with a specific format (as explained in the scoring/\*format\*.txt files). We provide examples of such files (scoring/example\*) and even helper conversion scripts if your system outputs results based on keyframes (scoring/convert\*). To score Scene Retrieval and Temporal Refinement results (refer to our [MMSys'15 paper](http://web.stanford.edu/~afaraujo/Araujo_et_al_MMSys_v14.pdf) for explanation of this terminology), respectively, do:

    > cd $mypath/videosearch/scoring
    > python evaluate_scene_retrieval.py example_scene_retrieval_results_file.txt light_dataset_public.txt 100
    > python evaluate_temporal_refinement.py example_temporal_refinement_results_file_frames.txt light_dataset_public.txt frames

## Citation
If you use the Stanford I2V dataset, please cite:

A. Araujo, J. Chaves, D. Chen, R. Angst and B. Girod. "Stanford I2V: A News Video Dataset for Query-by-Image Experiments", in Proc. ACM Multimedia Systems, 2015

Bibtex:

@inproceedings{AraujoMMSYS2015,

author = {Araujo, A. and Chaves, J. and Chen, D. and Angst, R. and Girod, B.},

booktitle = {Proc. ACM Multimedia Systems},

title = {{Stanford I2V: A News Video Dataset for Query-by-Image Experiments}},

year = {2015}

}
