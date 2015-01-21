# videosearch: Large-Scale Video Retrieval Using Images

Andre Araujo ([stanford.edu/~afaraujo](http://stanford.edu/~afaraujo), afaraujo@stanford.edu)

Image, Video and Multimedia Systems Group, Stanford University

This project currently contains code for 
- Keyframe extraction from videos
- SIFT descriptor extraction at large scale
- Scoring a retrieval system on the Stanford I2V dataset

## Quick start

Clone repository

    > cd $mypath
    > git clone https://github.com/andrefaraujo/videosearch.git
    > cd videosearch # the code repository

Creating executables

    > cd common/vlfeat-0.9.18/
    > make # Making vlfeat
    > sudo ln -s bin/glnxa64/libvl.so /usr/lib/libvl.so #sym-link it to your library
    > cd ../../indexer/extract_features/
    > make # Making programs for extracting, reading and writing features

At this point, you can test SIFT extraction on a single image by running:

    > ./test_extract

You can also test reading a ".siftb" file:

    > ./test_read

A larger scale extraction of features can be done with:

    > ./run_sift_extraction.sh # Look at the code of the shell script for more details

Test keyframe extraction

    > cd ../
    > python extract_keyframes.py test_video.mp4 test_video_out 1 scale=-1:480

Extracting keyframes and features from entire Stanford I2V dataset ([Dataset page](http://blackhole1.stanford.edu/vidsearch/dataset/stanfordi2v.html), [Download link](http://purl.stanford.edu/zx935qw7203))

    > cd ../news_videos/indexer/
    > python extract_database_keyframes.py # Look at script for more details and for changing parameters

Scoring the dataset

    > TODO

## Citation
If you use the Stanford I2V dataset, please cite:

A. Araujo, J. Chaves, D. Chen, R. Angst and B. Girod. "Stanford I2V: A News Video Dataset for Query-by-Image Experiments", in Proc. ACM Multimedia Systems (MMSys) 2015

@inproceedings{AraujoMMSYS2015,

author = {Araujo, A. and Chaves, J. and Chen, D. and Angst, R. and Girod, B.},

booktitle = {Proc. ACM Multimedia Systems},

title = {{Stanford I2V: A News Video Dataset for Query-by-Image Experiments}},

year = {2015}

}
