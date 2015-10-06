# Generate auxiliary file for scene re-ranking

import os, os.path, sys, glob, math, time, re, subprocess, shutil, errno

def print_usage():
    print "python %s <list_videos> <shot_thresh> <out_group_lists> " % os.path.basename(__file__);
    print "Example: python %s " % os.path.basename(__file__);
    print " test_video_lists.txt \ ";
    print " 0.8 \ ";
    print " test_group_lists_rerank.txt ";
    print " ";

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print_usage();
    else:
        # 0. Get arguments
        list_videos_file = sys.argv[1];
        with open(list_videos_file) as f:
            video_list = f.readlines();
            video_list = [line.strip() for line in video_list];
        shot_thresh = float(sys.argv[2]);
        out_group_lists = sys.argv[3];

        # Open output file
        fid_out = open(out_group_lists, 'w');

        # Shot extension
        new_extension = "_keyframes.shot_t%.1f" % shot_thresh

        # Loop over video paths, find shot files and get number of shots per video
        # Then, write to out file shot numbers per line for each video
        count_total_shots = 0;
        for video in video_list:
            # Get name for shot file associated to video
            (root, ext) = os.path.splitext(video);
            shot_file = "%s%s" % (root, new_extension);

            # Get total number of shots in this shot_file
            shots_this = file_len(shot_file);

            # Write to out file
            for shot in range(shots_this):
                fid_out.write("%d " % count_total_shots);
                count_total_shots += 1;
            fid_out.write("\n");

        # Close output file
        fid_out.close();

