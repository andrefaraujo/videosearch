import sys
from numpy import mean

## Tolerances used to extend the Ground-Truth segments
## to help assure all correct frames are included. 
## This is used because variations in the precise time segments, 
## due to the usage of different video players, was observed, 
## of up to 1 second, as reported in the MMSys'15 paper.
START_TOLERANCE = 1
END_TOLERANCE = 1

## A delay of 4 frames was observed for our version of ffmpeg frame extraction
## This number may be different for you, so please check for yourself
## Only used for "frame" mode results file format
FRAME_NUMBER_TIME_DELAY = 4



def timeStrToSeconds(time):
    time_split = time.split(':')
    if len(time_split) > 3:
        return -1
    val = 0
    for i in range(len(time_split)):
        val += int(time_split[-1 - i]) * (60**i)
    return val


def score_results_times(results_filename, ground_truth_filename):
    f = open(ground_truth_filename, 'r')
    ground_truth_lines = f.readlines()
    f.close()
    
    f = open(results_filename, 'r')

    jaccard_list = []

    while True:  ## For each Query
        file_pos = f.tell()
        query_line = f.readline()
        if query_line == '':
            break
        if not query_line.startswith("Query"):
            print "Error: Query line poorly formatted, at position %d" % file_pos
            exit()
        query_num = int(query_line.split()[-1])
        print "Query %d:" % query_num

        jaccard_this_query_list = []

        gt_line_split = ground_truth_lines[query_num].split()
        gt_clip_names = []
        for i in xrange(1, len(gt_line_split), 3):
            if gt_line_split[i] not in gt_clip_names:
                gt_clip_names.append(gt_line_split[i])
        number_gt_clips = len(gt_clip_names)


        while True:  ## For each Ground-Truth video for this query (in results file)
            file_pos = f.tell()
            result_line = f.readline()
            if result_line.startswith("Query"):
                f.seek(file_pos)
                break
            if result_line == '':
                break

            video_name = result_line.split(',')[0]

            result_reported_seconds = []
            for segment in result_line.split(',')[1:]:
                result_reported_seconds.extend(range(timeStrToSeconds(segment.split()[0]), timeStrToSeconds(segment.split()[1]) + 1))

            gt_reported_seconds = []
            for i in xrange(2, len(gt_line_split), 3):
                if gt_line_split[i-1] == video_name:
                    start_time = max(timeStrToSeconds(gt_line_split[i]) - START_TOLERANCE, 0)
                    end_time = timeStrToSeconds(gt_line_split[i+1]) + END_TOLERANCE
                    gt_reported_seconds.extend(range(start_time, end_time + 1))

            if len(gt_reported_seconds) == 0:
                print "Error: Problem finding ground-truth reported segments, Query %d, video %s" % (query_num, video_name)
                exit()

            result_reported_seconds = set(result_reported_seconds)
            gt_reported_seconds = set(gt_reported_seconds)

            num_intersection_frames = len(result_reported_seconds.intersection(gt_reported_seconds))
            num_union_frames = len(result_reported_seconds.union(gt_reported_seconds))

            jaccard_this_query_video = float(num_intersection_frames) / float(num_union_frames)

            ## Print Jaccard results
            print "  Video %s, Jaccard Index = %f" % (video_name, jaccard_this_query_video)
            jaccard_this_query_list.append(jaccard_this_query_video)


        ## Print Query Total Jaccard results
        mean_jaccard_this_query = sum(jaccard_this_query_list) / float(number_gt_clips)
        print "Total Results for Query %d: mean Jaccard Index = %f\n" % (query_num, mean_jaccard_this_query)
        jaccard_list.append(mean_jaccard_this_query)

    ## Print Total Retrieval Jaccard results
    print "\n\nTotal Results for All Queries: mean Jaccard Index = %f" % mean(jaccard_list)
    f.close()



def score_results_frames(results_filename, ground_truth_filename, short_list_size):
    f = open(ground_truth_filename, 'r')
    ground_truth_lines = f.readlines()
    f.close()
    
    f = open(results_filename, 'r')

    jaccard_list = []

    while True:  ## For each Query
        file_pos = f.tell()
        query_line = f.readline()
        if query_line == '':
            break
        if not query_line.startswith("Query"):
            print "Error: Query line poorly formatted, at position %d" % file_pos
            exit()
        query_num = int(query_line.split()[-1])
        print "Query %d:" % query_num

        jaccard_this_query_list = []

        gt_line_split = ground_truth_lines[query_num].split()
        gt_clip_names = []
        for i in xrange(1, len(gt_line_split), 3):
            if gt_line_split[i] not in gt_clip_names:
                gt_clip_names.append(gt_line_split[i])
        number_gt_clips = len(gt_clip_names)

        while True:  ## For each Ground-Truth video for this query (in results file)
            file_pos = f.tell()
            result_line = f.readline()
            if result_line.startswith("Query"):
                f.seek(file_pos)
                break
            if result_line == '':
                break

            video_name = result_line.split(',')[0]

            result_reported_seconds = []
            for keyframe in result_line.split(',')[1:]:
                result_reported_seconds.append(max(int(keyframe) - FRAME_NUMBER_TIME_DELAY, 0))

            if short_list_size != -1:
                result_reported_seconds = result_reported_seconds[0:short_list_size]

            gt_reported_seconds = []
            for i in xrange(2, len(gt_line_split), 3):
                if gt_line_split[i-1] == video_name:
                    start_time = max(timeStrToSeconds(gt_line_split[i]) - START_TOLERANCE, 0)
                    end_time = timeStrToSeconds(gt_line_split[i+1]) + END_TOLERANCE
                    gt_reported_seconds.extend(range(start_time, end_time + 1))

            if len(gt_reported_seconds) == 0:
                print "Error: Problem finding ground-truth reported segments, Query %d, video %s" % (query_num, video_name)
                exit()

            result_reported_seconds = set(result_reported_seconds)
            gt_reported_seconds = set(gt_reported_seconds)

            num_intersection_frames = len(result_reported_seconds.intersection(gt_reported_seconds))
            num_union_frames = len(result_reported_seconds.union(gt_reported_seconds))

            jaccard_this_query_video = float(num_intersection_frames) / float(num_union_frames)

            ## Print Jaccard results
            print "  Video %s, Jaccard Index = %f" % (video_name, jaccard_this_query_video)
            jaccard_this_query_list.append(jaccard_this_query_video)


        ## Print Query Total Jaccard results
        mean_jaccard_this_query = sum(jaccard_this_query_list) / float(number_gt_clips)
        print "Total Results for Query %d: mean Jaccard Index = %f\n" % (query_num, mean_jaccard_this_query)
        jaccard_list.append(mean_jaccard_this_query)

    ## Print Total Retrieval Jaccard results
    print "\n\nTotal Results for All Queries: mean Jaccard Index = %f" % mean(jaccard_list)
    f.close()






def printUsage():
    print "Usage: python " + sys.argv[0] + " results_filename ground_truth_filename mode [short_list_size]"
    print "results_filename: file following the temporal_refinement_results_file_format.txt rules"
    print "ground_truth_filename: light_dataset_public.txt or full_dataset_public.txt"
    print "mode: \'times\' or \'frames\', specifies the results file format used, doesn't change scoring details"
    print "[Optional] short_list_size: Number of keyframes to include in scoring.  Typical value is 50."
    print "                            Without this option, all keyframes reported in results_filename "
    print "                                will be used."
    print "                            This option should only be used with the 'frames' mode."
    print "                            In the MMSys'15 paper, we use a value of 50, so this option should"
    print "                                be set to 50 to get comparable results."
    print "                            Note: Use of this option assumes that keyframes listed in results_filename"
    print "                                have been listed in a decreasing order of confidence."
    print "                            Also, note that with our observed keyframe delay in ffmpeg, keyframes 1-4"
    print "                                all map to the time 0:00 and are not double counted.  Therefore if your"
    print "                                results report a list of keyframes short_list_size long, but that list"
    print "                                includes multiple keyframes in the range 1-4, then those keyframe will"
    print "                                not be double counted and the effective length of your list will be less"
    print "                                than short_list_size.  This is a very minor effect, but it can happen."



if __name__ == "__main__":
    if len(sys.argv) < 4:
        printUsage()
        exit()

    results_filename = sys.argv[1]
    ground_truth_filename = sys.argv[2]
    mode = sys.argv[3]

    short_list_size = -1
    if len(sys.argv) >= 5:
        short_list_size = int(sys.argv[4])

    if mode == "times":
        if short_list_size != -1:
            print "Error: Cannot use short_list_size optional argument with 'times' mode."
            exit()
        score_results_times(results_filename, ground_truth_filename)
    elif mode == "frames":
        score_results_frames(results_filename, ground_truth_filename, short_list_size)
    else:
        print "Invalid mode argument.  See usage below.\n"
        printUsage()
        exit()


