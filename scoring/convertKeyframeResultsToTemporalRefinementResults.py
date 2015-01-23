import sys

## Helper script to convert from simple keyframe-list formatted results to the temporal refinement results file format ('frames' mode)
## Not necessary to use this script, just a convenience if you want to print your system's results in a simple 
##     keyframe-list format and then use this script to convert your results to a format used by the evaluation scripts


def convertFileFormat_keyframeToTemporalRefinement(keyframe_results_filename, temporal_refinement_results_filename, ground_truth_filename):
    f = open(ground_truth_filename)
    ground_truth_lines = f.readlines()
    f.close()

    fin = open(keyframe_results_filename, 'r')
    fout = open(temporal_refinement_results_filename, 'w')

    while True:
        line = fin.readline()
        if line == '':
            break
        if not line.startswith("Query"):
            print "Error: query line is malformed: %s" % line
            exit()
        
        query_num = int(line.split()[-1])
        print "On Query %d" % query_num
        fout.write(line)

        results_list = []
        results_subdict = {}

        while True:
            file_pos = fin.tell()
            line = fin.readline()
            if line.startswith("Query"):
                fin.seek(file_pos)
                break
            if line == '':
                break

            video_name = line[0:line.rindex('_')] + ".mp4"

            if video_name not in ground_truth_lines[query_num]:
                continue

            if video_name not in results_list:
                results_list.append(video_name)
                results_subdict[video_name] = []

            keyframe_number = int(line[line.rindex('/')+1 : line.rindex('.')])
            if keyframe_number not in results_subdict[video_name]:
                results_subdict[video_name].append(keyframe_number)

        for video in results_list:
            fout.write("%s" % video)
            for keyframe_number in results_subdict[video]:
                fout.write(",%d" % keyframe_number)
            fout.write('\n')

    fin.close()
    fout.close()



def printUsage():
    print "Usage: python " + sys.argv[0] + " keyframe_results_filename temporal_refinement_results_filename ground_truth_filename"
    print "keyframe_results_filename: file following the keyframe_results_file_format.txt rules"
    print "temporal_refinement_results_filename: output filename for results in temporal refinement file format"
    print "ground_truth_filename: light_dataset.txt or full_dataset.txt"


if __name__ == "__main__":
    if len(sys.argv) < 4:
        printUsage()
        exit()

    keyframe_results_filename = sys.argv[1]
    temporal_refinement_results_filename = sys.argv[2]
    ground_truth_filename = sys.argv[3]

    convertFileFormat_keyframeToTemporalRefinement(keyframe_results_filename, temporal_refinement_results_filename, ground_truth_filename)
