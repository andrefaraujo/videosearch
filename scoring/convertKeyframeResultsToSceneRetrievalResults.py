import sys

## Helper script to convert from simple keyframe-list formatted results to the scene retrieval results file format
## Not necessary to use this script, just a convenience if you want to print your system's results in a simple 
##     keyframe-list format and then use this script to convert your results to a format used by the evaluation scripts


def convertFileFormat_keyframeToSceneRetrieval(keyframe_results_filename, scene_retrieval_results_filename, short_list_size):
    fin = open(keyframe_results_filename, 'r')
    fout = open(scene_retrieval_results_filename, 'w')

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

        while True:
            file_pos = fin.tell()
            line = fin.readline()
            if line.startswith("Query"):
                fin.seek(file_pos)
                break
            if line == '':
                break
            
            if len(results_list) <= short_list_size:
                video_name = line[0:line.rindex('_')] + ".mp4"
                if video_name not in results_list:
                    results_list.append(video_name)

        for video in results_list:
            fout.write("%s\n" % video)

    fin.close()
    fout.close()



def printUsage():
    print "Usage: python " + sys.argv[0] + " keyframe_results_filename scene_retrieval_results_filename short_list_size"
    print "keyframe_results_filename: file following the keyframe_results_file_format.txt rules"
    print "scene_retrieval_results_filename: output filename for results in scene retrieval file format"
    print "short_list_size: Number of scene results to keep per query.  Typical value is 100"



if __name__ == "__main__":
    if len(sys.argv) < 4:
        printUsage()
        exit()

    keyframe_results_filename = sys.argv[1]
    scene_retrieval_results_filename = sys.argv[2]
    short_list_size = int(sys.argv[3])

    convertFileFormat_keyframeToSceneRetrieval(keyframe_results_filename, scene_retrieval_results_filename, short_list_size)
