import sys
from numpy import mean


def score_results(results_filename, ground_truth_filename, short_list_size, FORCE):
    f = open(ground_truth_filename, 'r')
    ground_truth_lines = f.readlines()
    f.close()
    
    f = open(results_filename, 'r')

    AP_list = []
    pAt1_list = []
        
    while True:
        file_pos = f.tell()
        query_line = f.readline()
        if query_line == '':
            break
        if not query_line.startswith("Query"):
            print "Error: Query line poorly formatted, at position %d" % file_pos
            exit()
        query_num = int(query_line.split()[-1])
    
        ap_den = 0
        gt_line_split = ground_truth_lines[query_num].split()
        gt_videos_list = []
        for i in xrange(1, len(gt_line_split), 3):
            if gt_line_split[i] not in gt_videos_list:
                gt_videos_list.append(gt_line_split[i])
        ap_den = len(gt_videos_list)
        if ap_den == 0:
            print "Error: Ground Truth file error, query %d has no correct videos" % query_num
            exit()
    
        ap_aux = 0
        pAt1 = 0
    
        relevances = []
        num_results_used = 0
    
        already_scored_clips = []
    
        while num_results_used < short_list_size:
            file_pos = f.tell()
            result_line = f.readline()
            if result_line.startswith("Query"):
                f.seek(file_pos)
                break
            if result_line == '':
                break
            
            video_name = result_line.rstrip(' \n\r')
    
            if video_name in already_scored_clips and not FORCE:
                print "Error: duplicate video in ranked list for Query %d." % query_num
                print "Duplicate video name was %s" % video_name
                print "Ranked scene retrieval lists are expected to be without duplicates and at least short_list_size long."
                exit()
    
            already_scored_clips.append(video_name)
            num_results_used += 1
    
            if video_name in ground_truth_lines[query_num]:
                if num_results_used == 1:
                    pAt1 = 1
    
                relevances.append(1)
                ap_aux += float(sum(relevances)) / float(num_results_used)
    
            else:
                relevances.append(0)

        if num_results_used < short_list_size and not FORCE:
            print "Error: ranked list for Query %d includes less than short_list_size unique scenes." % query_num
            print "Ranked scene retrieval lists are expected to contain at least short_list_size number of unqiue scenes."
            exit()
    
        ap_this_query = float(ap_aux) / float(ap_den)

        # Print scores
        print "Query %d: AP = %f, p@1 = %d" % (query_num, ap_this_query, pAt1)
        AP_list.append(ap_this_query)
        pAt1_list.append(pAt1)
    
        while True:
            file_pos = f.tell()
            line = f.readline()
            if line.startswith("Query") or line == '':
                f.seek(file_pos)
                break
    
    f.close()

    print "Total Results: mAP = %f, mP@1 = %f" % (mean(AP_list), mean(pAt1_list))



def printUsage():
    print "Usage: python " + sys.argv[0] + " results_filename ground_truth_filename short_list_size [--force]"
    print "results_filename: file following the scene_retrieval_results_file_format.txt rules"
    print "ground_truth_filename: light_dataset.txt or full_dataset.txt"
    print "short_list_size: Number of scene results to score. Typical value is 100"
    print "[Optional] --force: Optional argument used to ignore requirements regarding duplicates and"
    print "                       results list length relative to short_list_size."



if __name__ == "__main__":
    if len(sys.argv) < 4:
        printUsage()
        exit()

    results_filename = sys.argv[1]
    ground_truth_filename = sys.argv[2]
    short_list_size = int(sys.argv[3])

    FORCE = False
    if len(sys.argv) > 4 and sys.argv[4] == "--force":
        FORCE = True
        

    score_results(results_filename, ground_truth_filename, short_list_size, FORCE)


