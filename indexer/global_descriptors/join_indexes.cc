/**********************************************
This program joins several GDIndex's

Note that several parameters of the GDIndex object are irrelevant
in this case, since all we're doing is loading and concatenating
indexes. So, many of them are set to default/0 values in the calls
to set_index_parameters() and set_query_parameters().
 **********************************************/

#include <iostream>
#include <fstream>
#include <thread>

#include "gdindex.h"

void usage(char** argv)
{
  printf("usage: %s [options] --indexes_list[-i] indexes_list_file_path --out_index[-o] out_index_file_path \n", argv[0]);
  printf("Example:\n");
  printf("./join_indexes -i test_frame_based_index_lists.txt -o test_frame_based_index_lists.sift_scfv_idx_k512 \n");
  printf("Options:\n");
  printf("--gdindex_parameters_path[-r] ARG: Path to folder of trained GDIndex parameters. (default: trained_parameters)\n");
  printf("--centroids[-c] ARG: Number of Gaussians/centroids to use in global descriptor (default: 512)\n");
  printf("--threads[-t] ARG: Number of threads to use (default: 1)\n");
  printf("--ld_mode[-l] ARG: Local descriptor mode to use. 0=SIFT; 1=SIFTGeo (default: 0)\n");
  printf("--verbose_level[-v] ARG (default: 1) \n");
}

void join_indexes(const vector<string>& index_files, 
                  const int start_ind, int end_ind,
                  GDIndex& gdindex) {
    uint number_files = index_files.size();
    if (end_ind == -1) {
        end_ind = number_files;
    }
    for (int count_i = start_ind; count_i < end_ind; count_i++) {
        string this_index = index_files.at(count_i);
        // Second argument is true because we want to read both
        // types of auxiliary information
        gdindex.read(this_index, true);
    }
}

int main(int argc, char** argv) {
    // Mandatory arguments
    string indexes_list_file_path = "";
    string out_index_file_path = "";
    
    // Default values for options
    int number_gaussians = 512;
    string gdindex_parameters_path = "trained_parameters";
    int number_threads = 1;
    uint ld_mode = 0;
    int verbose_level = 1;
    if (argc < 3) {
        cout << "Wrong usage!!!" << endl;
        cout << "***********************************" << endl;
        cout << "***********************************" << endl;
        cout << "***********************************" << endl;
        cout << "See usage below:" << endl;
        usage(argv);
        exit(EXIT_FAILURE);
    } else {
        for (int count_arg = 1; count_arg < argc; count_arg++) {
            if ((!strcmp(argv[count_arg], "--indexes_list")) || (!strcmp(argv[count_arg], "-i"))) {
                indexes_list_file_path = string(argv[count_arg + 1]);
                count_arg++;
            } else if ((!strcmp(argv[count_arg], "--out_index")) || (!strcmp(argv[count_arg], "-o"))) {
                out_index_file_path = string(argv[count_arg + 1]);
                count_arg++;
            } else if ((!strcmp(argv[count_arg], "--gdindex_parameters_path")) || (!strcmp(argv[count_arg], "-r"))) { 
                gdindex_parameters_path = string(argv[count_arg + 1]);
                count_arg++;
            } else if ((!strcmp(argv[count_arg], "--centroids")) || (!strcmp(argv[count_arg], "-c"))) {
                number_gaussians = atoi(argv[count_arg + 1]);
                count_arg++;
            } else if ((!strcmp(argv[count_arg], "--threads")) || (!strcmp(argv[count_arg], "-t"))) {
                number_threads = atoi(argv[count_arg + 1]);
                count_arg++;
            } else if ((!strcmp(argv[count_arg], "--ld_mode")) || (!strcmp(argv[count_arg], "-l"))) {
                ld_mode = static_cast<uint>(atoi(argv[count_arg + 1]));
                count_arg++;
            } else if ((!strcmp(argv[count_arg], "--verbose_level")) || (!strcmp(argv[count_arg], "-v"))) {
                verbose_level = atoi(argv[count_arg + 1]);
                count_arg++;
            } else {
                cout << "Incorrect argument!" << endl;
                cout << "See usage below:" << endl;
                usage(argv);
                exit(EXIT_FAILURE);
            }        
        }
    }

    // Check that all mandatory arguments were provided
    if (indexes_list_file_path == "") {
        cout << "indexes_list argument was not provided... see usage below:" << endl;
        usage(argv);
        exit(EXIT_FAILURE);
    }
    if (out_index_file_path == "") {
        cout << "out_index argument was not provided... see usage below:" << endl;
        usage(argv);
        exit(EXIT_FAILURE);
    }

    if (verbose_level) {
        cout << "Starting joining indexes using:" << endl;
        cout << "------>indexes_list_file_path = " << indexes_list_file_path << endl;
        cout << "------>out_index_file_path = " << out_index_file_path << endl;
        cout << "------>gdindex_parameters_path = " << gdindex_parameters_path << endl;
        cout << "------>number_gaussians = " << number_gaussians << endl;
        cout << "------>number_threads = " << number_threads << endl;
    }

    // Read index file paths
    vector<string> index_files;
    FILE* list_file = fopen(indexes_list_file_path.c_str(), "r");
    if (list_file == NULL) {
        fprintf(stderr, "Cannot open index list file: %s\n", indexes_list_file_path.c_str());
        exit(EXIT_FAILURE);
    }
    char text_buffer[4096];
    while (true) {
        int n_read = fscanf(list_file, "%s", text_buffer);
        if (n_read <= 0) break;
        index_files.push_back(string(text_buffer));
    }
    fclose(list_file);

    // Define some index parameters
    uint ld_length, ld_frame_length;
    string ld_extension, ld_name;
    if (ld_mode == GDIndex::SIFT_LOCAL_DESCRIPTOR) {
        ld_length = GDIndex::SIFT_LENGTH;
        ld_frame_length = GDIndex::SIFT_FRAME_LENGTH;
        ld_extension = GDIndex::SIFT_EXTENSION;
        ld_name = GDIndex::SIFT_NAME;
    } else if (ld_mode == GDIndex::SIFTGEO_LOCAL_DESCRIPTOR) {
        ld_length = GDIndex::SIFTGEO_LENGTH;
        ld_frame_length = GDIndex::SIFTGEO_FRAME_LENGTH;
        ld_extension = GDIndex::SIFTGEO_EXTENSION;
        ld_name = GDIndex::SIFTGEO_NAME;
    } else {
        cout << "Problem! ld_mode = " 
             << ld_mode
             << " is not recognized"
             << endl;
        exit(EXIT_FAILURE);
    }

    if (number_threads == 1) {
        // Join all indexes using one thread
        GDIndex gdindex;
        gdindex.set_index_parameters(ld_length, ld_frame_length, ld_extension, ld_name,
                                     GDIndex::LD_PCA_DIM, LD_PRE_PCA_POWER_DEFAULT, 
                                     number_gaussians,
                                     GD_POWER_DEFAULT,
                                     GD_INTRA_NORMALIZATION_DEFAULT,
                                     gdindex_parameters_path,
                                     verbose_level);
        gdindex.set_query_parameters(0, 0, 0, 0, 0, gdindex_parameters_path, verbose_level);

        // Join all indexes
        cout << "Joining " << index_files.size() << " indexes..." << endl;
        join_indexes(index_files, 0, -1, gdindex);
        cout << "done!" << endl;

        // Save index to file
        printf("Saving index to:\n%s\n", out_index_file_path.c_str());
        gdindex.write(out_index_file_path);
    } else {
        // In this case, use multiple threads to join indexes
        // Divide index_files' entries into each thread's ones
        int number_items = index_files.size();
        int computations_per_thread = static_cast<int>(floor(number_items/number_threads));
        int rest_computations = number_items % number_threads;
        uint number_threads_to_use = number_threads;

        if (computations_per_thread == 0) { 
            // Less work than thread, so each thread will do one item
            computations_per_thread = 1;
            rest_computations = 0;
            number_threads_to_use = number_items;
        }        

        // Create vector of index pointers, initialize them, assign them to 
        // subset of dataset
        vector<GDIndex*> partial_indexes(number_threads_to_use, NULL);
        vector<thread> join_threads;
        uint number_items_so_far = 0;
        for (uint count_thread = 0; count_thread < number_threads_to_use; count_thread++) {
            uint start_index = number_items_so_far;
            uint number_items_per_thread;
            // Casting is necessary here since we're comparing numbers that might be negative:
            if (static_cast<int>(count_thread) > rest_computations - 1) {
                number_items_per_thread = computations_per_thread;
            }
            else {
                number_items_per_thread = computations_per_thread + 1;
            }
            uint end_index = start_index + number_items_per_thread;

            // Initialize this index
            partial_indexes.at(count_thread) = new GDIndex();
            partial_indexes.at(count_thread)->set_index_parameters(ld_length, ld_frame_length, 
                                                                   ld_extension, ld_name,
                                                                   GDIndex::LD_PCA_DIM, LD_PRE_PCA_POWER_DEFAULT, 
                                                                   number_gaussians,
                                                                   GD_POWER_DEFAULT, 
                                                                   GD_INTRA_NORMALIZATION_DEFAULT,
                                                                   gdindex_parameters_path,
                                                                   verbose_level);
            partial_indexes.at(count_thread)->set_query_parameters(0, 0, 0, 0, 0, gdindex_parameters_path, 
                                                                   verbose_level);

            // Join indexes from desired range
            join_threads.push_back(thread(join_indexes,
                                          index_files, 
                                          start_index,
                                          end_index,
                                          ref(*partial_indexes.at(count_thread))));

            number_items_so_far += number_items_per_thread;
        }

        // Join threads, write each resulting index to a temp. file
        // and add these paths to vector of names
        vector<string> partial_indexes_paths;
        for (uint count_thread = 0; count_thread < number_threads_to_use; count_thread++) {
            join_threads.at(count_thread).join();
            string out_index_file_path_temp = out_index_file_path + to_string(count_thread);
            partial_indexes.at(count_thread)->write(out_index_file_path_temp);
            // Add to a vector of temp. indexes paths
            partial_indexes_paths.push_back(out_index_file_path_temp);
            // Clean up partial indexes
            if (partial_indexes.at(count_thread)) {
                delete partial_indexes.at(count_thread);
                partial_indexes.at(count_thread) = NULL;
            }
        }

        // Instantiate final index
        GDIndex gdindex_final;
        gdindex_final.set_index_parameters(ld_length, ld_frame_length, ld_extension, ld_name,
                                           GDIndex::LD_PCA_DIM, LD_PRE_PCA_POWER_DEFAULT, 
                                           number_gaussians,
                                           GD_POWER_DEFAULT, 
                                           GD_INTRA_NORMALIZATION_DEFAULT,
                                           gdindex_parameters_path,
                                           verbose_level);
        gdindex_final.set_query_parameters(0, 0, 0, 0, 0, gdindex_parameters_path, verbose_level);

        // Join partial indexes constructed by each thread
        cout << "Joining " << partial_indexes_paths.size() << " indexes from each thread..." << endl;
        join_indexes(partial_indexes_paths, 0, 
                     partial_indexes_paths.size(), gdindex_final);
        cout << "done!" << endl;
        
        // Delete temp files
        for (uint count_thread = 0; count_thread < number_threads_to_use; count_thread++) {
            if ( remove(partial_indexes_paths.at(count_thread).c_str()) != 0 ) {
                cout << "Error! Could not remove " 
                     << partial_indexes_paths.at(count_thread) << endl;
                exit(EXIT_FAILURE);
            }
        }

        // Write final index
        gdindex_final.write(out_index_file_path);
    }

    return EXIT_SUCCESS;
    
}



