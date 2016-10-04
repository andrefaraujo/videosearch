/**********************************************
This program builds a point-indexed database from local features.

Note that the query parameters of the GDIndex object are irrelevant
in this case, since all we're doing is generating global descriptors.
So, many of them are set to default/0 values in the call to
set_query_parameters().
 **********************************************/

#include <cassert>
#include <iostream>
#include <fstream>

#include "../../indexer/global_descriptors/gdindex.h"
#include "point_index_io.h"

// Function to get descriptor file name
string get_local_descriptor_filename(const string image_path,
                                     const string ext) {
    uint pos = image_path.find_last_of(".");
    return (image_path.substr(0, pos) + ext);
}

// Usage of this program
void usage(char** argv)
{
  printf("usage: %s [options] --image_list[-i] image_list_file_path --out_index[-o] index_file_path \n", argv[0]);
  printf("Example:\n");
  printf("./index_dataset -i test_db_all_keyframes.txt -o test_db_all_keyframes.fv_point_idx_k512\n");
  printf("Options:\n");
  printf("--gdindex_parameters_path[-r] ARG: path where GDIndex pretrained parameters are saved (default: trained_parameters)\n");
  printf("--centroids[-c] ARG: Number of centroids/Gaussians to use in global descriptor (default: 512)\n");
  printf("--ld_mode[-l] ARG: Local descriptor mode to use. 0=SIFT; 1=SIFTGeo (default: 0)\n");
  printf("--verbose_level[-v] ARG (default: 1) \n");
}

int main(int argc, char** argv)
{
    // Mandatory arguments
    string image_list_file_path = "";
    string out_file_path = "";
    
    // Default values for options
    uint number_gaussians = 512;
    uint ld_mode = 0;
    int verbose_level = 1;
    string gdindex_parameters_path = "trained_parameters";

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
            if ((!strcmp(argv[count_arg], "--image_list")) || (!strcmp(argv[count_arg], "-i"))) {
                image_list_file_path = string(argv[count_arg + 1]);
                count_arg++;
            } else if ((!strcmp(argv[count_arg], "--out_index")) || (!strcmp(argv[count_arg], "-o"))) {
                out_file_path = string(argv[count_arg + 1]);
                count_arg++;
            } else if ((!strcmp(argv[count_arg], "--gdindex_parameters_path")) || (!strcmp(argv[count_arg], "-r"))) { 
              gdindex_parameters_path = string(argv[count_arg + 1]);
              count_arg++;
            } else if ((!strcmp(argv[count_arg], "--centroids")) || (!strcmp(argv[count_arg], "-c"))) {
                number_gaussians = static_cast<uint>(atoi(argv[count_arg + 1]));
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
    if (image_list_file_path == "") {
        cout << "image_list argument was not provided... see usage below:" << endl;
        usage(argv);
        exit(EXIT_FAILURE);
    }
    if (out_file_path == "") {
        cout << "out_index argument was not provided... see usage below:" << endl;
        usage(argv);
        exit(EXIT_FAILURE);
    }

    if (verbose_level) {
        cout << "Starting index building using:" << endl;
        cout << "------>image_list_file_path = " << image_list_file_path << endl;
        cout << "------>out_file_path = " << out_file_path << endl;
        cout << "------>gdindex_parameters_path = " << gdindex_parameters_path << endl;
        cout << "------>number_gaussians = " << number_gaussians << endl;
        cout << "------>ld_mode = " << ld_mode << endl;
        cout << "------>verbose_level = " << verbose_level << endl;
    }

    // Initialize GDIndex
    GDIndex gdindex;
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
    gdindex.set_index_parameters(ld_length, ld_frame_length, ld_extension, ld_name,
                                 GDIndex::LD_PCA_DIM, LD_PRE_PCA_POWER_DEFAULT, 
                                 number_gaussians,
                                 GD_POWER_DEFAULT,
                                 false,
                                 false,
                                 gdindex_parameters_path,
                                 verbose_level);
    gdindex.set_query_parameters(0, 0, 0, 0, 0, gdindex_parameters_path, verbose_level);

    // Define feature files
    vector < string > feature_files;
    FILE* image_list_file = fopen(image_list_file_path.c_str(), "r");
    if (image_list_file == NULL) {
        fprintf(stderr, "Cannot open image list file: %s\n", image_list_file_path.c_str());
        exit(EXIT_FAILURE);
    }
    char text_buffer[2048];
    while (true) {
        int n_read = fscanf(image_list_file, "%s", text_buffer);
        if (n_read <= 0) break;
        string feature_file = get_local_descriptor_filename(text_buffer,
                                                            ld_extension);
        feature_files.push_back(feature_file);
    }
    fclose(image_list_file);
    if (verbose_level >= 3) cout << "Parsed " << feature_files.size() 
                                 << " feature files." << endl;

    // Generate database index
    vector < vector < uint > > vec_feat_assgns;
    vector < vector < float > > vec_feat_assgn_weights;
    vector < vector < vector < float > > > vec_feat_residuals;
    gdindex.generate_point_index(feature_files, verbose_level,
                                 vec_feat_assgns, vec_feat_assgn_weights,
                                 vec_feat_residuals);

	// Save index to file
	if (verbose_level)
        printf("Saving point index to:\n%s\n", out_file_path.c_str());

    write_fv_point_index(vec_feat_assgns, vec_feat_assgn_weights,
                         vec_feat_residuals, GDIndex::LD_PCA_DIM,
                         out_file_path);

	return EXIT_SUCCESS;
}

