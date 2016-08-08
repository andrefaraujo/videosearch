/**********************************************
This program builds a database index from local features.

Note that the query parameters of the GDIndex object are irrelevant
in this case, since all we're doing is generating global descriptors.
So, many of them are set to default/0 values in the call to
set_query_parameters().
 **********************************************/

#include <cassert>
#include <iostream>
#include <fstream>

#include "gdindex.h"

// Function to parse shot files
void parse_shot_file(string file_name, vector < uint >& output) {
    ifstream in_file(file_name.c_str());
    string line;
    output.clear();
    while (in_file >> line) {
        output.push_back(static_cast<uint>(atoi(line.c_str())));
    }
}

// Function to get descriptor file name
string get_local_descriptor_filename(const string image_path,
                                     const string ext) {
    uint pos = image_path.find_last_of(".");
    return (image_path.substr(0, pos) + ext);
}

// Usage of this program
void usage(char** argv)
{
  printf("usage: %s [options] --image_list[-i] image_list_file_path --out_index[-o] gdindex_file_path \n", argv[0]);
  printf("Example:\n");
  printf("./index_dataset -i test_db_all_keyframes.txt -o test_db_all_keyframes.scfv_idx_k512\n");
  printf("Options:\n");
  printf("--gdindex_parameters_path[-r] ARG: path where GDIndex pretrained parameters are saved (default: trained_parameters)\n");
  printf("--centroids[-c] ARG: Number of centroids/Gaussians to use in global descriptor (default: 512)\n");
  printf("--ld_mode[-l] ARG: Local descriptor mode to use. 0=SIFT; 1=SIFTGeo (default: 0)\n");
  printf("--gd_intra_normalization: Boolean that sets usage of intra-normalization mode for FVs (default: false)\n");
  printf("--shot_file_path[-s] ARG: Path to shot file. (default: not using it) \n");
  printf("--single_shot[-g]: Flag that sets the mode of using a single shot in entire database. This is equivalent to passing as shot_file_path a file that contains only one line with '0'. (default: not using it). This is used when aggregating over entire video clip. \n");
  printf("--shot_mode[-m] ARG: Mode to use when using shots: 0 = indep. keyframes per shot; 1 = feature aggregation in shot; (default = 0, used only if shot_file_path is passed or single_shot is set) \n");
  printf("--shot_keyf[-k] ARG: Max number of keyframes per shot to use, when using shot mode (it will try to use as many, but shots might be shorter than this number; in this case, it will use the number of frames in the shot). Set this to -1 if wanting to use all frames from shot (default: -1) \n");
  printf("--verbose_level[-v] ARG (default: 1) \n");
}

// Helper function
static bool file_exists(string filename)
{
  ifstream ifile(filename.c_str());
  return bool(ifile);
}

int main(int argc, char** argv)
{
    // Mandatory arguments
    string image_list_file_path = "";
    string out_file_path = "";
    
    // Default values for options
    uint number_gaussians = 512;
    uint ld_mode = 0;
    string shot_file_path = "";
    int shot_mode = 0;
    int shot_keyf = -1;
    int verbose_level = 1;
    string gdindex_parameters_path = "trained_parameters";
    bool single_shot = false;
    bool gd_intra_normalization = false;

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
            } else if ((!strcmp(argv[count_arg], "--shot_file_path")) || (!strcmp(argv[count_arg], "-s"))) {
                shot_file_path = string(argv[count_arg + 1]);
                count_arg++;
            } else if ((!strcmp(argv[count_arg], "--shot_mode")) || (!strcmp(argv[count_arg], "-m"))) {
                shot_mode = atoi(argv[count_arg + 1]);
                count_arg++;
            } else if ((!strcmp(argv[count_arg], "--shot_keyf")) || (!strcmp(argv[count_arg], "-k"))) {
                shot_keyf = atoi(argv[count_arg + 1]);
                count_arg++;
            } else if ((!strcmp(argv[count_arg], "--verbose_level")) || (!strcmp(argv[count_arg], "-v"))) {
                verbose_level = atoi(argv[count_arg + 1]);
                count_arg++;
            } else if ((!strcmp(argv[count_arg], "--single_shot")) || (!strcmp(argv[count_arg], "-g"))) {
                single_shot = true;
            } else if ((!strcmp(argv[count_arg], "--gd_intra_normalization"))) {
                gd_intra_normalization = true;
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
        cout << "------>gd_intra_normalization = " << gd_intra_normalization << endl;
        cout << "------>shot_file_path = " << shot_file_path << endl;
        cout << "------>single_shot = " << single_shot << endl;
        cout << "------>shot_mode = " << shot_mode << endl;
        cout << "------>shot_keyf = " << shot_keyf << endl;
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
                                 gd_intra_normalization,
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
    if ((shot_file_path == "") && (single_shot == false)) {
        // Using simple keyframe-based index
        gdindex.generate_index(feature_files, verbose_level);
    } else {
        // Read shot file
        vector<uint> shot_beg_frames;
        if (!single_shot) {
            assert(file_exists(shot_file_path));
            // Parse shots from shot_file_path
            parse_shot_file(shot_file_path, shot_beg_frames);
        } else {
            // In this case, we'll just use the entire database as a single shot.
            // So we just add it here
            shot_beg_frames.push_back(0);
        }
        if (verbose_level >= 3) cout << "Parsed " << shot_beg_frames.size() 
                                     << " shots." << endl;

        // Make sure last shot index is not greater than the number of frames in the
        // database; otherwise, it doesn't make sense!
        if (shot_beg_frames.size()) {
            assert(shot_beg_frames.at(shot_beg_frames.size() - 1) 
                   <= feature_files.size() - 1);
        }

        gdindex.generate_index_shot_based(feature_files, 
                                          shot_beg_frames, 
                                          shot_mode, 
                                          shot_keyf, 
                                          verbose_level);
        if (verbose_level >= 3) cout << "Index was generated." << endl;

        // In the indep-keyf mode using shots, we also need to store an additional file,
        // which will indicate which frame numbers (from the entire database) are the ones
        // in the database -- otherwise, there would be no way to know it
        if (shot_mode == GDIndex::SHOT_MODE_INDEP_KEYF) {
            if (verbose_level) 
                printf("Saving GDIndex frame list to:\n%s\n", (out_file_path + ".frmlist").c_str());
            gdindex.write_frame_list((out_file_path + ".frmlist"));        
        }
    }

    // Save index to file
    if (verbose_level) printf("Saving GDIndex to:\n%s\n", out_file_path.c_str());
    gdindex.write(out_file_path);

    return EXIT_SUCCESS;
}

