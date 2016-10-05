/**********************************************
This program binarizes a point-indexed index.
 **********************************************/

#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

#include "binarize_residuals.h"
#include "point_index_io.h"

// Usage of this program
void usage() {
    printf("usage: ./binarize_index --input[-i] input_path --output[-o] output_path \n");
    printf("For an example, see run_binarize_pi_descriptors.sh\n");
    printf("Options:\n");
    printf("--verbose_level[-v] ARG (default: 0) \n");
    printf("NOTE: residuals must be 32-dim. length! \n");
}

int main(int argc, char** argv) {
    // Mandatory arguments
    string input_path = "";
    string output_path = "";
    
    // Default values for options
    int verbose_level = 0;

    if (argc < 3) {
        cout << "Wrong usage!!!" << endl;
        cout << "***********************************" << endl;
        cout << "***********************************" << endl;
        cout << "***********************************" << endl;
        cout << "See usage below:" << endl;
        usage();
        exit(EXIT_FAILURE);
    } else {
        for (uint count_arg = 1; count_arg < argc; count_arg++) {
            if ((!strcmp(argv[count_arg], "--input")) || (!strcmp(argv[count_arg], "-i"))) {
                input_path = string(argv[count_arg + 1]);
                count_arg++;
            } else if ((!strcmp(argv[count_arg], "--output")) || (!strcmp(argv[count_arg], "-o"))) {
                output_path = string(argv[count_arg + 1]);
                count_arg++;
            } else if ((!strcmp(argv[count_arg], "--verbose_level")) || (!strcmp(argv[count_arg], "-v"))) {
                verbose_level = atoi(argv[count_arg + 1]);
                count_arg++;
            } else {
                cout << "Incorrect argument!" << endl;
                cout << "See usage below:" << endl;
                usage();
                exit(EXIT_FAILURE);
            }        
        }
    }

    // Check that all mandatory arguments were provided
    if (input_path == "") {
        cout << "input argument was not provided... see usage below:" << endl;
        usage();
        exit(EXIT_FAILURE);
    }
    if (output_path == "") {
        cout << "output argument was not provided... see usage below:" << endl;
        usage();
        exit(EXIT_FAILURE);
    }

    if (verbose_level) {
        cout << "Starting index building using:" << endl;
        cout << "------>input_path = " << input_path << endl;
        cout << "------>output_path = " << output_path << endl;
        cout << "------>verbose_level = " << verbose_level << endl;
    }

    // Index data that that will be read
    vector < vector < uint > > vec_feat_assgns;
    vector < vector < float > > vec_feat_assgn_weights;
    vector < vector < vector < float > > > vec_feat_residuals;  

    // Binarized residuals that will be output
    vector < vector < uint > > vec_feat_residuals_binarized;

    // Load input signatures
    if (verbose_level) cout << "Loading input index..." << endl;
    read_fv_point_index(input_path,
                        vec_feat_assgns,
                        vec_feat_assgn_weights,
                        vec_feat_residuals);
    if (verbose_level) cout << "done!" << endl;

    // Binarize residuals
    if (verbose_level) cout << "Binarizing residuals..." << endl;
    binarize_residuals(vec_feat_residuals, vec_feat_residuals_binarized);
    if (verbose_level) cout << "done!" << endl;

    // Write binarized index
    if (verbose_level) cout << "Writing binarized index..." << endl;
    write_bfv_point_index(vec_feat_assgns,
                          vec_feat_assgn_weights,
                          vec_feat_residuals_binarized,
                          output_path);
    if (verbose_level) cout << "done!" << endl;

    return EXIT_SUCCESS;
}

