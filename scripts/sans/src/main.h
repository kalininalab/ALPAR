#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <thread>
#include <regex>
#include <mutex>

#include "util.h"
#include "translator.h"
#include "cleanliness.h"

#ifdef useBF
    #include <bifrost/CompactedDBG.hpp>
    #include <bifrost/ColoredCDBG.hpp>

    #ifndef MAX_KMER_SIZE
        #define MAX_KMER_SIZE (maxK + 1)
    #endif
#endif



using namespace std;

// SANS ambages
// Symmetric Alignment-free phylogeNomic Splits
// phylogenomics with Abundance-filter, Multi-threading and Bootstrapping on Amino-acid or GEnomic Sequences
#define SANS_VERSION "2.4_01A"    // SANS ambages

/**
 * This is the entry point of the program.
 *
 * @param argc number of cmd args
 * @param argv cmd args
 * @return exit status
 */ 
int main(int argc, char* argv[]);

/** This function applies the specified filter to the given split list.
 * 
 * @param filter string specifying the type of filter
 * @param newick string with a file name to write a newick output to (or empty string)
 * @param map function that maps an integer to the original id, or null
 * @param split_list the list of splits to be filtered, e.g. graph::split_list
 * (@param support_values a hash map storing the absolut support values for each color set)
 * (@param bootstrap_no the number of bootstrap replicates for computing the per centage support)
 * @param verbose print progress
 * 
 */
void apply_filter(string filter, string newick, std::function<string(const uint64_t&)> map, multimap_<double, color_t>& split_list, bool verbose);
void apply_filter(string filter, string newick, std::function<string(const uint64_t&)> map, multimap_<double, color_t>& split_list, hash_map<color_t, uint32_t>* support_values, const uint32_t& bootstrap_no, bool verbose);
