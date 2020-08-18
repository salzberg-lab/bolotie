#include <iostream>
#include <cstring>

#include "arg_parse.hpp"

#include "find_best_path.hpp"
#include "build_probability_table.hpp"
#include "parents.hpp"
#include "aln.hpp"
#include "cons.hpp"

void print_help()
{
    std::cerr << "bolotie\n"
              << "Modes:\n"
              << "\tbuild - build a probability table for a set of clustered sequences\n"
              << "\tfind - find the best path through the probability table for a given set of sequences\n"
              << "\tdistmat - computes distance matrix for a given MSA\n"
              << "\tclean_clusters - removes outliers from a cluster based on distance to the consensus\n"
              << "\tparents - extract potential parents of recombinants based on seq. similarity\n"
              << "\taln - align sequences to the references and call variants\n"
              << "\tcons - generate consensus sequences based on the called variants\n";
}

int main(int argc, char ** argv)
{
    if (argc <= 1)
    {
        print_help();
    }
    else if (strcmp(argv[1], "build") == 0)
    {
        std::cerr << "building probability table\n";
        const int argc_build = argc - 1;
        char* argv_build[argc_build];
        memcpy(argv_build, argv + 1, argc_build * sizeof(char*));
        bolotie_build(argc_build, argv_build);
    }
    else if(strcmp(argv[1], "find") == 0)
    {
        std::cerr << "finding best path through the probability table\n";
        const int argc_find = argc - 1;
        char* argv_find[argc_find];
        memcpy(argv_find, argv + 1, argc_find * sizeof(char*));
        bolotie_find(argc_find, argv_find);
    }
    else if(strcmp(argv[1], "parents") == 0)
    {
        std::cerr << "output parent sequences of recombinant\n";
        const int argc_parents = argc - 1;
        char* argv_parents[argc_parents];
        memcpy(argv_parents, argv + 1, argc_parents * sizeof(char*));
        bolotie_parents(argc_parents, argv_parents);
    }
    else if(strcmp(argv[1], "aln") == 0)
    {
        std::cerr << "aligning sequences\n";
        const int argc_aln = argc - 1;
        char* argv_aln[argc_aln];
        memcpy(argv_aln, argv + 1, argc_aln * sizeof(char*));
        bolotie_aln(argc_aln, argv_aln);
    }
    else if(strcmp(argv[1], "cons") == 0)
    {
        std::cerr << "generating consensus sequences\n";
        const int argc_cons = argc - 1;
        char* argv_cons[argc_cons];
        memcpy(argv_cons, argv + 1, argc_cons * sizeof(char*));
        bolotie_cons(argc_cons, argv_cons);
    }
    else
    {
        std::cerr << "Unrecognized Mode: please consult the help manual\n";
        print_help();
    }

    return 0;
}
