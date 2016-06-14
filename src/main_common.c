#include <stdlib.h>
#include <stdio.h>

#include "main.h"

/***
  functions used by main but also used elsewhere
  code written by Ryan Liu
  moved into a separate file by Jonathan Stone to facilitate unit testing
  ***/

//***********************************************************************
// Function : Print Usage
// Caller   : main()
// Purpose  : print out the command line options of this program, and exit the program thereafter
// Input    : none
// Return   : none
// Display  : the complete list of command line options of this program
//***********************************************************************
void print_usage(void)
{
printf(
"\n\n Usage: swellix  [-#]  [-b] [-u] [-bhdm]  [-d [1~5 or a]]  [-h]  [-i INFILE]  [-k]  [-l]  [-m [c, i, n, s]]  [-mm]  [-noGU]  [-noSZG] \n"
"                     [-o OUTFILE]  [-p] [-s [1-3 or a]] [-M [MOTIF]]\n\n\n"

" Compact Option  Complete Option                        Description                                                           Default Value\n\n"
" -# | -hlx#m | --minimum_number_of_helices    (count) : Specify minimum number of helices to form a complete pair structure   [0]          \n\n"
" -b | -bundle|                                        : Activate the bundling mechanism to reduce similar structures          [FALSE]      \n\n"
" -u | -unbund|                                        : Unbundling of bundles for printout of complete solution set           [FALSE]      \n\n"
"      -bhdm  | --minimum_btwn_helix_distance  (size)  : Minimum between helix distance. used primarily in STMV                [0]          \n\n"
" -d | -dsply | --degree_of_detailed_display           : Print option. (only meaningful if compile cmd make dispOn is used)    [1]            \n"
"                                          all     | a : Print all, for detailed debugging purpose.                                           \n"
"                                          level 5 | 5 :                                                                                      \n"
"                                          level 4 | 4 :                                                                                      \n"
"                                          level 3 | 3 : Print component list, jump bush and strucutures                                      \n"
"                                          level 2 | 2 : Print minimum, only pair structures will be prints out                               \n"
"                                          level 1 | 1 : Print none, but statistics data and error, if any                                  \n\n"                                
" -s | -stat  |                                        : Sets which statistics are to be performed during computation          [1]            \n"
"                                          all     | a : Print all available statistics.                                                      \n"
"                                                  | 3 : Print the max free energy found in the computed structures.                          \n"
"                                                  | 2 : Print the max RNA distance found in the computed structures.                         \n"
"                                                  | 1 : Default action (currently do nothing except print number of structures)            \n\n"
" -M | -motif |                                        : Provided letter/dot&parenthesis pair, identify number of occurrences  [FALSE]        \n"
"                                                      : Give argument in the form \"GCUCUAACAGC&(((.....)))\"                                \n"
"                                                        with the '&' character dividing the sequence segment from the                        \n"
"                                                        structure to match with it.                                                        \n\n"
" -h | -help                                           : Prints this message and exits.                                                     \n\n"
" -i FILE                                              : Defines an input file, containing a sequence.                         [stdin]      \n\n"
" -k |        | -constraint                            : Introduce constraints from configuration file                         [FALSE]        \n"
"                                                      : Only give number of covariance pairs in command line, no detail                    \n\n"
" -l | -hlxlm | --minimum_length_of_helix      (size)  : Specify minimum number of nucleotides to form a complete helix        [1]          \n\n"
" -m | -mode  | --mode_of_algorithm                    : Specify the way to do calculation                                     [stru]         \n"
"                 component              | cmpnt | c   : Compute for component only                                                           \n"
"                 interval_look_up_table | intab | i   : For interval look-up table level                                                     \n"
"                 none                           | n   : Do nothing                                                                           \n"
"                 structure              | stru  | s   : Compute for all structures                                                         \n\n"                                
"      -mm    | --maximum_number_of_mismatches (count) : Specify Maximum number of mismatches inside any helix                 [0]          \n\n"
"      -noGU  | --no_G_U_allowed                       : Disallows G-U pairing                                                 [FALSE]      \n\n"
"      -noSZG | --no_sizing_allowed                    : Don't extend (swell) the size of the helix unless bulge presents      [FALSE]      \n\n"
" -o FILE                                              : Defines an output file, to overwrite with solutions.                  [stdout]     \n\n"
" -p | -hp#m  | --minimum_number_of_hairpins   (count) : Specify minimum number of hairpins to form a complete pair structure  [0]          \n\n");
 exit(1);
}  // end print_usage
