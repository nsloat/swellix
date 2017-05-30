#ifndef __SWELLIX_MAKE_BUNDLE_LIST2_H
#define __SWELLIX_MAKE_BUNDLE_LIST2_H

#include "main.h"

/*****************************************************************************************************************************
 *****************************************************************************************************************************
 *****************************************************************************************************************************
      Program : Swellix (Make Bundle List2 - created by Isaac) - h file
      Purpose : Ryan's RNA Sequence Folding Algorithm (modified Crumple)
      Input   : manual key-in or input file
      Return  : none
      Display : pair displays in multiple-bracket/parenthesis or alphabet layouts
              : non - only number of qualified pair structures, duplicates found and time consumed, and error messages
              : min - only non-level display and all qualified pair structures
              : prn - only min-level display and parenthesis-related display, for debugging purpose
              : brs - only prn-level display and brace-related display, for debugging purpose
              : all - every possible display, for detailed debugging purpose
      CAUTION : integer capability of 16 bit signed carries positive integers of 2^16/2 -1 = 32767
              :             "         32 bit          "         "        "       2^32/2 -1 = 2147483647
              :             "         64 bit          "         "        "       2^64/2 -1 = 9223372036854775807
       NOTE 1 : COMPONENT LIST                          STRUCTURE
              : seq             noSZG         SZG       noSZG                SZG
              : 14mer           21            35        118                  182
              : 28mer           96            156       124925
              : 42mer           227           374       197316084
              : tRNA-R1660.raw  835           1255
              : ecoli-5s.seq    2634          4126
              : stmv.seq        216207        358812
       NOTE 2 : helix find syntax ./find.py < /home/liu5227/code/seq/stmv.seq -s 9 | wc
 *****************************************************************************************************************************
 *****************************************************************************************************************************
 *****************************************************************************************************************************/

void initLabeledStructures(LabeledStructures* lab);
void resetLabeledStructures(LabeledStructures* lab, char* newtitle);
void freeLabeledStructures(LabeledStructures** lab);

void filename_to_indices(char* filename, int* position, int* length);
void run_sliding_windows(config* seq, global* crik);
int reset_dot_n_paren(config* seq);
int mark_paren(config* seq, knob* BNodeCursr);
int mark_mismatch(config* seq, knob* BNodeCursr);
int prt_partial_stru(config* seq, knob* BNodeCursr);
int prt_covari_n_V1_flag(config* seq, ToL* SNodeCursr);
int prt_hlix_info(config* seq, knob* BNodeCursr);
int prt_individual_partial_stru_on_bundle(config* seq, ToL* SNodeCursr);
int disp_dumi_node(config* seq, global* crik, int16_t bundleLB, int16_t bundleUB);
int disp_bundle_list(config* seq, global* crik);
void get_outer_pair_indices(char* structure, int len, int* start, int* end);
void remove_outer_pair(char* structure, int start, int end);
void make_bundles(config* seq, global* crik, LabeledStructures* lab);
void add_dumi_node(config* seq, global* crik, LabeledStructures* lab);
knob* init_b_node(int openOut, int openIn, int closeOut, int closeIn, int16_t** mismatch);
void add_s_b_node(global* crik, knob* refBNode, int16_t masterLB, int16_t masterUB, int8_t mms);
void add_b_node(global* crik, knob* refBNode, int16_t masterLB, int16_t masterUB);
int get_mismatches(char* structure, int openOut, int openIn, int closeOut, int closeIn, int16_t** mismatch);
int mark_top_notches_n_ltr(config* seq);
int mark_side_notches_n_ltr(config* seq, int16_t row);
int prt_ltr_n_num_inside_table(config* seq, global* crik, int16_t lukUpLB, int16_t row, int16_t bundleCnt);
int rev_order(config* seq, global* crik);
void remove_duplicates(char* filename);
int remove_duplicates_LL(ToL* start);
int nodes_equal(knob* first, knob* second);

#endif
