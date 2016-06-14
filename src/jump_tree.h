#ifndef __SWELLIX_MAKE_JUMP_TREE_H
#define __SWELLIX_MAKE_JUMP_TREE_H

#include "main.h"

/*****************************************************************************************************************************
 *****************************************************************************************************************************
 *****************************************************************************************************************************
      Program : Swellix (Make Jump Tree) - h file
      Purpose : Ryan's RNA Sequence Folding Algorithm (modified Crumple)
      Input   : manual key-in or input file
      Return  : none
      Display : pair displays in multiple-bracket/parenthesis or alphabet layouts
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
knob* copy_b_node_wo_links(knob* node);
knob* copy_b_node(knob* node);
void release_knob(knob* torelease);
void release_tol(ToL* torelease);
ToL* copy_tol(ToL* source);
ToL* combine(ToL* first, ToL* second);
int countAndPrintStructures(config* seq, global* crik);
int assign_new_hlix_link(global* crik, int16_t hlixBranchngIndx1);
int clear_old_hlix_link(global* crik);
int satisfy_constraint(config* seq, global* crik, local* todd);
int eval_prog_n_prt_stru(config* seq, global* crik);
int exit_curr_recur(config* seq, global* crik, local* todd);
int hook_dumi_node_on_stru(global* crik, local* todd, int16_t hlixBranchngIndx1);
local* init_todd(global* crik);
int        insert_beh_intrvl_normally(config* seq __unused, global* crik, local* todd);
int        bundle_is_available(config* seq, global* crik, local* todd);
int        is_duplicate_intrvl(global* crik, int16_t intrvlLB, int16_t intrvlUB);
int        is_intrvl_2b_rsto(global* crik, local* toddP, int8_t recurRoute);
int        time_to_quit(config* seq, local* todd);
int        is_there_bundle_somewhere_inside_this_hlix(config* seq, global* crik);
int        is_there_concern_on_ring_formation_of_linked_list(config* seq __unused, global* crik, local* todd);
int        is_there_no_worry_about_cmpnt_2_insert_behind(global* crik, knob* cmpntCursr);
int        is_there_prev_dumi_2_replace_curr_candidate_dumi_recur_lv_0(config* seq, global* crik, local* todd, int16_t cmpntOcupidID);
int        jump_stage_1_set_intrvl(config* seq, global* crik, local* todd, int16_t bundlePathFlag);
int        jump_stage_2_fit_hlix(config* seq, global* crik, local* toddP, int8_t recurRoute);
int        is_there_larger_dumi_2_replace_curr_candidate_cmpnt_with_concern_on_cmpnt_insert_behind(global* crik, knob* cmpntCursr, int16_t iCmpntTyp, int16_t cmpntOcupidID);
int        is_there_larger_dumi_2_replace_curr_candidate_cmpnt_wo_concern_on_cmpnt_insert_behind(config* seq, global* crik, knob* cmpntCursr, int16_t bundleSpan, int16_t iCmpntTyp);
int        remove_intrvl(config* seq, global* crik, local* toddP, int8_t recurRoute);
int        search_4_largest_rep_of_this_bundle_group(global* crik, local* todd, int16_t cmpntOcupidID);
int        swap_ins_n_beh_intrvl(config* seq __unused, global* crik, local* todd);
int        take_bundle_list_shortcut(config* seq, global* crik, local* todd, int16_t hlixBranchngIndx1);
int        take_cmpnt_list_normal_path(config* seq, global* crik, local* todd, int16_t hlixBranchngIndx1);

#endif // __SWELLIX_MAKE_JUMP_TREE_H
