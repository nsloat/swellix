#ifndef __SWELLIX_CMPNT_H
#define __SWELLIX_CMPNT_H

/*****************************************************************************************************************************
 *****************************************************************************************************************************
 *****************************************************************************************************************************
      Program : Swellix (Make Component List) - h file
      Purpose : Make Component List
      Input   : config* seq
              : global* crik
      Return  : none
      Display : none
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

int   assign_window_width(config* seq, global* crik, int callerType);
int   book_out(config* seq, global* crik);
int   check_bundle_max_width(config* seq, global* crik);
char* choose_brace_display_format(config* seq);
int   set_display_priority(config* seq, int dispCmpltFlag);
void  display_component_list(config* seq, global* crik, int dispCmpltFlag, int dispPrio, char* braceBound);
void  display_components(config* seq, global* crik, int dispCmpltFlag);
int   hook_new_node_on_component_list(config* seq __unused, global* crik, knob* newNode);
int   initialize_mismatch(config* seq, knob* newNode);
int   initialize_new_node(config* seq __unused, global* crik, knob* newNode, int16_t* coVariFlag);
int   base_pair_NOT_formed(config* seq, global* crik);
int   brace_is_formed(config* seq, global* crik, knob* newNode, int16_t* coVariFlag);
int   mismatch_by_cheating(config* seq, int16_t rowI, int16_t colJ);
int   mismatch_causes_lonely_pair(config* seq, global* crik, knob* newNode, int16_t cntOfMismatch);
int   mismatch_interrupt_covariance(config* seq, global* crik);
int   NOT_terminal_mismatch(config* seq, global* crik, int16_t cnt_of_mismatch);
int   is_s1(config* seq, int16_t rowI, int16_t colJ);
int   is_worth_exploring(global* crik, int8_t noBrsFlag);
int   recur_stage_1_brace_locatng(config* seq, global* crik);
int   recur_stage_2_brace_sizng(config* seq, global* crik);
int   recur_stage_3_brace_processng(config* seq, global* crik);
int   update_v1_covari_flag(config* seq, global* crik, int16_t* mustPairFlag);

#endif // __SWELLIX_CMPNT_H
