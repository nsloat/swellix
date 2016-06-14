#ifndef __SWELLIX_CLOSE_UP_H
#define __SWELLIX_CLOSE_UP_H

/*****************************************************************************************************************************
 *****************************************************************************************************************************
 *****************************************************************************************************************************
      Program : Swellix (Close Up) - h file
      Purpose : Free dynamically allocated memory spaces
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

int free_bundle_list(config* seq, global* crik);
int free_cmpnt_list(config* seq, global* crik);
int free_constraint(config* seq, global* crik);
int free_intrvl_luk_up_table(config* seq);
int free_parenthesis_look_up_table(config* seq);
int free_recycle_bin(config* seq);
int free_whole_intrvl(global* crik);

#endif // __SWELLIX_CLOSE_UP_H
