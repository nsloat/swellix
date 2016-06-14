#ifndef __SWELLIX_INITGEN_H
#define __SWELLIX_INITGEN_H

/*****************************************************************************************************************************
 *****************************************************************************************************************************
 *****************************************************************************************************************************
      Program : Swellix (Initialization - General) - h file
      Purpose : Initialize seq
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

int    initialize_commandline_argument_with_default_value(config* seq);
int    set_degree_of_detailed_disp(config* seq, int8_t arg, int argc, char** argv);
void   set_statMode(config* seq, int8_t arg, int argc, char** argv);
void   set_motif(config* seq, int8_t arg, int argc, char** argv);
int    set_mode_of_algo(config* seq, int8_t arg, int argc, char** argv);
void   cmd_line_err(char** argv, int arg);
int    init_seq_ltr(config* seq);

#endif // __SWELLIX_INITGEN_H
