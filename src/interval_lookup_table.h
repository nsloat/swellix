/*****************************************************************************************************************************
 *****************************************************************************************************************************
 *****************************************************************************************************************************
      Program : Swellix (Make Interval Look-up Table) - h file
      Purpose : Make interval look-up table
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

int format_the_intrvl_luk_up_table_prt_out(config* seq);
int fill_in_intrvl_luk_up_table_row_i(config* seq, global* crik, int16_t rowI);
int fill_in_intrvl_luk_up_table_cell_i_j(config* seq, global* crik, int16_t rowI, int16_t rowJ);
int init_i_j_cell(config* seq, int16_t rowI, int16_t colJ);
int lowerbound_fits_in(config* seq, global* crik, int16_t rowI, int16_t colJ, int16_t ocupidTypIndx);
int upperbound_fits_in(config* seq, global* crik, int16_t rowI, int16_t colJ, int16_t ocupidTypIndx);
int is_covarance_pair(config* seq, int16_t cmpntTyp);
int mark_reading_end_of_this_cell(config* seq, int16_t rowI, int16_t colJ, int16_t cmpntTyp);
int reach_end_of_loop(config* seq, global* crik, int16_t rowI, int16_t colJ, int16_t ocupidTypIndx);
