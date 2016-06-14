#ifndef __SWELLIX_INITCONST_H
#define __SWELLIX_INITCONST_H

/*****************************************************************************************************************************
 *****************************************************************************************************************************
 *****************************************************************************************************************************
      Program : Swellix (Initialize Constraints) - h file
      Purpose : Read in constraints from configuration files *.swlx, and initialize seq->coVari, seq->v1Pairng, seq->s1Pairng, seq->chemMod...
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

int    init_seq_constraint(config* seq);
int    init_covari_db(config* seq, FILE* fp);
int    init_s1_pairng_db(config* seq, FILE* fp);
int    init_v1_pairng_db(config* seq, FILE* fp);
int    init_max_pairng_dist(config* seq, FILE* fp);
int    init_min_pairng_dist(config* seq, FILE* fp);
int    read_in_covari_data(config* seq, FILE* fp, char* constraints, int i);
int    read_in_s1_pairng_data(config* seq, FILE* fp, char* constraints, int i);
int    read_in_v1_pairng_data(config* seq, FILE* fp, char* constraints, int i);
int    read_in_max_pairng_dist_data(config* seq, FILE* fp, char* constraints);
int    read_in_min_pairng_dist_data(config* seq, FILE* fp, char* constraints);
int    is_contradict_2_pairng_dist_constraints(config* seq, int i);
int    is_contradict_2_s1_pairng_constraints(config* seq, int i);
int    is_covarance_pair_i_defective(config* seq, int coVariI);
int    is_there_duplicate(config* seq, int coVariI, int coVariJ);
int    is_there_sudona_relationship(config* seq, int coVariI, int coVariJ);
int    sort_covari_pairs_in_ascend_order(config* seq);
int    chk_conflict_against_v1_pairng(config* seq);
int    display_constraints(config* seq);

int*** make_buckets(config* seq __unused, int numOfBaki);
int    toss_into_buckets(config* seq, int*** bucket);
int    display_full_buckets(config* seq, int*** bucket, int numOfBaki);
int    move_from_buckets_bak_2_covari(config* seq, int*** bucket, int numOfBaki);
int    shift_over(int*** bucket, int tag, int j);

#endif // __SWELLIX_INITCONST_H
