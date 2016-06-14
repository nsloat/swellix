/*****************************************************************************************************************************
 *****************************************************************************************************************************
 *****************************************************************************************************************************
      Program : Swellix (Initialize Constraints) - c file
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "main.h"
#include "init_constraint.h"

//*****************************************************************************
// Function : Initialize Covariance Database
// Caller   : initialize_sequence_by_getting_argument()
// Purpose  : set up the data related to covariance, and check any violation
// Input    : 
// Return   : 1: start over, 0: function ends properly
// Display  : inquiry for user input
//*****************************************************************************
int init_covari_db(config* seq, FILE* fp)
{
  char* constraints = calloc(CONSTRAINT_SIZE, sizeof(char));
  int   i = 0;
  int   coVariI, coVariJ, temp;
                                                                  // printf("chkpnt0\n");
  while( (fscanf(fp, "%s", constraints) != EOF) &&  (constraints[0] != '_') ){                               // verify '(' is present, regard it as hint for a new pair
    if(constraints[0] == '('){                                    // printf("chkpnt1\n");
      seq->coVari[i] = malloc(sizeof(int) * 2);
      read_in_covari_data(seq, fp, constraints, i);
      i++;
    }  // end if
  }    // end while

  seq->coVari = realloc(seq->coVari, i * sizeof(seq->coVari));
  seq->numCovari = i;

  for(coVariI = 0 ; coVariI < seq->numCovari ; coVariI++){

    if (seq->coVari[coVariI][0] > seq->coVari[coVariI][1]){             // [0] holds larger number than [1] does, so swap it!
      temp = seq->coVari[coVariI][0]; 
      seq->coVari[coVariI][0] = seq->coVari[coVariI][1];
      seq->coVari[coVariI][1] = temp;
    }  // end if 1

    if(is_covarance_pair_i_defective(seq, coVariI))
       exit(1);

    for(coVariJ = (coVariI + 1) ; coVariJ < seq->numCovari ; coVariJ++){
      if(is_there_duplicate(seq, coVariI, coVariJ))
        exit(1);
      if(!seq->sudona && is_there_sudona_relationship(seq, coVariI, coVariJ))
        exit(1);
    }  // end inner for
  }    // end outer for

  sort_covari_pairs_in_ascend_order(seq);

  free(constraints);

  return 0;
}  // end init_covari_db

//*****************************************************************************
// Function : Read In Covariance Data
// Caller   : init_covari_db()
// Purpose  : read in the covariance pair data from configuration file
//          : if there's error, pop it up, and exit the process
// Input    : 
// Return   : none
// Display  : none
//*****************************************************************************
int read_in_covari_data(config* seq, FILE* fp, char* constraints, int i)
{                                                                           disp(seq,DISP_ALL,"Enterng 'read_in_covari_data'\n");
  int j;
                                                                            // printf("chkpnt2\n");
  if(fscanf(fp, "%s\n", constraints) != '0' &&  atoi(constraints) > 0)      // check for a number
    seq->coVari[i][0] = atoi(constraints);
  else if(!strcmp("0", constraints))
    seq->coVari[i][0] = 0;
  else{                                                                     printf("Error in covariance configuration file\n");
                                                                            printf("%s is not a number\n", constraints);
    exit(1);
  }  // end if 1

  if((constraints[0] = fgetc(fp)) != ','){                                  // verify ',' is present
                                                                            printf("Error in covariance configuration file\n");
                                                                            printf("%c is not desinated symbol ','\n", constraints[0]);
    exit(1);
  }  // end if 2

  if(fscanf(fp, "%s\n", constraints) != '0' && atoi(constraints) > 0)       // check for a number
    seq->coVari[i][1] = atoi(constraints);
  else if(!strcmp("0", constraints))
    seq->coVari[i][1] = 0;
  else{                                                                     printf("Error in covariance configuration file\n");
                                                                            printf("%s is not a number\n", constraints);
    exit(1);
  }  // end if 3

  switch(is_contradict_2_pairng_dist_constraints(seq, i)){
    case 2:                                                                 fprintf(stderr, "\n\nContradict between covari (%d,%d) & Max pairng dist %d\n\n", seq->coVari[i][0], seq->coVari[i][1], seq->maxPairngDist);
     exit(1);
    case 1:                                                                 fprintf(stderr, "\n\nContradict between covari (%d,%d) & min pairng dist %d\n\n", seq->coVari[i][0], seq->coVari[i][1], seq->minPairngDist);
      exit(1);
  }  // end if 4

  if(is_contradict_2_s1_pairng_constraints(seq, i)){                        fprintf(stderr, "\n\nContradict between covari (%d,%d) & S1 pairng: ", seq->coVari[i][0], seq->coVari[i][1] );
    for(j = 0 ; j < seq->numS1Pairng ; j++)                                 fprintf(stderr, "%d ", seq->s1Pairng[j]);
                                                                            fprintf(stderr, "\n\n");
    exit(1);
  }  // end if 4

  if((constraints[0] = fgetc(fp)) != ')'){                                  // verify ')' is present
                                                                            printf("Error in covariance configuration file\n");
                                                                            printf("%c is not desinated symbol ')'\n", constraints[0]);
    exit(1);
  }  // end if 6

  return 0;
}  // end read_in_covari_data

//*****************************************************************************
// Function : Is Contradict to Pairing Distance Constraints?
// Caller   : read_in_covari_data()
// Purpose  : check if the distance of given covariance pair is too wide/narrow
//          : when comparing to the Max/min pairing distance
// Input    : 
// Return   : 2: larger than Max dist, 1: smaller than min dist, 0: no
// Display  : none
//*****************************************************************************
int is_contradict_2_pairng_dist_constraints(config* seq, int i)
{
  int16_t pairngDist = seq->coVari[i][1] - seq->coVari[i][0] - 1;
  if(pairngDist > seq->maxPairngDist) return 2;
  if(pairngDist < seq->minPairngDist) return 1;
  return 0;
}  // end is_contradict_2_pairng_dist_constraints

//*****************************************************************************
// Function : Is Contradict to S1 Pairing Constraints?
// Caller   : read_in_covari_data()
// Purpose  : check if the element(s) of the given covari pair is on the S1 pair list
// Input    : 
// Return   : 1: yes, 0: no
// Display  : none
//*****************************************************************************
int is_contradict_2_s1_pairng_constraints(config* seq, int i)
{
  int8_t j = seq->numS1Pairng;

  while(j){                                                        // printf("%d vs (%d,%d)\n", seq->s1Pairng[j-1], seq->coVari[i][0], seq->coVari[i][1]);
    if( (seq->s1Pairng[j-1] == seq->coVari[i][0]) || (seq->s1Pairng[j-1] == seq->coVari[i][1]) ) return 1;
    j--;
  }  // end while

  return 0;
}  // end is_contradict_2_pairng_dist_constraints

//*****************************************************************************
// Function : Is Covariance Pair I Defective?
// Caller   : init_covari_db()
// Purpose  : check if there's defects between the pair i, j
// Input    : 
// Return   : 1: yes, 0: no
// Display  : none
//*****************************************************************************
int is_covarance_pair_i_defective(config* seq, int coVariI)
{
  int   coVariIUB = seq->strLen - seq->minPairngDist - 1;
  int   coVariILB = seq->minPairngDist;

  if (seq->coVari[coVariI][0] == seq->coVari[coVariI][1]){                                              fprintf(stderr, "Error on pair #%d: (%d, %d)\n", (coVariI + 1), seq->coVari[coVariI][0], seq->coVari[coVariI][1]);
 												        fprintf(stderr, "Two numbers should be different.\n");
    return 1;
  } else if ( (seq->coVari[coVariI][1] - seq->coVari[coVariI][0]) <= seq->minPairngDist){               fprintf(stderr, "Error on pair #%d: (%d, %d)\n", (coVariI + 1), seq->coVari[coVariI][0], seq->coVari[coVariI][1]);
          	                                                                                        fprintf(stderr, "Two numbers should be separated by at least minimum loop size of %d.\n", seq->minPairngDist);
    return 1;
  } else if ( (seq->coVari[coVariI][0] < 0) || (seq->coVari[coVariI][0] >= coVariIUB) ){                fprintf(stderr, "Error on pair #%d: (%d, %d)\n", (coVariI + 1), seq->coVari[coVariI][0], seq->coVari[coVariI][1]);
         	                                                                                        fprintf(stderr, "The input number %d is improper. ", seq->coVari[coVariI][0]);
													fprintf(stderr, "It should be equal or larger than 0, and smaller than %d.\n", coVariIUB);
    return 1;
  } else if ( (seq->coVari[coVariI][1] <= coVariILB) || (seq->coVari[coVariI][1] >= seq->strLen) ){     fprintf(stderr, "Error on pair #%d: (%d, %d)\n", (coVariI + 1), seq->coVari[coVariI][0], seq->coVari[coVariI][1]);
         	                                                                                        fprintf(stderr, "The input number %d is improper.", seq->coVari[coVariI][1]);
													fprintf(stderr, "It should be equal or larger than %d, and smaller than %d.\n", coVariILB, seq->strLen);
    return 1;
  }  // end if 2

  return 0;
}  // end is_covarance_pair_i_defective

//*****************************************************************************
// Function : Is Duplicate of Previous?
// Caller   : init_covari_db()
// Purpose  : check if there's duplicates of number between current pair and previous pairs
// Input    : 
// Return   : 1: yes, 0: no
// Display  : none
//*****************************************************************************
int is_there_duplicate(config* seq, int coVariI, int coVariJ)
{
  if( (seq->coVari[coVariI][0] == seq->coVari[coVariJ][0]) || 
      (seq->coVari[coVariI][0] == seq->coVari[coVariJ][1]) ||
      (seq->coVari[coVariI][1] == seq->coVari[coVariJ][0]) || 
      (seq->coVari[coVariI][1] == seq->coVari[coVariJ][1])    ){
    fprintf(stderr, "There is duplicate between the pair #%d (%d, %d) and the pair #%d (%d, %d).\n\n", (coVariI + 1), seq->coVari[coVariI][0], seq->coVari[coVariI][1], (coVariJ + 1), seq->coVari[coVariJ][0], seq->coVari[coVariJ][1]);
    return 1;
  }  // end if

  return 0;
}  // end is_there_duplicate

//*****************************************************************************
// Function : Is Pseudoknot Relationship with Previous?
// Caller   : init_covari_db()
// Purpose  : check if there's duplicates of number between current pair and previous pairs
// Input    : 
// Return   : 1: yes, 0: no
// Display  : none
//*****************************************************************************
int is_there_sudona_relationship(config* seq, int coVariI, int coVariJ)
{
  int a, b, c, d;

  a = seq->coVari[coVariI][0] - seq->coVari[coVariJ][0];
  b = seq->coVari[coVariI][0] - seq->coVari[coVariJ][1];
  c = seq->coVari[coVariI][1] - seq->coVari[coVariJ][0];
  d = seq->coVari[coVariI][1] - seq->coVari[coVariJ][1];

  if( (a*b < 0 && c*d > 0) || (a*b > 0 && c*d < 0) ){
    fprintf(stderr, "There is pseudoknot relationship between the pair #%d (%d, %d) and the pair #%d (%d, %d).\n\n", (coVariI + 1), seq->coVari[coVariI][0], seq->coVari[coVariI][1], (coVariJ + 1), seq->coVari[coVariJ][0], seq->coVari[coVariJ][1]);
    return 1;
  }  // end if

  return 0;
}  // end is_there_sudona_relationship

//*****************************************************************************
// Function : Sort Covariance Pairs in Ascending Order
// Caller   : init_covari_db()
// Purpose  : reorder all the covariance pairs (i, j) in the ascending order based on 'i' value
//          : sorting algorithm: bucket sort
// Input    : seq
// Return   : none
// Display  : none
//*****************************************************************************
int sort_covari_pairs_in_ascend_order(config* seq)
{
  int    numBaki = seq->strLen / BAKI_SORT_DIVIDR + 1;
  int*** bucket = make_buckets(seq, numBaki);

  toss_into_buckets(seq, bucket);
  dispFullBaki(seq, bucket, numBaki);
  move_from_buckets_bak_2_covari(seq, bucket, numBaki);

  return 0;
}  // end sort_covari_pairs_in_ascend_order

//*****************************************************************************
// Function : Make Buckets
// Caller   : sort_covari_pairs_in_ascend_order()
// Purpose  : assign memory space to perform bucket sort
// Input    : seq
// Return   : bucket - n by 8 by 2 3D matrix containing the covariance data after this particular shift
// Display  : none
//*****************************************************************************
int*** make_buckets(config* seq __unused, int numBaki)
{
  int*** bucket = malloc(sizeof(*bucket) * numBaki);            // make several empty buckets
  int i, j;

  for(i = 0 ; i < numBaki ; i++){
    bucket[i] = malloc(sizeof(bucket[i]) * BAKI_SORT_DIVIDR);   // locate several seats inside each empty buckets

    for(j = 0 ; j < BAKI_SORT_DIVIDR ; j++){                    // divide each seat by two to fill in the covariance pair (i, j)
      bucket[i][j] = calloc(2, sizeof(int));                    disp(seq,DISP_ALL,"(%3d,%-3d) ", bucket[i][j][0], bucket[i][j][1]);
    }  // end inner for
                                                                disp(seq,DISP_ALL,"\n");
  }    // end outer for
                                                                disp(seq,DISP_ALL,"\n");
  return bucket;
}  // end make_buckets

//*****************************************************************************
// Function : Toss into Bucket
// Caller   : sort_covari_pairs_in_ascend_order()
// Purpose  : perform the bucket sort in ascending order on the seq->coVari
// Input    : seq
// Return   : none
// Display  : none
//*****************************************************************************
int toss_into_buckets(config* seq, int*** bucket)
{
  int i, j;
  int tag;                                                      // to specify the bucket it's going to be tossed into

  for(i = 0 ; i < seq->numCovari ; i++){
    tag = seq->coVari[i][0] / BAKI_SORT_DIVIDR;                 // for each covariance pair, determine which bucket to toss into
                                                                // based on the value i of (i, j)
    j = 0;
    while(bucket[tag][j][1]){                                   // if some covariance pair(s) already in the bucket, determine the order
      if      (bucket[tag][j][0] < seq->coVari[i][0]) j++;      // new pair is larger,  move new pair to next seat
      else if (bucket[tag][j][0] > seq->coVari[i][0]){          // new pair is smaller, move old pair to next seat
        shift_over(bucket, tag, j);
        bucket[tag][j][0] = seq->coVari[i][0];
        bucket[tag][j][1] = seq->coVari[i][1];
        break;
      } else {      	                                        fprintf(stderr, "Exception happened in 'sort_covari_pairs_in_ascend_order.\n");
        exit(1);
      }  // end if
    }    // end while
    bucket[tag][j][0] = seq->coVari[i][0];                      // fill the new covariance pair to the proper seat in the bucket
    bucket[tag][j][1] = seq->coVari[i][1];
  }  // end inner for

  return 0;
}  // end toss_into_buckets

//*****************************************************************************
// Function : Display Full Buckets
// Caller   : sort_covari_pairs_in_ascend_order()
// Purpose  : display buckets will covariance pairs in them, after the bucket sort process
// Input    : seq
// Return   : none
// Display  : buckets
//*****************************************************************************
int display_full_buckets(config* seq, int*** bucket, int numBaki)
{
  int i, j;

  for(i = 0 ; i < numBaki ; i++){                           
    for(j = 0 ; j < BAKI_SORT_DIVIDR ; j++){ disp(seq,DISP_ALL,"(%3d,%-3d) ", bucket[i][j][0], bucket[i][j][1]);
    }  // end inner for
                                             disp(seq,DISP_ALL,"\n");
  }    // end outer for
                                             disp(seq,DISP_ALL,"\n");
  if(bucket)   return 0;
  else if(seq) return 0;
  else         return 1;

}  // end display_full_buckets

//*****************************************************************************
// Function : Move from Buckets Back to Covariance
// Caller   : sort_covari_pairs_in_ascend_order()
// Purpose  : copy all covariance pairs from buckets to seq->coVari, and then free the space of buckets
// Input    : seq
// Return   : none
// Display  : none
//*****************************************************************************
int move_from_buckets_bak_2_covari(config* seq, int*** bucket, int numBaki)
{
  int    i, j, k = 0;

  for(i = 0 ; i < numBaki ; i++){                             // move all sorted covariance pairs from buckets back to seq->coVari
    for(j = 0 ; j < BAKI_SORT_DIVIDR ; j++){
      if( bucket[i][j][1] && (k < seq->numCovari) ){
        seq->coVari[k][0] = bucket[i][j][0];
        seq->coVari[k][1] = bucket[i][j][1];
        free(bucket[i][j]);
        k++;
      } else 
        free(bucket[i][j]);
    }  // end inner for
    
    free(bucket[i]);
  }    // end outer for

  free(bucket);

  return 0;
}  // end move_from_buckets_bak_2_covari

//*****************************************************************************
// Function : Shift Over
// Caller   : sort_covari_pairs_in_ascend_order()
// Purpose  : insert a certain covariance pair into a space which has been occupied
//          : doing this to maintain the ascending order of the over sorting process
// Input    : bucket - n by 4 by 2 3D matrix containing the covariance data before this particular shift
//          : tag    - index 1 of the bucket to be shifted
//          : i      - index 2 of the bucket to be shifted
//          :        - together bucket[tag][i] locates the cell to be shifted to next cell, which is bucket[tag][i+1]
// Return   : bucket - n by 4 by 2 3D matrix containing the covariance data after this particular shift
// Display  : none
//*****************************************************************************
int shift_over(int*** bucket, int tag, int j)
{
  int current_spot      = j;
  int shift_destination = j + 1;

  if (bucket[tag][shift_destination][1]){                             // if the spot is occupied
    shift_over(bucket, tag, shift_destination);                       // shift the pair in the next cell to next next cell
    bucket[tag][shift_destination][0] = bucket[tag][current_spot][0]; // shift current pair to the next cell
    bucket[tag][shift_destination][1] = bucket[tag][current_spot][1];
  } else {
    bucket[tag][shift_destination][0] = bucket[tag][current_spot][0]; // shift current pair to the next cell
    bucket[tag][shift_destination][1] = bucket[tag][current_spot][1];
  }  // end if

  return 0;
}  // end shift_over

//*****************************************************************************
// Function : Initialize S1 Pairing Database
// Caller   : initialize_sequence_by_getting_argument()
// Purpose  : set up the data related to S1 pairing, and check any violation
// Input    : seq
//          : fp
// Return   : 1: start over, 0: function ends properly
// Display  : inquiry for user input
//*****************************************************************************
int init_s1_pairng_db(config* seq, FILE* fp)
{
  char* constraints = calloc(CONSTRAINT_SIZE, sizeof(char));
  int64_t   i = 0;

  while(fscanf(fp, "%s", constraints) != EOF && (constraints[0] != '_')){
    if(constraints[0] == '('){                             // verify '(' is present, regard it as hint f$
      read_in_s1_pairng_data(seq, fp, constraints, i);     // printf("S1 = %d\n", seq->s1Pairng[i]);
      i++;                                                 // printf("i = %d\n", i);
    }  // end if
  }    // end while                                        // printf("s1 got %d stuff\n", i);

  seq->s1Pairng = realloc(seq->s1Pairng, i * sizeof(int));
  seq->numS1Pairng = i;
  free(constraints);

  return 0;
}  // end init_s1_pairng_db

//*****************************************************************************
// Function : Read In S1 Pairing Data
// Caller   : init_s1_pairng_db()
// Purpose  : read in the S1 pairing data from configuration file
//          : if there's error, pop it up, and exit the process
// Input    : 
// Return   : none
// Display  : none
//*****************************************************************************
int read_in_s1_pairng_data(config* seq, FILE* fp, char* constraints, int i)
{                                                                           disp(seq,DISP_ALL,"Enterng 'read_in_s1_pairng_data'\n");
  if(fscanf(fp, "%s\n", constraints) != '0' &&  atoi(constraints) > 0)      // check for a number
    seq->s1Pairng[i] = atoi(constraints) - 1;
  else if(!strcmp("0", constraints))
    seq->s1Pairng[i] = 0;
  else{                                                                     fprintf(stderr, "Error in S1 pairng configuration file\n");
    printf("%s is not a number\n", constraints);
    exit(1);
  }  // end if 1
  
  if(seq->s1Pairng[i] >= seq->strLen){                                      fprintf(stderr, "Error, S1 pairng %d shouldn't exceed sequence length %d\n", seq->s1Pairng[i], seq->strLen);
    exit(1);
  }  // end if 2

  if((constraints[0] = fgetc(fp)) != ')'){                                  // verify ')' is present
                                                                            printf("Error in S1 pairng configuration file\n");
                                                                            printf("%c is not desinated symbol ')'\n", constraints[0]);
    exit(1);
  }  // end if 4
                                                                            // printf("Constraint = %s\n", constraints);
  return 0;
}  // end read_in_s1_pairng_data

//*****************************************************************************
// Function : Initialize Maximum Pairing Distance
// Caller   : initialize_sequence_by_getting_argument()
// Purpose  : get the value of maximum pairing distance, and check any violation
// Input    : seq
//          : fp
// Return   : 1: start over, 0: function ends properly
// Display  : inquiry for user input
//*****************************************************************************
int init_max_pairng_dist(config* seq, FILE* fp)
{
  char* constraints = calloc(CONSTRAINT_SIZE, sizeof(char));

  while(fscanf(fp, "%s\n", constraints) != EOF && (constraints[0] != '_')){
    if(constraints[0] == '('){                           // verify '(' is present, regard it as hint f$
      read_in_max_pairng_dist_data(seq, fp, constraints);
      free(constraints);
      return 0;
    }  // end if
  }    // end while
                                                         // printf("maxPairDist = %d\n", seq->maxPairngDist);
  free(constraints);

  return 0;
}  // end init_max_pairng_dist

//*****************************************************************************
// Function : Read In Maximum Pairing Distance Data
// Caller   : init_max_pairng_dist()
// Purpose  : read in the max pairing distance value from configuration file
//          : if there's error, pop it up, and exit the process
// Input    : 
// Return   : none
// Display  : none
//*****************************************************************************
int read_in_max_pairng_dist_data(config* seq, FILE* fp, char* constraints)
{
  if(fscanf(fp, "%s\n", constraints) != '0' &&  atoi(constraints) > 0)      // check for a number
    seq->maxPairngDist = atoi(constraints);
  else if(!strcmp("0", constraints))
    seq->maxPairngDist = 0;
  else{                                                                     printf("Error in max pairng distanceconfiguration file\n");
                                                                            printf("%s is not a number\n", constraints);
    exit(1);
  }  // end if 1

  if((constraints[0] = fgetc(fp)) != ')'){                                  // verify ')' is present
                                                                            printf("Error in max pairng dist configuration file\n");
                                                                            printf("%c is not desinated symbol ')'\n", constraints[0]);
    exit(1);
  }  // end if 4
  return 0;
}  // end read_in_max_pairng_dist_data

//*****************************************************************************
// Function : Initialize Minimum Pairing Distance
// Caller   : initialize_sequence_by_getting_argument()
// Purpose  : get the value of minimum pairing distance, and check any violation
// Input    : seq
//          : fp
// Return   : 1: start over, 0: function ends properly
// Display  : inquiry for user input
//*****************************************************************************
int init_min_pairng_dist(config* seq, FILE* fp)
{
  char* constraints = calloc(CONSTRAINT_SIZE, sizeof(char));

  while(fscanf(fp, "%s\n", constraints) != EOF && (constraints[0] != '_')){
    if(constraints[0] == '('){                           // verify '(' is present, regard it as hint f$
      read_in_min_pairng_dist_data(seq, fp, constraints);
      free(constraints);
      return 0;
    }  // end if
  }    // end while
                                                         // printf("minPairDist = %d\n", seq->minPairngDist);
  free(constraints);

  return 0;
}  // end init_min_pairng_dist

//*****************************************************************************
// Function : Read In Minimum Pairing Distance Data
// Caller   : init_min_pairng_dist()
// Purpose  : read in the min pairing distance value from configuration file
//          : if there's error, pop it up, and exit the process
// Input    : 
// Return   : none
// Display  : none
//*****************************************************************************
int read_in_min_pairng_dist_data(config* seq, FILE* fp, char* constraints)
{
  if(fscanf(fp, "%s\n", constraints) != '0' &&  atoi(constraints) > 0)      // check for a number
    seq->minPairngDist = atoi(constraints);
  else if(!strcmp("0", constraints))
    seq->minPairngDist = 0;
  else{                                                                     printf("Error in min pairng distanceconfiguration file\n");
                                                                            printf("%s is not a number\n", constraints);
    exit(1);
  }  // end if 1

  if((constraints[0] = fgetc(fp)) != ')'){                                  // verify ')' is present
                                                                            printf("Error in min pairng dist configuration file\n");
                                                                            printf("%c is not desinated symbol ')'\n", constraints[0]);
    exit(1);
  }  // end if 4
  return 0;
}  // end read_in_min_pairng_dist_data

//*****************************************************************************
// Function : Initialize V1 Pairing Database
// Caller   : initialize_sequence_by_getting_argument()
// Purpose  : set up the data related to V1 pairing, and check any violation
// Input    : seq
//          : fp
// Return   : 1: start over, 0: function ends properly
// Display  : inquiry for user input
//*****************************************************************************
int init_v1_pairng_db(config* seq, FILE* fp)
{
  char* constraints = calloc(CONSTRAINT_SIZE, sizeof(char));
  int   i = 0;

  while(fscanf(fp, "%s", constraints) != EOF && (constraints[0] != '_')){
    if(constraints[0] == '('){                             // verify '(' is present, regard it as hint f$
      read_in_v1_pairng_data(seq, fp, constraints, i);     // printf("V1 = %d\n", seq->v1Pairng[i]);
      i++;
    }  // end if
  }    // end while
  
  seq->v1Pairng = realloc(seq->v1Pairng, i * sizeof(int));
  seq->numV1Pairng = i;
  chk_conflict_against_v1_pairng(seq);
  free(constraints);

  return 0;
}  // end init_v1_pairng_db

//*****************************************************************************
// Function : Read In V1 Pairing Data
// Caller   : init_v1_pairng_db()
// Purpose  : read in the V1 pairing data from configuration file
//          : if there's error, pop it up, and exit the process
// Input    : 
// Return   : none
// Display  : none
//*****************************************************************************
int read_in_v1_pairng_data(config* seq, FILE* fp, char* constraints, int i)
{                                                                           disp(seq,DISP_ALL,"Enterng 'read_in_v1_pairng_data'\n");
  if(fscanf(fp, "%s\n", constraints) != '0' &&  atoi(constraints) > 0)      // check for a number
    seq->v1Pairng[i] = atoi(constraints)-1;
  else if(!strcmp("0", constraints))
    seq->v1Pairng[i] = 0;
  else{                                                                     printf("Error in V1 pairng configuration file\n");
    printf("%s is not a number\n", constraints);
    exit(1);
  }  // end if 1

  if((constraints[0] = fgetc(fp)) != ')'){                                  // verify ')' is present
                                                                            printf("Error in V1 pairng configuration file\n");
                                                                            printf("%c is not desinated symbol ')'\n", constraints[0]);
    exit(1);
  }  // end if 4
                                                                            // printf("Constraint = %s\n", constraints);
  return 0;
}  // end read_in_v1_pairng_data

//*****************************************************************************
// Function : Check Conflict against V1 Pairings
// Caller   : init_v1_pairng_db()
// Purpose  : check if there is any contradiction between v1 pairings and other constraints
// Input    : seq
// Return   : none
// Display  : confliction message, then exit, if any
//*****************************************************************************
int chk_conflict_against_v1_pairng(config* seq)
{                                                                  disp(seq,DISP_ALL,"Enterng 'chk_conflict_against_v1_pairng'\n");
  int8_t v1Cursr, s1Cursr;

  for(v1Cursr = 0 ; v1Cursr < seq->numV1Pairng ; v1Cursr++){     // printf("v1 has %d, s1 has %d\n", seq->numV1Pairng, seq->numS1Pairng);
    for(s1Cursr = 0 ; s1Cursr < seq->numS1Pairng ; s1Cursr++){   // printf("(v1,s1) = (%d,%d)\n", seq->v1Pairng[v1Cursr], seq->s1Pairng[s1Cursr]);
      if(seq->v1Pairng[v1Cursr] == seq->s1Pairng[s1Cursr]){
        fprintf(stderr, "Error, conflict between V1 pairng %d & S1 pairng %d\n", seq->v1Pairng[v1Cursr], seq->s1Pairng[s1Cursr]);
        exit(1);
      }  // end if 1
    }    // end for 1
    for(s1Cursr = 0 ; s1Cursr < seq->numCovari ; s1Cursr++)
      if(seq->v1Pairng[v1Cursr] == seq->coVari[s1Cursr][0] || seq->v1Pairng[v1Cursr] == seq->coVari[s1Cursr][1]) {
        fprintf(stderr, "Error, duplication between V1 pairng %d & covariance (%d,%d)\n", seq->v1Pairng[v1Cursr], seq->coVari[s1Cursr][0], seq->coVari[s1Cursr][1]);
        exit(1);
      }  // end if 2
  }      // end outer for

  return 0;
}  // end chk_conflict_against_v1_pairng

//*****************************************************************************
// Function : Display Constraints
// Caller   : initialize_sequence_by_gettng_constraint()
// Purpose  : display the contents of all constraints, after verifying their legitimacy
//          : if there's error, pop it up, and exit the process
// Input    : 
// Return   : none
// Display  : contents of constraints
//*****************************************************************************
int display_constraints(config* seq)
{
  int16_t i;
                                          fprintf(seq->dispFile,"Covari (%d) = ", seq->numCovari);
  for(i = 0 ; i < seq->numCovari ; i++)   fprintf(seq->dispFile,"(%d,%d), ", seq->coVari[i][0], seq->coVari[i][1]);
                                          fprintf(seq->dispFile,"\nV1 Prng(%d) = ", seq->numV1Pairng);
  for(i = 0 ; i < seq->numV1Pairng ; i++) fprintf(seq->dispFile,"%d, ", seq->v1Pairng[i]);
                                          fprintf(seq->dispFile,"\nS1 Prng(%d) = ", seq->numS1Pairng);
  for(i = 0 ; i < seq->numS1Pairng ; i++) fprintf(seq->dispFile,"%d, ", seq->s1Pairng[i]);
                                          fprintf(seq->dispFile,"\nMax/Min pairng Dist = (%d,%d)\n", seq->maxPairngDist, seq->minPairngDist);
  return 0;
}  // end display_constraints

