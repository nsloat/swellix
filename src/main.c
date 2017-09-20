/*****************************************************************************************************************************************************
 *****************************************************************************************************************************************************
 *****************************************************************************************************************************************************
      Program : Swellix-main.c
      Purpose : Ryan's RNA Sequence Folding Algorithm (modified Crumple)
      Input   : manual key-in or input file
      Return  : none
      Display : pair displays in multiple-bracket/parenthesis or alphabet layouts
      CAUTION : integer capability of 16 bit signed carries positive integers of 2^16/2 -1 = 32767
              :             "         32 bit          "         "        "       2^32/2 -1 = 2147483647
              :             "         64 bit          "         "        "       2^64/2 -1 = 9223372036854775807
              :
       NOTE 1 :                        LENGTH   COMPONENT LIST          STRUCTURE                     EXE TIME (user)
              : seq                             noSZG         SZG       noSZG        SZG              noSZG                SZG
              :---------------------------------------------------------------------------------------------------------------------------------------
              : 14mer                  14       21            35        118          182              0.01sec              0.01sec
              : 2x14mer                28       96            156       124925       252307           0.14sec              0.24sec
              : alfalfa_AMV4_843-881   39       220           318       50781503     114864454        18sec                38sec
              : 3x14mer                42       227           374       197316084    548588270        1min 10sec           3min
              : gND7-506_12-55_seq     44       354                     5370612992                   57min 30sec
              : Mic_A_50               50       420           678       74616503570  217425780862     9hr 43min 57.64sec   27hr 13min 49.13sec
              : Mic_A_72               72
              : tRNA-R1660             74       835           1255
              : tRMA-R1660 (-l 4)      74       48            76        2622         9946             0.02sec              0.05sec
              : tRMA-R1660 (-l 3)      74                     201                    955889                                0.4sec
              : tRMA-R1660 (-l 2)      74       337           547       592166437    3219896561       4min 30sec           17min 30sec
              : ecoli-5s               120      2634          4126
              : ecoli-5s   (-l 2)      120      950           1492      
              : stmv.seq               1058     216207        358812
              :
       NOTE 2 : helix find syntax ./find.py < /home/liu5227/code/seq/stmv.seq -s 9 | wc    
 *****************************************************************************************************************************************************
 *****************************************************************************************************************************************************
 *****************************************************************************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "main.h"
#include "init_general.h"
#include "init_constraint.h"
#include "paren_lookup_table.h"
#include "component_list.h"
#include "interval_lookup_table.h"
#include "bundle_list.h"
#include "jump_tree.h"
#include "unit_tests.h"
#include "close_up.h"
#include "statistics.h"
#ifdef _MPI
#include <mpi.h>
#include "mpi_jump_tree.h"
#endif

int rank, wsize;

//*****************************************************************************
// Function : Main
// Caller   : none
// Purpose  : main
// Input    : argv - the array with all the command contents
//          : argc - the total number of of cells with data in argv[]
// Return   : none
// Note     : for comand line options, see 'prtUsage()'
//          : 'constraint' refers to the experimental constraints
//          : MPI_Wtime() and printf here are for load balancing evaluation
//          : two ways to get input - 
//          : 1. initialize_sequence_by_getting_argument() from command line
//          : 2. initialize_sequence_by_getting_constraints() from conf file
//*****************************************************************************

int main(int argc, char** argv) {

  int start, end, span, rem;
#ifdef _MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&wsize);
#else
  rank = 0; wsize = 1;
#endif
 
  config* seq = initialize_sequence_by_getting_argument(argc, argv);
  initialize_sequence_by_settng_memory(seq);
  initialize_sequence_by_gettng_constraint(seq);

  span = (seq->strLen)/(wsize);
  rem = seq->strLen % wsize;
  start = rank*span;
  end = (rank == wsize-1) ? start+span+rem : start+span;

  if(rank == 0) {
    printf("%s\n",seq->ltr);
  }
  init_vrna(seq->ltr);
  float energy = get_mfe_structure(seq->ltr, seq->mfe);
  const int INF = 10000000;
  seq->minenergy = INF;

  global* crik = initialize_crik(seq);

  if(seq->algoMode > MODE_NONE) {
    if(rank==0)disp(seq,DISP_LV1,"started paren lookup table...\n");
    make_parenthesis_look_up_table(seq);
    if(rank==0)disp(seq,DISP_LV1,"completed paren lookup table\nstarted component list\n");
    make_component_list(seq, crik);
    if(rank==0)disp(seq,DISP_LV1,"completed component list\n");
    if(seq->algoMode > MODE_CMPNT) {
      if(rank==0)disp(seq,DISP_LV1,"started interval lookup table...\n");
      make_interval_look_up_table(seq, crik);
      if(rank==0)disp(seq,DISP_LV1,"completed interval lookup table\n");
      if (seq->bundle) {
        if(rank==0)disp(seq,DISP_LV1,"started bundle list...\n");
    	make_bundle_list(seq, crik);
        if(rank==0)disp(seq,DISP_LV1,"completed bundle list\n");
      }
      if(seq->algoMode > MODE_INTAB) {
        if(rank==0)disp(seq,DISP_LV1,"started jump tree...\n");
          make_jump_tree(seq, crik, start, end);
        if(rank==0)disp(seq,DISP_LV1,"completed\n");
      }
    }    // end if 2
  }      // end if 1

#ifdef _MPI
  uint64_t sumStructures = 0, sumUBStructures = 0, maxDist = 0, motifCount = 0;
  float minEng = 0.0;
  MPI_Reduce(&crik->numStru, &sumStructures, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&crik->numUnbundledStru, &sumUBStructures, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&seq->minenergy, &minEng, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&seq->maxdist, &maxDist, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&seq->motifCount, &motifCount, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  if(rank==0) {
    seq->minenergy = minEng;
    seq->maxdist = maxDist;
    seq->motifCount = motifCount;
    crik->numUnbundledStru = sumUBStructures;
    crik->numStru = sumStructures;
  }
#endif  

  if(rank==0) {
    print_results(seq, crik);
  }
  close_up(seq, crik);

#ifdef _MPI
  MPI_Finalize();
#endif

  return 0;
}  // end main

//*****************************************************************************
// Function : Initialize Sequence by Getting Arguments
// Caller   : main()
// Purpose  : intepret command line options and tune the processing path based on it the structure seq is generated here
// Input    : argv - the array with all the command contents
//          : argc - the total number of of cells with data in argv[]
// Return   : seq
// Display  : error message, if any
//*****************************************************************************
config* initialize_sequence_by_getting_argument(int argc, char** argv) {

  config* seq = malloc(sizeof(config));
  int8_t arg;                                                                                                 // the index counter refering to the address of the argv[]
                                                                                                              // used to constantly check with argc
  initialize_commandline_argument_with_default_value(seq);

  for (arg = 1; arg < argc; arg++) {
    if        (!strcmp("-#", argv[arg]) || !strcmp("-hlx#m", argv[arg]) || !strcmp("--minimum_number_of_helices", argv[arg]))    
      if(++arg == argc || !atoi(argv[arg])) cmd_line_err(argv, arg);
      else          	                    seq->minNumHlix = atoi(argv[arg]);                                // min helices to form a complete structure
    else if   (!strcmp("-b", argv[arg])     || !strcmp("-bundle", argv[arg]))                                 // activate the bundling mechanism -- similar structures reduction referencing (SSRR)
                                            seq->bundle = TRUE;
    else if   (!strcmp("-u", argv[arg])     || !strcmp("-unbund", argv[arg]))				      // activate the unbundling process to verify that the bundling mechanism 
													      // accurately computes all possible structures. Only use for debugging
                                            seq->unbundle = TRUE;					      // since this defeats the purpose of bundling otherwise.
    else if   (!strcmp("-s", argv[arg])     || !strcmp("-stat", argv[arg]))
                                            set_statMode(seq, arg++, argc, argv);
    else if   (!strcmp("-M", argv[arg])     || !strcmp("-motif", argv[arg])) {
                                            set_motif(seq, arg++, argc, argv);
                                            seq->motif = TRUE;
    }
    else if   (!strcmp("-bhdm", argv[arg])  || !strcmp("--minimum_between_helix_distance", argv[arg]))    
      if(++arg == argc || !atoi(argv[arg])) cmd_line_err(argv, arg);
      else          	                    seq->minBtwnHlixDist = atoi(argv[arg]);                           // minimum between helix distance. used primarily in STMV
    else if   (!strcmp("-d", argv[arg])     || !strcmp("-dsply", argv[arg]) || !strcmp("--degree_of_detailed_display", argv[arg])){
      set_degree_of_detailed_disp(seq, arg++, argc, argv);
    } else if (!strcmp("-h", argv[arg])     || !strcmp("-help",  argv[arg]) )
      print_usage();
    else if (!strcmp("-i",   argv[arg]))                                                                      // take in the sequence data
      seq->srcFile = fopen(argv[++arg], "r");
    else if   (!strcmp("-k", argv[arg])     || !strcmp("-constraint", argv[arg]))                             // introduce the constraints: V1, S1, Covariance, etc
      seq->constraintActive = 1;
    else if   (!strcmp("-l", argv[arg])     || !strcmp("-hlxlm", argv[arg]) || !strcmp("--minimum_length_of_helix", argv[arg]))    
      if(++arg == argc || !atoi(argv[arg])) cmd_line_err(argv, arg);
      else          	                    seq->minLenOfHlix = atoi(argv[arg]);                              // min nucleotides to form a complete helix
    else if   (!strcmp("-m", argv[arg])     || !strcmp("-mode", argv[arg]) || !strcmp("--mode_of_algorithm", argv[arg])){
      set_mode_of_algo(seq, arg++, argc, argv);
    } else if (!strcmp("-mm", argv[arg])    || !strcmp("--maximum_number_of_mismatches", argv[arg]))    
      if(++arg == argc || !atoi(argv[arg])) cmd_line_err(argv, arg);
      else          	                    seq->maxNumMismatch = atoi(argv[arg]);                            // maximum count of mismatches per helix 
    else if   (!strcmp("-noGU", argv[arg])  || !strcmp("--no_G_U_allowed", argv[arg]))                        // disallows G-U pairing
      seq->noGU = TRUE;
    else if   (!strcmp("-noSZG", argv[arg]) || !strcmp("--no_sizing_allowed", argv[arg]))                     // disallows swelling helix unless necessary
      seq->noSZG = TRUE;
    else if   (!strcmp("-o",      argv[arg++])) {                                                               // designate the printing destination
      char *ofile = calloc(1, strlen(argv[arg]) + 15);
      memcpy(ofile, argv[arg], strlen(argv[arg]));
      seq->dispFile = fopen(ofile, "w");
      free(ofile);
    } else if (!strcmp("-p", argv[arg])     || !strcmp("-hp#m", argv[arg]) || !strcmp("--minimum_number_of_hairpins", argv[arg]))    
      if(++arg == argc || !atoi(argv[arg])) cmd_line_err(argv, arg);
      else          	                    seq->minNumHP = atoi(argv[arg]);                                  // min hairpins per structure
    else {
      printf("\n Unknown option %s", argv[arg]);
      print_usage();
    }  // end if
  }    // end for

  if(seq->maxNumMismatch >= seq->minLenOfHlix) {
    fprintf(stderr, "Mismatch too huge. Minimum helix length '-l' "
                    "has to be larger than maximum mismatch counts '-mm'\n");
    exit(1);
  }  // end if

  return seq;
}  // end initialize_sequence_by_getting_argument

//*****************************************************************************
// Function : Initialize Sequence by Setting Memories
// Caller   : main()
// Purpose  : evaluate and initialize seq substructures memory allocations
// Input    : seq
// Return   : none
// Display  : none
//*****************************************************************************
int initialize_sequence_by_settng_memory(config* seq) {

  if(!seq->dispFile) seq->dispFile = stdout;                         disp(seq,DISP_ALL,"################################################################################################  Enterng 'init_seq'\n");
  init_seq_ltr(seq);                                                 // initialize seq->ltr (get the sequence data, and store in seq->ltr)
  seq->strLen    = strlen(seq->ltr);                                 // count the number of characters in this string
  seq->dotNParen = calloc((seq->strLen + 1), sizeof(char));          // assign the size of seq->dotNParen

  seq->insIntrvlMinSize = seq->minPairngDist + seq->minLenOfHlix * 2;
  seq->behIntrvlMinSize = seq->minPairngDist + seq->minLenOfHlix * 2;
  seq->maxdist = 0;
  seq->minenergy = 0;
   
  return 0;
}  // end initialize_sequence_by_settng_memory

//*****************************************************************************
// Function : Initialize Sequence by Getting Constraints
// Caller   : main()
// Purpose  : intepret constraints from sequence file, to set up constraint database
// Input    : seq
// Return   : none
// Display  : none
//*****************************************************************************
int initialize_sequence_by_gettng_constraint(config* seq) {

  disp(seq,DISP_ALL,"Entering 'init_seq_constraint'\n");
  if(!seq->constraintActive) return 0;

  char*  constraints = calloc(1000, sizeof(char));
  char*  keywdCatch  = calloc(1000, sizeof(char));
  FILE*  fp          = seq->srcFile;
  
  seq->s1Pairng = calloc(256, sizeof(int));                // set aside the room for S1 pairing data
  seq->v1Pairng = calloc(256, sizeof(int));                // set aside the room for V1 pairing data
  seq->coVari = calloc(256, sizeof(seq->coVari));                 // set aside the room for covariance data
  
  while(fscanf(fp, "%s", keywdCatch) != EOF){             // if conf file is not provided, constraints part is bypasse3d

    if      (!strcmp("[MAX", keywdCatch)          &&      // max pairing distance
            (fscanf(fp, "%s", keywdCatch) != EOF) &&      // ||
	    !strcmp("PAIRING", keywdCatch)        &&      // ||
	    (fscanf(fp, "%s", keywdCatch) != EOF) &&      // ||
            !strcmp("DISTANCE]", keywdCatch)         )    // ||
      init_max_pairng_dist(seq, fp);                      // \/
    else if (!strcmp("[MIN", keywdCatch)          &&      // min pairing distance
	    (fscanf(fp, "%s", keywdCatch) != EOF) &&      // ||
	    !strcmp("PAIRING", keywdCatch)        &&      // ||
	    (fscanf(fp, "%s", keywdCatch) != EOF) &&      // ||
	    !strcmp("DISTANCE]", keywdCatch)         )    // ||
      init_min_pairng_dist(seq, fp);                      // \/
    else if (!strcmp("[S1", keywdCatch)           &&      // S1 pairing
	    (fscanf(fp, "%s", keywdCatch) != EOF) &&      // ||
	    !strcmp("PAIRING]", keywdCatch)          )    // ||
      init_s1_pairng_db(seq, fp);                         // \/
    else if (!strcmp("[COVARIANCE]", keywdCatch))         // covariance
      init_covari_db(seq, fp);                            // \/
    else if (!strcmp("[V1", keywdCatch)           &&      // V1 pairing
	    (fscanf(fp, "%s", keywdCatch) != EOF) &&      // ||
	    !strcmp("PAIRING]", keywdCatch)          )    // ||
      init_v1_pairng_db(seq, fp);                         // \/
  }  // end while

  seq->maxPairngDist = (seq->maxPairngDist > seq->strLen) ? seq->strLen : seq->maxPairngDist;
	    
  dispConstraints(seq);
  free(keywdCatch);
  free(constraints);

  return 0;
}  // end initialize_sequence_by_gettng_constraint

//*****************************************************************************
// Function : Initialize Crick
// Caller   : main()
// Purpose  : initialize 'global* crik'
// Input    : seq
// Return   : crik
// Display  : none
//*****************************************************************************
global* initialize_crik(config* seq) {

  disp(seq,DISP_ALL,"################################################################################################  Enterng 'initialize_crik'\n");

  global* crik = malloc(sizeof(global));                                                 // Crick the Book-keeper, our handy right-hand man #1
  crik->linkedmms         = 0;
  crik->rstoCounter       = 0;
  crik->rstoErrCounter    = 0;
  crik->skippedStru	      = 0;
  crik->lvlOfRecur        = 0;
  crik->numCmpnt          = 0;
  crik->numHP             = 1;
  crik->numHlix           = 1;
  crik->numStru           = 0;
  crik->numBundles	  = 0;
  crik->numUnbundledStru  = 0;
  crik->specialRstoFlag   = 0;
  crik->intrvlCntr        = 0;                                                                 // it's like batch number, used to distinguish one batch from another
  crik->intrvlCntr        = 0;
  crik->opnBrsStop        = seq->strLen - ((seq->minLenOfHlix - 1) << 1) - seq->minPairngDist; // farthest location the opnBrsOutIndx may reach
  crik->numCmpntTyp       = seq->strLen;
  crik->numCmpntTypOcupid = 0;
  crik->test1ErrTotal     = 0;
  crik->test2ErrTotal     = 0;
  crik->test3ErrTotal     = 0;
  crik->test4ErrTotal     = 0;
  crik->interval          = NULL;
  crik->hlixInStru        = NULL;

  if(seq->constraintActive)
    crik->struMustPairFlag = calloc((seq->numCovari + seq->numV1Pairng), sizeof(int16_t));     // used to track the covariance pairs the current structure contains  
  crik->mustPairLength = seq->numCovari + seq->numV1Pairng;
  crik->numStru = 0;

  return crik;
}  // end initialize_crik

//*****************************************************************************
// Function : Make Parenthesis Look-up Table
// Caller   : main()
// Purpose  : tabulate the base pair forming relationships to speed up the parentheses matching up (A-U, C-G, G-U) process
// Input    : seq
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
int make_parenthesis_look_up_table(config* seq) {
  
  disp(seq,DISP_ALL,"################################");
  disp(seq,DISP_ALL,"################################");
  disp(seq,DISP_ALL,"################################");
  disp(seq,DISP_ALL,"  Entering 'make_paren_luk_up'\n");

  int16_t rowI;

  seq->parenLukUpTable = malloc(sizeof(*seq->parenLukUpTable) * seq->strLen);   // set up the backbone for the rows of the table
                                                                                disp(seq,DISP_LV3,"Parenthesis\n");
                                                                                disp(seq,DISP_LV3,"Look-up Table: \n");
										disp(seq,DISP_LV3,"                 %s\n", seq->ltr);
  for(rowI = 0 ; rowI < seq->strLen ; rowI++){                                                  
    fill_in_parenthesis_look_up_table_row_i(seq, rowI);                                // fill the table, row by row
  }  // end for
                                                                                disp(seq,DISP_LV3,"\n");
//                                                                                fprintf(seq->dispFile,"\nString Length = %d\n\n", seq->strLen);
  return 0;
}  // end make_parenthesis_look_up_table

//*****************************************************************************
// Function : Make Component List
// Caller   : main()
// Purpose  : gather and pile up the helices as components to speed up structure matching by making component (look-up) list
// Input    : seq
//          : crik
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
int make_component_list(config* seq, global* crik) {
  
  disp(seq,DISP_ALL,"################################################################################################  Entering 'make_component_list'\n");
 
  crik->cmpntList = calloc(seq->strLen, sizeof(edgeList));       // the backbone of 'edge list' to store all the components arbitrarily assgin some number of rooms to accomodate unknown number of components trim out the extra empty unused rooms later after all components are found
  crik->cmpntListOcupidTyp = calloc(seq->strLen, sizeof(int));

  recur_stage_1_brace_locatng(seq, crik);

  book_out(seq, crik);                                          // final summary and stats about component list

  return 0;
}  // end make_component_list

//*****************************************************************************
// Function : Make Interval Look-up Table
// Caller   : main()
// Purpose  : tabulate the availability of components based on the needs of interval ends to facilitate the interval fitting process
// Input    : 
// Return   : 
// Display  : 
// Note     : 1st dimension - rowI : row    i, lower bound of interval with fix size of string length (seq->strLen)
//          : 2nd dimension - colJ : column j, upper              "                        "
//          : 3rd dimension - hgtK : height k, number of cmpnt types which may fit into the interval
//*****************************************************************************
int make_interval_look_up_table(config* seq, global* crik) {

  disp(seq,DISP_ALL,"################################################################################################  Entering 'make_intrvl_luk_up'\n");

  int16_t rowI;

  seq->intrvlLukUpTable = malloc(sizeof(seq->intrvlLukUpTable) * seq->strLen);     // set up the backbone for the rows of the table
                                                                                   disp(seq,DISP_LV3,"\n\nInterval Look-up Table :\n\n");
                                                                                   disp(seq,DISP_ALL,"  Interval Boundary    Cmpnt Type\n\n");
  if(DISP) format_the_intrvl_luk_up_table_prt_out(seq);                            // if display isn't desired for speed-up, this for loop may be skipped
    
  for(rowI = 0 ; rowI < seq->strLen ; rowI++)
    fill_in_intrvl_luk_up_table_row_i(seq, crik, rowI);                            // fill the table, row by row
                                                                                   disp(seq,DISP_LV4,"\n");
  return 0;
}  // end make_interval_look_up_table

//*****************************************************************************
// Function : Make Bundle List
// Caller   : main()
// Purpose  : make mini jump tree to ease the burden of make_jump_tree() when there's no worry about mult-branch loop inside a interval of interest
// Input    : seq
//          : crik
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
int make_bundle_list(config* seq, global* crik) {
  disp(seq,DISP_ALL,"################################################################################################  Entering 'make_bundle_st'\n");

  initialize_eden(seq, crik);
  run_sliding_windows(seq, crik);
  dispEden(seq, crik);
  return 0;
}  // end make_bundle_list

void make_bundled_cmpntList(config* seq, global* crik) {
  int j;
  int k;

  // remove components that are less than max bundle len
  for (j = 0; j < seq->strLen; j++) {
    knob* cmpntCursr = crik->cmpntList[j].knob;
    knob* prev = crik->cmpntList[j].knob;
    k = 0;
    while (cmpntCursr) {
      int len = cmpntCursr->closeBrsOutIndx - cmpntCursr->opnBrsOutIndx + 1;
      int maxBundleLen = (seq->minPairngDist * 2) + (seq->minLenOfHlix * 6) - 1;
      if (len <= maxBundleLen) {
        if (k == 0) { // first element needs to be deleted
          knob* toDelete = cmpntCursr;
          crik->cmpntList[j].knob = cmpntCursr->cmpntListNext;
          toDelete->cmpntListNext = NULL;
          cmpntCursr = crik->cmpntList[j].knob;
          free(toDelete);
          crik->numCmpnt--;
          crik->cmpntList[j].cnt--;
          continue;
        } else { // this current element needs to be deleted, in the middle of the list
          knob* toDelete = cmpntCursr;
          prev->cmpntListNext = cmpntCursr->cmpntListNext;
          toDelete->cmpntListNext = NULL;
          cmpntCursr = prev->cmpntListNext;
          free(toDelete);
          crik->numCmpnt--;
          crik->cmpntList[j].cnt--;
          continue;
        }
      }

      prev = cmpntCursr;
      cmpntCursr = cmpntCursr->cmpntListNext;
      k++;
    }
  }

  for (j = 0; j < seq->strLen; j++) {
    for (k = seq->strLen-1; k > j; k--) {
      if (crik->eden[j][k].dumiNode) {
        crik->eden[j][k].dumiNode->cmpntListNext = crik->cmpntList[j].knob;
        crik->cmpntList[j].knob = crik->eden[j][k].dumiNode;
        crik->cmpntList[j].cnt++;
        crik->numCmpnt++;
      }
    }
  }
  display_components(seq, crik, 1);

  // reset the component list occupied type list
  crik->numCmpntTypOcupid = 0;
  for (j = 0; j < seq->strLen; j++) {
    if (crik->cmpntList[j].knob) {
      crik->cmpntListOcupidTyp[crik->numCmpntTypOcupid] = j;
      crik->numCmpntTypOcupid++;
    }
  }
}

//*****************************************************************************
// Function : Make Jump Tree
// Caller   : main()
// Purpose  : make use of bundle list, component list and interval look-up table to form the structures, consume over 90% (or more) of total execution time
// Input    : seq
//          : crik
// Return   : none
// Display  : structures, if display priority is lv 2 or lower
// Note     : 1st dimension - rowI : row    i, lower bound of interval with fix size of string length (seq->strLen)
//          : 2nd dimension - colJ : column j, upper              "                        "
//          : 3rd dimension - hgtK : height k, number of cmpnt types which may fit into the interval
//*****************************************************************************
int make_jump_tree(config* seq, global* crik, int start, int end) {
    disp(seq,DISP_ALL,"################################################################################################  Entering 'make_jump_tree'\n");
    disp(seq,DISP_ALL,"Recursion level starts at zero from here\n");

  local* todd = init_todd(crik);
  int16_t    i;
  int16_t    hlixBranchngIndx1 = 0; // all branches of hlix are rooted (numbered) from make_jump_tree()

  crik->hlixInStru = NULL;
  crik->intrvlCntr = 0; // this reset is necessary because it's just used by make_bundle_list()

  // Combine bundles into component list if necessary
  if(seq->bundle)
    make_bundled_cmpntList(seq, crik);

  int keepgoing;

#ifdef _MPI

  knob* cmpnts[crik->numCmpnt];
  int counter = 0;
  for(i = 0; i < crik->numCmpntTypOcupid; i++) {
    cmpnts[counter] = crik->cmpntList[crik->cmpntListOcupidTyp[i]].knob;
    cmpnts[counter]->newCLindex = counter;
    counter++;
    while(cmpnts[counter-1]->cmpntListNext != NULL) {
      cmpnts[counter] = cmpnts[counter-1]->cmpntListNext;
      cmpnts[counter]->newCLindex = counter;
      counter++;
    }
  }

  crik->numCmpnt = counter;

  int cmpntTracker[counter];
  for(i = 0; i < counter; i++) cmpntTracker[i] = cmpnts[i]->newCLindex;
  
  crik->mpiCList = cmpnts;
  make_jump_tree_parallel(seq, crik, todd, cmpnts, cmpntTracker);

#else

  // scan thru the whole cmpnt list
  for(i = 0; i < crik->numCmpntTypOcupid; i++) {
    // grab one component on the component type indicated by crik->cmpntListOcupidTyp[i]
    todd->cmpntLLCursr = crik->cmpntList[crik->cmpntListOcupidTyp[i]].knob;
    keepgoing = 1;

    // check if it's worth it to go thru the rest of the components
    if(!time_to_quit(seq, todd)) {
      while(todd->cmpntLLCursr) {
        disp(seq,DISP_ALL,"at 'make_jump_tree' while loop, lvl 0 hlix selected\n");
        disp(seq,DISP_ALL,"{%2d-%-2d{    }%2d-%-2d}\n", todd->cmpntLLCursr->opnBrsOutIndx, 
                                                        todd->cmpntLLCursr->opnBrsInnIndx, 
                                                        todd->cmpntLLCursr->closeBrsInnIndx, 
                                                        todd->cmpntLLCursr->closeBrsOutIndx);
        hlixBranchngIndx1++;
        clear_old_hlix_link(crik);

        // from here, a complex series of recursion starts here    <----- CORE OF SWELLIX
        take_cmpnt_list_normal_path(seq, crik, todd, hlixBranchngIndx1);

        todd->RSTO = NULL;

        todd->intrvlIns->opnBrsInnIndx = -1;
        todd->intrvlIns->closeBrsInnIndx = -1;
        todd->intrvlIns->lvlOfRecur = -1;
        todd->intrvlIns->intrvlInsFormdFlag = 0;
        todd->intrvlIns->jumpTreeNext = NULL;
        todd->intrvlIns->intrvlCntr = 0;

        todd->intrvlBeh->opnBrsInnIndx = -1;
        todd->intrvlBeh->closeBrsInnIndx = -1;
        todd->intrvlBeh->lvlOfRecur = -1;
        todd->intrvlBeh->intrvlInsFormdFlag = 0;
        todd->intrvlBeh->jumpTreeNext = NULL;
        todd->intrvlBeh->intrvlCntr = 0;

        todd->cmpntLLCursr = todd->cmpntLLCursr->cmpntListNext;
      }  // end while
    }
  }      // end for

#endif

  exit_curr_recur(seq,crik, todd);

  return 0;
}  // end make_jump_tree

//*****************************************************************************
// Function : Print Statistics
// Caller   : main()
// Purpose  : print final summary
// Input    : seq
//          : crik
// Return   : none
// Display  : number of structures and other necessary statistics
//*****************************************************************************
void print_results(config* seq, global* crik)
{                                                             disp(seq,DISP_ALL,"Entering 'print_results'\n");
  if(seq->algoMode > MODE_INTAB){
    if (DISP)  fprintf(seq->dispFile, "\n\nComponents: %ld\n", crik->numCmpnt);

    if (seq->bundle) {
      if (DISP) fprintf(seq->dispFile, "Bundles: %ld\nBundled structures: %ld\nTotal structures: %ld\n",crik->numBundles, crik->numStru, crik->numUnbundledStru);
      else fprintf(seq->dispFile, "Bundled structures: %ld\nTotal structures: %ld\n", crik->numStru, crik->numUnbundledStru);
    } else {
      fprintf(seq->dispFile, "Total structures: %ld\n", crik->numStru);
    }

    print_stats(seq);    

    if (DISP) {  // if compiled with dispOn, print these statistics
      fprintf(seq->dispFile, "skipped structures: %d\n", crik->skippedStru);
    }
  }  // end outer if
}  // end print_results

//*****************************************************************************
// Function : close Up
// Caller   : main()
// Purpose  : clear up all the memory allocations
// Input    : seq
//          : crik
// Return   : none
// Display  : none
//*****************************************************************************
int close_up(config* seq, global* crik)
{                                    disp(seq,DISP_ALL,"################################################################################################  Entering 'close_up''\n");
  if(seq->srcFile != stdin)
    fclose(seq->srcFile);            // input  file closing
  if(seq->dispFile != stdout)
    fclose(seq->dispFile);           // output file closing

  free_intrvl_luk_up_table(seq);     // clean up interval look-up table
  free_parenthesis_look_up_table(seq);      // clean up parenthesis look-up table
  free_cmpnt_list(seq, crik);         // clean up component list 'crik->cmpntList'
  free_bundle_list(seq, crik);        // clean up 2D bundle list 'crik->eden'

  free_whole_intrvl(crik);           // clean up interval 'crik->interval'
  free_constraint(seq, crik);        // clean up constraints 'seq->coVari', 'seq->v1Pairng', 'seq->s1Pairng', 'seq->chemMod', etc
  free(crik);
  free_recycle_bin(seq);             // some not-yet cleaned todd stuff
  free(seq->ltr);
  free(seq->motifSeq);
  free(seq->motifStruc);
  free(seq->dotNParen);
  free(seq->mfe);
  free(seq);

  cleanup_stats();

  return 0;
}  // end close_up

//*****************************************************************************
// Function : Display
// Caller   : any
// Purpose  : displays most of the general information, such as status, events, warnings, results etc
// Input    : seq      -- seq->dotNParen only will probably be used
//          : dispPrio -- Display Priority, indicates priority of display
//          : msg      -- "blah blah, whatever contents meant to be displayed, or something like %d, %s, etc"
//          : ...      -- arguments which holds the value/contents to be displayed
// Return   : none
// Display  : printing of data &/or words
//*****************************************************************************
void display(config* seq, int16_t dispPrio, char* msg, ...)
{
  if  (dispPrio < seq->dispMinPrio){return;}    // the priority is not high enough to be shown
  else{                                  // all other displays go this way
    va_list args;
    va_start(args, msg);
    vfprintf(seq->dispFile, msg, args);  // print outs, where the arguments are stored in the variable 'args'
    va_end(args);
  }  // end if
}  // end display

//*****************************************************************************
// Function : Display Structures
// Caller   : make_jump_tree() & prepare_4_another_jump_stage_1()
// Purpose  : display all the hlix in the present structure in the dots and parentheses form
// Input    : seq    - sequence
//          : crik   - crick
// Return   : none
// Display  : dots and parentheses
//*****************************************************************************
void display_structures(config* seq, global* crik, int8_t locationMark)
{                                        disp(seq,DISP_ALL,"Enterng disp_stru from source %d\n", locationMark);

  knob*  hlixCursr = crik->hlixInStru;

  // fix any duplicate helices in hlixCursr
  knob* prev = NULL;
  while (hlixCursr) {
    if (hlixCursr == hlixCursr->jumpTreeNext || (hlixCursr->jumpTreeNext == prev && crik->hlixInStru->jumpTreeNext != NULL)) {
      hlixCursr->jumpTreeNext = NULL;
      continue;
    }
    prev = hlixCursr;
    hlixCursr = hlixCursr->jumpTreeNext;
  }

  hlixCursr = crik->hlixInStru;

  int16_t  flagUB    = seq->numCovari + seq->numV1Pairng;
  int8_t   ASCIIHlixNotation = 65;
  int8_t   asciiHlixNotation = 97;
  int      i;

  int monkeyfist = seq->strLen;
  for(i = 0 ; i < monkeyfist; i++) seq->dotNParen[i] = '.';

  while(hlixCursr){
    if (hlixCursr->bundleFlag){
      knob* repHelix = crik->eden[hlixCursr->opnBrsOutIndx][hlixCursr->closeBrsOutIndx].repNode;
      knob* it = repHelix;
      while(it->bundleListNext) {
        it->jumpTreeNext = it->bundleListNext;
        it = it->bundleListNext;
      }
      it->jumpTreeNext = hlixCursr->jumpTreeNext;
      hlixCursr = repHelix;
    }    

    for(i = hlixCursr->opnBrsOutIndx ; i <= hlixCursr->opnBrsInnIndx ; i++) seq->dotNParen[i] = '(';
    for(i = hlixCursr->closeBrsInnIndx ; i <= hlixCursr->closeBrsOutIndx ; i++) seq->dotNParen[i] = ')';

    if(seq->maxNumMismatch)
      for(i = 0 ; i < seq->maxNumMismatch ; i++)
        if(hlixCursr->mismatchFlag[i][1] >= 0) {
          seq->dotNParen[hlixCursr->mismatchFlag[i][0]] = '_';
          seq->dotNParen[hlixCursr->mismatchFlag[i][1]] = '_';
        }  // end inner if

    ASCIIHlixNotation++;
    asciiHlixNotation++;
    hlixCursr = hlixCursr->jumpTreeNext;
  }  // end while

  update_stats(seq);

  if(seq->dispMinPrio < DISP_LV1){
    fprintf(seq->dispFile, "%s", seq->dotNParen);
    if(flagUB){   
      printf("[ ");
      for(i = 0 ; i < flagUB ; i++){                
        fprintf(seq->dispFile, "%d ", crik->struMustPairFlag[i]);
      }  // end for
      fprintf(seq->dispFile, "] ");
    }  // end if
    hlixCursr = crik->hlixInStru;
    fprintf(seq->dispFile, "\n");
    testParenBal(seq, crik);
    testLupSze(seq, crik);
    testHlixSze(seq, crik);
    locationMark++;    // purely dummy operation
  }  // end if
}  // end display_structures


//*****************************************************************************
// Function : Display Linked List
// Caller   : any
// Purpose  : display the contents of linked lists: Restore Interval, Interval, Structures, toddP Components, todd Components
// Input    : seq  -- seq->dotNParen only will probably be used
//          : crik -- crick
// Return   : none
// Display  : printing of critical features of linked list
//*****************************************************************************
void display_linked_list(config* seq, global* crik, local* todd, local* toddP)
{                                              disp(seq,DISP_ALL,"Entering 'display_linked_list'\n");
  knob* intrvlCursr = crik->interval;
  knob* hlixCursr   = crik->hlixInStru;
  knob* toddCursr   = NULL;
  knob* toddPCursr  = NULL;
  knob* rstoPCursr  = NULL;
  knob* rstoCursr   = NULL;
  int8_t  i;
                                               disp(seq,DISP_LV4,"\n--------------------------------------------");
					       disp(seq,DISP_LV4,"LVL of RECUR = < %d >\n", crik->lvlOfRecur);
  if(todd){
  				               disp(seq,DISP_LV4,"todd   lvRcr = < %d >\n", todd->lvlOfRecur);
                                               disp(seq,DISP_LV4,"                                              ");
					       disp(seq,DISP_LV4,"intrvlCntr =  %d \n", todd->intrvlCntr);
                                               disp(seq,DISP_ALL,"todd intrvl bound = %d-%d\n", todd->intrvlLB, todd->intrvlUB);
  }  // end if 1

  if(toddP){
                                               disp(seq,DISP_ALL,"toddP intrvl bound = %d-%d\n", toddP->intrvlLB, toddP->intrvlUB);
    if(toddP->cmpntLLCursr) toddPCursr = toddP->cmpntLLCursr;
    rstoPCursr = toddP->RSTO;
                                               disp(seq,DISP_LV4,  "Rsto P  Interval = ");     // crik->RSTO P
    while(rstoPCursr){                                                                         // ||
                                               disp(seq,DISP_LV4,"[%3d[ %d-%d-%d ]%-3d] ",     // ||
						    rstoPCursr->opnBrsInnIndx,                 // ||
						    rstoPCursr->hlixBranchngIndx1,             // ||
						    rstoPCursr->intrvlCntr,                    // ||
						    rstoPCursr->lvlOfRecur,                    // ||
						    rstoPCursr->closeBrsInnIndx);              // ||
      rstoPCursr = rstoPCursr->jumpTreeNext;                                                   // ||
    }  // end while 1                                                                          // ||
                                               disp(seq,DISP_LV4,"\n");                        // \/
  }    // end if 2

  if(todd) {
    toddCursr = todd->cmpntLLCursr;
    rstoCursr = todd->RSTO;
                                               disp(seq,DISP_LV4,  "Restore Interval = ");     // crik->RSTO

    while(rstoCursr){                                                                          // ||
                                               disp(seq,DISP_LV4,"[%3d[ %d-%d-%d ]%-3d] ",     // ||
						    rstoCursr->opnBrsInnIndx,                  // ||
						    rstoCursr->hlixBranchngIndx1,              // ||
						    rstoCursr->intrvlCntr,                     // ||
						    rstoCursr->lvlOfRecur,                     // ||
						    rstoCursr->closeBrsInnIndx);               // ||
      if(rstoCursr == rstoCursr->jumpTreeNext) rstoCursr->jumpTreeNext = NULL;
      rstoCursr = rstoCursr->jumpTreeNext;                                                     // ||
    }  // end while 1                                                                          // ||
                                               disp(seq,DISP_LV4,"\n");                        // \/
  }    // end outer if
                                               disp(seq,DISP_LV4,  "Interval         = ");     // crik->interval
  while(intrvlCursr){                                                                          // ||
                                               disp(seq,DISP_LV4,"[%3d[ %d-%d-%d ]%-3d] ",     // ||
						    intrvlCursr->opnBrsInnIndx,                // ||
						    intrvlCursr->hlixBranchngIndx1,            // ||
						    intrvlCursr->intrvlCntr,                   // ||
						    intrvlCursr->lvlOfRecur,                   // ||
						    intrvlCursr->closeBrsInnIndx);             // ||
    intrvlCursr = intrvlCursr->jumpTreeNext;                                                   // ||
  }    // end while 1                                                                          // \/

                                               disp(seq,DISP_LV4,"\nStructure    = ");         // crik->hlixInStru
  while(hlixCursr){                                                                            // ||
                                               disp(seq,DISP_LV4,"{%3d-%-3d{ %2d-%d-%d }%3d-%-3d} ",
						    hlixCursr->opnBrsOutIndx,                  // ||
						    hlixCursr->opnBrsInnIndx,                  // ||
						    hlixCursr->hlixBranchngIndx1,              // ||
						    hlixCursr->intrvlCntr,                     // ||
						    hlixCursr->lvlOfRecur,                     // ||
						    hlixCursr->closeBrsInnIndx,                // ||
						    hlixCursr->closeBrsOutIndx);               // ||
    if (hlixCursr == hlixCursr->jumpTreeNext)                                                  // ||
      break;                                                                                   // ||
    hlixCursr = hlixCursr->jumpTreeNext;	                                               // ||
                                                                                               // ||
  }  // end while 2                                                                            // \/

                                               disp(seq,DISP_LV4,"     ( ");                   // crik->struMustPairFlag
  for(i = 0 ; i < seq->numCovari ; i++){                                                       // ||
    disp(seq,DISP_LV4,"%d ", crik->struMustPairFlag[i]);                                       // ||
  }  // end for                                                                                // ||
                                               disp(seq,DISP_LV4,") ");                        // \/

   
                                               disp(seq,DISP_LV4,"\ntodd         = ");         // todd->cmpntLLCursr
  if(toddCursr){                                                                               // ||
    while(toddCursr){                                                                          // ||
                                               disp(seq,DISP_LV4,"[%3d-%-3d[   ]%3d-%-3d] ",   // ||
						    toddCursr->opnBrsOutIndx,                  // ||
						    toddCursr->opnBrsInnIndx,                  // ||
						    toddCursr->closeBrsInnIndx,                // ||
						    toddCursr->closeBrsOutIndx);               // ||
      if (toddCursr == toddCursr->jumpTreeNext)                                                // ||
        break;                                                                                 // ||
      toddCursr = toddCursr->jumpTreeNext;                                                     // ||
                                                                                               // ||
    }  // end while 3                                                                          // ||
                                               disp(seq,DISP_LV4,"\n    ins: [%d[ %d ]%d]",    // ||
						    todd->intrvlIns->opnBrsInnIndx,            // ||
						    todd->intrvlIns->lvlOfRecur,               // ||
						    todd->intrvlIns->closeBrsInnIndx );        // ||
                                               disp(seq,DISP_LV4,"    beh: [%d[ %d ]%d]\n",    // ||
						    todd->intrvlBeh->opnBrsInnIndx,            // ||
						    todd->intrvlBeh->lvlOfRecur,               // ||
						    todd->intrvlBeh->closeBrsInnIndx );        // ||
  } else {                                     disp(seq,DISP_LV4,"\n");                        // ||
  }  // end if                                                                                 // \/

                                               disp(seq,DISP_LV4,"toddP        = ");           // toddP->cmpntLLCursr
  if(toddPCursr){                                                                              // ||
    while(toddPCursr){                                                                         // ||
                                               disp(seq,DISP_LV4,"[%3d-%-3d[   ]%3d-%-3d] ",   // ||
						    toddPCursr->opnBrsOutIndx,                 // ||
						    toddPCursr->opnBrsInnIndx,                 // ||
						    toddPCursr->closeBrsInnIndx,               // ||
						    toddPCursr->closeBrsOutIndx  );            // ||
      if (toddPCursr == toddPCursr->jumpTreeNext)                                              // ||
        break;                                                                                 // ||
      toddPCursr = toddPCursr->jumpTreeNext;                                                   // ||
    }  // end while 3                                                                          // ||
                                               disp(seq,DISP_LV4,"\n    ins: [%d[ %d ]%d]",    // ||
						    toddP->intrvlIns->opnBrsInnIndx,           // ||
						    toddP->intrvlIns->lvlOfRecur,              // ||
						    toddP->intrvlIns->closeBrsInnIndx );       // ||
                                               disp(seq,DISP_LV4,"    beh: [%d[ %d ]%d]\n",    // ||
						    toddP->intrvlBeh->opnBrsInnIndx,           // ||
						    toddP->intrvlBeh->lvlOfRecur,              // ||
						    toddP->intrvlBeh->closeBrsInnIndx );       // ||
  }    // end if 2                                                                             // ||
                                               disp(seq,DISP_LV4,"\n");                        // \/
}  // end display_linked_list
