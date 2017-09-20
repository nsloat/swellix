/*****************************************************************************************************************************
 *****************************************************************************************************************************
 *****************************************************************************************************************************
      Program : Swellix (Initialization - General) - c file
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "main.h"
#include "init_general.h"

extern int rank;

//*****************************************************************************
// Function : Initialize Command Line Argument Default Values
// Caller   : initialize_sequence_by_getting_argument()
// Purpose  : assign the default values to some variables, in case those variables aren't assigned via command line options
// Input    : 
// Return   : none
// Display  : none
//*****************************************************************************
int initialize_commandline_argument_with_default_value(config* seq)
{
  seq->algoMode          = MODE_STRU;                         // see main.h for this value definition                                                               
  seq->bundle            = FALSE;                             // by default, bundle feature isn't activated, which means SSRR (similar structure reduction referencing isn't activated)
  seq->unbundle		     = FALSE;			      // by default, unbundling is off. Only activate when bundling is also activated and only for verifying that bundling is 
							      // accurately computing all relevant structures
  seq->motif             = FALSE;
  seq->motifCount        = 0;
  seq->motifSeq          = NULL;
  seq->motifStruc        = NULL;
  seq->constraintActive  = 0;                                 // by default, there's no constraint input
  seq->dispFile          = stdout;                            // display file. By default, there's none
  seq->dispMinPrio       = DISP_LV1;                          // set up default display priority here
  seq->maxNumBlg         = DEFAULT_MAX_NUMBER_OF_BULGE;       // see main.h for this value definition
  seq->maxNumMismatch    = DEFAULT_MAX_NUMBER_OF_MISMATCH;    // see main.h for this value definition
  seq->maxPairngDist     = DEFAULT_MAX_PAIRING_DISTANCE;      // just an arbitrary number, to be updated by seq->strLen
  seq->minBtwnHlixDist   = 0;                                 // by default, two side-by-side helices may be right against each other
  seq->minLenOfHlix      = DEFAULT_MIN_LENGTH_OF_HELIX;       // see main.h for this value definition
  seq->minNumHlix        = DEFAULT_MIN_NUMBER_OF_HELIX;       // see main.h for this value definition
  seq->minNumHP          = DEFAULT_MIN_NUMBER_OF_HP;          // see main.h for this value definition
  seq->minPairngDist     = DEFAULT_HAIR_PIN_LOOP_SIZE;        // see main.h for this value definition
  seq->noGU              = FALSE;                             // by default, GU pair is acceptable
  seq->noSZG             = FALSE;                             // by default, there's no expansion (swelling) of helix
  seq->numCovari         = 0;                                 // by default, there's no covariance
  seq->numS1Pairng       = 0;                                 // by default, there's no S1 Pairing nucleotide
  seq->numV1Pairng       = 0;                                 // by default, there's no V1 Pairing nucleotide
  seq->srcFile           = stdin;                             // source  file. By default, there's none
  seq->sudona            = FALSE;                             // by default, there's no pseudoknot formation allowed
  seq->statMode          = STAT_DEFAULT;

  return 0;
}  // end initialize_commandline_argument_with_default_value

//*****************************************************************************
// Function : Set Degree of Detailed Display
// Caller   : initialize_sequence_by_getting_argument()
// Purpose  : determine how much detail will be displayed during the present execution time
// Input    : 
// Return   : none
// Display  : none
//*****************************************************************************
int set_degree_of_detailed_disp(config* seq, int8_t arg, int argc, char** argv)
{
  if(++arg == argc)
    cmd_line_err(argv, arg);
  else if (!strcmp("all",     argv[arg]) || !strcmp("a", argv[arg]))
    seq->dispMinPrio = DISP_ALL;                                        // maximum display option
  else if (!strcmp("level 5", argv[arg]) || !strcmp("5", argv[arg]))
    seq->dispMinPrio = DISP_LV5;                                
  else if (!strcmp("level 4", argv[arg]) || !strcmp("4", argv[arg]))
    seq->dispMinPrio = DISP_LV4;                                
  else if (!strcmp("level 3", argv[arg]) || !strcmp("3", argv[arg]))
    seq->dispMinPrio = DISP_LV3;                                        // print cmpnt list, jump bush and stru
  else if (!strcmp("level 2", argv[arg]) || !strcmp("2", argv[arg]))
    seq->dispMinPrio = DISP_LV2;                                        // print structures
  else if (!strcmp("level 1", argv[arg]) || !strcmp("1", argv[arg]))
    seq->dispMinPrio = DISP_LV1;                                        // print statistics data and error, if any
  else {
    fprintf(stderr, "\n\nUnknown option '%s %s'\n", argv[arg -1 ], argv[arg]);
    print_usage();
  }  // end inner if 

  return 0;
}  // end set_degree_of_detailed_disp

//*****************************************************************************
// Function : Set Mode of Statistical Output
// Caller   : initialize_sequence_by_getting_argument()
// Purpose  : determine which statistical parameters to calculate during runtime
// Input    : config* seq to set global values, commandline parameters arg, argc, argv
// Return   : none
// Display  : none
//*****************************************************************************
void set_statMode(config* seq, int8_t arg, int argc, char** argv) 
{
  if(++arg == argc)
    cmd_line_err(argv, arg);
  else if(!strcmp("1", argv[arg]))
    seq->statMode = STAT_DEFAULT;
  else if(!strcmp("2", argv[arg]))
    seq->statMode = STAT_MAX_DISTANCE;
  else if(!strcmp("3", argv[arg]))
    seq->statMode = STAT_MAX_FREE_ENERGY;
  else if(!strcmp("all", argv[arg]) || !strcmp("a", argv[arg]))
    seq->statMode = STAT_ALL;
  else {
    fprintf(stderr, "\n\nUnknown option '%s %s'\n", argv[arg-1], argv[arg]);
    print_usage();
  }
}

void set_motif(config* seq, int8_t arg, int argc, char** argv)
{
  if(++arg == argc)
    cmd_line_err(argv, arg);
  else {
    if(rank == 0) fprintf(seq->dispFile,"\nmotif read from cmdline    : %s\n", argv[arg]);
    char* tmp;
    char* tok = strtok(argv[arg], "&");
    
    unsigned int len1 = strlen(tok);
    seq->motifStruc = calloc(len1+1, sizeof(char));
    strcpy(seq->motifStruc, tok);
    tok = strtok(NULL, "&");
    char* tok1 = strtok(seq->motifStruc, "x");
    if(len1 != strlen(tok1)) { 
      tmp = tok1; 
      tok1 = strtok(NULL, "x"); 
      sprintf(seq->motifStruc, "%sx%s", tmp, tok1);
      len1 = strlen(seq->motifStruc);
    }

    if(tok == NULL) { 
      fprintf(stderr,"\nERROR:\nAn error occurred with the formatting of the command-line motif.\nCheck for quotation marks around the motif and ensure that a '&' separates the\nletters from the dot-parenthesis structure.\nRun Swellix with -h to see the motif reference input.\n\n");
      free(seq->motifStruc);
      fclose(seq->srcFile);
      free(seq);
      exit(1);
    }
    unsigned int len2 = strlen(tok);
    seq->motifSeq = calloc(len2+1, sizeof(char));
    strcpy(seq->motifSeq, tok);
    tok1 = strtok(seq->motifSeq, "x");
    if(len2 != strlen(tok1)) { 
      tmp = tok1; 
      tok1 = strtok(NULL, "x"); 
      sprintf(seq->motifSeq, "%sx%s", tmp, tok1);
      len2 = strlen(seq->motifSeq);
    }

    if(len1 != len2) {
      fprintf(stderr, "\nERROR:\nThe motif pieces are not the same length. Ensure that the length of the RNA sequence\nmatches the length of the matching structure you want to find.\nRun Swellix with -h to see the motif reference input.\n\n");
      free(seq->motifSeq);
      free(seq->motifStruc);
      fclose(seq->srcFile);
      free(seq);
      exit(1);//      print_usage();
    }
    if(rank==0) printf("motif stored: xxxxxxx%sxxxxxxx\n              xxxxxxx%sxxxxxxx\n", seq->motifSeq, seq->motifStruc);
  }
}

//*****************************************************************************
// Function : Set Mode of Algorithm
// Caller   : initialize_sequence_by_getting_argument()
// Purpose  : determine which stage (how far) will be reached during the present execution time
// Input    : 
// Return   : none
// Display  : none
//*****************************************************************************
int set_mode_of_algo(config* seq, int8_t arg, int argc, char** argv)
{
  if(++arg == argc)
    cmd_line_err(argv, arg);
  else if (!strcmp("cmpnt",   argv[arg]) || !strcmp("c", argv[arg]) || !strcmp("component", argv[arg]))
    seq->algoMode = MODE_CMPNT;                                                                                       // compute for components only
  else if (!strcmp("intab",   argv[arg]) || !strcmp("i", argv[arg]) || !strcmp("interval_look_up_table", argv[arg]))
    seq->algoMode = MODE_INTAB;                                                                                       // compute for interval look-up table level
  else if (!strcmp("none",    argv[arg]) || !strcmp("n", argv[arg]))
    seq->algoMode = MODE_NONE;                                                                                        // compute for nothing
  else if (!strcmp("stru",    argv[arg]) || !strcmp("s", argv[arg]) || !strcmp("structure", argv[arg]))
    seq->algoMode = MODE_STRU;                                                                                        // compute for all structures (by-pass bush)
  else {
    fprintf(stderr, "\n\nUnknown option '--mode %s'\n", argv[arg]);
    print_usage();
  }  // end inner if 

  return 0;
}  // end set_mode_of_algo

//*****************************************************************************
// Function : Command Line Error
// Caller   : initialize_sequence_by_getting_argument()
// Purpose  : report the error and print out the instruction table during the command input session
// Input    : argv - the array with all the command contents
//          : arg  - the index counter refering to the address of the argv[]
// Return   : none
// Display  : error message concerning the sympton
//*****************************************************************************
void cmd_line_err(char** argv, int arg)
{
  printf("\n\nError: The option '%s' needs to go with some parameter(s). See table below for more information\n\n", argv[arg - 1]);
  print_usage();
}  // end cmd_line_err

//*****************************************************************************
// Function : Initialize Sequence -> Letters
// Caller   : initialize_sequence_by_settng_memory()
// Purpose  : evaluate and initialize the memory allocations for 'seq->ltr', which stores the sequence in English letters 
// Input    : seq
// Return   : none
// Display  : request for sequence, if not given via file
//*****************************************************************************
int init_seq_ltr(config* seq)
{                                                                           disp(seq,DISP_ALL,"Entering 'init_seq_ltr'\n");
  seq->ltr = calloc(1000000, sizeof(char));                                 // temporarily assign arbitrary 1000000 Bytes to put letters

  if(!seq->srcFile){                                                        // inquire user to input sequence, in case the sequence file isn't given
    seq->srcFile = stdin;
    printf("No sequence given. Please input sequence here: ");
  }  // end if

  fscanf(seq->srcFile, "%s\n",seq->ltr);                                      
  seq->ltr = realloc(seq->ltr, (strlen(seq->ltr) + 1) * sizeof(char));      // resize the seq->ltr to just enough to fit the string
  seq->mfe = calloc(strlen(seq->ltr)+1, sizeof(char));
  seq->recycleBin = NULL;

  return 0;
}  // end init_seq_ltr__is

