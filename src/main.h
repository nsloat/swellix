#ifndef __SWELLIX_MAIN_H
#define __SWELLIX_MAIN_H

/*****************************************************************************************************************************
 *****************************************************************************************************************************
 *****************************************************************************************************************************
      Program : Swellix-main.h
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

#include <stdint.h>
#include <stdio.h>


long long g_x1;
long long g_x2;

#define TRUE                           1
#define FALSE                          0
#define RST_INTRVL                     12    // restore interval
#define INS_INTRVL                     11    // inside interval
#define BEH_INTRVL                     10    // behind interval
#define DEFAULT_HAIR_PIN_LOOP_SIZE     3
#define DEFAULT_MIN_LENGTH_OF_HELIX    1
#define DEFAULT_MAX_PAIRING_DISTANCE   16384 // just an arbitrary huge number (2^14), to be updated by seq->strLen during initialization process
#define DEFAULT_MAX_NUMBER_OF_BULGE    0
#define DEFAULT_MAX_NUMBER_OF_MISMATCH 0
#define DEFAULT_MIN_NUMBER_OF_HELIX    0
#define DEFAULT_MIN_NUMBER_OF_HP       0
#define CONSTRAINT_SIZE                64    // just an arbitrary number (2^6) large enought to accomodate largest possble number of constraints, used in init_constraint.c
#define BAKI_SORT_DIVIDR               8     // bucket sort divider                        : used for bucket sort, to determine the size of each bucket

#ifdef  _display
#define disp(...)                      display(__VA_ARGS__)
#define dispLL(...)                    display_linked_list(__VA_ARGS__)
#define dispConstraints(...)           display_constraints(__VA_ARGS__)
#define dispCmpnt(...)                 display_components(__VA_ARGS__)
#define dispEden(...)                  display_eden(__VA_ARGS__)
#define dispStru(...)                  display_structures(__VA_ARGS__)
#define dispFullBaki(...)              display_full_buckets(__VA_ARGS__)
#define DISP                           TRUE 
#else   
#define disp(...)
#define dispLL(...)
#define dispConstraints(...)        
#define dispCmpnt(...)      
#define dispEden(...)
#define dispStru(...)                  display_structures(__VA_ARGS__)
#define dispFullBaki(...)              
#define DISP                           FALSE
#endif

#ifdef  _report
#define disp_stru(...)                 display_structures(__VA_ARGS__)
#else
#define disp_stru(...)
#endif

#ifdef  _test
#define testParenBal(...)              test_parentheses_balance(__VA_ARGS__)
#define testLupSze(...)                test_loop_size(__VA_ARGS__)
#define testHlixOvrlap(...)            test_helix_overlap(__VA_ARGS__)
#define testHlixSze(...)               test_helix_size(__VA_ARGS__)
#else
#define testParenBal(...)
#define testLupSze(...)
#define testHlixOvrlap(...)
#define testHlixSze(...)
#endif

#ifdef  _greeting
#define greeting(...)                  printf("Hello Formosa!\n");    
#else
#define greeting(...)                  
#endif

//#if __GNUC_PREREQ == (2, 7)               // used to check if the user's compiler is gcc, if yes, then __unused will be in the form of '__attribute__((__unused__))' to ignore unused 'config* seq', but if the compiler isn't gcc, this one is unnessasary
//#define __unused                       __attribute__((__unused__))
//#else
#define __unused
//#endif 

#ifdef _PARAM
#define paramfile  _PARAM
#else
#define paramfile  NULL
#endif

enum dispPrio { // Define: indicate the complexity of display (display priority), used by seq->dispMinPrio
  DISP_ALL,
  DISP_LV5,
  DISP_LV4,
  DISP_LV3,
  DISP_LV2,
  DISP_LV1
};

enum callerFlags {
   FIRST_SESSION,
   SECND_SESSION,
   WHILE_SESSION,
   RSTO_SESSION
};

enum algoMode { // Define: indicate the mode and goal of calculation (algorithm mode), used by seq->algoMode
  MODE_NONE,        //             "    none
  MODE_CMPNT,       //          compute for components only
  MODE_INTAB,       //             "    until interval look-up table (by-pass jump bush)
  MODE_STRU         //             "    for all structures (by-pass jump bush)
};

enum statMode {
  STAT_DEFAULT,
  STAT_MAX_DISTANCE,
  STAT_MAX_FREE_ENERGY,
  STAT_ALL
};

extern char* NODE_TYPE_TO_STRING[];

typedef struct LabeledStructures {
  int titlesize;    // the size of the title string
  int buffsize;     // the size of the structures buffer allocated
  char* title;      // the #x#.lab file name string
  char* structures; // newline-delimited list of output from the subopt+label process
} LabeledStructures;

typedef struct knob       knob;
typedef struct interval   interval;
typedef struct local      local;
typedef struct ToL        ToL;

struct interval { // Define: the linked list, direct adoption from crumple type 'interval_struct'
  int16_t   cmpntTyp;    // component type
  interval* prev;        // mark the linked list card before the present one
  interval* next;        //            "              after         "
};

struct knob { // Define: the linked list, direct adoption from crumple type 'interval_struct'
  int16_t** mismatchFlag;       // used for cmpntList, to store which nucleotide contains the mismatch
  int16_t*  mustPairFlag;       // used for cmpntList, to store what kind of covari/V1 pairing nucleotides this particular cmpnt contains 
  int64_t   intrvlCntr;         // act like batch number, used to distinguish one batch from another, to facilitate proper new/old helix addition/removal. Now it's also used by make_bundle_list() as serial number of a structure/stem
  int64_t   bundleCntr;         // act like bundle serial number
  int64_t   parentIntrvlCntr;   // record the intrvlCntr of its parent interval, so that when parent hlix is gone, there'll be no more interval restore
  int16_t   closeBrsInnIndx;    // close  "   inner   "              : close       "             :      "          close      trailing   "   
  int16_t   closeBrsOutIndx;    // close  "   outer   "              : close       "             :      "          close      trailing   "   
  int16_t   hlixBranchngIndx1;  // distinguish one branch of interval from another, which comes from other interval branches from the recur lvl 0
  int16_t   hlixBranchngIndx2;  //                     "                      "                        "                       "                1
  int16_t   opnBrsInnIndx;      // open   "   inner   "              : open brace outer location : location of the open brace leading letter
  int16_t   opnBrsOutIndx;      // open brace outer index            : open brace outer location : location of the open brace leading letter
  int16_t   outsideIntrvlLB;    // outside interval lower bound      : for the interval housing this knob
  int16_t   outsideIntrvlUB;    // outside interval upper bound      : for the interval housing this knob
  int8_t    intrvlTypFlag;      // interval type flag                : 0 for behind interval, and 1 for inside interval
  int8_t    intrvlInsFormdFlag; // inside interval formed            : to inform behind intrvl to raise intrvlBothFormdFlag, if behind interval may be formed, too
  int16_t   lvlOfRecur;         // record the level of recursion at which this helix was made
  int8_t    rstoOnQFlag;        // restore on queue flag             : whichever interval holding this is rsto intrvl on Q
  int8_t    specialRstoFlag;    // special restore flag              : accompany with specialRstoIntrvl, mark the presence of this interval for special treatment during interval restore
  int8_t    bundleFlag;         // bundle flag                       : used to remind 'display_structure' that this one is from bundle list and therefore handle specially
  int       newCLindex;         // new component list index          : new, array-based implementation of component list. Better suited for MPI processing.
  knob*     cmpntListNext;      // component list next               : for use in component list 'crik->cmpntList'
  knob*     bundleListNext;     // bundle list next                  : for use in bundle list (eden) for simple (pure nested, not side-by-side) recursion
  knob*     jumpTreeNext;       // jump tree next                    : for use in jump tree to link the intervals
};

typedef struct {
  knob*       knob;             // sub-type                          : linked list growing out of the backbone of the edge list containing the components it contains under this component type
  int16_t     cnt;              // count                             : total numbers of components in that component type
} edgeList;

struct ToL { // Define: tree of life
  knob*       branchNode;       // branch node                       : One branch is one structure. One branch is composed of many branch nodes.
  int16_t*    mustPairFlag;     // must pair flag                    : the combined result of all the flags of the helices of this structure
  int16_t     numHlix;          // number of helices                 : number of helices in this particular structure
  int16_t     bundleCntr;       // bundle counter                    : serial number of this structure/stem
  ToL*        next;
  int8_t      numOfMismatches;  // number of mismatches in the outer helix
};         

typedef struct {
  ToL*    stemNode;             // stem node                         : One tree is one bundle. One tree is composed of many stem nodes, from each of which one complete structure extends. 
                                //                                      Precisely, each node contains one structure, containing many helices.
  knob*   dumiNode;             // dummy node                        : used by 'make_jump_tree' to display the span of a particular bundle, to facilitate SSR (similar structure reduction) operation
  knob*   repNode;              // representative node               : instead of printing each bundled structure, identify one node to represent a whole bundle. 
                                //                                      This representative is to be chosen based on which structure in the bundle has the most base pairs.
  int16_t numStru;              // number of structures              : height of the tree, namely, number of nodes. which indicates the number of structures on this tree
  int8_t  dumiUsdFlag4RecurLv0; // dummy used flag for recur level 0 : record all dummies which have been used, so that all the later smaller dummies which may be covered by previous larger dummies may be skipped from printing
} thuong;                       // beloved


typedef struct { // Define: the main type of specimen (sequence) to be tested
  FILE*      srcFile;           // source file                       : file pointer to the exterior source file of sequence, if Src_Of_SeqFlag = FILE_IN
  FILE*      dispFile;          // display file                      :       "         "        file storing all possible pairing structures
  char*      ltr;               // letter                            : where the letters of A, C, G, U are stored
  char*      dotNParen;         // dot & parentheses                 : where the symbols '[__[..]___]...' or 'AA.BBB....bb...aa' are stored
  int16_t*** intrvlLukUpTable;  // interval look-up table            : used to store all components which fits in the interval
  int16_t**  coVari;            // covariance                        : store the data of covariance pairs
  int16_t*   s1Pairng;          // S1 Pairing                        : store the data of S1 pairing nucleotides, must not pair nucleotides
  int16_t*   v1Pairng;          // V1 Pairing                        :        "          V1      "         "   , must pair nucleotides
  int8_t**   parenLukUpTable;   // parenthesis look-up table         : used to store all pairing comparison for every single pair of nucleotides
  int16_t    strLen;            // string length                     : true string length of the sequence
  int16_t    minPairngDist;     // minimum pairing distance          : minimum space between open and close brace, which is also the MIN PAIRNG DISTANCE in the conf file. It was named "in hairpin loop size"
  int16_t    maxPairngDist;     // maximum pairing distance          : maximum space    "               "                            MAX    "      "
  int8_t     bundle;            // bundle                            : bundle flag to activate the bundling mechanism -- similar structures reduction referencing (SSRR) 
  int8_t     unbundle;		// unbundle			     : unbundling flag -- for debugging/verifying combinatorial completeness of bundling mechanism.
  int8_t     statMode;          // statistics mode                   : the state of this flag determines which statistical calculations are performed
  int8_t     maxdist;
  float      minenergy;
  char*      mfe;
  int8_t     motif;
  uint64_t   motifCount;
  char*      motifSeq;
  char*      motifStruc;
  int8_t     maxNumBlg;         // maximum number of bulge           : maximum count of bulges per helix
  int8_t     maxNumMismatch;    // maximum number of mismatch        : maximum count of mismatches per helix
  int8_t     minNumHlix;        // minimum number of helix           : minimum count of helices in one structure to be eligible to print8_tout
  int8_t     minNumHP;          // minimum number of hairpins        : minimum count of hairpins in one structure to be eligible to print8_tout
  int8_t     minLenOfHlix;      // minimum length of helix           : minimum count of nucleotides in one helix
  int8_t     minBtwnHlixDist;   // minimum between helix distance    : used primarily in STMV, since there has to be a min length of 2 nt in order to turn from one helix to next
  int8_t     noGU;              // no G-U allowed                    : disallow G-U pair to be formed in a sequence
  int8_t     noSZG;             // no sizing allowed                 : don't expand (swell) the size of the helix unless necessary (i.e. bulge)
  int8_t     sudona;            // pseudoknot                        : allow the existence of pseudoknot configuration
  int8_t     dispMinPrio;       // display minimum priority          : minimum priority level required to be displayed on 'dispFile'
  int8_t     algoMode;          // mode of algorithm                 : depending on this flag, Swellix may perform some parts of the algorithm but not others.
  int8_t     insIntrvlMinSize;  // inside interval minimum size      : used mainly for making intervals inside
  int8_t     behIntrvlMinSize;  // behind interval minimum size      : used mainly for making intervals behind
  int16_t    numCovari;         // number of covariance              : accompany the 2x2 array seq->covari, value '0' as default
  int16_t    numS1Pairng;       // number of S1 pairings             : accompany the vector seq->s1Pairng, value '0' as default
  int16_t    numV1Pairng;       // number of V1 pairings             :            "         seq->v1Pairng                "                     
  int16_t    constraintActive;  // constraints active                : determined by the presence of -k option: with '-k' - 1: active; w/o '-k' -  0: non-active
  local*     recycleBin;        // recycle bin                       : used to troubleshoot the bug associate with covariance
} config;

typedef struct { // Define: the parameter book-keeper of global nature (only one in the whole run)
  thuong**     eden;                // eden                              : data type is mini_jump_tree n by n 2D plane with many trees, where n is sequence length, to store the bundle list, as short cut substitute of component list
  edgeList*    cmpntList;           // component list                    : store the list of all components, which is based on the structure of edge list of the graph theory
  knob*        interval;            // interval                          : used to run the recursion inside and behind
  knob*        hlixInStru;          // helices in structure              : final helices for print out
  int16_t*     struMustPairFlag;    // structure must-pair flag          : used to track all the covariance pairs V1 pairing pairs present in the structure
  int          mustPairLength;      // must-pair flag array length
  int16_t*     cmpntListOcupidTyp;  // component list occupied type      : record which types are occupied, so searching for loop may be faster for interval seeking
  int64_t      intrvlCntr;          // interval counter                  : act like batch number, used to distinguish one batch from another, to facilitate proper new/old helix addition/removal
  int64_t      numCmpnt;            // number of component               : keep track of the size of 'crik->cmpntList', used in making components
  uint64_t     numStru;             // number of structures              : the total count of the structures (or if bundling is on, total structures with bundles)
  uint64_t     numUnbundledStru;    // number of total unbundled structures : the total count of all the possible structures
  uint64_t     numBundles; 	    // number of bundles				 : total number of bundles
  int16_t      lvlOfRecur;          // level of recursion                : depth of roots into recursion loops from ground level
  int16_t      closeBrsInnIndx;     // close  "   inner   "              : close       "            
  int16_t      closeBrsOutIndx;     // close  "   outer   "              : close       "            
  int16_t      closeBrsWidth;       // close  "   width                  : size of close brace
  int16_t      closeParenIndx;      // close     "        "              : location of close      "
  int16_t      closeBrsCursr;       // close  "     "                    :             "               close                 "                    close brace 
  int16_t      numCmpntTyp;         // number of component types         : keep track of the size of 'crik->cmpntList', used in making components
  int16_t      numCmpntTypOcupid;   // number of component types used    : no specific purpose, may be used as a counter or something
  int16_t      numHlix;             // number of helices                 : keep track of number of helices in a structure, used in 'mix_n_match_4_stru' & 'hash table 4 stru'
  int16_t      numHP;               // number of hairpin                 : in attemp to build 4th dimension on top of the 3D hash table: CG location and mass and ip
  int16_t      opnBrsCursr;         // open brace cursor                 : used to scan thru the whole open brace to probe each nucleotide inside open  brace
  int16_t      opnBrsInnIndx;       // open   "   inner   "              : open brace outer location
  int16_t      opnBrsOutIndx;       // open brace outer index            : open brace outer location
  int16_t      opnBrsStop;
  int16_t      opnBrsWidth;         // open   "   width                  : size of open  brace
  int16_t      opnParenIndx;        // open parenthesis index            : location of open parenthesis cadidate 
  int16_t      test1ErrTotal;       // test 1 error total                : sum parentheses balance test error total occurences
  int16_t      test2ErrTotal;       // test 2 error total                : sum loop size test error total occurences
  int16_t      test3ErrTotal;       // test 3 error total                : sum helix overlap test error total occurences
  int16_t      test4ErrTotal;       // test 4 error total                : sum ????? test error total occurences
  int8_t       specialRstoFlag;     // special restore flag              : accompany with specialRstoIntrvl, mark the presence of this i
  uint64_t       rstoCounter;
  uint64_t     rstoErrCounter;
  int16_t      skippedStru;         // debugging purposes				 : number of structures skipped during bundling
//  int8_t       struLinked;          // boolean for whether or not the previous helix is linked to current helix
  int8_t       linkedmms;           // number of linked mismatches
  knob**       mpiCList;            // array-ified component list pointers for easier cross-PE communication
} global;

struct local { // Define: the parameter book-keeper of local nature (one per recursion)
  knob*      cmpntLLCursr;       // component linked list cursor
  knob*      intrvlIns;          // interval inside
  knob*      intrvlBeh;          // interval behind
  knob*      RSTO;
  int64_t    intrvlCntr;
  int16_t    intrvlUB;           // interval upper bound
  int16_t    intrvlLB;           // interval lower bound
  int16_t    lukUpCmpntTypUB;    // look-up component type upper bound: used in jump_stage_2_fit_hlix
  int16_t    lukUpCmpntTypLB;    // look-up component type lower bound:          "            "
  int8_t     intrvlInsFormdFlag; // inside interval formed            : to inform behind intrvl to raise intrvlBothFormdFlag, if behind interval may be formed, too
  int16_t    lvlOfRecur;         // level of recursion
  local*     next;
};

int          close_up(config* seq, global* crik);
void         display(config* seq, int16_t dispPrio, char* msg, ...);
void         display_eden(config* seq, global* crik);
void         display_linked_list(config* seq, global* crik, local* todd, local* toddP);
void         display_progress_periodically(global* crik, int8_t locationMark);
void         display_structures(config* seq, global* crik, int8_t locationMark);
global*      initialize_crik(config* seq);
int          initialize_eden(config* seq, global* crik);
config*      initialize_sequence_by_getting_argument(int argc, char** argv);
int          initialize_sequence_by_gettng_constraint(config* seq);
int          initialize_sequence_by_settng_memory(config* seq);
int          make_component_list(config* seq, global* crik);
int          make_bundle_list(config* seq, global* crik);
int          make_interval_look_up_table(config* seq, global* crik);
int          make_jump_tree(config* seq, global* crik, int start, int end);
int          make_parenthesis_look_up_table(config* seq);
void         print_results(config* seq, global* crik);
void         print_usage(void);
int          set_outer_interval(config* seq, global* crik);

#endif // __SWELLIX_MAIN_H

