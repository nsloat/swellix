/*****************************************************************************************************************************************************
 *****************************************************************************************************************************************************
 *****************************************************************************************************************************************************
      Program : Swellix (Make Bundle List2 - created by Isaac) - c file
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
#include <unistd.h>
#include <ctype.h>
#include <time.h>
#include "main.h"
#include "bundle_list.h"

//*****************************************************************************
// Function : Initialize Eden
// Caller   : make_bundle_list()
// Purpose  : initialize the memory allocation of the eden ground
// Input    : seq
//          : crik
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
int initialize_eden(config* seq, global* crik)
{
	disp(seq,DISP_ALL,"Enterng 'initialize_eden'");
	int16_t i;

	crik->eden = calloc(seq->strLen, sizeof(thuong*));      // the main plain of 'edge list' to store all the bundle fruits
                                                          // arbitrarily assgin some number of rooms to accomodate unknown number of bundle plans
	for(i = 0 ; i < seq->strLen ; i++)
		crik->eden[i] = calloc(seq->strLen, sizeof(thuong));

	return 0;
}  // end initialize_eden

//*****************************************************************************
// Function : Get Paren Indices
// Caller   : make_bundles()
// Purpose  : Gets the index values for mismatches (if any) in the given structure
// Input    : structure, len, openOut, openIn, closeOut, closeIn
//            openOut - open parenthesis, outer index
//            openIn  - "     "           inner index
//            closeOut - close parenthesis, outer index
//            closeIn  - "      "           inner index
// 			  mismatchFlags - marks the location of mismatches
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
int get_mismatches(char* structure, int openOut, int openIn, int closeOut, int closeIn, int16_t** mismatch) {
	int i;
	// Find any mismatches
	int j = 0;
	int toReturn = 0;
	for (i = openOut; i <= openIn; i++) {
		if (i < 0) {
			mismatch[j][0] = i;
			j++;
		} else if (structure[i] == '.') {
	//		fprintf(seq->dispFile, "j = %d\n", j);
			mismatch[j][0] = i;
			j++;

			if (i == openOut+1 || i == openIn-1) {
				toReturn = 1;
			}
		}
	}

	j = 0;
	for (i = closeOut; i >= closeIn; i--) {
		if (i >= (int)strlen(structure)) {
			mismatch[j][1] = i;
			j++;
		} else if (structure[i] == '.') {
			mismatch[j][1] = i;
			j++;
		}
	}

	return toReturn;
}

//*****************************************************************************
// Function : Remove outer pair
// Caller   : make_bundles()
// Purpose  : Removes the outer pair of parentheses
// Input    : structure, start, end
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
void remove_outer_pair(char* structure, int start, int end) {
	structure[start] = '.';
	structure[end] = '.';
}

//*****************************************************************************
// Function : Filename to Indices
// Caller   : add_dumi_node()
// Purpose  : Parses the indices from a given filename: POSxLEN.lab
// Input    : filename, position, length
// Return   : none
// Display  : none
//*****************************************************************************
void add_dumi_node(config* seq, global* crik, char* filename) {
	int begin, end;
	filename_to_indices(filename, &begin, &end);

	knob* dumi;
	int16_t hlixLen = (end - begin + 1 - seq->minPairngDist) >> 1;

	crik->eden[begin][end].dumiNode = calloc(1, sizeof(knob));

	dumi                    = crik->eden[begin][end].dumiNode;
	dumi->opnBrsOutIndx     = begin;
	dumi->opnBrsInnIndx     = begin + hlixLen;
	dumi->closeBrsInnIndx   = end - hlixLen;
	dumi->closeBrsOutIndx   = end;
	dumi->lvlOfRecur        = -1;
	dumi->hlixBranchngIndx1 = -1;
	dumi->intrvlCntr        = -1;
	dumi->bundleCntr        = -1;
	dumi->outsideIntrvlLB   = -1;
	dumi->outsideIntrvlUB   = -1;
	dumi->specialRstoFlag   =  0;
	dumi->rstoOnQFlag       =  0;
	dumi->bundleFlag        =  1;
	dumi->bundleListNext    = NULL;
	dumi->jumpTreeNext      = NULL;
}

void filename_to_indices(char* filename, int* begin, int* end) {
	char pos[50], len[50], suffix[4];
	sscanf(filename, "%[^x]x%[^.].%s", pos, len, suffix);
	*end = atoi(pos);
	*begin = atoi(len);
	*begin = *end - *begin + 1;
}

void remove_duplicates(char* filename) {
	char tmp[100];
	strcpy(tmp, "sort -u ");
	strcat(tmp, filename);
	strcat(tmp, " > tmp.txt");
	system(tmp);
	strcpy(tmp, "rm ");
	strcat(tmp, filename);
	system(tmp);
	strcpy(tmp, "mv tmp.txt ");
	strcat(tmp, filename);
	system(tmp);
}

/* Function to test whether two nodes are equal */
int nodes_equal(knob* first, knob* second) {

	// Compare combined helices as well
	knob* firstCursr = first;
	int openIn1 = first->opnBrsInnIndx;
	int openOut1 = first->opnBrsOutIndx;
	int closeIn1 = first->closeBrsInnIndx;
	int closeOut1 = first->closeBrsOutIndx;
	int formed = 0;
	firstCursr = firstCursr->bundleListNext;
	while(firstCursr) {
		if (openOut1 - firstCursr->opnBrsInnIndx == 1 && firstCursr->closeBrsInnIndx - closeOut1 == 1) {
			openOut1 = firstCursr->opnBrsOutIndx;
			closeOut1 = firstCursr->closeBrsOutIndx;
			formed = 1;

		//	printf("formed: %d|%d-%d|%d\n", openOut1, openIn1, closeIn1, closeOut1);

		} else if (firstCursr->opnBrsOutIndx - openIn1 == 1 && closeIn1 - firstCursr->closeBrsOutIndx == 1) {
			openIn1 = firstCursr->opnBrsInnIndx;
			closeIn1 = firstCursr->closeBrsInnIndx;
			formed = 1;

	//		printf("formed: %d|%d-%d|%d\n", openOut1, openIn1, closeIn1, closeOut1);

		} else if (!formed) {
			openIn1 = firstCursr->opnBrsInnIndx;
			openOut1 = firstCursr->opnBrsOutIndx;
			closeIn1 = firstCursr->closeBrsInnIndx;
			closeOut1 = firstCursr->closeBrsOutIndx;
		}

		firstCursr = firstCursr->bundleListNext;
	}

	knob* secondCursr = second;
	int openIn2 = second->opnBrsInnIndx;
	int openOut2 = second->opnBrsOutIndx;
	int closeIn2 = second->closeBrsInnIndx;
	int closeOut2 = second->closeBrsOutIndx;
	formed = 0;
	secondCursr = secondCursr->bundleListNext;
	while(secondCursr) {
		if (openOut2 - secondCursr->opnBrsInnIndx == 1 && secondCursr->closeBrsInnIndx - closeOut2 == 1) {
			openOut2 = secondCursr->opnBrsOutIndx;
			closeOut2 = secondCursr->closeBrsOutIndx;
			formed = 1;

		//	printf("formed: %d|%d-%d|%d\n", openOut2, openIn2, closeIn2, closeOut2);

		} else if (!formed) {
			openIn2 = secondCursr->opnBrsInnIndx;
			openOut2 = secondCursr->opnBrsOutIndx;
			closeIn2 = secondCursr->closeBrsInnIndx;
			closeOut2 = secondCursr->closeBrsOutIndx;
		}

		secondCursr = secondCursr->bundleListNext;
	}

	int escaped1 = 0;
	int escaped2 = 0;
	while (first && second) {
		if ((first->opnBrsInnIndx == second->opnBrsInnIndx && first->opnBrsOutIndx == second->opnBrsOutIndx && first->closeBrsInnIndx == second->closeBrsInnIndx && first->closeBrsOutIndx == second->closeBrsOutIndx)
				|| (openOut1 == openOut2 && openIn1 == openIn2 && closeIn1 == closeIn2 && closeOut1 == closeOut2)) {
			if (!escaped1 && openOut1 == openOut2 && openIn1 == openIn2 && closeIn1 == closeIn2 && closeOut1 == closeOut2) {
				while (first && first->opnBrsOutIndx >= openOut1 && first->closeBrsOutIndx <= closeOut1 && first->opnBrsInnIndx <= openIn1
						&& first->closeBrsInnIndx >= closeIn1) {
					first = first->bundleListNext;
				}
				escaped1 = 1;
			} else {
				first = first->bundleListNext;
			}

			if (!escaped2 && openOut1 == openOut2 && openIn1 == openIn2 && closeIn1 == closeIn2 && closeOut1 == closeOut2) {
				while (second && second->opnBrsOutIndx >= openOut2 && second->closeBrsOutIndx <= closeOut2 && second->opnBrsInnIndx <= openIn2
									&& second->closeBrsInnIndx >= closeIn2) {
					second = second->bundleListNext;
				}
			} else {
				second = second->bundleListNext;
			}

			if ((first && !second) || (!first && second)) {
				return 0;
			}
		} else {
			return 0;
		}
	}

	return 1;
}

/* Function to remove duplicates from a unsorted linked list */
int remove_duplicates_LL(ToL* start)
{
  ToL *ptr1, *ptr2, *dup;
  ptr1 = start;
  int count = 0;
  /* Pick elements one by one */
  while(ptr1 != NULL && ptr1->next != NULL)
  {
     ptr2 = ptr1;

     /* Compare the picked element with rest of the elements */
     while(ptr2->next != NULL)
     {
       /* If duplicate then delete it */
       if(nodes_equal(ptr1->branchNode,ptr2->next->branchNode))
       {
          /* sequence of steps is important here */
          dup = ptr2->next;
          ptr2->next = ptr2->next->next;
          free(dup);
          count++;
       }
       else /* This is tricky */
       {
          ptr2 = ptr2->next;
       }
     }
     ptr1 = ptr1->next;
  }

  return count;
}

//*****************************************************************************
// Function : Run Sliding Windows
// Caller   : make_bundle_component_list()
// Purpose  : Calls sliding windows and converts the results into nodes in eden
// Input    : seq, crik
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
void run_sliding_windows(config* seq, global* crik) {
	
	int pid = getpid();	     // Used for creating unique directories of structures so that many instances of Swellix can be run 
				     // at once without interfering with the others' data.

	char idstring[50];
	char bundleDir[256];

	int pidlength = sprintf(idstring, "%d", pid);

	sprintf(idstring,"%s%d%s","/sliding_windows/", pid, "-seq.conf");

	FILE* outfile;                       // Write a .conf file for sliding windows to use
	char outputFilename[1024];
	if (getcwd(outputFilename, sizeof(outputFilename)) != NULL)
		strcat(outputFilename, idstring);
	else {
		fprintf(seq->dispFile, "Error getting current working directory while trying to write file.");
		exit(1);
	}

	outfile = fopen(outputFilename, "w");

	if (outfile == NULL) {
	  fprintf(stderr, "Can't open output file %s!\n",
	          outputFilename);
	  exit(1);
	}

	char mods[seq->strLen+1];
	mods[seq->strLen] = '\0';

	int i;
	for (i = 0; i < seq->strLen; i++)
		mods[i] = '.';
	for (i = 0; i < seq->numS1Pairng; i++) {
		mods[seq->s1Pairng[i]] = 'X';
	}

	int16_t window, tmm, asymmetry;
	window = (seq->minPairngDist * 2) + (seq->minLenOfHlix * 6) - 1;
	if (window > seq->strLen) window = seq->strLen;
	tmm = seq->maxNumMismatch;
	asymmetry = 0;

	fprintf(outfile, "SEQ = \"%s\"\nMODS = \"%s\"\nWINDOW = %d\nLENGTH = %d\nTERMINAL_MISMATCHES = %d\nASYMMETRY = %d\nMISMATCHES = %d\n", seq->ltr, mods, window, seq->minLenOfHlix, tmm, asymmetry, seq->maxNumMismatch);

	fclose(outfile);

	// Use system calls to run sliding windows
	char cwd[1024];
	if (getcwd(cwd, sizeof(cwd)) != NULL) {
		if (chdir(strcat(cwd, "/sliding_windows")) == -1) {
			fprintf(seq->dispFile, "Error changing directory to %s!\n", cwd);
			exit(1);
		}
	} else {
		fprintf(seq->dispFile, "Error getting current working directory while trying to change directory to sliding_windows.");
		exit(1);
	}
	char command[128];
	sprintf(command, "%s %d%s %d %s", "./create_structures.sh", pid, "-seq.conf", pid, _BUNDLE);
	system(command);
	
	sprintf(bundleDir, "%s%s%d", _BUNDLE, "labeled/", pid);
	if (chdir(bundleDir) == -1) {
		fprintf(seq->dispFile, "Error changing to bundleDir: %s!\n", bundleDir);
		exit(1);
	}

	system("ls > list.txt"); // get list of labeled files to read from


	FILE* infile;
	infile = fopen(strcat(bundleDir, "/list.txt"), "r");

	if (infile == NULL) {
		fprintf(stderr, "Can't open file list.txt!\n");
		exit(1);
	}

	bundleDir[strlen(bundleDir)-9-pidlength] = '\0';
//	if (chdir(cwd) == -1) {
//		fprintf(seq->dispFile, "Error changing directory to %s!\n", cwd);
//		exit(1);
//	}
	char filename[50];
	while (fscanf(infile, "%s", filename) != EOF) {
		if (isalpha(filename[0])) // not a labeled file
			continue;

		// add a new node to represent a new bundle
		add_dumi_node(seq, crik, filename);

		// remove duplicates from labeled file
		remove_duplicates(filename);

		// is a labeled file, so open it and convert text to structures
		make_bundles(seq, crik, filename);
	}
	strcat(cwd, "/");
	if (chdir(cwd) == -1) {
		fprintf(seq->dispFile, "Error changing to cwd: %s!\n", cwd);
		exit(1);
	}
	
	sprintf(command, "%s%d%s", "rm ", pid, "-seq.conf");
	system(command);
	fclose(infile);

//	sprintf(command, "rm -r %s%d", bundleDir, pid);
//	system(command);

	// remove duplicates from linked list
	int j;
	for (i = 0; i < seq->strLen-seq->minLenOfHlix; i++) {
		for (j = seq->minLenOfHlix; j < seq->strLen; j++) {
			crik->eden[i][j].numStru -= remove_duplicates_LL(crik->eden[i][j].stemNode);
		}
	}
}

//*****************************************************************************
// Function : Add Dumi Node
// Caller   : run_sliding_windows()
// Purpose  : Adds a new dumi node (bundle) to the tree
// Input    : seq, crik, filename
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
knob* init_b_node(int openOut, int openIn, int closeOut, int closeIn, int16_t** mismatch) {
	knob* newBNode = calloc(1, sizeof(knob)); // new branch node
   	 newBNode->mismatchFlag = mismatch;
	newBNode->opnBrsInnIndx = openIn;
	newBNode->opnBrsOutIndx = openOut;
	newBNode->closeBrsInnIndx = closeIn;
	newBNode->closeBrsOutIndx = closeOut;

	return newBNode;
}

//*****************************************************************************
// Function : Add Stem and Branch Node
// Caller   : set_outer_interval()
// Purpose  : add one more structure node to the linked list, where there probably has been at least one node. Also add the first helix node for that structure node
// Input    : crik
//          : refNode: a component node taken from a sliding_windows file that will be copied into the new branch node
//          : masterLB
//          : masterUB
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
void add_s_b_node(global* crik, knob* refBNode, int16_t masterLB, int16_t masterUB, int8_t mms)
{
  knob* newBNode = calloc(1, sizeof(knob));          // new branch node
  ToL*    newSNode = calloc(1, sizeof(ToL));             // new stem node

  crik->intrvlCntr++;
  newBNode->bundleFlag = 0;
  newBNode->mustPairFlag      = refBNode->mustPairFlag;   // refer directly to the reference component, since it's fix, won't change
  newBNode->mismatchFlag      = refBNode->mismatchFlag;   // refer directly to the reference component, since it's fix, won't change
  newBNode->opnBrsOutIndx     = refBNode->opnBrsOutIndx;
  newBNode->opnBrsInnIndx     = refBNode->opnBrsInnIndx;
  newBNode->closeBrsInnIndx     = refBNode->closeBrsInnIndx;
  newBNode->closeBrsOutIndx     = refBNode->closeBrsOutIndx;
  newBNode->lvlOfRecur        = 100;
  newBNode->bundleCntr        = crik->intrvlCntr;
  newBNode->intrvlCntr        = -1;
  newBNode->hlixBranchngIndx1 = -1;
  newBNode->bundleListNext     = NULL;                     // bundleListNext holds the original link after the jumpTreeNext lost it during the make_jump_tree() operations
  newBNode->jumpTreeNext      = NULL;

  newSNode->branchNode = newBNode;
  newSNode->numHlix  = 1;
  newSNode->numOfMismatches = mms;

  crik->eden[masterLB][masterUB].numStru++;
  newSNode->next = crik->eden[masterLB][masterUB].stemNode;
  crik->eden[masterLB][masterUB].stemNode = newSNode;
  newSNode->bundleCntr = crik->intrvlCntr;
}  // end add_another_s_b_node

//*****************************************************************************
// Function : Add Stem and Branch Node
// Caller   : set_outer_interval()
// Purpose  : add one more structure node to the linked list, where there probably has been at least one node. Also add the first helix node for that structure node
// Input    : crik
//          : refNode: a component node taken from a sliding_windows file that will be copied into the new branch node
//          : masterLB
//          : masterUB
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
void add_b_node(global* crik, knob* refBNode, int16_t masterLB, int16_t masterUB)
{
  knob* newBNode = calloc(1, sizeof(knob));          // new branch node
  newBNode->bundleFlag = 0;
  newBNode->mustPairFlag      = refBNode->mustPairFlag;   // refer directly to the reference component, since it's fix, won't change
  newBNode->mismatchFlag      = refBNode->mismatchFlag;   // refer directly to the reference component, since it's fix, won't change
  newBNode->opnBrsOutIndx     = refBNode->opnBrsOutIndx;
  newBNode->opnBrsInnIndx     = refBNode->opnBrsInnIndx;
  newBNode->closeBrsInnIndx     = refBNode->closeBrsInnIndx;
  newBNode->closeBrsOutIndx     = refBNode->closeBrsOutIndx;
  newBNode->lvlOfRecur        = crik->eden[masterLB][masterUB].stemNode->branchNode->lvlOfRecur++;
  newBNode->bundleCntr        = crik->eden[masterLB][masterUB].stemNode->branchNode->bundleCntr;
  newBNode->intrvlCntr        = -1;
  newBNode->hlixBranchngIndx1 = -1;
  newBNode->bundleListNext     = crik->eden[masterLB][masterUB].stemNode->branchNode; // bundleListNext holds the original link after the jumpTreeNext lost it during the make_jump_tree() operations
  newBNode->jumpTreeNext      = crik->eden[masterLB][masterUB].stemNode->branchNode; // link it to the current branch node

  crik->eden[masterLB][masterUB].stemNode->branchNode = newBNode;
  crik->eden[masterLB][masterUB].stemNode->numHlix++;
}  // end add_another_s_b_node

//*****************************************************************************
// Function : Make Bundles
// Caller   : run_sliding_windows()
// Purpose  : Opens the labeled file and converts the structures in the file
//            into bundles in swellix.
// Input    : seq, crik
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
void make_bundles(config* seq, global* crik, char* filename) {
  FILE* infile;

  char cwd[1024];
  if (getcwd(cwd, sizeof(cwd)) != NULL) {
    infile = fopen(filename, "r");
  } else {
    fprintf(seq->dispFile, "Error getting current working directory while trying to open labeled files.");
    exit(1);
  }

  if (infile == NULL) {
    fprintf(stderr, "Can't open file %s!\n", filename);
    exit(1);
  }


  int begin, end;
  char line[1024];
  int maxBasePairs = 0;
  int8_t repFlag;

  while (fscanf(infile, "%s", line) != EOF) {
    char structure[999999]; // arbitrary length. probably won't be longer than this. if it is, won't complete in reasonable time anyway
    char s1[20], s2[20], s3[20], remainder[50];
    sscanf(line, "%[^,],%[^,],%[^,],%[^,],%s", s1, s2, structure, s3, remainder);
		
    repFlag = FALSE;
    int helices = atoi(s3);
    if(helices > maxBasePairs) { 
      repFlag = TRUE;
      maxBasePairs = helices;
    }
    int helices_added = 0;
    // create nodes based on the helix info in the file
    int i, LB, UB;
    int must_connect = 0;
    int prev_openIn = 0;
    int prev_closeIn = 0;
    int mm_counts = 0;
    int linked = 0;

    for (i = 0; i < helices; i++) {
      int openOut, openIn, closeOut, closeIn;
      int16_t** mismatch; // initialize mismatch flag
      mismatch = calloc(1, seq->maxNumMismatch*sizeof(int16_t*));
      int j;
      for (j = 0; j < seq->maxNumMismatch; j++) {
        mismatch[j] = calloc(1, 2*sizeof(int16_t));
        mismatch[j][0] = -1;
        mismatch[j][1] = -1;
      }

      char a1[20], a2[20], a3[20], a4[20];
      sscanf(remainder, "%[^/]/%[^|]|%[^/]/%[^,],%s", a1, a2, a3, a4, remainder);
      filename_to_indices(filename, &begin, &end);
      openOut = atoi(a1) + begin;
      openIn = atoi(a2) + begin;
      closeIn = atoi(a3) + begin;
      closeOut = atoi(a4) + begin;

	//		fprintf(seq->dispFile, "%s %d - %d | %d - %d\n", structure, openOut, openIn, closeIn, closeOut);
      int window = (seq->minPairngDist * 2) + (seq->minLenOfHlix * 6) - 1;
      if (closeOut >= seq->strLen || openOut < 0 || closeOut-openOut+1 > window) {
	//			fprintf(seq->dispFile, "This one skipped. Out of bounds.\n");
        break; // due to terminal mismatches, length is longer than sequence. skip
      }

      if (must_connect) {
        if (openOut - prev_openIn != 1 || prev_closeIn - closeOut != 1) {
          if (helices_added > 0) {
            crik->eden[LB][UB].stemNode = crik->eden[LB][UB].stemNode->next;
            crik->eden[LB][UB].numStru--;
            if (!crik->eden[LB][UB].stemNode) {
              crik->eden[LB][UB].dumiNode = NULL;
            }
          }
          crik->skippedStru++;
          break;
        }

        must_connect = 0;
      }

      int result = get_mismatches(structure, openOut-begin, openIn-begin, closeOut-begin, closeIn-begin, mismatch);

      if (result) { // lonely pair issue
        if (i == helices-1 || (mismatch[seq->maxNumMismatch-1][0]+begin-openOut+1 < seq->minLenOfHlix && !must_connect)) {
          if (helices_added > 0) {
            crik->eden[LB][UB].stemNode = crik->eden[LB][UB].stemNode->next;
            crik->eden[LB][UB].numStru--;
            if (!crik->eden[LB][UB].stemNode) {
              crik->eden[LB][UB].dumiNode = NULL;
            }
          }
          crik->skippedStru++;
          break;
        } 

        must_connect = 1;
      }

      int skip = 0;
      for (j = 0; j < seq->maxNumMismatch; j++) {
        if (mismatch[j][1] >= 0) {
          mismatch[j][0] += begin; // calibrate the mismatch indices so they match the actual sequence
          mismatch[j][1] += begin;

          if (seq->parenLukUpTable[mismatch[j][0]][mismatch[j][1]]) { // this mismatch can actually pair. skip
			//			fprintf(seq->dispFile, "mismatch[j][0]: %d, mismatch[j][1]: %d\n", mismatch[j][0], mismatch[j][1]);
            skip = 1;
          }
        }
      }

      if (skip) {
		//		fprintf(seq->dispFile, "This one skipped.\n");
        if (helices_added > 0) {
          crik->eden[LB][UB].stemNode = crik->eden[LB][UB].stemNode->next;
          crik->eden[LB][UB].numStru--;
          if (!crik->eden[LB][UB].stemNode) {
            crik->eden[LB][UB].dumiNode = NULL;
          }
        }
        crik->skippedStru++;
        break; // skip because of pairing
      }

      int q;
      for (q = 0; q < seq->maxNumMismatch; q++) {
        if (mismatch[q][0] != -1) {
          mm_counts++;
        }
      }

      if (linked && openOut - prev_openIn == 1 && prev_closeIn - closeOut == 1) {
        crik->eden[LB][UB].stemNode->numOfMismatches = mm_counts;
      } else {
        linked = 0;
      }

      prev_openIn = openIn;
      prev_closeIn = closeIn;

      knob* newBnode = init_b_node(openOut, openIn, closeOut, closeIn, mismatch);

      if (helices_added == 0) {
        linked = 1;

        LB = openOut;
        UB = closeOut;
        if (!crik->eden[openOut][closeOut].dumiNode) { // this dumi node needs to be added
          knob* dumi;
          int16_t hlixLen = (closeOut - openOut + 1 - seq->minPairngDist) >> 1;
          crik->eden[openOut][closeOut].dumiNode = calloc(1, sizeof(knob));

          dumi                    = crik->eden[openOut][closeOut].dumiNode;
          dumi->opnBrsOutIndx     = openOut;
          dumi->opnBrsInnIndx     = openOut + hlixLen;
          dumi->closeBrsInnIndx   = closeOut - hlixLen;
          dumi->closeBrsOutIndx   = closeOut;
          dumi->lvlOfRecur        = -1;
          dumi->hlixBranchngIndx1 = -1;
          dumi->intrvlCntr        = -1;
          dumi->bundleCntr        = -1;
          dumi->outsideIntrvlLB   = -1;
          dumi->outsideIntrvlUB   = -1;
          dumi->specialRstoFlag   =  0;
          dumi->rstoOnQFlag       =  0;
          dumi->bundleFlag        =  1;
          dumi->bundleListNext    = NULL;
          dumi->jumpTreeNext      = NULL;
        }

        int mms = 0;
        int l;
        for (l = 0; l < seq->maxNumMismatch; l++) {
          if (mismatch[l][0] != -1) {
            mms++;
          }
        }

        add_s_b_node(crik, newBnode, openOut, closeOut, mms); // add the new stem and branch nodes
        helices_added++;
      } else {
//				fprintf(seq->dispFile, "Adding new b node.\n");
        add_b_node(crik, newBnode, LB, UB); // add the new branch to the existing stem
        helices_added++;
      }
      free(newBnode);
      free(mismatch);
    }
    if(repFlag) crik->eden[LB][UB].repNode = crik->eden[LB][UB].stemNode->branchNode; 
  }

  if (!crik->eden[begin][end].stemNode) {
    crik->eden[begin][end].dumiNode = NULL;
    crik->eden[begin][end].numStru = 0;
  }

  fclose(infile);
}

//*****************************************************************************
// Function : Reset Dots and Parentheses
// Caller   : prt_individual_partial_stru_on_bundle()
// Purpose  : erase all parentheses and mismatch markings and reset all with dots
// Input    : seq
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
int reset_dot_n_paren(config* seq)
{
  int16_t i;
  for(i = 0 ; i < seq->strLen ; i++) seq->dotNParen[i] = '.';
  return 0;
}  // end reset_dot_n_paren

//*****************************************************************************
// Function : Mark Parentheses
// Caller   : prt_individual_partial_stru_on_bundle()
// Purpose  : mark the parentheses notations ( or ) on the paired-up nucleotides of the structure
// Input    : seq
//          : BNodeCursr
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
int mark_paren(config* seq, knob* BNodeCursr)
{
  int16_t i;
  for(i = BNodeCursr->opnBrsOutIndx ; i <= BNodeCursr->opnBrsInnIndx ; i++) seq->dotNParen[i] = '(';
  for(i = BNodeCursr->closeBrsInnIndx ; i <= BNodeCursr->closeBrsOutIndx ; i++) seq->dotNParen[i] = ')';
  return 0;
}  // end mark_paren

//*****************************************************************************
// Function : Mark Mismatches
// Caller   : prt_individual_partial_stru_on_bundle()
// Purpose  : mark the mismatch notations _ on the mismatch nucleotides of the structure to truthfully reflect the real situation on the structure
// Input    : seq
//          : BNodeCursr
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
int mark_mismatch(config* seq, knob* BNodeCursr)
{
  int16_t m;

  if(seq->maxNumMismatch)
    for(m = 0 ; m < seq->maxNumMismatch ; m++)
      if(BNodeCursr->mismatchFlag[m][1] >= 0){
	seq->dotNParen[BNodeCursr->mismatchFlag[m][0]] = '_';
	seq->dotNParen[BNodeCursr->mismatchFlag[m][1]] = '_';
      }  // end inner if
  return 0;
}  // end mark_mismatch

//*****************************************************************************
// Function : Print Parital Structure
// Caller   : prt_individual_partial_stru_on_bundle()
// Purpose  : print out the id and the dots and parentheses of a partial structure of a bundle
// Input    : seq
//          : BNodeCursr
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
int prt_partial_stru(config* seq, knob* BNodeCursr)
{
  if(!seq || !BNodeCursr) return 0;  // purely useless but simply to shut up the compler warning
  disp(seq,DISP_LV3,"#%-4d - ", BNodeCursr->bundleCntr);
  disp(seq,DISP_LV3,"     %s ", seq->dotNParen);
  return 0;
}  // end prt_partial_stru

//*****************************************************************************
// Function : Print Covariance and V1 Flags
// Caller   : prt_individual_partial_stru_on_bundle()
// Purpose  : mark the covariance and V1 pairing status of that bundle partial structure
// Input    : seq
//          : SNodeCursr
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
int prt_covari_n_V1_flag(config* seq, ToL* SNodeCursr)
{
  int16_t mustPairSize = seq->numCovari + seq->numV1Pairng;
  int16_t i;

  if(mustPairSize && SNodeCursr){
    disp(seq,DISP_LV3,"[ ");
    for(i = 0 ; i < mustPairSize ; i++) disp(seq,DISP_LV3,"%d ", SNodeCursr->mustPairFlag[i]);
    disp(seq,DISP_LV3,"] ");
  }  // end if

  return 0;
}  // end prt_covari_n_V1_flag

//*****************************************************************************
// Function : Print Helix Information
// Caller   : prt_individual_partial_stru_on_bundle()
// Purpose  : list the detailed information about that helix
// Input    : seq
//          : BNodeCursr
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
int prt_hlix_info(config* seq, knob* BNodeCursr)
{
  while(BNodeCursr && seq){
    disp(seq,DISP_LV3,"{%3d %-3d{ %d %d }%3d %-3d} ", BNodeCursr->opnBrsOutIndx, BNodeCursr->opnBrsInnIndx, BNodeCursr->lvlOfRecur, BNodeCursr->bundleCntr, BNodeCursr->closeBrsInnIndx, BNodeCursr->closeBrsOutIndx);
    BNodeCursr = BNodeCursr->jumpTreeNext;
  }  // end while
  return 0;
}  // end prt_hlix_info

//*****************************************************************************
// Function : Print Individual Parital Structure on Bundle
// Caller   : disp_detail_bundle_list(), disp_SSR_ref_list()
// Purpose  : print out one partial structure of a bundle
// Input    : seq
//          : SNodeCursr
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
int prt_individual_partial_stru_on_bundle(config* seq, ToL* SNodeCursr)
{
  knob*  BNodeCursr;

  if(SNodeCursr->numHlix){
    reset_dot_n_paren(seq);
    BNodeCursr = SNodeCursr->branchNode;
    while(BNodeCursr){
      mark_paren(seq, BNodeCursr);
      mark_mismatch(seq, BNodeCursr);
      BNodeCursr = BNodeCursr->jumpTreeNext;
    }  // end while
    BNodeCursr = SNodeCursr->branchNode;
    prt_partial_stru(seq, BNodeCursr);
    prt_covari_n_V1_flag(seq, SNodeCursr);
    prt_hlix_info(seq, BNodeCursr);
    disp(seq,DISP_LV3,"\n");
    SNodeCursr = SNodeCursr->next;
  }    // end if
  return 0;
}  // end prt_individual_partial_stru_on_bundle

//*****************************************************************************
// Function : Display Dummy Node
// Caller   : disp_SSR_ref_list()
// Purpose  : display the node to be actually used by 'make_jump_tree'
// Input    : seq
//          : crik
//          : i2
//          : j2
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
int disp_dumi_node(config* seq, global* crik, int16_t bundleLB, int16_t bundleUB)
{                                                              disp(seq,DISP_ALL,"Enterng 'disp_dumi_node'");
  knob* dumi = crik->eden[bundleLB][bundleUB].dumiNode;
  int16_t i;

  if(!dumi) return 0;
  for(i = 0 ; i < seq->strLen ; i++)                           seq->dotNParen[i] = '.';
  for(i = dumi->opnBrsOutIndx ; i < dumi->closeBrsOutIndx ; i++) seq->dotNParen[i] = '#';
                                                               seq->dotNParen[dumi->closeBrsOutIndx] = '/';
                                                               disp(seq,DISP_LV3,"             %s\n", seq->dotNParen);
  return 0;
}  // end disp_dumi_node

//*****************************************************************************
// Function : Display Similar Structure Reduction Reference List
// Caller   : display_eden()
// Purpose  : format very similar with 'detail_bundle_list', but contains all bundles of and within that interval span
// Input    : seq
//          : crik
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
int disp_bundle_list(config* seq, global* crik)
{
  ToL*     SNodeCursr;
  int16_t  i, j;

  disp(seq,DISP_LV3,"\n\n\n<<<<  SSRR - Similar Structure Reduction Reference List  >>>>\n\n");
  for(i = 0 ; i < seq->strLen ; i++){                            // bundle interval open
    for(j = (i + seq->minPairngDist) ; j < seq->strLen ; j++){  // bundle interval close
      if(crik->eden[i][j].dumiNode){
        disp(seq,DISP_LV3," [%3d-%-3d] - %s\n", i, j, seq->ltr);  // repeat the original sequence as bundle division line
        disp_dumi_node(seq, crik, i, j);
        SNodeCursr = crik->eden[i][j].stemNode;               // individual structure of a bundle
        while(SNodeCursr){
          prt_individual_partial_stru_on_bundle(seq, SNodeCursr);
          SNodeCursr = SNodeCursr->next;
        }
      }
    }
  }

  disp(seq,DISP_LV3,"\n");

  return 0;
}  // end disp_bundle_list

//*****************************************************************************
// Function : Mark Top Notches and Letters
// Caller   : disp_2D_look_up_table()
// Purpose  : mark a '*' every 10 nucleotides as scale and also the letters accross the top banner
// Input    : seq
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
int mark_top_notches_n_ltr(config* seq)
{
  int16_t i;
  for(i = 0 ; i < seq->strLen ; i++)       // || make top notches
    if (i % 10 == 0){                      // ||
      disp(seq,DISP_LV3,"* ");             // ||
    } else {                               // ||
      disp(seq,DISP_LV3,"  ");             // ||
    }  // end if                           // ||
  disp(seq,DISP_LV3,"\n   ");              // \/
  for(i = 0 ; i < seq->strLen ; i++)       // || top line letters
    disp(seq,DISP_LV3,"%c ",seq->ltr[i]);  // \/
  return 0;
}  // end mark_top_notches_n_ltr

//*****************************************************************************
// Function : Mark Top Notches and Letters
// Caller   : disp_2D_look_up_table()
// Purpose  : mark a '*' every 10 nucleotides as scale and also the letters accross the top banner
// Input    : seq
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
int mark_side_notches_n_ltr(config* seq, int16_t row)
{
  if(row % 10 == 0 && seq){
    disp(seq,DISP_LV3,"\n*%c ", seq->ltr[row]);  // mark letter and notch
  } else {
    disp(seq,DISP_LV3,"\n %c ", seq->ltr[row]);  // mark letter only
  }  // end if

  return 0;
}  // end mark_side_notches_n_ltr

//*****************************************************************************
// Function : Print Letters and Numbers inside Table
// Caller   : disp_2D_look_up_table()
// Purpose  : fill in the look-up table by the number of partial structures available at that cell (i, j)
// Input    : seq
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
int prt_ltr_n_num_inside_table(config* seq, global* crik, int16_t lukUpLB, int16_t row, int16_t bundleCnt)
{
  int16_t bundleUB = (seq->minPairngDist << 1) + (seq->minLenOfHlix * 6);
  int16_t lukUpUB;
  int16_t j;

  if(lukUpLB < seq->strLen){
    disp(seq,DISP_LV3,"%c   ", seq->ltr[row]);
    lukUpUB = (row + bundleUB < seq->strLen) ? (row + bundleUB) : seq->strLen;
    for(j = lukUpLB ; j < lukUpUB ; j++){
      if (crik->eden[row][j].numStru){
	disp(seq,DISP_LV3,"%-2d", crik->eden[row][j].numStru);
      } else {
	disp(seq,DISP_LV3,"  ");
      }  // end if
      bundleCnt = crik->eden[row][j].numStru ? (bundleCnt + 1) : bundleCnt;
    }  // end inner for
  }    // end if
  return bundleCnt;
}  // end prt_ltr_n_num_inside_table

//*****************************************************************************
// Function : Display 2D Look-up Table
// Caller   : display_eden()
// Purpose  : display the (i, j) = (interval LB, interval UB) in the bird view perspective
//          : number in each cell indicates the number of structures included in that interval (i, j)
// Input    : seq
//          : crik
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
int disp_2D_look_up_table(config* seq, global* crik)
{                                         disp(seq,DISP_ALL,"Enterng 'disp_2D_look_up_table'");
  int16_t bundleLB = seq->minPairngDist;
  int16_t lukUpLB;
  int16_t bundleCnt = 0;                  // total number of bundles/intervals
  int16_t i, j;

  disp(seq,DISP_LV3,"\n\n   ");
  mark_top_notches_n_ltr(seq);
  for(i = 0 ; i < seq->strLen ; i++){
    mark_side_notches_n_ltr(seq, i);
    lukUpLB = i + bundleLB;
    for(j = 0 ; j < (lukUpLB - 2) ; j++)  // || leading spaces of each row (there's no code designed for trailing spaces)
      disp(seq,DISP_LV3,"  ");            // \/
    bundleCnt = prt_ltr_n_num_inside_table(seq, crik, lukUpLB, i, bundleCnt);
  }      // end outer for
  disp(seq,DISP_LV3,"\n\n");

  return bundleCnt;
}  // end disp_2D_look_up_table

//*****************************************************************************
// Function : Display Eden
// Caller   : make_bundle_list()
// Purpose  : display the complete picture of eden
// Input    : seq
//          : crik
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
void display_eden(config* seq, global* crik)
{
  disp(seq,DISP_ALL,"Entering 'display eden'");
  int16_t bundleCnt = disp_2D_look_up_table(seq, crik);
  crik->numBundles = bundleCnt;
  disp(seq,DISP_LV3,"Total Bundles = %d. Structures in Bundle = %ld\n", bundleCnt, crik->intrvlCntr);
  disp_bundle_list(seq, crik);
}  // end display_eden

//*****************************************************************************
// Function : Reverse Order
// Caller   : make_bundle_list()
// Purpose  : rearrange/reverse the order of each structure such that in each cell, the order of linked list of the partial structures will be the same as the jump_tree crik->hlixInStru and the structures on the cell will have the ascending order
// Input    : crik
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
int rev_order(config* seq, global* crik)
{
  int16_t i, j;
  ToL*    SCursr;
  ToL*    newStem;
  ToL*    tempS;

  for(i = 0 ; i < seq->strLen ; i++){
    for(j = 0 ; j < seq->strLen ; j++){
      newStem = NULL;
      tempS   = NULL;
      if(crik->eden[i][j].numStru){
	SCursr = crik->eden[i][j].stemNode;
	while(SCursr){
	  tempS = SCursr->next;
	  SCursr->next = newStem;
	  newStem = SCursr;
	  SCursr = tempS;
	}  // end while
	crik->eden[i][j].stemNode = newStem;
      }    // end if
    }  // end inner for
  }    // end outer for
  return 0;
}  // end rev_order
