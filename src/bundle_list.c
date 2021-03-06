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
#include "subopt.h"

#ifdef _MPI
#include "mpi.h"
#endif

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
int get_mismatches(char* structure, int openOut, int openIn, int closeOut, int closeIn, int16_t** mismatch)
{
  int i;
  // Find any mismatches
  int j = 0;
  int toReturn = 0;
  for (i = openOut; i <= openIn; i++) {
    if (i < 0) {
      mismatch[j][0] = i;
      j++;
    } else if (structure[i] == '.') {
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
void remove_outer_pair(char* structure, int start, int end)
{
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
void add_dumi_node(config* seq, global* crik, LabeledStructures* lab)
{
	int begin, end;
	filename_to_indices(lab->title, &begin, &end);
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

void filename_to_indices(char* filename, int* begin, int* end)
{
	char pos[50], len[50], suffix[4];
	sscanf(filename, "%[^x]x%[^.].%s", pos, len, suffix);
	*end = atoi(pos);
	*begin = atoi(len);
	*begin = *end - *begin + 1;
}

/* Function to test whether two nodes are equal */
int nodes_equal(knob* first, knob* second)
{

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
    }
    else if (firstCursr->opnBrsOutIndx - openIn1 == 1 && closeIn1 - firstCursr->closeBrsOutIndx == 1) {
      openIn1 = firstCursr->opnBrsInnIndx;
      closeIn1 = firstCursr->closeBrsInnIndx;
      formed = 1;
    }
    else if (!formed) {
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
    }
    else if (!formed) {
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
    if ((first->opnBrsInnIndx == second->opnBrsInnIndx && 
         first->opnBrsOutIndx == second->opnBrsOutIndx && 
         first->closeBrsInnIndx == second->closeBrsInnIndx && 
         first->closeBrsOutIndx == second->closeBrsOutIndx) || 
        (openOut1 == openOut2 && 
         openIn1 == openIn2 && 
         closeIn1 == closeIn2 && 
         closeOut1 == closeOut2)) {

      if (!escaped1 && 
          openOut1 == openOut2 && 
          openIn1 == openIn2 && 
          closeIn1 == closeIn2 && 
          closeOut1 == closeOut2) {

        while (first && 
               first->opnBrsOutIndx >= openOut1 && 
               first->closeBrsOutIndx <= closeOut1 && 
               first->opnBrsInnIndx <= openIn1 && 
               first->closeBrsInnIndx >= closeIn1) {
          first = first->bundleListNext;
        }
        escaped1 = 1;
      } else {
        first = first->bundleListNext;
      }

      if (!escaped2 && 
          openOut1 == openOut2 && 
          openIn1 == openIn2 && 
          closeIn1 == closeIn2 && 
          closeOut1 == closeOut2) {

        while (second && 
               second->opnBrsOutIndx >= openOut2 && 
               second->closeBrsOutIndx <= closeOut2 && 
               second->opnBrsInnIndx <= openIn2 && 
               second->closeBrsInnIndx >= closeIn2) {
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
void run_sliding_windows(config* seq, global* crik)
{
  extern int rank, wsize;
  int istart, iend, span, rem;

  span = (seq->strLen)/(wsize);
  rem = seq->strLen % wsize;
  istart = rank*span;
  iend = (rank == wsize-1) ? istart+span+rem : istart+span-1;

  // Get ready to run sliding windows
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

  LabeledStructures* labs = calloc(seq->strLen*100, sizeof(LabeledStructures*));
  int labsSize = 0;
  int labsMax = seq->strLen*100;

  int index;
  char* subSeq = calloc(window+1, sizeof(char));
  char* subMod = calloc(window+1, sizeof(char));
  int start, stroffset, substrLen;
  for(index = istart; index < iend+1; index++) {
    start = index-window-1;
    stroffset = start+1 > 0 ? start+1 : 0;
    substrLen = index-stroffset > window ? window : index-stroffset;
    strncpy(subSeq, seq->ltr+stroffset, substrLen);
    strncpy(subMod, mods+stroffset, substrLen);
    slide_those_windows(subSeq, subMod, start, mods, window, tmm, asymmetry, seq, labs, &labsSize, &labsMax);
  }

#ifdef _MPI
  // Create MPI datatype to pass LabeledStructures structs across PEs
  MPI_Datatype mpi_labeled_structures;
  MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_UB}; //due to the dynamically allocated arrays in the structs
                                                      //MPI needs to take the integer fields plus however much more
                                                      //memory the struct takes up.
  int blocklens[3] = {1, 1, 1};
  LabeledStructures setup;
  MPI_Aint displacements[3];
  displacements[0] = (MPI_Aint)&setup.titlesize - (MPI_Aint)&setup;
  displacements[1] = (MPI_Aint)&setup.buffsize - (MPI_Aint)&setup;
  displacements[2] = (MPI_Aint)sizeof(setup);
  MPI_Type_create_struct(3, blocklens, displacements, types, &mpi_labeled_structures);
  MPI_Type_commit(&mpi_labeled_structures);

  /*
    The way this parallelization works is that all PEs do their own bundling on a pseudo-evenly distributed
    chunk of the input RNA sequence. After each PE finishes its own portion of the sliding windows step, said PE
    needs to transmit its data to all other PEs so they can construct their own appropriate data locally. For this,
    MPI_Allgatherv will be used to handle the (probable) case that not all PEs finish with the same number of 
    LabeledStructures instances populated.
  */

  // We need to sum up the total number of LabeledStructures instances that got filled with data.
  int totalLabs = 0;
  MPI_Allreduce(&labsSize, &totalLabs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // Now, gather the local counts of LabeledStructures instances from each PE an store them 
  // in an array for Allgatherv.
  int counts[wsize];
  MPI_Allgather(&labsSize, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);

  // Construct an array to store the displacement of each PEs data in the final complete array as per Allgatherv.
  int* displs = calloc(wsize, sizeof(int));
  for(index = 1; index < wsize; index++)
    displs[index] = displs[index-1] + counts[index-1];

  // Allocate space for the complete array of LabeledStructures instances, using the sum totalLabs from above to
  // allocate just enough memory as needed.
  LabeledStructures* completeLabs = (LabeledStructures*)calloc(totalLabs, sizeof(LabeledStructures));
  
  // Now, from each PEs local labs array, gather the respective elements and construct a complete array in each
  // PE according to our count and displacement arrays from above.
  MPI_Allgatherv(labs, labsSize, mpi_labeled_structures, 
                 completeLabs, counts, displs, mpi_labeled_structures, 
                 MPI_COMM_WORLD);

  // Since the structs we sent contain dynamically allocated character arrays, we need to now allocate memory for
  // these cstrings in each PE's memory. In addition, on each PE we need to copy over the actual cstring values to
  // the appropriate places in the completeLabs array.
  for(index = 0; index < totalLabs; index++) {
    completeLabs[index].title = (char*)calloc(completeLabs[index].titlesize, sizeof(char));
    completeLabs[index].structures = (char*)calloc(completeLabs[index].buffsize, sizeof(char));
    if(index >= displs[rank] && index < displs[rank]+counts[rank]) {
      strcpy(completeLabs[index].title, labs[index-displs[rank]].title);
      if(completeLabs[index].buffsize != labs[index-displs[rank]].buffsize) {
        completeLabs[index].buffsize = labs[index-displs[rank]].buffsize;
        completeLabs[index].structures = (char*)realloc(completeLabs[index].structures, completeLabs[index].buffsize);
      }
      strcpy(completeLabs[index].structures, labs[index-displs[rank]].structures);
    }
  }

  // Finally, to complete the transmission of the data between all of the PEs, each PE needs to send its own local
  // information to the other PEs in their completeLabs struct arrays.
  LabeledStructures* ptr;
  for(index = 0; index < totalLabs; index++) {
    ptr = &completeLabs[index];
    if(index < displs[rank] || index >= displs[rank] + counts[rank]) {
      initLabeledStructures(ptr);
    }
    int root, i = 0;
    while(i < wsize) {
      // Use the displacements and counts to figure out which PE should be the source for this index in the
      // completeLabs array.
      if(index >= displs[i] && index < displs[i]+counts[i]) root = i;
      i++;
    }

    MPI_Bcast(ptr->title, ptr->titlesize, MPI_CHAR, root, MPI_COMM_WORLD);
    MPI_Bcast(&ptr->buffsize, 1, MPI_INT, root, MPI_COMM_WORLD);

    if(rank != root && ptr->buffsize > 4096) {
        ptr->structures = (char*)realloc(ptr->structures, ptr->buffsize);
    }
    MPI_Bcast(ptr->structures, ptr->buffsize, MPI_CHAR, root, MPI_COMM_WORLD);
  }

  // Now that the information has been distributed across all PEs, let each PE construct its bundles locally using
  // the full array of LabeledStructures instances.
  for(index = 0; index < totalLabs; index++) {
    add_dumi_node(seq, crik, &completeLabs[index]);
    make_bundles(seq, crik, &completeLabs[index]);
    free(completeLabs[index].title);
    free(completeLabs[index].structures);

    free(labs[index].title);
    free(labs[index].structures);
  }
  
  free(completeLabs);
  MPI_Type_free(&mpi_labeled_structures);

#else

  // The default (serial) behavior: the process should have populated labs with all possible LabeledStructures
  // so just make the bundles as above with the MPI-accumulated completeLabs array.
  for(index = 0; index < labsSize; index++) {
    add_dumi_node(seq, crik, &labs[index]);
    make_bundles(seq, crik, &labs[index]);
    free(labs[index].title);
    free(labs[index].structures);
  }

#endif

  free(labs);
  free(subSeq);
  free(subMod);

  // remove duplicates from linked list
  int j;
  for (i = 0; i < seq->strLen-seq->minLenOfHlix; i++) {
    for (j = seq->minLenOfHlix; j < seq->strLen; j++) {
      crik->eden[i][j].numStru -= remove_duplicates_LL(crik->eden[i][j].stemNode);
    }
  }
}

void initLabeledStructures(LabeledStructures *lab)
{
  lab->titlesize = 64;
  lab->title = (char*)calloc(lab->titlesize, sizeof(char));
  lab->buffsize = 4096;  // default buffsize to 4kB arbitrarily
  lab->structures = (char*)calloc(lab->buffsize, sizeof(char)); 
}

void resetLabeledStructures(LabeledStructures* lab, char* newtitle) 
{
  lab->title = memset(lab->title, '\0', lab->titlesize);
  lab->structures = memset(lab->structures, '\0', lab->buffsize);
  if(newtitle != NULL) strcat(lab->title, newtitle);
}

void freeLabeledStructures(LabeledStructures **lab)
{
  free((*lab)->title);
  free((*lab)->structures);

  free(*lab);
  *lab = NULL;
}

//*****************************************************************************
// Function : Add Dumi Node
// Caller   : run_sliding_windows()
// Purpose  : Adds a new dumi node (bundle) to the tree
// Input    : seq, crik, filename
// Return   : none
// Display  : error message, if necessary
//*****************************************************************************
knob* init_b_node(int openOut, int openIn, int closeOut, int closeIn, int16_t** mismatch)
{
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
void make_bundles(config* seq, global* crik, LabeledStructures* lab)
{

  int begin, end;
  char* tok;
  int maxBasePairs = 0;
  int8_t repFlag;
  tok = strtok(lab->structures, "\n");

  while (tok != NULL) {
    char structure[999999]; // arbitrary length. probably won't be longer than this. if it is, won't complete in reasonable time anyway
    char s1[20], s2[20], s3[20], remainder[50];
    sscanf(tok, "%[^,],%[^,],%[^,],%[^,],%s", s1, s2, structure, s3, remainder);
		
    repFlag = FALSE;
    int helices = atoi(s3);
    int numPairs = 0, index = 0;
    while(structure[index] != '\0') {
      if(structure[index++] == '(') numPairs++;
    }
    if(numPairs > maxBasePairs) {
      repFlag = TRUE;
      maxBasePairs = numPairs;
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
      filename_to_indices(lab->title, &begin, &end);
      openOut = atoi(a1) + begin;
      openIn = atoi(a2) + begin;
      closeIn = atoi(a3) + begin;
      closeOut = atoi(a4) + begin;

      int window = (seq->minPairngDist * 2) + (seq->minLenOfHlix * 6) - 1;
      if (closeOut >= seq->strLen || openOut < 0 || closeOut-openOut+1 > window) {
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

          skip = seq->parenLukUpTable[mismatch[j][0]][mismatch[j][1]]; // this mismatch can actually pair. skip
        }
      }

      if (skip) {
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
        add_b_node(crik, newBnode, LB, UB); // add the new branch to the existing stem
        helices_added++;
      }
      free(newBnode);
      free(mismatch);
    }
    if(repFlag) crik->eden[LB][UB].repNode = crik->eden[LB][UB].stemNode->branchNode; 

    tok = strtok(NULL, "\n");
  }

  if (!crik->eden[begin][end].stemNode) {
    crik->eden[begin][end].dumiNode = NULL;
    crik->eden[begin][end].numStru = 0;
  }
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
    for(m = 0 ; m < seq->maxNumMismatch ; m++) {
      if(BNodeCursr->mismatchFlag[m][1] >= 0){
	seq->dotNParen[BNodeCursr->mismatchFlag[m][0]] = '_';
	seq->dotNParen[BNodeCursr->mismatchFlag[m][1]] = '_';
      }  // end inner if
    }
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
