/*****************************************************************************************************************************
 *****************************************************************************************************************************
 *****************************************************************************************************************************
 Program : Swellix (Make Jump Tree) - c file
 Purpose : Ryan's RNA Sequence Folding Algorithm (modified Crumple)
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
#include "close_up.h"
#include "jump_tree.h"
#include "unit_tests.h"
#include "statistics.h"
#include "mpi_jump_tree.h"

char* CALLER_FLAG_TO_STRING[] = { "FIRST_SESSION", "SECND_SESSION",
		"WHILE_SESSION", "RSTO_SESSION" };
char* NODE_TYPE_TO_STRING[] = { "HEAD", "CURR", "TAIL" };

//Counts the number of calls made to jump_stage_1(). Used to stop recursion after a set number of calls
//when running the GNU Debugger.
int numOfCalls = 0;
int numCalls = 0;
int jst2Count = 0;
int exitCount = 0;
int rstoCount = 0;

//*****************************************************************************
// Function : Assign New Helices Links
// Caller   : take_cmpnt_list_normal_path()
// Purpose  : connect crik->hlixInStru to another components
// Input    : seq  
//          : crik 
//          : todd
// Return   : none
// Display  : none
//*****************************************************************************
int assign_new_hlix_link(global* crik, int16_t hlixBranchngIndx1) {
  knob* hlixCursr = crik->hlixInStru;

  crik->intrvlCntr++;
  hlixCursr->hlixBranchngIndx1 = hlixBranchngIndx1;
  hlixCursr->intrvlCntr = crik->intrvlCntr;
  hlixCursr->lvlOfRecur = 0;
  crik->numHlix = 1;
  return 0;
}  // end assign_new_hlix_link

//*****************************************************************************
// Function : Clear Old Helices Links
// Caller   : make_jump_tree()
// Purpose  : clear all the helix links of ->jumpTreeNext, to prevent infinite loops
// Input    : seq  
//          : crik 
// Return   : none
// Display  : none
//*****************************************************************************
int clear_old_hlix_link(global* crik) {
  knob* hlixCursr = crik->hlixInStru;
  knob* temp;

  while (hlixCursr) {
    temp = hlixCursr->jumpTreeNext;
    hlixCursr->lvlOfRecur = -1;
    hlixCursr->jumpTreeNext = NULL;
    crik->hlixInStru = temp;
    hlixCursr = temp;
  }  // end while
  hlixCursr = crik->hlixInStru;

  return 0;
}  // end clear_old_hlix_link

//*****************************************************************************
// Function : Does Satisfy Constraints?
// Caller   : take_cmpnt_list_normal_path()
// Purpose  : update the constraint status and check if the the any of the constraint requirements are still satisfied
// Input    : seq  
//          : crik 
//          : todd
// Return   : none
// Display  : none
//*****************************************************************************
int satisfy_constraint(config* seq, global* crik, local* todd) {
  int16_t i;
  int16_t flagUB = seq->numCovari + seq->numV1Pairng;
  int8_t failFlag = 0;

  for (i = 0; i < flagUB; i++) {
    if (todd->cmpntLLCursr->mustPairFlag[i]) { // it's the first helix of the whole structure, therefore whatever covariances this helix got, will be the covari the whole structure has no more, no less than that
      crik->struMustPairFlag[i] = 1;
    } else {
      crik->struMustPairFlag[i] = 0;
      failFlag = 1;
    }  // end if
  }    // end for
  return !failFlag; // if any covariance pair is missing, the current structure won't count
}  // end satisfy_constraint

//*****************************************************************************
// Function : Count Structures
// Caller   : eval_prog_n_prt_stru and jump_stage_2_fit_hlix
// Purpose  : counts the total number of unbundled structures in a structure containing bundles
// Input    : crik
// Return   : none
// Display  : none
//*****************************************************************************
int countAndPrintStructures(config* seq, global* crik) {
  knob* hlixCursr = crik->hlixInStru;
  int num = 1;
  while (hlixCursr) {
    if (hlixCursr->bundleFlag) {
      if (hlixCursr->jumpTreeNext && !hlixCursr->jumpTreeNext->bundleFlag && hlixCursr->opnBrsOutIndx - hlixCursr->jumpTreeNext->opnBrsInnIndx == 1
          && hlixCursr->jumpTreeNext->closeBrsInnIndx - hlixCursr->closeBrsOutIndx == 1) {
        int structures = crik->eden[hlixCursr->opnBrsOutIndx][hlixCursr->closeBrsOutIndx].numStru;
        int mms = 0;
        int i;
        for (i = 0; i < seq->maxNumMismatch; i++) {
          if (hlixCursr->jumpTreeNext->mismatchFlag[i][0] != -1) {
            mms++;
          }
        }

        ToL* tolCursr = crik->eden[hlixCursr->opnBrsOutIndx][hlixCursr->closeBrsOutIndx].stemNode;
        while (tolCursr) {
          if (mms + tolCursr->numOfMismatches > seq->maxNumMismatch) {
            structures--;
          }
          tolCursr = tolCursr->next;
        }

        num *= structures;
      } else {
        num *= crik->eden[hlixCursr->opnBrsOutIndx][hlixCursr->closeBrsOutIndx].numStru;
      }
    }
    hlixCursr = hlixCursr->jumpTreeNext;
  }

  crik->numUnbundledStru += num;
//	printf("#structs: %d ", num);

  if(seq->unbundle) {	// UNBUNDLING
    ToL* results;
    ToL* first;
    ToL* next;
    hlixCursr = crik->hlixInStru;
    if (!hlixCursr->bundleFlag) {
      first = calloc(1, sizeof(ToL));
      first->branchNode = copy_b_node_wo_links(hlixCursr);
      first->next = NULL;
    } else {
      first = crik->eden[hlixCursr->opnBrsOutIndx][hlixCursr->closeBrsOutIndx].stemNode;
    }
	
    if (hlixCursr->jumpTreeNext) {
      if (hlixCursr->jumpTreeNext->bundleFlag) {
        next = crik->eden[hlixCursr->jumpTreeNext->opnBrsOutIndx][hlixCursr->jumpTreeNext->closeBrsOutIndx].stemNode;
      } else {
        next = calloc(1, sizeof(ToL));
        next->branchNode = copy_b_node_wo_links(hlixCursr->jumpTreeNext);
        next->next = NULL;
      }
    } else {
      next = calloc(1, sizeof(ToL));
      next->branchNode = NULL;
      next->next = NULL;
    }
	
    results = combine(first, next);
    if(!hlixCursr->bundleFlag) release_tol(first);
    if(!hlixCursr->jumpTreeNext || (hlixCursr->jumpTreeNext && !hlixCursr->jumpTreeNext->bundleFlag)) release_tol(next);

    hlixCursr = hlixCursr->jumpTreeNext;
    if (hlixCursr) {
      hlixCursr = hlixCursr->jumpTreeNext;
	
      while (hlixCursr) {
        if (!hlixCursr->bundleFlag) {
          next = calloc(1, sizeof(ToL));
          next->branchNode = copy_b_node_wo_links(hlixCursr);
          next->next = NULL;
        } else {
          next = crik->eden[hlixCursr->opnBrsOutIndx][hlixCursr->closeBrsOutIndx].stemNode;
        }
        ToL* temp = results;
        results = combine(results, next);
        release_tol(temp);
        if(!hlixCursr->bundleFlag) release_tol(next);
	
        hlixCursr = hlixCursr->jumpTreeNext;
      }
    }

    hlixCursr = crik->hlixInStru;
    ToL* resultsCursr = results;	
	
    // print out all structures in results
    while (results) {
      crik->hlixInStru = results->branchNode;
      disp_stru(seq, crik, 11);                  // for report  (official use)
      update_max_energy(seq);
      update_max_distance(seq);
      numCalls++;
      dispStru(seq,crik,11);                     // for display (debug use)
      results = results->next;
    }
	
    // return hlixInStru back to original helix
    crik->hlixInStru = hlixCursr;
    release_tol(resultsCursr);
  } else {

    disp_stru(seq, crik, 11);                  // for report  (official use)
    dispStru(seq,crik,11);                     // for display (debug use)
  }

  return num;
}

knob* copy_b_node_wo_links(knob* node) {
  knob* newNode = calloc(1, sizeof(knob));
  if (!node)
    newNode = NULL;
  else {
    newNode->bundleCntr = node->bundleCntr;
    newNode->bundleFlag = node->bundleFlag;
    newNode->bundleListNext = node->bundleListNext;
    newNode->closeBrsInnIndx = node->closeBrsInnIndx;
    newNode->closeBrsOutIndx = node->closeBrsOutIndx;
    newNode->opnBrsInnIndx = node->opnBrsInnIndx;
    newNode->opnBrsOutIndx = node->opnBrsOutIndx;
    newNode->cmpntListNext = node->cmpntListNext;
    newNode->hlixBranchngIndx1 = node->hlixBranchngIndx1;
    newNode->intrvlInsFormdFlag = node->intrvlInsFormdFlag;
    newNode->mismatchFlag = node->mismatchFlag;
    newNode->mustPairFlag = node->mustPairFlag;
    newNode->jumpTreeNext = NULL;
  }

  return newNode;
}

knob* copy_b_node(knob* node) {
  knob* newNode = NULL;
  if (!node)
    newNode = NULL;
  else {
    newNode = calloc(1, sizeof(knob));
    newNode->bundleCntr = node->bundleCntr;
    newNode->bundleFlag = node->bundleFlag;
    newNode->bundleListNext = node->bundleListNext;
    newNode->closeBrsInnIndx = node->closeBrsInnIndx;
    newNode->closeBrsOutIndx = node->closeBrsOutIndx;
    newNode->opnBrsInnIndx = node->opnBrsInnIndx;
    newNode->opnBrsOutIndx = node->opnBrsOutIndx;
    newNode->cmpntListNext = node->cmpntListNext;
    newNode->hlixBranchngIndx1 = node->hlixBranchngIndx1;
    newNode->intrvlInsFormdFlag = node->intrvlInsFormdFlag;
    newNode->mismatchFlag = node->mismatchFlag;
    newNode->mustPairFlag = node->mustPairFlag;
    newNode->bundleListNext = node->bundleListNext;
    newNode->jumpTreeNext = NULL;
    // make a deep copy of next node
    if (node->jumpTreeNext) {
      knob* newNext = copy_b_node(node->jumpTreeNext);
      newNode->jumpTreeNext = newNext;
    } else {
      newNode->jumpTreeNext = NULL;
    }
  }
	
  return newNode;
}

ToL* copy_tol(ToL* source) {	//deep copy a ToL object
  ToL* newToL = NULL;

  if(source) {
    newToL = calloc(1, sizeof(ToL));
    newToL->branchNode = NULL;
    if(source->branchNode)		//deep copy the branch nodes, so that the newly created structure can be comprehensively freed
      newToL->branchNode = copy_b_node(source->branchNode);
    newToL->mustPairFlag = NULL;
    newToL->numHlix = source->numHlix;
    newToL->bundleCntr = source->bundleCntr;
    newToL->numOfMismatches = source->numOfMismatches;
    newToL->next = copy_tol(source->next);		//recursively copy the next ToL object in the linked list
  }

  return newToL;	//the new copy of source
}

void release_knob(knob* obj) {	//free the memory used by a knob object. 
  if(obj) {
    if(obj->jumpTreeNext) {
      release_knob(obj->jumpTreeNext);
      obj->jumpTreeNext = NULL;
    }
    free(obj);
    obj = NULL;
  }
}



void release_tol(ToL* obj) {
  if(obj) {
    release_knob(obj->branchNode);
    release_tol(obj->next);
    free(obj);
    obj = NULL;
  }
}

ToL* combine(ToL* first, ToL* second) {
  ToL* combined = NULL;

  while (first) {
    knob* newB = first->branchNode;
    ToL* secondCursr = second;
    while (secondCursr) {
      ToL* newS = calloc(1, sizeof(ToL));

      knob* copyB = copy_b_node(newB);
      knob* it = copyB;
      while (it->jumpTreeNext) {
        it = it->jumpTreeNext;
      }
			
      it->jumpTreeNext = copy_b_node(secondCursr->branchNode);

      newS->branchNode = copyB;
      newS->next = combined;
      combined = newS;

      secondCursr = secondCursr->next;
    }

    first = first->next;
  }

  return combined;
}

//*****************************************************************************
// Function : Evaluate Progress and Print Structure
// Caller   : take_bundle_list_shortcut()
// Purpose  : check if the current structure is qualified to increment structure count and print out
// Input    : seq
//          : crik
// Return   : none
// Display  : none
//*****************************************************************************
int eval_prog_n_prt_stru(config* seq, global* crik) {
  if ((crik->numHlix >= seq->minNumHlix) && // This segment of code is visited to see if 'one_helix-one_hairpin structure' is qualified to be a complete structure
      (crik->numHP >= seq->minNumHP)) {
    crik->numStru++;

    if (seq->bundle) { // count unbundled structures
      countAndPrintStructures(seq, crik);
      crik->linkedmms = 0;
    } else {
      disp_stru(seq, crik, 11);              // for report  (official use)
      dispStru(seq,crik,11);                    // for display (debug use)
      crik->linkedmms = 0;
    }
  }  // end if
  return 0;
}  // end eval_prog_n_prt_stru

//*****************************************************************************
// Function : Exit Current Recursion
// Caller   : make_jump_tree() & jump_stage_2_fit_hlix()
// Purpose  : clean up all memories inside todd
// Input    : seq
//          : crik
//          : todd
// Return   : none
// Display  : none
//*****************************************************************************
int exit_curr_recur(config* seq, global* crik, local* todd) {
  disp(seq,DISP_ALL,"Enterng 'exit_curr_recur'\n");
  knob* temp;
  knob* hlixCursr = crik->hlixInStru;
  int16_t flagUB = seq->numCovari + seq->numV1Pairng;
  int16_t i;
  exitCount++;

  /* Bug fix note from Nathan Sloat - July 16, 2015
   * 
   * I found an edge case here involving a new recursion level being entered. Most of the time, when a new interval is chosen, there will be at least one component which
   * can be fit into it. This is just the general case for some sequence with no imposed constraints. However, the situation that I came across was in the 23rd 
   * hlixBranchngIndx and it involved an interval in which no components in the component list could fit. This meant that no new helices were tacked on to the beginning 
   * of the crik->hlixInStru linked list. This wouldn't have been a problem, if it weren't for the fact that exit_curr_recur still executed normally, and specifically, 
   * this next conditional statement evaluated as true. The only logic check here was that hlixCursr != NULL. This meant that if there was any structure currently being 
   * stored in the crik->hlixInStru list, this next block would execute. Inside this block, there is a piece of code to remove the most recently added helix -- assuming 
   * that one is added at every recursion level. Herein lies my problem, as in this case a new level of recursion was entered with no new helices added to the structure. 
   * So, the helix from one recursion level back is removed which in turn causes a multitude of computational problems later on. These include memory leaks and circularly 
   * linked lists. I added a single logic check to make sure that the front of the linked list is actually a new helix from the current recursion level's todd object, 
   * and the problem was fixed.
   */
  if (hlixCursr && hlixCursr->intrvlCntr == todd->intrvlCntr) { // clear the helix of this present recursion level of before moving one level back (up)
    for (i = 0; i < flagUB; i++) {                        // || REMOVE HELIX
      if (hlixCursr->mustPairFlag[i]) {
        disp(seq,DISP_ALL,"hlix removed: {%d-%d{ %d }%d-%d}\n", hlixCursr->opnBrsOutIndx, hlixCursr->opnBrsInnIndx, hlixCursr->lvlOfRecur, hlixCursr->closeBrsInnIndx, hlixCursr->closeBrsOutIndx);
        crik->struMustPairFlag[i] = 0;                    // ||
      }  // end inner if                                  // ||
    }    // end for                                       // ||
                                                          // ||
    temp = hlixCursr->jumpTreeNext;                       // || remove linked list node
    hlixCursr->specialRstoFlag = 0;                       // || ||
    hlixCursr->intrvlCntr = -1;                           // || ||
    hlixCursr->jumpTreeNext = NULL;                       // || ||
    crik->hlixInStru = temp;                              // || \/
    crik->numHlix--;                                      // \/

    if (todd->RSTO
        && !is_there_concern_on_ring_formation_of_linked_list(seq, crik,
        todd)) {
      disp(seq,DISP_ALL,"            special On Q %d-%d\n",todd->RSTO->opnBrsInnIndx, todd->RSTO->closeBrsInnIndx);
      crik->interval = todd->RSTO;
      rstoCount++;
      disp(seq,DISP_ALL,"Interval restore todd on queue!\n");
    }  // end inner if
  }    // end outer if

  disp(seq,DISP_ALL,"     ins  intrvl: [%2d,%-2d]\n", todd->intrvlIns->opnBrsInnIndx, todd->intrvlIns->closeBrsInnIndx);
  free(todd->intrvlIns);
	
//	if(!hlixCursr && todd->intrvlBeh != todd->RSTO) { free(todd->RSTO); todd->RSTO = NULL; }

  disp(seq,DISP_ALL,"     beh  intrvl: [%2d,%-2d]\n", todd->intrvlBeh->opnBrsInnIndx, todd->intrvlBeh->closeBrsInnIndx);
  free(todd->intrvlBeh);
  disp(seq,DISP_ALL,"Exiting 'jump_stage_2', recur lv: %d->%d\n", todd->lvlOfRecur, (todd->lvlOfRecur - 1));
  free(todd);
  crik->lvlOfRecur--;

  return 0;
}  // end exit_curr_recur

//*****************************************************************************
// Function : Initialize Todd
// Caller   : make_jump_tree()
// Purpose  : assign beginning values to todd
// Input    : crik
// Return   : todd
// Display  : none
//*****************************************************************************
local* init_todd(global* crik) {
  local* todd = calloc(1, sizeof(local));
  knob* insCursr;
  knob* behCursr;
  knob* intrvlCursr = crik->interval;

  todd->intrvlIns = calloc(1, sizeof(knob));
  insCursr = todd->intrvlIns;
  insCursr->opnBrsInnIndx = -1; // init as -1 instead of 0, since 0 presents the location 0 on the sequence
  insCursr->closeBrsInnIndx = -1;
  insCursr->lvlOfRecur = -1;
  insCursr->intrvlInsFormdFlag = 0;
  insCursr->jumpTreeNext = NULL;
  insCursr->intrvlCntr = 0;
  todd->intrvlBeh = calloc(1, sizeof(knob));
  behCursr = todd->intrvlBeh;
  behCursr->opnBrsInnIndx = -1; // init as -1 instead of 0, since 0 presents the location 0 on the sequence
  behCursr->closeBrsInnIndx = -1;
  behCursr->lvlOfRecur = -1;
  behCursr->intrvlInsFormdFlag = 0;
  behCursr->jumpTreeNext = NULL;
  behCursr->intrvlCntr = 0;
  todd->intrvlInsFormdFlag = 0;
  todd->intrvlCntr = -1;
  todd->RSTO = NULL;

  if (intrvlCursr) {
    todd->intrvlLB = intrvlCursr->opnBrsInnIndx;
    todd->intrvlUB = intrvlCursr->closeBrsInnIndx;
  } else {
    todd->intrvlLB = -1; // init as -1 instead of 0, since 0 presents the location 0 on the sequence
    todd->intrvlUB = -1;
  }  // end if 2

  todd->lvlOfRecur = 0;

  return todd;
}  // end init_todd

//*****************************************************************************
// Function : Insert Behind Interval Normally
// Caller   : make_jump_beh_intrvl()
// Purpose  : make interval with reasonable size behind the helix NORMALLY
// Input    : crik - crick
//          : todd  
// Return   : none
// Display  : none
//*****************************************************************************
int insert_beh_intrvl_normally(config* seq __unused, global* crik, local* todd) {
  disp(seq,DISP_ALL,"enterng \"insert_beh_intrvl_normally\"\n");
  knob* behCursr = todd->intrvlBeh;
  knob* hlixCursr = crik->hlixInStru;

  behCursr->lvlOfRecur = hlixCursr->lvlOfRecur;
  behCursr->intrvlCntr = hlixCursr->intrvlCntr;
  behCursr->hlixBranchngIndx1 = hlixCursr->hlixBranchngIndx1;
  behCursr->parentIntrvlCntr = hlixCursr->intrvlCntr;
  behCursr->intrvlTypFlag = BEH_INTRVL;                // as a jump_beh 'flag'
  behCursr->jumpTreeNext = crik->interval;
  disp(seq,DISP_LV4,"behInterval formed on normally:");
  crik->interval = behCursr;
  disp(seq,DISP_LV4," (%d,%d), recur lv = %d\n", behCursr->opnBrsInnIndx, behCursr->closeBrsInnIndx, behCursr->lvlOfRecur);

  if (crik->interval == crik->interval->jumpTreeNext) {
    crik->interval->jumpTreeNext = NULL;
    printf("CIRCUS\n");
  }

  return 0;
}  // end insert_beh_intrvl_normally

//*****************************************************************************
// Function : Is Bundle Activated and Available?
// Caller   : make_jump_tree()
// Purpose  : check if bundle feature is activated and if so, if there is duminode available to use
// Input    : seq  
//          : crik 
//          : todd
// Return   : 0 (no) or 1 (yes)
// Display  : none
//*****************************************************************************
int bundle_is_available(config* seq, global* crik, local* todd) {
  return seq->bundle && 
         crik->eden[todd->cmpntLLCursr->opnBrsOutIndx][todd->cmpntLLCursr->closeBrsOutIndx].dumiNode;
}

//*****************************************************************************
// Function : Is Duplicate Interval?
// Caller   : make_jump_ins_intrvl() & make_jump_beh_intrvl()
// Purpose  : check if there's already another identical interval in front
// Input    : seq      - sequence
//          : intrvlLB - interval lower bound
//          : intrvlUB - interval upper bound
// Return   : 0 (no) or 1 (yes)
// Display  : none
//*****************************************************************************
int is_duplicate_intrvl(global* crik, int16_t intrvlLB, int16_t intrvlUB) {
  if (crik->interval) {
    if ((intrvlLB == crik->interval->opnBrsInnIndx) && 
        (intrvlUB == crik->interval->closeBrsInnIndx))
      return 1;
  }  // end if

  return 0;
}  // end is_duplicate_intrvl

//*****************************************************************************
// Function : Is Interval to Be Restored?
// Caller   : remove_intrvl()
// Purpose  : verify if this interval is to be restored
// Input    : crik   - crick
// Return   : 1: yes, 0: no
// Display  : none
//*****************************************************************************
int is_intrvl_2b_rsto(global* crik, local* toddP, int8_t recurRoute) {
  knob* rstoCursr;

  if (recurRoute == WHILE_SESSION) {
    rstoCursr = toddP->RSTO;

    while (rstoCursr) {
      if ((rstoCursr->opnBrsInnIndx == crik->interval->opnBrsInnIndx) && 
          (rstoCursr->closeBrsInnIndx == crik->interval->closeBrsInnIndx))
        return 0; // duplicate interval, shouldn't be added to toddP->RSTO
      else
        rstoCursr = rstoCursr->jumpTreeNext; // scan thru all the rsto list to make sure there's no duplicate
    }  // end while

    return 1;      // in WHILE_SESSION, and no duplicate. So, to be restored
  } else {
    return 0;             // not even in WHILE_SESSION, out of consideration
  }    // end outer if
}  // end is_intrvl_2b_rsto

//*****************************************************************************
// Function : Is It Time to Quit
// Caller   : make_jump_tree()
// Purpose  : check if it's still worth it to go thru the rest of the components
//          : since there may not be enough room for minimum required helices or hairpins 
//          : or since the some covariance pair has been missed
// Input    : seq  
//          : todd
// Return   : 1: time to wrap up; 0: not yet, we should keep looking for more structures
// Display  : none
//*****************************************************************************
int time_to_quit(config* seq, local* todd) {
  disp(seq,DISP_ALL,"enterng 'time_to_quit");
  int8_t numMustPair = seq->numCovari + seq->numV1Pairng;
  int8_t numV1Pairng = seq->numV1Pairng;

  while (numMustPair > seq->numCovari) {  // if there's v1 pairing constraints
    if (seq->v1Pairng[numV1Pairng - 1] < todd->cmpntLLCursr->opnBrsOutIndx) {
      disp(seq,DISP_ALL,"passed v1 constraints %d, we gotta wrap up\n", seq->v1Pairng[numV1Pairng - 1]);
      return 1;                                                      // ||
    } else
			// ||
    numMustPair--;                                          // ||
    numV1Pairng--;                                              // ||
  } // end while                                                              // \/

  while (numMustPair) { // if there's covariance and/or v1 pairing constraints
    if (seq->coVari[numMustPair - 1][0] < todd->cmpntLLCursr->opnBrsOutIndx) {
      disp(seq,DISP_ALL,"passed covari LB (%d,%d), we gotta wrap up\n\n", seq->coVari[numMustPair - 1][0], seq->coVari[numMustPair - 1][1]);
      return 1;                                                      // ||
    } else
      numMustPair--;                                          // ||
  } // end while                                                              // \/

  if (seq->minNumHlix) {              // if there's hlix min length constraint
    if (seq->strLen - todd->cmpntLLCursr->opnBrsOutIndx - seq->minPairngDist - 2 * seq->minLenOfHlix * seq->minNumHlix < 0) {           // ||
      disp(seq,DISP_LV1,"comp. start = %d\nminPairingDist = %d\nMinHelixLen = %d\nmin#ofHelices = %d\n",todd->cmpntLLCursr->opnBrsOutIndx, seq->minPairngDist, seq->minLenOfHlix, seq->minNumHlix);
      disp(seq,DISP_LV1,"process stopped prematurally at component type of %d\n", todd->cmpntLLCursr->opnBrsOutIndx);
      return 1;                                                      // ||
    } else {                                                           // ||
      return 0;                               // || we should keep looking
    } // end if                                                               // \/
  }  // end if

  return 0;                                 // no need to answer this question
}  // end time_to_quit

//*****************************************************************************
// Function : Is There Concern on Formation of Ring Linked List?
// Caller   : jump_stage_3_verify_hlix(), remove_intrvl(), remove_hlix()
// Purpose  : link crik->interval to the end of toddP->RSTO
// Input    : seq    - sequence
//          : crik   - crick
//          : todd
// Return   : 1 (yes), 0 (no)
// Display  : none
//*****************************************************************************
int is_there_concern_on_ring_formation_of_linked_list(config* seq __unused, global* crik, local* todd) {
  knob* rstoCursr = todd->RSTO;
  knob* intrvlCursr = crik->interval;

  while (intrvlCursr) {
    rstoCursr = todd->RSTO;
    while (rstoCursr) {
      if (rstoCursr->opnBrsInnIndx == crik->interval->opnBrsInnIndx || rstoCursr == crik->interval) {
        disp(seq,DISP_LV3,"Bad restoring process detected: Candidate node to be restored: [ %d[ ]%d ]\n", crik->interval->opnBrsInnIndx, crik->interval->closeBrsInnIndx);
        return 1;
      }  // end if
      rstoCursr = rstoCursr->jumpTreeNext;
    }    // end inner while
    intrvlCursr = intrvlCursr->jumpTreeNext;
  }      // end outer while

  return 0;
}  // end is_there_concern_on_ring_formtion_of_linked_list

//*****************************************************************************
// Function : Jump Stage 1: Set Interval
// Caller   : make_jump_tree(), jump_stage_3_verify_hlix
// Purpose  : set up the platform, ie. upper/lower bound, in which the helix may be placed
// Input    : seq 
//          : crik
//          : todd
//          : toddP
//          : bundlePathFlag - 1: if it's called under shortcut bundle list path
//          :                  0: if it's called under normal component list path
// Return   : none
// Display  : error messages, if necessary
// Note     : for speed up purpose, this function is rather long
//*****************************************************************************
int jump_stage_1_set_intrvl(config* seq, global* crik, local* todd, int16_t bundlePathFlag) {

  //Counts the number of calls for debugging purposes can be deleted safely if necessary (see top of file)
  numOfCalls++;

  disp(seq,DISP_ALL,"Entering 'jump_stage_1', ");
  g_x1++;
  dispLL(seq,crik,todd,NULL);

  knob* insCursr = todd->intrvlIns;
  knob* hlixCursr = crik->hlixInStru;
  int16_t tempLB;
  int16_t tempUB;

  bundlePathFlag++;
  // || MAKE JUMP INSIDE INTERVAL
  todd->intrvlIns->intrvlInsFormdFlag = 0;                              // ||
  todd->intrvlInsFormdFlag = 0;                                          // ||
  tempLB = hlixCursr->opnBrsInnIndx + 1;                                 // ||
  tempUB = (hlixCursr->opnBrsInnIndx - hlixCursr->opnBrsOutIndx) > (seq->minLenOfHlix - 1) ? // || if true, that means it conflict against l*p + k rule, which will bring duplicate
			hlixCursr->closeBrsInnIndx - 1 : // || sacrefice one nt space to prevent duplicate
			hlixCursr->closeBrsInnIndx - 1; // || no worry about duplicate, since l*p + k isn't violated
  if (((tempUB - tempLB + 1) < seq->insIntrvlMinSize) ||                // ||
      (seq->intrvlLukUpTable[tempLB][tempUB][0] < 0) || // || chk intrvl look-up table to see if this interval is on the table, if not, there's no cmpnt to fit
      (seq->intrvlLukUpTable[tempLB][tempUB][1] < 0) || is_duplicate_intrvl(crik, tempLB, tempUB)) {
    disp(seq,DISP_ALL,"No ins intvl made\n");
  } else {                                                               // ||
    insCursr->opnBrsInnIndx = tempLB;                                  // ||
    insCursr->closeBrsInnIndx = tempUB;                                // ||
    insCursr->hlixBranchngIndx1 = hlixCursr->hlixBranchngIndx1;        // ||
    insCursr->parentIntrvlCntr = hlixCursr->intrvlCntr;                // ||
    todd->intrvlInsFormdFlag = 1; // || as reference 4 'make_jump_beh_intrvl'
    insCursr->intrvlInsFormdFlag = 1;                                  // ||
    insCursr->lvlOfRecur = hlixCursr->lvlOfRecur;                      // ||
    insCursr->intrvlCntr = hlixCursr->intrvlCntr;                      // ||
    insCursr->intrvlTypFlag = INS_INTRVL;         // || as a jump_ins 'flag'
    insCursr->jumpTreeNext = crik->interval;
    disp(seq,DISP_LV4,"Intrvl inside formed:");
    crik->interval = insCursr;
    disp(seq,DISP_LV4," (%d,%d), recur lv = %d\n", insCursr->opnBrsInnIndx, insCursr->closeBrsInnIndx, insCursr->lvlOfRecur);
  } // end inner if                                                                            // \/

  tempLB = hlixCursr->closeBrsOutIndx + seq->minBtwnHlixDist + 1;        // ||
  tempUB = (todd->intrvlUB > 0) ? todd->intrvlUB : (seq->strLen - 1);    // ||

  if (((tempUB - tempLB + 1) < seq->insIntrvlMinSize) ||
      (seq->intrvlLukUpTable[tempLB][tempUB][0] < 0) || // || chk intrvl look-up table to see if this interval is on the table, if not, there's no cmpnt to fit
      (seq->intrvlLukUpTable[tempLB][tempUB][1] < 0) || is_duplicate_intrvl(crik, tempLB, tempUB)) {
    disp(seq,DISP_ALL,"No beh intvl made, either introvl too narrow, intrvl doesn't not fit or is duplicate\n");
  } else if ((crik->interval && // || This is used to fix the issue of identical crik->interval node mixed with todd->beh node
              crik->interval->jumpTreeNext && todd->intrvlBeh->opnBrsInnIndx
              == crik->interval->jumpTreeNext->opnBrsInnIndx) || crik->interval == todd->intrvlBeh) {
    disp(seq,DISP_LV3,"Error, ring formation of linked list detected\n");
    printf("infinite linked list\n");
  } else {                                                               // ||
    todd->intrvlBeh->opnBrsInnIndx = tempLB;                           // ||
    todd->intrvlBeh->closeBrsInnIndx = tempUB;                         // ||
    if (todd->intrvlInsFormdFlag) {                                    // ||
      swap_ins_n_beh_intrvl(seq, crik, todd);
      dispLL(seq,crik,todd,0);
    } else {                                                           // ||
      insert_beh_intrvl_normally(seq, crik, todd);
      dispLL(seq,crik,todd,0);
    } // end inner if 2                                                                         // ||
  } // end outer if                                                                           // \/

  while (crik->interval && crik->hlixInStru) {
    if (crik->hlixInStru->intrvlCntr == crik->interval->intrvlCntr)
      jump_stage_2_fit_hlix(seq, crik, todd, FIRST_SESSION); // <----------RECURSION this way (1st session)
    else
      break;
  }

  while (crik->interval) // exhaust all intervals present in crik->interval (while session)
    jump_stage_2_fit_hlix(seq, crik, todd, WHILE_SESSION); // <---------- RECURSION this way

  dispLL(seq,crik,todd,NULL);
  return 0;
}  // end jump_stage_1_set_intrvl

//*****************************************************************************
// Function : Jump Stage 2: Fit Helix
// Caller   : jump_stage_1_set_intrvl()
// Purpose  : search for proper helix to fit into the interval
//          : used when a new component type is to be visited (not the first helix of the structure, cuz the first one is covered by 'make_jump_tree() )
// Input    : seq    - sequence
//          : crik   - crick
//          : toddP  - todd's parent
// Return   : none
// Display  : error message, if necessary
// Note     : new recursion levels generally start here, assigned in init_todd()
//          : stage 2 scans thru component types & the components under the same type, also check the LCM & combinatory types of duplicates
//          : stage 3 scans verify if a given helix may be fit into the structure, which has been integrated into stage 2
//          : l*p + k rule is applied here
// Note     : for speed up purpose, this function is rather long
//*****************************************************************************
int jump_stage_2_fit_hlix(config* seq, global* crik, local* toddP, int8_t recurRoute) {
  jst2Count++;
  disp(seq,DISP_ALL,"Enterng 'jump_stage_2', recur lv: %d->%d\n", crik->interval->lvlOfRecur, (crik->interval->lvlOfRecur + 1) ); disp(seq,DISP_ALL,"location source = %s\n", CALLER_FLAG_TO_STRING[recurRoute]);
  local* todd = NULL;
  todd = init_todd(crik);
  knob* intrvlCursr = crik->interval;
  knob* intrvlMsgr = crik->interval; // to carry on interval information, and keep moving, since the top interval is gonna be removed in 'remove_intrvl()'
  knob* cmpntCursr;
  knob* hlixCursr;
  int16_t* tempInt;
  int16_t i = 0; // i has to be initialized cuz it's going to be used in while loop
  int16_t j, k;
  int16_t flagUB = seq->numCovari + seq->numV1Pairng;
  int8_t V1FailFlag = 0;

									   // || select interval-look-up table
  if (toddP) {
    disp(seq,DISP_ALL,"crik->interval = (%d,%d)\n", intrvlCursr->opnBrsInnIndx, intrvlCursr->closeBrsInnIndx);
    disp(seq,DISP_ALL,"intrvlLukUpTab = %d\n", seq->intrvlLukUpTable[intrvlCursr->opnBrsInnIndx][intrvlCursr->closeBrsInnIndx][0]);

    todd->lukUpCmpntTypLB = seq->intrvlLukUpTable[intrvlCursr->opnBrsInnIndx][intrvlCursr->closeBrsInnIndx][0];
    todd->lukUpCmpntTypUB = seq->intrvlLukUpTable[intrvlCursr->opnBrsInnIndx][intrvlCursr->closeBrsInnIndx][1];   } // end if

  if (intrvlCursr) {               // set todd intrvl upper and lower bound
    todd->intrvlLB = intrvlCursr->opnBrsInnIndx;
    todd->intrvlUB = intrvlCursr->closeBrsInnIndx;
  } else {
    todd->intrvlLB = -1; // init as -1 instead of 0, since 0 presents the location 0 on the sequence
    todd->intrvlUB = -1;
  } // end if 2 

  todd->lvlOfRecur = toddP ? toddP->lvlOfRecur + 1 : 0;
  disp(seq,DISP_ALL,"Todd made, recur lvl = %d\n", todd->lvlOfRecur);
  crik->lvlOfRecur++;

  dispLL(seq,crik,0,toddP); 
  disp(seq,DISP_ALL,"intrvl: (%d,%d) -> (LB UB): (%d,%d)\n", crik->interval->opnBrsInnIndx, 
                                                             crik->interval->closeBrsInnIndx, 
                                                             todd->lukUpCmpntTypLB, 
                                                             todd->lukUpCmpntTypUB);
  dispLL(seq,crik,0,toddP);

  remove_intrvl(seq, crik, toddP, recurRoute);
  disp(seq,DISP_ALL,"This 2 b hooked bak later\n");
  crik->intrvlCntr++; // it's like batch number, used to distinguish one batch from another
                      // to facilitate new helix insert and old helix removal

  todd->intrvlCntr = crik->intrvlCntr;
  while ((i < crik->numCmpntTypOcupid) && (crik->cmpntListOcupidTyp[i] < todd->lukUpCmpntTypLB))
  i++; // the lowerbound of the looked-up component is at the left of the target lowerbound

  if (intrvlMsgr->intrvlTypFlag == BEH_INTRVL)
    crik->numHP++; // hairpin increases only if the new interval is behind interval 
                   // (inside interval only nest in further at the same hairpin)

  while ((i < crik->numCmpntTypOcupid) && (crik->cmpntListOcupidTyp[i] <= todd->lukUpCmpntTypUB)) { 
    // scan thru the whole cmpnt list of interested types. previously it scan only up to 
    // todd->lukUpCmpntTypUB + 1 due to l*p + k rule, we have to add one up to accomidate 
    // the missing bit (see tempUB making of make_jump_ins_intrvl in jump_stage_1_set_intrvl()

    cmpntCursr = crik->cmpntList[crik->cmpntListOcupidTyp[i]].knob;
    disp(seq,DISP_ALL,"Jump Stage 2 'While' loop cmpnt:\n");
    hlixCursr = crik->hlixInStru;
    dispLL(seq,crik,0,toddP);

    while (cmpntCursr) { 
      // SCAN A COMPONENT TYPE : scan thru the components of the same component type
      hlixCursr = crik->hlixInStru;

      todd->cmpntLLCursr = cmpntCursr;
      dispLL(seq,crik,todd,toddP);

      g_x2++;
      disp(seq,DISP_ALL,"cmpnt cadidate: {%2d-%-2d{    }%2d-%-2d}\n", cmpntCursr->opnBrsOutIndx, 
                                                                      cmpntCursr->opnBrsInnIndx, 
                                                                      cmpntCursr->closeBrsInnIndx, 
                                                                      cmpntCursr->closeBrsOutIndx);
 
     if (cmpntCursr->closeBrsOutIndx <= todd->intrvlUB) { 
        // HELIX FITS! helix size fittable to the interval, and going to be placed in

        cmpntCursr->intrvlCntr = todd->intrvlCntr;
        dispLL(seq,crik,todd,toddP);
        if (hlixCursr) {
          disp(seq,DISP_ALL,"            (lvlOfRecur, branchngIndx, intrvlCntr)\n"); 
          disp(seq,DISP_ALL,"hlixInStru: (%10d, %12d, %10d)\n", hlixCursr->lvlOfRecur, 
                                                                hlixCursr->hlixBranchngIndx1, 
                                                                hlixCursr->intrvlCntr); 
          disp(seq,DISP_ALL,"intrvlMsgr: (%10d, %12d, %10d)\n", intrvlMsgr->lvlOfRecur, 
                                                                intrvlMsgr->hlixBranchngIndx1, 
                                                                cmpntCursr->intrvlCntr);

          if ((hlixCursr->lvlOfRecur == intrvlMsgr->lvlOfRecur + 1) && 
              (hlixCursr->hlixBranchngIndx1 == intrvlMsgr->hlixBranchngIndx1) &&
              (hlixCursr->intrvlCntr == todd->intrvlCntr)) { 
            // REPLACE OLD HELIX BY NEW ONE : it takes three parameters to ensure that 
            // the present helix to be erased is indeed from the same level

            int maxMismatches = 0;
            int numMismatches = 0;
            if (!cmpntCursr->bundleFlag) {
              maxMismatches = 1;
              int q;
              for (q = 0; q < seq->maxNumMismatch; q++) {
                if (cmpntCursr->mismatchFlag[q][0] == -1) {
                  maxMismatches = 0;
                } else {
                  numMismatches++;
                }
              }
            }

            if (hlixCursr && cmpntCursr && hlixCursr->jumpTreeNext && 
                (hlixCursr->jumpTreeNext->opnBrsInnIndx - cmpntCursr->opnBrsOutIndx == -1) && 
                (hlixCursr->jumpTreeNext->closeBrsInnIndx - cmpntCursr->closeBrsOutIndx == 1)) {
              // l*p + k rule in the case of add new hlix w/o removing old one
              int q;
              if (!hlixCursr->jumpTreeNext->bundleFlag) {
                for (q = 0; q < seq->maxNumMismatch; q++) {
                  if (hlixCursr->jumpTreeNext->mismatchFlag[q][0] != -1) {
                    numMismatches++;
                  }
                }
              }

              crik->linkedmms += numMismatches;
              if (crik->linkedmms > seq->maxNumMismatch) {
                crik->linkedmms = 0;
                break;
              }

              if (seq->maxNumMismatch > 0) {
                if (((hlixCursr->jumpTreeNext->opnBrsInnIndx - 
                      hlixCursr->jumpTreeNext->opnBrsOutIndx != (seq->minLenOfHlix - 1)) && 
                     ((hlixCursr->jumpTreeNext->mismatchFlag[seq->maxNumMismatch-1][0] != 
                       hlixCursr->jumpTreeNext->opnBrsInnIndx || maxMismatches) || 
                      ((cmpntCursr->opnBrsInnIndx - hlixCursr->jumpTreeNext->opnBrsOutIndx + 1) >= 
                       (seq->minLenOfHlix * 3) && !cmpntCursr->bundleFlag && 
                       !hlixCursr->jumpTreeNext->bundleFlag)) && 
                     (!cmpntCursr->bundleFlag && (cmpntCursr->mismatchFlag[0][0] != 
                      cmpntCursr->opnBrsOutIndx || cmpntCursr->opnBrsInnIndx - cmpntCursr->opnBrsOutIndx + 1 != 
                      seq->minLenOfHlix)))) { 
                  // if this check point is passed, the following will be 'REPLACE OLD HELIX BY NEW ONE'
                  crik->linkedmms = 0;
                  break;
                }
              } else if(hlixCursr->jumpTreeNext->opnBrsInnIndx - 
                        hlixCursr->jumpTreeNext->opnBrsOutIndx != (seq->minLenOfHlix - 1)) {
                crik->linkedmms = 0;
                break;
              }
            } else {
              crik->linkedmms = 0;
            }
            dispLL(seq,crik,todd,toddP);
            cmpntCursr->intrvlCntr = todd->intrvlCntr;
            dispLL(seq,crik,todd,toddP);

            if (seq->numCovari) { 
              // update constraint flags (currently, only V1 and covariance are included)
              flagUB = seq->numCovari + seq->numV1Pairng;
              tempInt = calloc(flagUB, sizeof(int16_t));

              for (j = 0; j < flagUB; j++) {
                if (!cmpntCursr->mustPairFlag[j] && hlixCursr->mustPairFlag[j]) { 
                  // replace old hlix will lose a covariance pair or v1 pairing in helix
                  crik->struMustPairFlag[j] = 0;
                } // end inner if
                dispLL(seq,crik,todd,toddP);

                // if there's no issue of losing covariance pair, then the new flag set should contain all flags
                tempInt[j] = (cmpntCursr->mustPairFlag[j] || crik->struMustPairFlag[j]); 
              } // end for

              free(crik->struMustPairFlag);
              crik->struMustPairFlag = tempInt; // replace old flags set by new one
              crik->mustPairLength = flagUB;
            } // end outer if

            // unhook the reference to old helix, and hook reference to the new helix. 
            // erase the record of the old helix, so that it may be used again somewhere else
            cmpntCursr->jumpTreeNext = hlixCursr->jumpTreeNext;
            hlixCursr->jumpTreeNext = NULL;
            hlixCursr->outsideIntrvlLB = -1;
            hlixCursr->outsideIntrvlUB = -1;
            crik->hlixInStru = cmpntCursr;
            hlixCursr = crik->hlixInStru;

            // write in the new helix
            hlixCursr->outsideIntrvlLB = todd->intrvlLB;
            hlixCursr->outsideIntrvlUB = todd->intrvlUB;
            hlixCursr->lvlOfRecur = intrvlMsgr->lvlOfRecur + 1;
            hlixCursr->hlixBranchngIndx1 = intrvlMsgr->hlixBranchngIndx1;

            int rstoFlag = 0;
            if (todd->RSTO && !is_there_concern_on_ring_formation_of_linked_list(seq, crik, todd)) {
              disp(seq,DISP_ALL,"            special On Q %d-%d\n", todd->RSTO->opnBrsInnIndx,
                                                                    todd->RSTO->closeBrsInnIndx);
              crik->interval = todd->RSTO;
              rstoFlag = 1;
              rstoCount++;
              disp(seq,DISP_ALL,"Interval restore todd on queue!\n");
            } // end outer if

            dispLL(seq,crik,todd,toddP);

// || make sure there'r enough helices (equal or larger than the min requirement)// || make sure there'r enough hairpins (equal or larger than the min requirement)
            if ((crik->numHlix >= seq->minNumHlix) && (crik->numHP >= seq->minNumHP)) { 
              for (k = 0; k < flagUB; k++)
// || make sure all the must pair nucleotides are paired, otherwise this structure isn't qualified, and nothing will happen
                if (!crik->struMustPairFlag[k]) { 
                  V1FailFlag = 1;
                  break; // || no point to continue looping, since this structure is disqualified for sure
                } // end if

              dispLL(seq,crik,todd,toddP);

              if (!V1FailFlag) { // make sure all V1 nucleotides are paired
// if this point is reached, that means this structure contains enough helices and hairpins, and all the V1 pairs are paird
                crik->numStru++; 
                if (seq->bundle) { // count unbundled structures
                  countAndPrintStructures(seq, crik);
                  crik->linkedmms = 0;
                } else {
                  disp_stru(seq, crik, 11); // for report  (official use)
                  dispStru(seq,crik,11); // for display (debug use)
                  crik->linkedmms = 0;
                }
              } // end inner 1 if
            } // end outer if

            todd->intrvlInsFormdFlag = 0;
#ifdef _MPI
            if(numOfCalls % CHUNK_SIZE == (CHUNK_SIZE - 1)) {
              int to = is_work_needed();
              if(to != -1) {
                send_work(crik, todd, to);
              }
            }
#else
            jump_stage_1_set_intrvl(seq, crik, todd, 0); // regular path, w/o inside interval restoration
#endif
          } else {     // ADD A NEW HELIX, NO HELIX IS REMOVED
            cmpntCursr->intrvlCntr = todd->intrvlCntr;
            dispLL(seq,crik,todd,toddP);

            for (j = 0; j < flagUB; j++) // || update constraint flags, since no helix is removed, no need to uncheck any box, but just add new checks
              if (cmpntCursr->mustPairFlag[j]) 
                crik->struMustPairFlag[j] = 1; // || covari flag got '1'   <==================== Must Pair (covari or V1)
            cmpntCursr->jumpTreeNext = hlixCursr;      // ||
            crik->numHlix++;                           // ||
            crik->hlixInStru = cmpntCursr;             // ||
            hlixCursr = crik->hlixInStru;              // ||
            hlixCursr->outsideIntrvlLB = todd->intrvlLB; // || write new one
            hlixCursr->outsideIntrvlUB = todd->intrvlUB; // ||
            hlixCursr->lvlOfRecur = intrvlMsgr->lvlOfRecur + 1;             // ||
            hlixCursr->hlixBranchngIndx1 = intrvlMsgr->hlixBranchngIndx1;     // ||

            int maxMismatches = 0;
            int numMismatches = 0;
            if (!hlixCursr->bundleFlag) {
              maxMismatches = 1;
              int q;
              for (q = 0; q < seq->maxNumMismatch; q++) {
                if (hlixCursr->mismatchFlag[q][0] == -1) {
                  maxMismatches = 0;
                } else {
                  numMismatches++;
                }
              }
            }

            int skip = 0;
            if((hlixCursr->jumpTreeNext->opnBrsInnIndx - hlixCursr->opnBrsOutIndx == -1) && // || l*p + k rule in the case of add new hlix w/o removing old one
               (hlixCursr->jumpTreeNext->closeBrsInnIndx - hlixCursr->closeBrsOutIndx == 1) && 
               !hlixCursr->jumpTreeNext->bundleFlag) {
              int q;
              for (q = 0; q < seq->maxNumMismatch; q++) {
                if(hlixCursr->jumpTreeNext->mismatchFlag[q][0] != -1) {
                  numMismatches++;
                }
              }
              crik->linkedmms += numMismatches;
              if (crik->linkedmms > seq->maxNumMismatch) {
                skip = 1;
              }
            } else {
              crik->linkedmms = 0;
            }

            if ((hlixCursr->jumpTreeNext->opnBrsInnIndx - hlixCursr->opnBrsOutIndx == -1) && // || l*p + k rule in the case of add new hlix w/o removing old one
                (hlixCursr->jumpTreeNext->closeBrsInnIndx - hlixCursr->closeBrsOutIndx == 1) && 
                seq->maxNumMismatch > 0 && 
                (hlixCursr->jumpTreeNext->opnBrsInnIndx - hlixCursr->jumpTreeNext->opnBrsOutIndx != 
                 (seq->minLenOfHlix - 1)) && 
                ((hlixCursr->jumpTreeNext->mismatchFlag[seq->maxNumMismatch - 1][0] != 
                  hlixCursr->jumpTreeNext->opnBrsInnIndx || maxMismatches) || 
                 ((hlixCursr->opnBrsInnIndx	- hlixCursr->jumpTreeNext->opnBrsOutIndx + 1) >= 
                  (seq->minLenOfHlix * 3) && !hlixCursr->bundleFlag && !hlixCursr->jumpTreeNext->bundleFlag)) &&
                (!hlixCursr->bundleFlag && (hlixCursr->mismatchFlag[0][0] != hlixCursr->opnBrsOutIndx || 
                                            hlixCursr->opnBrsInnIndx - hlixCursr->opnBrsOutIndx + 1	!= 
                                            seq->minLenOfHlix))) { // ||
              crik->linkedmms = 0;

            } else if((hlixCursr->jumpTreeNext->opnBrsInnIndx - hlixCursr->opnBrsOutIndx == -1) && // || l*p + k rule in the case of add new hlix w/o removing old one
                      (hlixCursr->jumpTreeNext->closeBrsInnIndx - hlixCursr->closeBrsOutIndx == 1) && 
                      (hlixCursr->jumpTreeNext->opnBrsInnIndx - hlixCursr->jumpTreeNext->opnBrsOutIndx != 
                       (seq->minLenOfHlix - 1)) && 
                      seq->maxNumMismatch == 0) {
              crik->linkedmms = 0;
            } else if (skip) {
              crik->linkedmms = 0;
            } else {
// || make sure there'r enough helices (equal or larger than the min requirement)// || make sure there'r enough hairpins (equal or larger than the min requirement)
              if ((crik->numHlix >= seq->minNumHlix) && (crik->numHP >= seq->minNumHP)) { 
                for (k = 0; k < flagUB; k++)
                  if (!crik->struMustPairFlag[k])
                    V1FailFlag = 1; // || make sure all the must pair nucleotides are paired, otherwise this structure isn't qualified, and nothing will happen

                if (!V1FailFlag) { // || if this point is reached, that means this structure contains enough helices and hairpins, and all the V1 pairs are paird
                  crik->numStru++;
                  if(seq->bundle) { // count unbundled structures
                    countAndPrintStructures(seq, crik);
                    crik->linkedmms = 0;
                  } else {
                    disp_stru(seq, crik, 11); // for report  (official use)
                    dispStru(seq,crik,11); // for display (debug use)
                    crik->linkedmms = 0;
                  }
                } // end inner 1 if
              } // end outer if
              todd->intrvlInsFormdFlag = 0;

#ifdef _MPI
              if(numOfCalls % CHUNK_SIZE == (CHUNK_SIZE - 1)) {
                int to = is_work_needed();
                if(to != -1) {
                  send_work(crik, todd, to);
                }
              }            
#else
              jump_stage_1_set_intrvl(seq, crik, todd, 0); // regular path, w/o inside interval restoration
#endif
            }  // end if (used by l*p + k feature)
          }  // end if B2

        } else { // || NO CHANGE OF HELIX : to prevent duplicate (currently there's no helix in structure)
          cmpntCursr = cmpntCursr->cmpntListNext;
          continue;
        } // end if B1 1

      } else {
        disp(seq,DISP_ALL,"hlix tail exceeds intrvl UB: %d >= %d\n", cmpntCursr->closeBrsOutIndx, 
                                                                     todd->intrvlUB); 
        disp(seq,DISP_ALL,"chk pnt 8       \n"); 
        dispLL(seq,crik,todd,toddP); 
        disp(seq,DISP_ALL,"Head bak to jump_stage_2 recur lv: %d->%d\n", todd->lvlOfRecur, todd->lvlOfRecur);
        break; // should be break, not continue, since we are sure there's no more helix in this component whose tail will fit in
      }      // end if gnd lvl 1
      dispLL(seq,crik,todd,toddP); 
      disp(seq,DISP_ALL,"Exiting 'jump_stage_3', recur lv: %d->%d\n", intrvlMsgr->lvlOfRecur, 
                                                                      intrvlMsgr->lvlOfRecur);
      cmpntCursr = cmpntCursr->cmpntListNext; // when reaching this point, either the futher potential recursions have been exhausted
    }      // end inner while
    i++;

  } // end while 2                                      // all legitimate cmpnt types have been exhausted, which is the end of this hairpin

  if (intrvlMsgr->intrvlTypFlag == BEH_INTRVL)
    crik->numHP--; // that concluded this round of hairpin quest, since the legitimate cmpnt types have been exhausted

  exit_curr_recur(seq, crik, todd);
  return 0;
}  // end jump_stage_2_fit_hlix

//*****************************************************************************
// Function : Is There Larger Dummy to Replace Current Candidate Component with Concern on Component Insert behind
// Caller   : is_there_prev_dumi_2_replace_curr_candidate_dumi_recur_lv_0()
// Purpose  : carefully search around for larger dummy to replace the current candidate component, and use the smallest dummy whose 5' end leg (component type) is larger than the 3' end leg (closeBrsOutIndx) of current candidate component
// Input    : crik
//          : cmpntCursr
//          : iCmpontTyp
//          : cmpntOcupidID
// Return   : 1 (yes), 0 (no)
// Display  : none
//*****************************************************************************
int is_there_larger_dumi_2_replace_curr_candidate_cmpnt_with_concern_on_cmpnt_insert_behind(
		global* crik, knob* cmpntCursr, int16_t iCmpntTyp, int16_t cmpntOcupidID) {
  int16_t j;
  int16_t k = cmpntOcupidID;

  while ((k < crik->numCmpntTypOcupid - 1) && // find the smallest component whose leading edge (5' leg) is larger than the trailing edge (3' leg) of the current candidate component
         (cmpntCursr->closeBrsOutIndx >= crik->cmpntListOcupidTyp[k]))
    k++; // search for the component type which is closest to trailing edge (closeBrsOutIndx) of current compoent, but larger
  for (j = cmpntCursr->closeBrsOutIndx; j < crik->cmpntListOcupidTyp[k]; j++) // use the component type 'crik->cmpntListOcupidTyp[k]' as the upperbound to increment j
    if (crik->eden[iCmpntTyp][j].dumiUsdFlag4RecurLv0)
      return 1; // if flag is on, there must be a dummy, but there's dummy doesn't mean it has been used. So, we have to check this flag, not dumiNode

  return 0;
} // end is_there_larger_dumi_2_replace_curr_candidate_cmpnt_with_concern_on_cmpnt_insert_behind

//*****************************************************************************
// Function : Is There Larger Dummy to Replace Current Candidate Component without Concern on Component Insert behind
// Caller   : is_there_prev_dumi_2_replace_curr_candidate_dumi_recur_lv_0()
// Purpose  : freely search around for larger dummy to replace the current candidate component, as long as don't make the component span to large, since no dummy can be too large
// Input    : seq
//          : crik
//          : cmpntCursr
//          : bundleSpan
//          : iCmpontTyp
// Return   : 1 (yes), 0 (no)
// Display  : none
//*****************************************************************************
int is_there_larger_dumi_2_replace_curr_candidate_cmpnt_wo_concern_on_cmpnt_insert_behind(
		config* seq, global* crik, knob* cmpntCursr, int16_t bundleSpan, int16_t iCmpntTyp) {

  int16_t jUB = (iCmpntTyp + bundleSpan < seq->strLen) ? (iCmpntTyp + bundleSpan) : seq->strLen; // no worry, no component can insert behind, so just go ahead and expand 'j' to its allowed limit, which is the bundle span
  int16_t j;

  for (j = cmpntCursr->closeBrsOutIndx; j < jUB; j++) // increment closing parenthesis (3' leg) to see if if there's any's any larger bundle to cover the current candidate component
    if (crik->eden[iCmpntTyp][j].dumiUsdFlag4RecurLv0)
      return 1; // if flag is on, there must be a dummy, but there's dummy doesn't mean it has been used. So, we have to check this flag, not dumiNode

  return 0;
} // end is_there_larger_dumi_2_replace_curr_candidate_cmpnt_wo_concern_on_cmpnt_insert_behind

//*****************************************************************************
// Function : Remove Interval
// Caller   : jump_stage_2_fit_hlix()
// Purpose  : detach the top interval node (not delete)
// Input    : seq    - sequence
//          : crik   - crick
// Return   : none
// Display  : none
//*****************************************************************************
int remove_intrvl(config* seq, global* crik, local* toddP, int8_t recurRoute) {
  disp(seq,DISP_ALL,"Enterng 'remove_intrvl'\n");
  knob* temp;
  knob* rstoCursr = toddP->RSTO;

  if(is_intrvl_2b_rsto(crik, toddP, recurRoute) && 
      !is_there_concern_on_ring_formation_of_linked_list(seq, crik, toddP)) {
    disp(seq,DISP_ALL,"intrvl [%d-%d] left behind, 2b rsto\n", crik->interval->opnBrsInnIndx, crik->interval->closeBrsInnIndx);
    if (rstoCursr) {
      while (rstoCursr->jumpTreeNext)
        rstoCursr = rstoCursr->jumpTreeNext; // || RESTORE!! scan to locate the tail of toddP->RSTO
      rstoCursr->jumpTreeNext = crik->interval; // || hook the whole crik->interval on!
    } else if (crik->interval) {                                       // ||
      toddP->RSTO = crik->interval; // || toddP->RSTO is empty, therefore no need for scan
    } else {                                                           // ||
      toddP->RSTO = NULL;                                            // ||
    } // end if                                                           // ||
    crik->interval = crik->interval->jumpTreeNext; // \/ left behind the first node of crik->interval
  } else {
    disp(seq,DISP_ALL,"intrvl [%d-%d] removed\n", crik->interval->opnBrsInnIndx, crik->interval->closeBrsInnIndx );
    temp = crik->interval->jumpTreeNext; // || NO RESTORE  just go ahead and remove this interval
    crik->interval = temp;                                             // ||
  } // end if                                                             // \/

  return 0;
}  // end remove_intrvl

//*****************************************************************************
// Function : Search for Largest Representative of This Bundle Group
// Caller   : take_bundle_list_shortcut()
// Purpose  : use the following strategy to reduce the similar structure to max
//          : 1) if the current component trailing edge (closeBrsOutIndx) equal or greater than the leading edge (cmpntListOcupidTyp) of the very last component of the whole component list, then we don't need to worry about anything, but just go ahead and explore the largest component as representative of the whole group
//          : 2) if the current component trailing edge (closeBrsOutIndx) is smaller than the leading edge (cmpntListOcupidTyp) of the very last component of the whole component list, then we have to:
//          :    a) search for the first component type whose leading edge greater than the  current component trailing edge
//          :    b) check for the largest dummy node whose trailing edge is still smaller than the leading edge of the component of a)
//          :    c) use this dummy node as the largest representative
// Input    : crik
//          : todd  
//          : cmpntOcupidID
// Return   : none
// Display  : none
//*****************************************************************************
int search_4_largest_rep_of_this_bundle_group(global* crik, local* todd, int16_t cmpntOcupidID) {
  int16_t i = cmpntOcupidID;

  if(todd->cmpntLLCursr->closeBrsOutIndx >= crik->cmpntListOcupidTyp[crik->numCmpntTypOcupid - 1]) { // when if is true, it means there's no other component available to insert on the right of this current component referred by 'todd->cmpntLLCursr'
    while (todd->cmpntLLCursr->cmpntListNext && // make sure there's larger component available at cmpntListNext
           crik->eden[todd->cmpntLLCursr->cmpntListNext->opnBrsOutIndx][todd->cmpntLLCursr->cmpntListNext->closeBrsOutIndx].dumiNode) // make sure this larger component isn't too large to exceed bundle list upper bound
      todd->cmpntLLCursr = todd->cmpntLLCursr->cmpntListNext; // if the above case is true, there's no need to print out each of the component indvidually, but can just pick the longest one as representative under the same 'bundle' in a sense
  } else {
    while (todd->cmpntLLCursr->closeBrsOutIndx >= crik->cmpntListOcupidTyp[i])
      i++; // search for the component type which is closest to trailing edge (closeBrsOutIndx) of current compoent, but larger
    while ((todd->cmpntLLCursr->cmpntListNext) && // make sure there's larger component available at cmpntListNext
           (todd->cmpntLLCursr->cmpntListNext->closeBrsOutIndx < crik->cmpntListOcupidTyp[i]) && 
           (crik->eden[todd->cmpntLLCursr->cmpntListNext->opnBrsOutIndx][todd->cmpntLLCursr->cmpntListNext->closeBrsOutIndx].dumiNode)) // make sure this larger component isn't too large to exceed bundle list upper bound
      todd->cmpntLLCursr = todd->cmpntLLCursr->cmpntListNext; // if the above case is true, there's no need to print out each of the component indvidually, but can just pick the longest one as representative under the same 'bundle' in a sense
  }  // end if

  return 0;
}  // end search_4_largest_rep_of_this_bundle_group

//*****************************************************************************
// Function : Swap Inside and Behind Interval
// Caller   : make_jump_beh_intrvl()
// Purpose  : to restore interval (inside interval) once before behind interval is ever investigated
// Input    : seq
//          : crik
//          : todd  
// Return   : none
// Display  : none
//*****************************************************************************
int swap_ins_n_beh_intrvl(config* seq __unused, global* crik, local* todd) {
  disp(seq,DISP_ALL,"enterng 'swap_ins_n_beh_intrvl'\n");
  if (!crik->interval)
    return 0;

  knob* behCursr = todd->intrvlBeh;
  knob* hlixCursr = crik->hlixInStru;

  behCursr->lvlOfRecur = hlixCursr->lvlOfRecur;
  behCursr->intrvlCntr = hlixCursr->intrvlCntr;
  behCursr->hlixBranchngIndx1 = hlixCursr->hlixBranchngIndx1;
  behCursr->parentIntrvlCntr = hlixCursr->intrvlCntr;
  behCursr->intrvlTypFlag = BEH_INTRVL;                // as a jump_beh 'flag'
  behCursr->jumpTreeNext = crik->interval->jumpTreeNext;
  disp(seq,DISP_ALL,"Swap of ins vs beh\n");
  crik->interval->jumpTreeNext = behCursr;
  disp(seq,DISP_LV4," (%d,%d), recur lv = %d\n", behCursr->opnBrsInnIndx, behCursr->closeBrsInnIndx, behCursr->lvlOfRecur);

  if (crik->interval == crik->interval->jumpTreeNext)
    fprintf(stderr, "BAD BAD BAD!!! CIRCUS!!!\n");

  return 0;
}  // end swap_ins_n_beh_intrvl

//*****************************************************************************
// Function : Take Bundle List - Shortcut
// Caller   : make_jump_tree()
// Purpose  : Seek structures from 'bundle list', so that a whole series of partial structures inside the same interval may be replaced by a dummy node
// Input    : seq  
//          : crik 
//          : todd
//          : hlixBranchngIndx1
// Return   : none
// Display  : none
//*****************************************************************************
//int take_bundle_list_shortcut(config* seq, global* crik, local* todd,
//		int16_t hlixBranchngIndx1) {
//  if(is_there_prev_dumi_2_replace_curr_candidate_dumi_recur_lv_0(seq, crik, todd, cmpntOcupidID)) return 0; // if there's other dummy large enough to cover this one, this one should be skipped // EDIT: no need for this since each bundle is independent
//	search_4_largest_rep_of_this_bundle_group(crik, todd, cmpntOcupidID); // no need for this either
//	hook_dumi_node_on_stru(crik, todd, hlixBranchngIndx1);
//	eval_prog_n_prt_stru(seq, crik);
//	jump_stage_1_set_intrvl(seq, crik, todd, 0);
//	crik->numHlix--;
//	return 0;
//}  // end take_bundle_list_shortcut

//*****************************************************************************
// Function : Take Component List - Normal Path
// Caller   : make_jump_tree()
// Purpose  : Seek structures from 'component list', which takes extra work to combine into structures, instead of from the 'bundle list', which has completed the combination jobs
// Input    : seq  
//          : crik 
//          : todd
//          : toddP
//          : hlixBranchngIndx1
// Return   : none
// Display  : none
//*****************************************************************************
int take_cmpnt_list_normal_path(config* seq, global* crik, local* todd, int16_t hlixBranchngIndx1) {
  disp(seq,DISP_ALL,"Enterng 'take_cmpnt_list_normal_path'\n");
  crik->hlixInStru = todd->cmpntLLCursr;
  assign_new_hlix_link(crik, hlixBranchngIndx1);
  dispLL(seq,crik,todd,0);
  if (satisfy_constraint(seq, crik, todd))
    eval_prog_n_prt_stru(seq, crik);
  jump_stage_1_set_intrvl(seq, crik, todd, 0);
  return 0;
}  // end take_cmpnt_list_normal_path
