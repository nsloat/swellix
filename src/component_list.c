/*****************************************************************************************************************************
 *****************************************************************************************************************************
 *****************************************************************************************************************************
      Program : Swellix (Make Component List) - c file
      Purpose : Make Component List
      Input   : config* seq
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

#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "component_list.h"

//*****************************************************************************
// Function : Assign Window Width
// Caller   : recur_stage_2_brace_sizng()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//*****************************************************************************
int assign_window_width(config* seq, global* crik, int callerType)
{
                                                                                                     disp(seq,DISP_ALL,"Entering 'assign_window_width'\n");
  int16_t windowWidth;                                                                               // farthest possible extension for window sizing  
  int     refVal;

  if (seq->noSZG){
    windowWidth = seq->minLenOfHlix;
  } else {
    switch(callerType){
    case 0:                                                                                          // evaluate open window width largest value
      windowWidth = (crik->closeBrsOutIndx - seq->minPairngDist - crik->opnBrsOutIndx + 1) / 2 + seq->maxNumBlg;
      if(windowWidth <  seq->minLenOfHlix)     windowWidth = seq->minLenOfHlix;                      // min length of helix must be kept
      if(windowWidth >= seq->minLenOfHlix * 2) windowWidth = seq->minLenOfHlix * 2 - 1;              // to prevent duplicate, size double the min length is aboslutely redundant w/o doubt - the 2n-1 rule of thumb

      break;
    case 1:                                                                                          // evaluate close window width smallest value
      windowWidth = crik->opnBrsWidth - seq->maxNumBlg;
      if (windowWidth < seq->minLenOfHlix) windowWidth = seq->minLenOfHlix;
      break;
    case 2:                                                                                          // evaluate close window width largest value
      windowWidth = crik->opnBrsWidth + seq->maxNumBlg;
      refVal      = crik->closeBrsOutIndx - seq->minPairngDist - crik->opnBrsWidth - crik->opnBrsOutIndx + 1;
      if (windowWidth > refVal) windowWidth = refVal;
      break;
    }  // end switch
  }    // end if

  return windowWidth;
}  // end assign_window_width

//*****************************************************************************
// Function : Book Out
// Caller   : make_component_list()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//*****************************************************************************
int book_out(config* seq, global* crik)
{
  int i;
  disp(seq,DISP_ALL,"Entering 'book_out'\n");
  dispCmpnt(seq, crik, 1);
  crik->cmpntListOcupidTyp = realloc(crik->cmpntListOcupidTyp, crik->numCmpntTypOcupid * sizeof(int16_t));

  if (DISP) {
      fprintf(seq->dispFile,"\nNumber of Components = %ld\n\n", crik->numCmpnt);
      printf("Occupied Compnent List:");
      for(i = 0 ; i < crik->numCmpntTypOcupid ; i++) {
          printf("%d, ", crik->cmpntListOcupidTyp[i]);
      }  // end for
      printf("\n");
  }

  return 0;
}  // end book_out
  
//*****************************************************************************
// Function : Determine Brace Coordinate Display Format
// Caller   : display_components()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//*****************************************************************************
char* choose_brace_display_format(config* seq)
{
  if      (seq->strLen <  100) return "    {%2d~%-2d{     }%2d~%-2d}";
  else if (seq->strLen < 1000) return "    {%3d~%-3d{     }%3d~%-3d}";
  else                         return "    {%4d~%-4d{     }%4d~%-4d}";

  return 0;
}  // end choose_brace_display_format

//*****************************************************************************
// Function : Determine Display Priority
// Caller   : display_components()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//*****************************************************************************
int set_display_priority(config* seq, int dispCmpltFlag)
{
  int i;
  if (dispCmpltFlag){                          // display only at the final stage in 'book_out' function
                                               disp(seq,DISP_LV3, "                     Cmpnt");
    for(i = 0 ; i < (seq->strLen - 1) ; i++){
 					       disp(seq,DISP_LV3," ");
    }  // end for
					       disp(seq,DISP_LV3, "Coord                         Covariance Pair Matching Status\n");
					       disp(seq,DISP_LV3, "Cmpnt List:  Type   0 %s------------------------------------------------------------------------------------------\n", seq->ltr);
    return DISP_LV3;
  } else {                                     // display each time 'brace_is_formed' is true
    return DISP_LV4;   
  }  // end if
  return 0;
}  // end set_display_priority

//*****************************************************************************
// Function : Display Component List
// Caller   : display_components()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//*****************************************************************************
void display_component_list(config* seq, global* crik, int dispCmpltFlag, int dispPrio, char* braceBound)
{                                                                               disp(seq,DISP_ALL,"Entering display_component_listc\n");
  knob* LLCursr;
  int16_t i, j, k, m;
  int16_t flagUB = seq->numCovari + seq->numV1Pairng;
  int8_t  mismatchFlag = 0;                                                     // if true, that means the spot should mark as mismatch '_', not parenthesis

  for(i = 0 ; i < crik->numCmpntTyp ; i++){                                     // display the dots and the parentheses
    LLCursr = crik->cmpntList[i].knob;
    
    if (i && LLCursr){                                                          disp(seq,dispPrio,"            Type %3d ", i);
                                    					        disp(seq,dispPrio, "%s------------------------------------------------------------------------------------------\n", seq->ltr);
    }  // end if

    while(LLCursr){
                                                                                disp(seq,dispPrio,"                     ");
      for(j = 0                          ; j <  LLCursr->opnBrsOutIndx ; j++){  disp(seq,dispPrio,".");         }
      for(j = LLCursr->opnBrsOutIndx     ; j <= LLCursr->opnBrsInnIndx ; j++){
	for(m = 0 ; m < seq->maxNumMismatch ; m++){                             // very inefficient way of watching for mismatch location
	  if (LLCursr->mismatchFlag[m][0] == j){                                disp(seq,dispPrio,"_"); // display mismatch spots by '_'
	    mismatchFlag = 1;
	  }  // end if
	}    // end for	    
	if (mismatchFlag){
	  mismatchFlag = 0;
	  continue;
	} else {
                                                                disp(seq,dispPrio,"("); // display '(' when it's sure that mismatch isn't at this point
	}  // end if
      }      // end for
      for(j = LLCursr->opnBrsInnIndx + 1 ; j <  LLCursr->closeBrsInnIndx ; j++){  disp(seq,dispPrio,".");     }
      for(j = LLCursr->closeBrsInnIndx     ; j <= LLCursr->closeBrsOutIndx ; j++){
	for(m = 0 ; m < seq->maxNumMismatch ; m++){                             // very inefficient way of watching for mismatch location
	  if (LLCursr->mismatchFlag[m][1] == j){                                disp(seq,dispPrio,"_"); // display mismatch spots by '_'
	    mismatchFlag = 1;
	  }  // end if
	}    // end for	    
	if (mismatchFlag){
	  mismatchFlag = 0;
	  continue;
	} else {                                                                disp(seq,dispPrio,")"); // display ')' when it's sure that mismatch isn't at this point
	}  // end if
      }      // end for
      for(j = LLCursr->closeBrsOutIndx + 1 ; j <  seq->strLen            ; j++){  disp(seq,dispPrio,".");         }
                                                                                disp(seq,dispPrio, braceBound,
										     LLCursr->opnBrsOutIndx,
										     LLCursr->opnBrsInnIndx, 
										     LLCursr->closeBrsInnIndx, 
										     LLCursr->closeBrsOutIndx  );
    if(seq->numCovari + seq->numV1Pairng){  		         	        disp(seq,dispPrio,"     -     [ ");
	for(k = 0 ; k < flagUB ; k++){                                          disp(seq,dispPrio,"%d ",LLCursr->mustPairFlag[k]);
	}  // end for
                                                                                disp(seq,dispPrio,"]");
      } else {                                                                  disp(seq,dispPrio,"                ");
      }    // end if

      if(seq->maxNumMismatch){				         		disp(seq,dispPrio,"     -     [ ");
	for(k = 0 ; k < seq->maxNumMismatch ; k++){                             disp(seq,dispPrio,"%d,%d ",LLCursr->mismatchFlag[k][0], LLCursr->mismatchFlag[k][1]);
	}  // end for
                                                                                disp(seq,dispPrio,"]");
      }    // end if  
         	                                                                disp(seq,dispPrio,"\n");
	

										disp(seq,DISP_ALL,"Flag = %d\n", dispCmpltFlag);
      LLCursr = LLCursr->cmpntListNext;
    }  // end while
  }    // end for

  if(braceBound) dispPrio = dispCmpltFlag;                                      // completely useless dummy operation, which won't execute under report or profile mode
}  // end display_component_list

//*****************************************************************************
// Function : Display Components
// Caller   : recur_stage_3_brace_processng()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//*****************************************************************************
void display_components(config* seq, global* crik, int dispCmpltFlag)
{
  int   dispPrio   = set_display_priority(seq, dispCmpltFlag);
  char* braceBound = choose_brace_display_format(seq);

  display_component_list(seq, crik, dispCmpltFlag, dispPrio, braceBound);
}  // end display_components

//*****************************************************************************
// Function : Initialize New Node
// Caller   : recur_stage_3_brace_processng()
// Purpose  : insert the new node to the top of the linked list of component list
// Input    : 
// Return   : none
// Display  : none
//*****************************************************************************
int hook_new_node_on_component_list(config* seq __unused, global* crik, knob* newNode)
{                                                                              disp(seq,DISP_ALL,"Enterng 'hook_new_node_on_component_list'\n");

	if(crik->cmpntList[crik->opnBrsOutIndx].knob){
      newNode->cmpntListNext = crik->cmpntList[crik->opnBrsOutIndx].knob;      // for not-first card inserting
      crik->cmpntList[crik->opnBrsOutIndx].knob = newNode;
    } else {
      crik->cmpntList[crik->opnBrsOutIndx].knob = newNode;                    // for first card inserting
      crik->cmpntListOcupidTyp[crik->numCmpntTypOcupid] = crik->opnBrsOutIndx;  // take note that this cmpnt type is used (occupied)
      crik->numCmpntTypOcupid++;
    }  // end if

    crik->cmpntList[crik->opnBrsOutIndx].cnt++;
  return 0;
}  // end hook_new_node_on_component_list

//*****************************************************************************
// Function : Initialize Mismatch
// Caller   : recur_stage_3_brace_processng()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//*****************************************************************************
int initialize_mismatch(config* seq, knob* newNode)
{                                                                 disp(seq,DISP_ALL,"Enterng 'initialize_mismatch'\n");
  int i;

  if(seq->maxNumMismatch){
    newNode->mismatchFlag = malloc(sizeof(newNode->mismatchFlag) * seq->maxNumMismatch);
    for(i = 0 ; i < seq->maxNumMismatch ; i++){
      newNode->mismatchFlag[i] = malloc(sizeof(int16_t) * 2);     disp(seq,DISP_ALL,"(%d,%d) ", newNode->mismatchFlag[i][0], newNode->mismatchFlag[i][1]);
      newNode->mismatchFlag[i][0] = -1;                           // shouldn't initialize as 0, since 0 maybe mis-interpreted as the location of '0', 
      newNode->mismatchFlag[i][1] = -1;                           // and mis-place a '_' at the location of 0
    }  // end for
                                                                  disp(seq,DISP_ALL, "\n");
  }    // end if

  return 0;
}  // end initialize_mismatch

//*****************************************************************************
// Function : Is Base Pair Not Formed?
// Caller   : brace_is_formed()
// Purpose  : look up the parenthesis look-up table to check if a given pair satisfied WC or GU pairing
// Input    : seq
//          : crik
// Return   : 1 - not paired
//          : 0 - paired
// Display  : 
//*****************************************************************************
int base_pair_NOT_formed(config* seq, global* crik)
{ 
  return !seq->parenLukUpTable[crik->opnParenIndx][crik->closeParenIndx];
}

//*****************************************************************************
// Function : Is Brace Formed without Similar Helix Reduction
// Caller   : recur_stage_3_brace_processng()
// Purpose  : 
// Input    : 
// Return   : 1:yes, 0:no
// Display  : 
//*****************************************************************************
int brace_is_formed(config* seq, global* crik, knob* newNode, int16_t* mustPairFlag)
{                                                                              disp(seq,DISP_ALL,"Entering 'brace_is_formed'\n");
  int8_t noBrsFlag     = 0;
  int8_t cntOfMismatch = 0;

  crik->opnBrsCursr = 0;                                                       // open brace cursor : used to probe each nucleotide
  crik->closeBrsCursr = 0;                                                       // close  "      "   :        "        "        "

  while(is_worth_exploring(crik, noBrsFlag)){
    crik->opnParenIndx = crik->opnBrsOutIndx + crik->opnBrsCursr;
    crik->closeParenIndx = crik->closeBrsOutIndx - crik->closeBrsCursr;
								
    if(base_pair_NOT_formed(seq, crik)){                                    disp(seq,DISP_ALL,"(%d,%d): (%c,%c) pair not matched, bye...\n", crik->opnParenIndx, crik->closeParenIndx, seq->ltr[crik->opnParenIndx], seq->ltr[crik->closeParenIndx]);
      if(mismatch_causes_lonely_pair(seq, crik, newNode, cntOfMismatch)) return 2;
      if(is_s1(seq, crik->opnParenIndx, crik->closeParenIndx)) return 0;
      if(mismatch_by_cheating(seq, crik->opnParenIndx, crik->closeParenIndx)) return 2;   // mismatch not allowed at either ends of a helix
      if (cntOfMismatch < seq->maxNumMismatch){                                // exploring mismatch options
	if(mismatch_interrupt_covariance(seq,crik)) return 2;
        newNode->mismatchFlag[cntOfMismatch][0] = crik->opnParenIndx;
	newNode->mismatchFlag[cntOfMismatch][1] = crik->closeParenIndx;
	cntOfMismatch++;                                                       disp(seq,DISP_ALL,"1 mismatch created: (%d,%d)\n", crik->opnParenIndx, crik->closeParenIndx);
      } else
	return 0;                                                              // too many mismatches, so just report 0!
    } else {                                                                   disp(seq,DISP_ALL,"(%d,%d): (%c,%c) pair matched, keep going...\n", crik->opnParenIndx, crik->closeParenIndx, seq->ltr[crik->opnParenIndx], seq->ltr[crik->closeParenIndx]  );
      update_v1_covari_flag(seq, crik, mustPairFlag);                          // this base pair formed, therefore it's time to update the must pair flags
    }  // end outer if                                                        
    
    crik->opnBrsCursr++;
    crik->closeBrsCursr++;
  }    // end while

                                                                               disp(seq,DISP_LV4,"      !!!!!!!Brace Formed: (%d~%d) (%d~%d)!!!!!!!\n", crik->opnBrsOutIndx, crik->opnBrsInnIndx, crik->closeBrsInnIndx, crik->closeBrsOutIndx);
  return 1;                                                                    // brace forming is confirmed here, if the program can pass all tests above and reach this point
}  // end brace_is_formed

//*****************************************************************************
// Function : Initialize New Node
// Caller   : recur_stage_3_brace_processng()
// Purpose  : set memory space and initial values
// Input    : 
// Return   : none
// Display  : none
//*****************************************************************************
int initialize_new_node(config* seq __unused, global* crik, knob* newNode, int16_t* mustPairFlag)
{                                                    disp(seq,DISP_ALL,"Enterng 'initialize_new_node'\n");
  newNode->opnBrsOutIndx     = crik->opnBrsOutIndx;
  newNode->opnBrsInnIndx     = crik->opnBrsInnIndx;
  newNode->closeBrsInnIndx     = crik->closeBrsInnIndx;
  newNode->closeBrsOutIndx     = crik->closeBrsOutIndx;
  newNode->mustPairFlag      = mustPairFlag;
  newNode->lvlOfRecur        = -1;
  newNode->hlixBranchngIndx1 = -1;
  newNode->intrvlCntr        = -1;
  newNode->bundleCntr        = -1;
  newNode->outsideIntrvlLB   = -1;
  newNode->outsideIntrvlUB   = -1;
  newNode->specialRstoFlag   =  0;
  newNode->rstoOnQFlag       =  0;
  newNode->bundleFlag        =  0;
  newNode->bundleListNext     = NULL;                 disp(seq,DISP_LV4,"new card inserted: ");
  newNode->jumpTreeNext      = NULL;                 disp(seq,DISP_LV4,"{%2d~%-2d{     }%2d~%-2d}\n", newNode->opnBrsOutIndx, newNode->opnBrsInnIndx, newNode->closeBrsInnIndx, newNode->closeBrsOutIndx);
  crik->numCmpnt++;

  return 0;
}  // end initialize_new_node

//*****************************************************************************
// Function : Is Cheating?
// Caller   : brace_is_formed()
// Purpose  : check if there's cheating, which means, a pair is marked mismatch due to S1 or CM, but actually they are qualified to pair as WC or GU
// Input    : seq
//          : rowI
//          : colJ
// Return   : 1 (yes), 0 (no)
// Display  : none
//*****************************************************************************
int mismatch_by_cheating(config* seq, int16_t rowI, int16_t colJ)
{
  char rowChar = seq->ltr[rowI];
  char colChar = seq->ltr[colJ];

  if(seq->noGU == FALSE){                                // check regular case by chkng the type 'char'
    if(   ((rowChar == 'A') && (colChar == 'U'))         // ||
       || ((rowChar == 'C') && (colChar == 'G'))         // ||
       || ((rowChar == 'G') && (colChar == 'C'))         // ||
       || ((rowChar == 'G') && (colChar == 'U'))         // ||
       || ((rowChar == 'U') && (colChar == 'A'))         // ||
       || ((rowChar == 'U') && (colChar == 'G')) ){      disp(seq,DISP_ALL,"Cheating potential captured!\n");
      return 1;                                          // \/ cheating!
    }  // end inner if 1
  }else{                                                 // check no GU pair case by chkng the type 'char'
    if(   ((rowChar == 'A') && (colChar == 'U'))         // ||
       || ((rowChar == 'C') && (colChar == 'G'))         // ||
       || ((rowChar == 'G') && (colChar == 'C'))         // ||
       || ((rowChar == 'U') && (colChar == 'A')) ){      disp(seq,DISP_ALL,"Cheating potential captured!\n");
      return 1;                                          // \/ cheating!
    }  // end inner if 2
  }    // end outer if

  return 0;                                              // save, ok to put as mismatch, no cheating concern
}  // end mismatch_by_cheating

//*****************************************************************************
// Function : Is Lonely Pair?
// Caller   : brace_is_formed()
// Purpose  : check if certain mismatch nucleotide will be isolated such that it looks like lonely pair
// Input    : seq
//          : crik
// Return   : 1 (yes), 0 (no)
// Display  : 
//*****************************************************************************
int mismatch_causes_lonely_pair(config* seq, global* crik, knob* newNode, int16_t cntOfMismatch)
{
  if(crik->opnBrsCursr == 1) return 1;                                                  // causing 5' end lonely pair
  if(cntOfMismatch)
    if(newNode->mismatchFlag[cntOfMismatch - 1][0] == crik->opnParenIndx - 2) return 1; // causing inner lonely pair
  if(crik->closeBrsCursr == crik->closeBrsWidth - 2 && seq->parenLukUpTable[crik->opnParenIndx + 1][crik->closeParenIndx - 1]) return 1; // causing 3' end lonely pair

  return 0;
}  // end mismatch_causes_lonely_pair

//*****************************************************************************
// Function : Is Mismatch on Covariance Spots?
// Caller   : brace_is_formed()
// Purpose  : check if mismatch nucleotides happen to be covariance nucleotides, which is not allowed
// Input    : seq
//          : crik
// Return   : 1 (yes), 0 (no)
// Display  : 
//*****************************************************************************
int mismatch_interrupt_covariance(config* seq, global* crik)
{
  int8_t i;

  for(i = 0 ; i < seq->numCovari ; i++)
    if( (seq->coVari[i][0] == crik->opnParenIndx) ||
	(seq->coVari[i][0] == crik->closeParenIndx) ||
	(seq->coVari[i][1] == crik->opnParenIndx) ||
	(seq->coVari[i][1] == crik->closeParenIndx)   )  return 1;  // oops, mismatch happened on covariance spot
  
  return 0;                                                       // safe, mismatch spot has nothing to do with covariance
}  // end mismatch_interrupt_covariance

//*****************************************************************************
// Function : Is Non-terminal Mismatch?
// Caller   : brace_is_formed()
// Purpose  : check if mismatch nucleotides happen at nonterminal spots of the nucleotides, which is not allowed
// Input    : crik
//          : cntOfMismatch
// Return   : 1 (yes), 0 (no)
// Display  : 
//*****************************************************************************
int NOT_terminal_mismatch(config* seq, global* crik, int16_t cntOfMismatch)
{
  return !( (crik->opnBrsCursr == 0) || (cntOfMismatch == 1 && crik->opnBrsCursr == 1) || (cntOfMismatch == 2 && crik->opnBrsCursr == 2) || (crik->opnBrsCursr == crik->opnBrsWidth - 1) || (crik->closeBrsCursr == crik->closeBrsWidth - 1) || (crik->opnBrsCursr == crik->opnBrsWidth - 2 && !seq->parenLukUpTable[crik->opnParenIndx + 1][crik->closeParenIndx - 1]) || (crik->closeBrsCursr == crik->closeBrsWidth - 2 && !seq->parenLukUpTable[crik->opnParenIndx + 1][crik->closeParenIndx - 1]));
}  // end NOT_terminal_mismatch

//*****************************************************************************
// Function : Is S1?
// Caller   : brace_is_formed()
// Purpose  : check if there's cheating, which means, a pair is marked mismatch due to S1 or CM, but actually they are qualified to pair as WC or GU
// Input    : seq
//          : rowI
//          : colJ
// Return   : 1 (yes), 0 (no)
// Display  : 
//*****************************************************************************
int is_s1(config* seq, int16_t rowI, int16_t colJ)
{
  int16_t i;

  for(i = 0 ; i < seq->numS1Pairng ; i++)
    if( (seq->s1Pairng[i] == rowI) ||
	(seq->s1Pairng[i] == colJ)   ) return 1;

  return 0;
}  // is_s1

//*****************************************************************************
// Function : Is Terminal Mismatch?
// Caller   : brace_is_formed()
// Purpose  : check if the mismatch happened at the 5' or 3' very end or hairpin loop of the helix
// Input    : crik
// Return   : 1 (yes), 0 (no)
// Display  : 
//*****************************************************************************
int is_terminal_mismatch(global* crik)
{
  return !crik->opnBrsCursr                           ||  // terminal mismatch at nt #0
         (crik->opnBrsCursr == crik->opnBrsWidth - 1) ||  // terminal mismatch at the edge of the 5' side of the hairpin loop
         (crik->closeBrsCursr == crik->closeBrsWidth - 1);    // terminal mismatch at the edge of the 3' side of the hairpin loop
}  // end is_terminal_mismatch

//*****************************************************************************
// Function : Is Worth Exploring?
// Caller   : brace_is_formed()
// Purpose  : check if the cursr already exceed the exploration limit or if this brace already be prounced hopeless (noBrsFlag)
// Input    : crik
//          : noBrsFlag
// Return   : 1 (yes), 0 (no)
// Display  : 
//*****************************************************************************
int is_worth_exploring(global* crik, int8_t noBrsFlag)
{
  return (crik->opnBrsCursr < crik->opnBrsWidth) && (crik->closeBrsCursr < crik->closeBrsWidth) && (!noBrsFlag);
}

//*****************************************************************************
// Function : Recursion Stage 1: Brace Locating
// Caller   : make_component_list()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//*****************************************************************************
int recur_stage_1_brace_locatng(config* seq, global* crik)
{                                                                                                        disp(seq,DISP_ALL,"Entering 'recur_stage_1'\n");
  int16_t bakUpLB;

  for(crik->opnBrsOutIndx = 0 ; crik->opnBrsOutIndx < crik->opnBrsStop ; crik->opnBrsOutIndx++){                                                     
    bakUpLB = crik->opnBrsOutIndx + (seq->minLenOfHlix - 1) * 2 + seq->minPairngDist;                    // closest location closeOutChar moves to opnOutChar
    for(crik->closeBrsOutIndx = seq->strLen - 1 ; crik->closeBrsOutIndx > bakUpLB ; crik->closeBrsOutIndx--){  disp(seq,DISP_ALL,"------------------------------------------------------------In brace loop (%d,%d):\n", crik->opnBrsOutIndx, crik->closeBrsOutIndx);
      recur_stage_2_brace_sizng(seq, crik);
    }  // end inner for
  }    // end outer for

  return 0;
}  // end recur_stage_1_brace_locatng

//*****************************************************************************
// Function : Recursion Stage 2: Brace Sizing
// Caller   : recur_stage_1_brace_locatng()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//*****************************************************************************
int recur_stage_2_brace_sizng(config* seq, global* crik)
{                                                                                                          disp(seq,DISP_ALL,"Entering 'recur_stage_2'\n");
  int16_t opnBrsWidthUB;
  int16_t closeBrsWidthLB;
  int16_t closeBrsWidthUB;
  int     sizngBypassFlag = 0;                                                                           // skip the rest of hopeless sizing attempts
 
  opnBrsWidthUB = assign_window_width(seq, crik, 0);

  for(crik->opnBrsWidth = seq->minLenOfHlix ; crik->opnBrsWidth <= opnBrsWidthUB ; crik->opnBrsWidth++){ // open brace sizing
    crik->opnBrsInnIndx = crik->opnBrsOutIndx + crik->opnBrsWidth - 1;                                   // location of open brace inner edge
    closeBrsWidthLB = assign_window_width(seq, crik, 1);
    closeBrsWidthUB = assign_window_width(seq, crik, 2);
    for(crik->closeBrsWidth = closeBrsWidthLB; crik->closeBrsWidth <= closeBrsWidthUB ; crik->closeBrsWidth++){    // close brace sizing
      crik->closeBrsInnIndx = crik->closeBrsOutIndx - crik->closeBrsWidth + 1;                                 // location of close brace inner edge

      sizngBypassFlag = recur_stage_3_brace_processng(seq, crik);                                        disp(seq,DISP_ALL,"Test on : (%d~%d) (%d~%d)\n", crik->opnBrsOutIndx, crik->opnBrsInnIndx, crik->closeBrsInnIndx, crik->closeBrsOutIndx);
      if(sizngBypassFlag) return 0;                                                                      // skip all the rest of the check
    }    // end inner for
  }      // end outer for

  return 0;
}  // end recur_stage_2_brace_sizng

//*****************************************************************************
// Function : Recursion Stage 3: Brace Processing
// Caller   : recur_stage_2_brace_sizng()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//*****************************************************************************
int recur_stage_3_brace_processng(config* seq, global* crik)
{                                                                          disp(seq,DISP_ALL,"Entering 'recur_stage_3_brace_processng'\n");
  knob*  newNode      = calloc(1, sizeof(knob));
  int16_t* mustPairFlag = seq->constraintActive ? calloc((seq->numCovari + seq->numV1Pairng), sizeof(int16_t)) : NULL;
  int16_t  i;

  initialize_mismatch(seq, newNode);

  if( (crik->opnBrsWidth < 0) || (crik->closeBrsWidth < 0) ) {               disp(seq,DISP_LV1,"Error! negative brace width detected!!!\n\n");
    return 0;
  }  // end if 1

  int result = brace_is_formed(seq, crik, newNode, mustPairFlag);
  if (result == 2) {
	  for(i = 0 ; i < seq->maxNumMismatch ; i++)
	        free(newNode->mismatchFlag[i]);                                      // clean constraints
	  free(newNode->mismatchFlag);
	  free(newNode);
	  free(mustPairFlag);
	  return 0;
  } else if (result == 1){           // new component qualified to form!
    initialize_new_node(seq, crik, newNode, mustPairFlag);
    hook_new_node_on_component_list(seq, crik, newNode);
    return 0;
  } else {                                                                 // component disqualified to form :(
    for(i = 0 ; i < seq->maxNumMismatch ; i++)
      free(newNode->mismatchFlag[i]);                                      // clean constraints
    free(newNode->mismatchFlag);                                          
    free(newNode);                                                                
    free(mustPairFlag);
  } // end if 2

  return 1;                                                                // 'sizngBypassFlag' will be assigned 1, which will skip all rest of sizng process
}  // end recur_stage_3_brace_processng

//*****************************************************************************
// Function : Update V1 and Covariance Flags
// Caller   : brace_is_formed()
// Purpose  : update the contents of mustPairFlag
// Input    : seq
//          : crik
//          : mustPairFlag
// Return   : 
// Display  : 
//*****************************************************************************
int update_v1_covari_flag(config* seq, global* crik, int16_t* mustPairFlag)
{                                                                   disp(seq,DISP_ALL,"Enterng 'update_v1_covari_flag'\n");
  int8_t mustPairFlagCntr = seq->numCovari + seq->numV1Pairng;
  int8_t v1FlagCntr = seq->numV1Pairng;

  while(mustPairFlagCntr > seq->numCovari){                       // || V1 PAIRING : update flag status
    if( (seq->v1Pairng[v1FlagCntr-1] == crik->opnParenIndx) ||      // || unlike covari, the v1 may be matched by either opn paren or close paren
	(seq->v1Pairng[v1FlagCntr-1] == crik->closeParenIndx)    )    // ||
      mustPairFlag[mustPairFlagCntr-1] = 1;                         // ||
    mustPairFlagCntr--;                                             // ||
    v1FlagCntr--;                                                   // ||
  }  // end while 1                                                 // \/

  while(mustPairFlagCntr > 0){                                      // || COVARIANCE PAIRING : update flag status
    if(seq->coVari[mustPairFlagCntr-1][0] == crik->opnParenIndx)    // ||
      mustPairFlag[mustPairFlagCntr-1] = 1;                         // ||
    mustPairFlagCntr--;                                             // ||
  }  // end while 2                                                 // \/
  
  return 0;
}  // end update_v1_covari_flag

