/*****************************************************************************************************************************
 *****************************************************************************************************************************
 *****************************************************************************************************************************
      Program : Swellix (Close Up) - c file
      Purpose : Free dynamically allocated memory spaces
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
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "close_up.h"

//***********************************************************************
// Function : Free Bundle List
// Caller   : close_up()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//***********************************************************************
int free_bundle_list(config* seq, global* crik)
{
  if(seq->algoMode < MODE_INTAB || !seq->bundle) return 0;

  ToL*    tempS;
  knob* tempB;
  int16_t i, j;

  for(i = 0 ; i < seq->strLen ; i++){                                    // || clean the whole 2D eden plan
    for(j = 0 ; j < seq->strLen ; j++){                                  // ||
      if(crik->eden[i][j].dumiNode){
        if(crik->eden[i][j].dumiNode->cmpntListNext) printf("dumiNode at (%d,%d) points to cmpntListNext\n",i,j);//free(crik->eden[i][j].dumiNode->cmpntListNext);
        free(crik->eden[i][j].dumiNode);     // || || clean all the stem node (the whole tree of life) on a particular eden spot
        crik->eden[i][j].dumiNode = NULL;
      }
      while(crik->eden[i][j].stemNode){                                  // || ||
        while(crik->eden[i][j].stemNode->branchNode){                    // || || || clean all branches on a particular stem node
          tempB = crik->eden[i][j].stemNode->branchNode->bundleListNext;  // || || || bundleListNext holds the original link after the jumpTreeNext lost it during the make_jump_tree() operations, therefore it's more reliable during close_up process
          free(crik->eden[i][j].stemNode->branchNode);                   // || || || 
          crik->eden[i][j].stemNode->branchNode = NULL;
          crik->eden[i][j].stemNode->branchNode = tempB;                 // || || ||
        }  // end inner while                                            // || || \/
        if(seq->constraintActive)                                        // || ||
          free(crik->eden[i][j].stemNode->mustPairFlag);                 // || ||
        tempS = crik->eden[i][j].stemNode->next;                         // || ||
        free(crik->eden[i][j].stemNode);                                 // || ||
        crik->eden[i][j].stemNode = NULL;
        crik->eden[i][j].stemNode = tempS;                               // || ||
      }    // end outer while                                            // || \/
    }      // end inner for                                              // ||
    free(crik->eden[i]);                                                 // ||
  }        // end outer for                                              // ||
  free(crik->eden);                                                      // \/

  return 0;
}  // end free_bundle_list

//***********************************************************************
// Function : Free Component List
// Caller   : close_up()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//***********************************************************************
int free_cmpnt_list(config* seq, global* crik)
{
  if(seq->algoMode < MODE_CMPNT) return 0;

  knob* cmpntCursr;
  knob* temp;
  int     i, j;

  for(i = 0 ; i < crik->numCmpntTypOcupid; i++){                                          // || clean the entire edge list
    cmpntCursr = crik->cmpntList[crik->cmpntListOcupidTyp[i]].knob;
//    cmpntCursr = crik->cmpntList[i].knob;                                           // ||
    while(cmpntCursr){                                                               // || || clean entire branch on a particular backbone point
      temp = cmpntCursr->cmpntListNext;                                               // || ||
      for(j = 0 ; j < seq->maxNumMismatch ; j++) free(cmpntCursr->mismatchFlag[j]);  // || || || clean constraints
      //free(cmpntCursr->mismatchFlag);                                                // || || ||
      free(cmpntCursr->mustPairFlag);                                                // || || \/
      if(!cmpntCursr->bundleFlag) { free(cmpntCursr);                                                              // || || 
      }	else {
        cmpntCursr->cmpntListNext = NULL;
      }
      cmpntCursr = temp;                                                             // || ||
    }  // end while                                                                  // || \/
  }    // end for                                                                    // || 
  free(crik->cmpntListOcupidTyp);                                                     // || 
  free(crik->cmpntList);                                                              // \/

  return 0;
}  // end free_cmpnt_list

//***********************************************************************
// Function : Free Constraints
// Caller   : close_up()
// Purpose  : clean up constraints 'seq->coVari', 'seq->v1Pairng', 'seq->s1Pairng', 'seq->chemMod', etc
// Input    : 
// Return   : 
// Display  : 
//***********************************************************************
int free_constraint(config* seq, global* crik)
{
  if(!seq->constraintActive) return 0;       // if the configuration file is not used to introduce constraint, then this function is useless

  int i = seq->numCovari;

  while(i){
    free(seq->coVari[i - 1]);
    i--;
  }  // end while

  free(seq->coVari);
  free(crik->struMustPairFlag);
  free(seq->s1Pairng);
  free(seq->v1Pairng);

  return 0;
}  // end free_covari

//***********************************************************************
// Function : Free Interval Look-up Tables
// Caller   : close_up()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//***********************************************************************
int free_intrvl_luk_up_table(config* seq)
{
  if(seq->algoMode < MODE_INTAB) return 0;

  int16_t i, j, jLB;

  for(i = 0 ; i < seq->strLen ; i++){                         // clean look-up table, row by row
    jLB = i + seq->minPairngDist + seq->minLenOfHlix * 2 - 3;
    if(seq->algoMode > MODE_CMPNT){                      
      for(j = jLB ; j < seq->strLen ; j++){              
        free(seq->intrvlLukUpTable[i][j]);               
      }  // end inner for                                
      free(seq->intrvlLukUpTable[i]);                    
    }    // end if 1                                     
  }      // end outer for                               
  if(seq->algoMode > MODE_CMPNT) 
    free(seq->intrvlLukUpTable);                         

  return 0;
}  // end free_intrvl_luk_up_table

//***********************************************************************
// Function : Free Parenthesis Look-up Tables
// Caller   : close_up()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//***********************************************************************
int free_parenthesis_look_up_table(config* seq)
{
  if(seq->algoMode < MODE_CMPNT) return 0;

  int16_t i;
  for(i = 0 ; i < seq->strLen ; i++){                         // clean look-up table, row by row
    free(seq->parenLukUpTable[i]);                       
  }  // end for
  free(seq->parenLukUpTable);                            
  return 0;
}  // end free_parenthesis_look_up_table

//***********************************************************************
// Function : Free Recycle Bin
// Caller   : close_up()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//***********************************************************************
int free_recycle_bin(config* seq)
{
  local* temp;
  local* toddCursr;

  if(seq->numCovari){
    toddCursr = seq->recycleBin;
    while(toddCursr){
//      free(toddCursr->intrvl2BRsto);
      free(toddCursr->intrvlIns);
      free(toddCursr->intrvlBeh);
      temp = toddCursr;
      toddCursr = toddCursr->next;
      free(temp);
    }  // end while
  }  // end if

  return 0;
}  // end free_recycle_bin

//***********************************************************************
// Function : Free Whole Interval
// Caller   : close_up()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//***********************************************************************
int free_whole_intrvl(global* crik)
{
  knob* temp = crik->interval;

  while(temp){
    crik->interval = temp->jumpTreeNext;
    free(temp);
    temp = crik->interval;
  }  // end while
  return 0;
}  // end free_whole_intrvl

