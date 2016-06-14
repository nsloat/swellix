/*****************************************************************************************************************************
 *****************************************************************************************************************************
 *****************************************************************************************************************************
      Program : Swellix (Unit Test) - c file
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
#include "main.h"
#include "unit_tests.h"

//***********************************************************************
// Function : Test - Parentheses Balance
// Caller   : display_structures()
// Purpose  : make sure number of open parentheses matches close ones
// Input    : seq  - sequence
//          : crik - crick
// Return   : none
// Display  : Error message, if necessary
//***********************************************************************
int test_parentheses_balance(config* seq, global* crik)
{                                            disp(seq,DISP_ALL,"Enterng Test 1: 'testParenBal'\n");
  int16_t parenUnbalanceIndx = 0;
  int16_t i;
  
  for(i = 0 ; i < seq->strLen ; i++){
    if(seq->dotNParen[i] == '(') parenUnbalanceIndx++;
    if(seq->dotNParen[i] == ')') parenUnbalanceIndx--;
  }  // end for

  if(parenUnbalanceIndx != 0){
    fprintf(stdout, "$$$$$$$$$$$$$$$$$$$$$$$$$  Unbalanced Parentheses found: #%ld %s\n", crik->numStru, seq->dotNParen);
    crik->test1ErrTotal++;
    return 0;
  }  // end if

  return 0;
}  // end test_parentheses_balance
				      
//***********************************************************************
// Function : Test - Loop Size
// Caller   : display_structures()
// Purpose  : make sure number of nucleotides inside any helix is large enough
// Input    : seq  - sequence
//          : crik - crick
// Return   : none
// Display  : Error message, if necessary
//***********************************************************************
int test_loop_size(config* seq, global* crik)
{                                            disp(seq,DISP_ALL,"Enterng Test 2: 'testLupSze'\n");
  int16_t lupSze     = 0;
  int16_t opnBrsFlag = 0;
  int16_t inLupFlag  = 0;
  int16_t i;
  
  for(i = 0 ; i < seq->strLen ; i++){
    if (seq->dotNParen[i] == '('){
      opnBrsFlag = 1;
      inLupFlag  = 0;
      lupSze     = 0;
    } else if ( ((seq->dotNParen[i] == '.') || (seq->dotNParen[i] == '#') || (seq->dotNParen[i] == '/')) && opnBrsFlag ){
      opnBrsFlag = 0;
      inLupFlag  = 1;
      lupSze     = 1;
    } else if ( ((seq->dotNParen[i] == '.') || (seq->dotNParen[i] == '#') || (seq->dotNParen[i] == '/')) && inLupFlag ) {  // . is normal dot, # is bundle mark
      lupSze++;
    } else if ( (seq->dotNParen[i] == ')') && inLupFlag ) {
      if(lupSze < seq->minPairngDist){
	fprintf(stdout, "$$$$$$$$$$$$$$$$$$$$$$$$$  Insufficient loop size 1 found: #%ld %s\n", crik->numStru, seq->dotNParen);
	crik->test2ErrTotal++;
	return 0;
      }  // end inner if
      opnBrsFlag = 0;
      inLupFlag  = 0;
      lupSze     = 0;
    } else if ( (seq->dotNParen[i] == ')') && (seq->dotNParen[i-1] == '(') ) {
      fprintf(stdout, "$$$$$$$$$$$$$$$$$$$$$$$$$  Insufficient loop size 2 found: #%ld %s\n", crik->numStru, seq->dotNParen);
      crik->test2ErrTotal++;
      return 0;
    }  // end if
  }    // end for
      
  return 0;
}  // end test_parentheses_balance
				      
//***********************************************************************
// Function : Test - Helix Overlap
// Caller   : jump_stage_3_verify_hlix()
// Purpose  : detect overlap situation between the new helix and any other helix
// Input    : seq  - sequence
//          : crik - crick
// Return   : none
// Display  : Error message, if necessary
//***********************************************************************
int test_helix_overlap(config* seq __unused, global* crik, local* todd, int8_t locationMark)
{                                                                                                 disp(seq,DISP_ALL,"Enterng Test 3: 'test_helix_overlap'\n");
  knob* newHlix   = todd->cmpntLLCursr;
  knob* hlixCursr = locationMark ? crik->hlixInStru : crik->hlixInStru->jumpTreeNext;
  int16_t d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15, d16;

  while(hlixCursr){
    d1  = newHlix->opnBrsOutIndx - hlixCursr->opnBrsOutIndx;
    d2  = newHlix->opnBrsOutIndx - hlixCursr->opnBrsInnIndx;
    d3  = newHlix->opnBrsOutIndx - hlixCursr->closeBrsInnIndx;
    d4  = newHlix->opnBrsOutIndx - hlixCursr->closeBrsOutIndx;
    d5  = newHlix->opnBrsInnIndx - hlixCursr->opnBrsOutIndx;
    d6  = newHlix->opnBrsInnIndx - hlixCursr->opnBrsInnIndx;
    d7  = newHlix->opnBrsInnIndx - hlixCursr->closeBrsInnIndx;
    d8  = newHlix->opnBrsInnIndx - hlixCursr->closeBrsOutIndx;
    d9  = newHlix->closeBrsInnIndx - hlixCursr->opnBrsOutIndx;
    d10 = newHlix->closeBrsInnIndx - hlixCursr->opnBrsInnIndx;
    d11 = newHlix->closeBrsInnIndx - hlixCursr->closeBrsInnIndx;
    d12 = newHlix->closeBrsInnIndx - hlixCursr->closeBrsOutIndx;
    d13 = newHlix->closeBrsOutIndx - hlixCursr->opnBrsOutIndx;
    d14 = newHlix->closeBrsOutIndx - hlixCursr->opnBrsInnIndx;
    d15 = newHlix->closeBrsOutIndx - hlixCursr->closeBrsInnIndx;
    d16 = newHlix->closeBrsOutIndx - hlixCursr->closeBrsOutIndx;

    if (d1 && d2 && d3 && d4 && d5 && d6 && d7 && d8 && d9 && d10 && d11 && d12 && d13 && d14 && d15 && d16){ 
    } else {                                                                                        // identical index present in two different hlix
      fprintf(stdout, "$$$$$$$$$$$$$$$$$$$$$$$$$  hlix overlap type 1 detected after #%ld\n", crik->numStru);
      disp(seq,DISP_LV1,"list: %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
	   d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15, d16);
      crik->test3ErrTotal++;
      return 1;
    }  // end if 1
    if(d2  * d5  < 0){                                                                              // overlap of 2 open  braces
      fprintf(stdout, "$$$$$$$$$$$$$$$$$$$$$$$$$  hlix overlap type 2 detected after #%ld\n", crik->numStru);
      disp(seq,DISP_LV1,"list: %d %d\n", d2, d5);
      crik->test3ErrTotal++;
      return 1;
    }  // end if 2
    if(d4  * d7  < 0){                                                                              // overlap of open braces vs close
      fprintf(stdout, "$$$$$$$$$$$$$$$$$$$$$$$$$  hlix overlap type 3 detected after #%ld\n", crik->numStru);
      disp(seq,DISP_LV1,"list: %d %d\n", d4, d7);
      crik->test3ErrTotal++;
      return 1;
    }  // end if 3
    if(d10 * d13 < 0){                                                                              // overlap of open braces vs close
      fprintf(stdout, "$$$$$$$$$$$$$$$$$$$$$$$$$  hlix overlap type 4 detected after #%ld\n", crik->numStru);
      disp(seq,DISP_LV1,"list: %d %d\n", d10, d13);
      crik->test3ErrTotal++;
      return 1;
    }  // end if 4
    if(d12 * d15 < 0){                                                                              // overlap of close braces
      fprintf(stdout, "$$$$$$$$$$$$$$$$$$$$$$$$$  hlix overlap type 5 detected after #%ld\n", crik->numStru);
      disp(seq,DISP_LV1,"list: %d %d\n", d12, d15);
      crik->test3ErrTotal++;
      return 1;
    }  // end if 5
    if(d1  * d4 * d13 * d16 < 0){                                                                   // pseudoknot
      fprintf(stdout, "$$$$$$$$$$$$$$$$$$$$$$$$$  hlix overlap type 6 detected after #%ld\n", crik->numStru);
      disp(seq,DISP_LV1,"list: %d %d %d %d\n", d1, d4, d13, d16);
      crik->test3ErrTotal++;
      return 1;
    }  // end if 6

    hlixCursr = hlixCursr->jumpTreeNext;
  }  // end while
                                                                                                    disp(seq,DISP_LV4,"test 3 passed\n"); 
  return 0;
}  // end test_helix_overlap

//***********************************************************************
// Function : Test - Helix Size
// Caller   : display_structures()
// Purpose  : detect helix of insufficient size
// Input    : seq  - sequence
//          : crik - crick
// Return   : none
// Display  : Error message, if necessary
//***********************************************************************
int test_helix_size(config* seq, global* crik)
{                                                                                      disp(seq,DISP_ALL,"Enterng Test 4: 'test_helix_size'\n");
  knob* hlixCursr = crik->hlixInStru;

  while(hlixCursr){
    if(hlixCursr->opnBrsInnIndx - hlixCursr->opnBrsOutIndx + 1 < seq->minLenOfHlix){   // chk the number of nucleotides owned by opn brace
      fprintf(stdout, "$$$$$$$$$$$$$$$$$$$$$$$$$  open  brace not large enough after #%ld\n", crik->numStru);
      crik->test4ErrTotal++;
    }  // end if 1

    if(hlixCursr->closeBrsOutIndx - hlixCursr->closeBrsInnIndx + 1 < seq->minLenOfHlix){   // chk the number of nucleotides owned by close brace    
      fprintf(stdout, "$$$$$$$$$$$$$$$$$$$$$$$$$  close brace not large enough after #%ld\n", crik->numStru);
      crik->test4ErrTotal++;
    }  // end if 2

    hlixCursr = hlixCursr->jumpTreeNext;
  }  // end while
                                                                                       disp(seq,DISP_LV4,"test 4 passed\n"); 
  return 0;
}  // end test_helix_size
