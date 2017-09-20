/*****************************************************************************************************************************
 *****************************************************************************************************************************
 *****************************************************************************************************************************
      Program : Swellix (Make Parenthesis Look-up Table) - c file
      Purpose : Make Parenthesis Look-up Table
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
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "paren_lookup_table.h"

//*****************************************************************************
// Function : Fill In Parenthesis Look-up Table Row I
// Caller   : make_parenthesis_look_up_table()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//*****************************************************************************
int fill_in_parenthesis_look_up_table_row_i(config* seq, int16_t rowI)
{
  int16_t pairngLB = rowI + seq->minPairngDist + 1;                                         // determine the length of empty space
  int16_t pairngUB = ( (seq->maxPairngDist + rowI + 2) > seq->strLen) ? seq->strLen : (seq->maxPairngDist + rowI + 2);
  int16_t colJ;
  int16_t blankFillngIndx;

  seq->parenLukUpTable[rowI] = malloc(sizeof(*seq->parenLukUpTable[rowI]) * seq->strLen);   // set up all columns for the i'th row
  disp(seq,DISP_LV3,"               %d%c", rowI, seq->ltr[rowI]);
  for(blankFillngIndx = 0 ; blankFillngIndx < pairngLB ; blankFillngIndx++){
    disp(seq,DISP_LV3," "); // fill in empty space
  }  // end for 1

  for(colJ = pairngLB ; colJ < pairngUB ; colJ++){                                          // fill in 0: can't pair, or 1: can pair
    seq->parenLukUpTable[rowI][colJ] = paren_is_formed(seq, rowI, colJ);
    disp(seq,DISP_LV3,"%d", seq->parenLukUpTable[rowI][colJ]);
  }  // end for 2

  while(colJ < seq->strLen){                                                                // fill in the trailing zeros
    seq->parenLukUpTable[rowI][colJ] = 0;
    disp(seq,DISP_LV3,"0");
    colJ++;
  }  // end for 3

  disp(seq,DISP_LV3,"\n");
  return 0;
}  // end fill_in_parenthesis_look_up_table_row_i

//*****************************************************************************
// Function : Is Parenthesis Formed?
// Caller   : fill_in_parenthesis_look_up_table_row_i()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//*****************************************************************************
int paren_is_formed(config* seq, int16_t rowI, int16_t colJ)
{
  char opnParenChar = seq->ltr[rowI];
  char closeParenChar = seq->ltr[colJ];
  char opnCoVariMark;
  char closeCoVariMark;
  int  a, b, c, d;

  if(seq->numCovari){                                                                 // check basic covariation violation by chkng the type 'int'
    int i = seq->numCovari;                                                           // ||
                                                                                      // ||
    while(i){                                                                         // ||
      opnCoVariMark = seq->coVari[i-1][0];                                            // ||
      closeCoVariMark = seq->coVari[i-1][1];                                            // ||
                                                                                      // ||
      if      (   ( rowI == opnCoVariMark) && (colJ == closeCoVariMark)  ) return 1;    // || the paren pair matches covari pair
      else if (   ((rowI == opnCoVariMark) && (colJ != closeCoVariMark))                // || heads match, tails don't
	       || ((rowI != opnCoVariMark) && (colJ == closeCoVariMark))                // || heads don't match, tails does
   	       || (rowI == closeCoVariMark)                                             // || head matches tail
 	       || (colJ == opnCoVariMark)                              ) return 0;    // \/ tail matches head
                                                                                   
    if(seq->sudona == FALSE){                                                         // check axiliary covariation violation (about pseudoknot) by chkng the type 'int'
      a = rowI - opnCoVariMark;                                                       // ||
      b = rowI - closeCoVariMark;                                                       // ||
      c = colJ - opnCoVariMark;                                                       // ||
      d = colJ - closeCoVariMark;                                                       // ||
      if(   ((a*b > 0) && (c*d < 0))                                                  // || rowI is outside coVariMark window, while colJ is inside coVariMark window
         || ((a*b < 0) && (c*d > 0)) ) return 0;                                      // || rowI is inside coVariMark window, while colJ is outside coVariMark window
    }  // end if                                                                      // \/
      i--;                                                               
    }  // end while                                                      
  }    // end if                                                         

  if(seq->numS1Pairng){                                                               // check basic S1 pairing violation by chkng the type 'int'
    int i = seq->numS1Pairng - 1;                                                     // ||
                                                                                      // ||
    while(i >= 0){                                                                    // ||
      if( (rowI == seq->s1Pairng[i]) || (colJ == seq->s1Pairng[i]) ) return 0;        // || either rowI or colJ is on the S1 pairng list, which is prohibited to form a pair
      i--;                                                                            // ||
    }  // end while                                                                   // ||
  }    // end if                                                                      // \/
    

  if(seq->noGU == FALSE){                                                             // check regular case by chkng the type 'char'
    if(   ((opnParenChar == 'A') && (closeParenChar == 'U'))                            // ||
       || ((opnParenChar == 'C') && (closeParenChar == 'G'))                            // ||
       || ((opnParenChar == 'G') && (closeParenChar == 'C'))                            // ||
       || ((opnParenChar == 'G') && (closeParenChar == 'U'))                            // ||
       || ((opnParenChar == 'U') && (closeParenChar == 'A'))                            // ||
       || ((opnParenChar == 'U') && (closeParenChar == 'G')) ) return 1;                // \/
  }else{                                                                              // check no GU pair case by chkng the type 'char'
    if(   ((opnParenChar == 'A') && (closeParenChar == 'U'))                            // ||
       || ((opnParenChar == 'C') && (closeParenChar == 'G'))                            // ||
       || ((opnParenChar == 'G') && (closeParenChar == 'C'))                            // ||
       || ((opnParenChar == 'U') && (closeParenChar == 'A')) ) return 1;                // \/
  }  // end if

  return 0;
}  // end isParenFormed

