/*****************************************************************************************************************************
 *****************************************************************************************************************************
 *****************************************************************************************************************************
      Program : Swellix (Make Interval Look-up Table) - c file
      Purpose : Make interval look-up table
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
#include "interval_lookup_table.h"

//*****************************************************************************
// Function : Format the Interval Look-up Table Print Out
// Caller   : make_interval_look_up_table()
// Purpose  : set up the frame work for the table
// Input    : 
// Return   : 
// Display  : 
//*****************************************************************************
int format_the_intrvl_luk_up_table_prt_out(config* seq)
{
  int16_t i;
                                         disp(seq,DISP_LV3,"      ");
  for(i = 0 ; i < seq->strLen ; i++){    disp(seq,DISP_LV3,"%2d    ", i);
  }  // end for 1
                                         disp(seq,DISP_LV3,"\n");
                                         disp(seq,DISP_LV3,"       ");
  for(i = 0 ; i < seq->strLen ; i++){    disp(seq,DISP_LV3,"%c     ", seq->ltr[i]);
  }  // end for 1
                                         disp(seq,DISP_LV3,"\n");

return 0;
}  // end format_the_intrvl_luk_up_table_prt_out

//*****************************************************************************
// Function : Fill In Interval Look-up Table Row I
// Caller   : make_interval_look_up_table()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//*****************************************************************************
int fill_in_intrvl_luk_up_table_row_i(config* seq, global* crik, int16_t rowI)
{                                                                                           disp(seq,DISP_ALL,"entering fill_in_intervl_luk_up_table_row_i'\n");
  int16_t jLB = rowI + seq->minPairngDist + seq->minLenOfHlix * 2 - 3;                      // determine the length of empty space between open & close helix
  int16_t colJ;
  int16_t i;
  int16_t emptySpaceUB = seq->minPairngDist + seq->minLenOfHlix * 2 - 3;

  seq->intrvlLukUpTable[rowI] = malloc(sizeof(int16_t**) * seq->strLen);  // set up all columns for the row i
  disp(seq,DISP_ALL,"<< #%-3d %c >>", rowI, seq->ltr[rowI]);
  disp(seq,DISP_LV3,"%-2d %c", rowI, seq->ltr[rowI]);

  if(jLB >= seq->strLen){
    disp(seq,DISP_ALL,"           -- interval too narrow\n");
    disp(seq,DISP_LV3,"\n");
    return 0;
  }  
  for(i = 0 ; i < emptySpaceUB ; i++){
    disp(seq,DISP_LV3,"      ");
  }
  for(i = 0 ; i < rowI ; i++){
    disp(seq,DISP_LV3,"      ");
  }  // end for 1

  for(colJ = jLB ; colJ < seq->strLen ; colJ++){                                            // fill in the cell with array of intervals
    fill_in_intrvl_luk_up_table_cell_i_j(seq, crik, rowI, colJ);
  }  // end for 2
  disp(seq,DISP_LV3,"\n");
  return 0;
}  // end fill_in_parenthesis_look_up_table_row_i__mplut

//*****************************************************************************
// Function : Fill In Interval Look-up Table Cell (I, J)
// Caller   : fill_in_intrvl_luk_up_table_row_i()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//*****************************************************************************
int fill_in_intrvl_luk_up_table_cell_i_j(config* seq, global* crik, int16_t rowI, int16_t colJ)
{                                                                                          disp(seq,DISP_ALL,"entering fill_in_intervl_luk_up_table_cell'\n");
  int16_t ocupidTypIndx;
  int16_t fittableTypUB = 0;                                                               // keep track of largest cmpntTyp whose UB really fit into i,j
  int16_t cmpntTyp;                                                                        // just a component type

  init_i_j_cell(seq, rowI, colJ);

  for(ocupidTypIndx = 0 ; ocupidTypIndx < crik->numCmpntTypOcupid ; ocupidTypIndx++){
    disp(seq,DISP_ALL,"\n(i,j, ocupidTyp, MAX) = (%d,%d,%d,%d)\n", rowI, colJ, ocupidTypIndx, crik->numCmpntTypOcupid );

    switch(lowerbound_fits_in(seq, crik, rowI, colJ, ocupidTypIndx)){                          // check LB
      case 0:                                                                              // 0: cmpnt lower bound larger than cursrUB
	      return 0;                                                                          
      case 1:                                                                              // 1: cmpnt typ is too small (falls behind), keep looping and see how it goes
	      continue;                            
      case 2:                                                                              // 2: cmpnt typ falls in the interval
	      break;                                                                        
      case 3:                                                                              // 3: fittableTypUB is ideal as upper bound
        seq->intrvlLukUpTable[rowI][colJ][1] = crik->cmpntListOcupidTyp[fittableTypUB];
	      reach_end_of_loop(seq, crik, rowI, colJ, ocupidTypIndx);
	      return 0;
    }  // end switch

    if(upperbound_fits_in(seq, crik, rowI, colJ, ocupidTypIndx)){                              // check UB
      fittableTypUB = ocupidTypIndx;
      cmpntTyp      = crik->cmpntListOcupidTyp[ocupidTypIndx];

      if(is_covarance_pair(seq, cmpntTyp)){                                                   // that means crik->cmpntListOcupidTyp[ocupidTypIndx] is on the covariance pair
	      mark_reading_end_of_this_cell(seq,rowI,colJ, cmpntTyp);                            // therefore, seq->intrvlLukUpTable[rowI][colJ][1] = crik->cmpntListOcupidTyp[ocupidTypIndx]
                                                                                           // so that during the jump tree process, it will not pick a particular component of next component type
                                                                                           // which will eliminate the present component with covariance pair on it
	      return 0;
      }  // end inner if
    }    // end outer if
  }      // end for

  reach_end_of_loop(seq, crik, rowI, colJ, ocupidTypIndx);

  return 0;
}  // end fill_in_intrvl_luk_up_table_cell_i_j

//*****************************************************************************
// Function : Initialize Cell (I, J)
// Caller   : fill_in_intrvl_luk_up_table_cell_i_j()
// Purpose  : set memory space & initial value for the cell (I, J) on the interval look-up table
// Input    : 
// Return   : 
// Display  : 
//*****************************************************************************
int init_i_j_cell(config* seq, int16_t rowI, int16_t colJ)
{                                                                                     disp(seq,DISP_ALL,"entering init_i_j_cell'\n");
  seq->intrvlLukUpTable[rowI][colJ] = calloc(2, sizeof(int16_t*));                     // assign 2 spaces in cell (i, j),
  seq->intrvlLukUpTable[rowI][colJ][0] = -1;                                          // one marks the reading beginning on the cmpntListOcupidTyp
  seq->intrvlLukUpTable[rowI][colJ][1] = -1;                                          //            "          ending                "
                                                                                      // both initialized with -1 is to indicate that
                                                                                      // no useful information is in there yet
  
  disp(seq,DISP_ALL,"\n  (#%-2d, #%-2d) (%c, %c) ", rowI, colJ, seq->ltr[rowI], seq->ltr[colJ]);

  return 0;
}  // end init_i_j_cell

//*****************************************************************************
// Function : Does Lower Bound Fit in?
// Caller   : fill_in_intrvl_luk_up_table_cell_i_j()
// Purpose  : 
// Input    : 
// Return   : 
// Display  : 
//*****************************************************************************
int lowerbound_fits_in(config* seq, global* crik, int16_t rowI, int16_t colJ, int16_t ocupidTypIndx)
{                                                                                   disp(seq,DISP_ALL,"entering lowerbound_fits_in'\n");
  int16_t cursrUB = colJ - seq->minPairngDist - seq->minLenOfHlix * 2 + 1;            

  if(crik->cmpntListOcupidTyp[ocupidTypIndx] < rowI){
    disp(seq,DISP_ALL,"cmpntTyp %d not big enough\n", crik->cmpntListOcupidTyp[ocupidTypIndx]);
    return 1;                                                                       // the cmpnt type must be rowI < cmpntTyp < colJ to fall into the ballpark
  }    // end if 1                                                                  // if not, then move on to check next larger one to see if it's large enough

  if(crik->cmpntListOcupidTyp[ocupidTypIndx] > cursrUB){
    disp(seq,DISP_ALL,"cmpnt type %d too big\n", crik->cmpntListOcupidTyp[ocupidTypIndx]);
    if(seq->intrvlLukUpTable[rowI][colJ][0] == -1){
      disp(seq,DISP_ALL,"Searched thru the whole list and didn't find any pair\n");
      disp(seq,DISP_LV3,"      ");
      return 0;
    } else if (seq->intrvlLukUpTable[rowI][colJ][0] != -1){
      return 3;                                                                     // cursrUB reached, therefore there's no point to keep searching.
                                                                                    // so, fittableTypUB is to be used as ocupidTyp entry for cmpnt UB
    } else {
      disp(seq,DISP_LV1,"Exception happen in 'lowerbound_fits_in'");
    }  // end inner if
    return 0;                                                                       // the cmpnt type must be cmpntTyp < colJ
  }    // end if 2
  disp(seq,DISP_ALL,"LB fit in, exiting 'lowerbound_fits_in'\n");
  return 2;                                                                         // return 1 means LB does fit into the interval
                                                                                    // now we'r gonna chk UB
} // end lowerbound_fits_in

//*****************************************************************************
// Function : Does Upper Bound Fit in?
// Caller   : fill_in_intrvl_luk_up_table_cell_i_j()
// Purpose  : 
// Input    : 
// Return   : 1 (yes), 0 (no)
// Display  : 
//*****************************************************************************
int upperbound_fits_in(config* seq, global* crik, int16_t rowI, int16_t colJ, int16_t ocupidTypIndx)
{                                                                                             disp(seq,DISP_ALL,"entering upperbound_fits_in'\n");
  if(crik->cmpntList[crik->cmpntListOcupidTyp[ocupidTypIndx]].knob->closeBrsOutIndx <= colJ){   // outer edge of helix is within the upper bound.
    if(seq->intrvlLukUpTable[rowI][colJ][0] == -1){                                           // read beginning value hasn't been assigned
      seq->intrvlLukUpTable[rowI][colJ][0] = crik->cmpntListOcupidTyp[ocupidTypIndx];
      disp(seq,DISP_ALL,"\nread beginning of (%d, %d) is found to be %d\n", rowI, colJ, crik->cmpntListOcupidTyp[ocupidTypIndx]);
    }  // end inner if
    return 1;
  } else {
    disp(seq,DISP_ALL,"try next round, keep trying!\n");
    return 0;
  }    // end outer if
}  // end upperbound_fits_in

//*****************************************************************************
// Function : Is Covariance Pair?
// Caller   : fill_in_intrvl_luk_up_table_cell_i_j()
// Purpose  : check if a specific component type is on the covariance list
// Input    : 
// Return   : 1 (yes), 0 (no)
// Display  : 
//*****************************************************************************
int is_covarance_pair(config* seq, int16_t cmpntTyp)
{                                                                                              disp(seq,DISP_ALL,"enterng 'is_covarance_pair'\n");
  int8_t i;

  for(i = 0 ; i < seq->numCovari ; i++){
    if      (seq->coVari[i][0] == cmpntTyp) return 1;
    else if (seq->coVari[i][1] == cmpntTyp) fprintf(stderr, "Exception in is_covarance_pair\n");  // if there's match with cmptTyp, it should only happen at the head, not the tail
    else                                    return 0;
  }  // end for

  return 0;
}  // end is_covarance_pair

//*****************************************************************************
// Function : Mark Reading End of This Cell
// Caller   : fill_in_intrvl_luk_up_table_cell_i_j()
// Purpose  : end the investigation on this cell (i, j), and print out that cell %d-%d
// Input    : 
// Return   : 
// Display  : (i,j) cell look-up data
//*****************************************************************************
int mark_reading_end_of_this_cell(config* seq, int16_t rowI, int16_t colJ, int16_t cmpntTyp)
{                                                    disp(seq,DISP_ALL,"enterng 'mark_reading_end_of_this_cell'\n");
  seq->intrvlLukUpTable[rowI][colJ][1] = cmpntTyp;
  disp(seq,DISP_ALL,"read end of (%d,%d) is found to be, due to covari cut-off %d\n",  rowI, colJ, seq->intrvlLukUpTable[rowI][colJ][1]);
  disp(seq,DISP_ALL,"cmpnt type from %3d to %3d\n", seq->intrvlLukUpTable[rowI][colJ][0], seq->intrvlLukUpTable[rowI][colJ][1]);
  disp(seq,DISP_LV3," %2d-%-2d", seq->intrvlLukUpTable[rowI][colJ][0], seq->intrvlLukUpTable[rowI][colJ][1]);
  return 0;
}  // end mark_reading_end_of_this_cell

//*****************************************************************************
// Function : Reach End of Loop
// Caller   : fill_in_intrvl_luk_up_table_cell_i_j()
// Purpose  : make a proper ending of the loop
// Input    : seq
//          : crik
//          : rowI
//          : colJ
//          : ocupidTypIndx
// Return   : 
// Display  : 
//*****************************************************************************
int reach_end_of_loop(config* seq, global* crik, int16_t rowI, int16_t colJ, int16_t ocupidTypIndx)
{                                                                                     disp(seq,DISP_ALL,"entering reach_end_of_loop'\n");
  if ( (seq->intrvlLukUpTable[rowI][colJ][0] != -1) &&
       (seq->intrvlLukUpTable[rowI][colJ][1] == -1)){
    disp(seq,DISP_ALL,"read ending  of (%d, %d) is found to be %d\n",  rowI, colJ, crik->cmpntListOcupidTyp[ocupidTypIndx-1]);
    seq->intrvlLukUpTable[rowI][colJ][1] = crik->cmpntListOcupidTyp[ocupidTypIndx - 1];
    disp(seq,DISP_ALL,"cmpnt type from %3d to %3d\n", seq->intrvlLukUpTable[rowI][colJ][0], seq->intrvlLukUpTable[rowI][colJ][1]);
    disp(seq,DISP_LV3," %2d-%-2d", seq->intrvlLukUpTable[rowI][colJ][0], seq->intrvlLukUpTable[rowI][colJ][1]);
  } else if ( (seq->intrvlLukUpTable[rowI][colJ][0] != -1) &&
              (seq->intrvlLukUpTable[rowI][colJ][1] != -1)){
    disp(seq,DISP_ALL,"cmpnt type from %3d to %3d\n", seq->intrvlLukUpTable[rowI][colJ][0], seq->intrvlLukUpTable[rowI][colJ][1]  );
    disp(seq,DISP_LV3," %2d-%-2d", seq->intrvlLukUpTable[rowI][colJ][0], seq->intrvlLukUpTable[rowI][colJ][1]  );
  } else {
    disp(seq,DISP_ALL,"Searched thru the whole list and didn't find any pair\n");
    disp(seq,DISP_LV3,"      ");
  }  // end if

  return 0;
}  // end reach_end_of_loop
