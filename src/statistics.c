#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "main.h"
#include "statistics.h"
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/read_epars.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/eval.h>



//extern vrna_md_t md;
//extern vrna_param_t* P;
//extern vrna_fold_compound_t* vc;
//extern vrna_fold_compound_t* evalc;
static vrna_md_t md;
static vrna_param_t* P;
static vrna_fold_compound_t* vc;
static vrna_fold_compound_t* evalc;
float energy;
int distance;


/********************************************************************************************************************
 ********************************************************************************************************************
 ********************************************************************************************************************

	This file contains functions used to analyse structures/sequences. Initially consisting of implementations
	from ViennaRNA, the statistics available to use may grow in variety.

 ********************************************************************************************************************
 ********************************************************************************************************************
 *******************************************************************************************************************/


void init_vrna(char* sequence) {
  vrna_md_set_default(&md);
  if(paramfile != NULL) read_parameter_file(paramfile);
  P = vrna_params(&md);
  evalc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);
}

//******************************************************************************
//
//
//
//
//
//
//******************************************************************************
float get_mfe_structure(char* sequence, char* structure) {
//  if(md == NULL) { printf("error in get_mfe_structure; model uninitialized\n"); return -10000000; }
  if(vc == NULL) //vrna_fold_compound_free(vc);
    vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_MFE);
  return vrna_mfe(vc, structure);
}

void update_stats(config* seq) {
  switch (seq->statMode) {
    case STAT_DEFAULT:
      break;
    case STAT_MAX_DISTANCE:
      update_max_distance(seq);
      break;
    case STAT_MAX_FREE_ENERGY:
      update_max_energy(seq);
      break;
    case STAT_ALL:
      update_max_energy(seq);
      update_max_distance(seq);
      break;
    default:
      break;
  }
  if(seq->motif) update_motif_count(seq);
}

//******************************************************************************
//
//
//
//
//
//
//******************************************************************************
void update_max_energy(config* seq) {
  if(evalc == NULL) //vrna_fold_compound_free(evalc);
    evalc = vrna_fold_compound(seq->ltr, &md, VRNA_OPTION_EVAL_ONLY);
  energy = vrna_eval_structure(evalc, seq->dotNParen);
  if(energy > seq->maxenergy) {
    seq->maxenergy = energy;
  }
}

//******************************************************************************
//
//
//
//
//
//
//******************************************************************************
void update_max_distance(config* seq) {
  distance = vrna_bp_distance(seq->dotNParen, seq->mfe);
  if(distance > seq->maxdist) seq->maxdist = distance;
}

void update_motif_count(config* seq) {
//printf("updating motif count\n");
  char* s = seq->motifStruc;               // motif structure dot-parenthesis form
  char* l = seq->motifSeq;                 // motif sequence "letter" form
  char dnc = 'x';                          // "do not care" character
  char wcp = 'N';                          // Watson-Crick pair character
  int mlen = strlen(s);                    // length of motif string
  char* L = seq->ltr;                      // target sequence "letter" form
  char* S = seq->dotNParen;                // target sequence dot-parenthesis form
  int slen = seq->strLen;                  // target sequence length
  int match;                               // flag variable to signal a match
  int i = 0, j, k, m;                      // accessory indexing variables
//  printf("S: %s\n bundle: %i", S, seq->);

  while(i <= slen-mlen) {
    match = 1;
    k = 0;
    m = mlen-1;
    while((s[k] != dnc) && (k < mlen)) { // if this loop exits normally, then either the 5' side of a bulge
                                         // motif matched or the whole hairpin motif matched.
//printf("S: %s\ns: %s\nL: %s\nl: %s\nentered 5' match section\n", S, s, L, l);
      if((s[k] != S[i+k]) || ((l[k] != L[i+k]) && (l[k] != wcp))) { match = 0; break; }

      k++;
    }
    if(k < mlen && match) {
//printf("s: %s\nl: %s\nentered 3' match section\n", s, l);
    // if this condition is still TRUE, then the matched motif was the 5' side 
    // of a variable-length motif and we need to look for the 3' side.
      j = slen-1;
      int kmem = k;
      while(j >= i+mlen-1) { // check the whole 3' side of the sequence for a match.
        match = 1;
        m = mlen-1;
        k = kmem;
        while(s[m] != dnc) { // until hitting the "don't care" symbol, try to match the 3' side of the motif.
          if((s[m] != S[j-mlen+m+1]) || ((l[m] != L[j-mlen+m+1]) && (l[m] != wcp))) { match = 0; break; }
          m--;
        }
        if(match) {
          int opn = 0, cls = 0;
          while(i+k < j-mlen+m+1) { // check the "don't care" section for balanced pairing to ensure that
                                    // the motif matches are actually paired together accordingly.
            opn = ((S[i+k] == '(') || (S[j-mlen+m+1] == '(')) ? opn+1 : opn;
            cls = ((S[i+k] == ')') || (S[j-mlen+m+1] == ')')) ? cls+1 : cls;
            k++;
            m--;
          }
          if(opn == cls) { seq->motifCount++; return; } 
        }
        j--;
      }
      i++; continue; // an edge case could present where the above logic counts a valid motif during the final
                     // iteration of the j loop. this edge case would count a motif above and exit the loop with
                     // match == 1, then the following conditional statement would increment the motif count again.
    }
    if(match) { seq->motifCount++; return; }
    i++;
  }
}

void cleanup_stats() {
  if(vc != NULL) vrna_fold_compound_free(vc);
  if(evalc != NULL) vrna_fold_compound_free(evalc);
  if(P != NULL) free(P);
}

void print_stats(config* seq) {
  switch (seq->statMode) {
    case STAT_DEFAULT:
      break;
    case STAT_MAX_DISTANCE:
      fprintf(seq->dispFile, "Max RNAdistance: %d\n", seq->maxdist);
      break;
    case STAT_MAX_FREE_ENERGY:
      fprintf(seq->dispFile, "Max Free Energy: %.4f\n", seq->maxenergy);
      break;
    case STAT_ALL:
      fprintf(seq->dispFile, "Max RNAdistance: %d\n", seq->maxdist);
      fprintf(seq->dispFile, "Max Free Energy: %.4f\n", seq->maxenergy);
      break;
    default:
      break;
  }
  if(seq->motif) fprintf(seq->dispFile, "Motif match count: %d\n", seq->motifCount);
}
