#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include "subopt.h"
#include "bundle_list.h"

#if 0
{
#ifdef _MPI

#include <mpi.h>

#endif // _MPI
}
#endif

options OPTIONS;
char saveto[128];
char cwd[256];
char configID[64];
int START;
char* SEQ;
char* MODS;
int WINDOW;
int LENGTH;
int TERMINAL_MISMATCHES;
int ASYMMETRY;
int MISMATCHES;
LabeledStructures* mylab; //pointer to the labeledstructure object to use
LabeledStructures* mylabs; //target array of objects to build up
int* mylabsSize;
int* mylabsMax;
config* seq;

int slide_those_windows(char* subSeq, 
                        char* subMod, 
                        int startindex, 
                        char* mods, 
                        int window, 
                        int tmms, 
                        int asymm, 
                        config* pseq, 
                        LabeledStructures* labs,
                        int* labsSize,
                        int* labsMax) {

//printf("new slide windows call\n");

    // These parameters are a direct result of copying the label.c code into this file, so they could be more efficiently integrated
    // and possibly omitted upon further development. As for now, a working version of this change to the algorithm is of utmost importance.
    //
    START = startindex; //indicates the index within the sequence marking the beginning of the window
    SEQ = calloc(strlen(pseq->ltr)+1, sizeof(char));
    strcpy(SEQ, pseq->ltr);//SEQ = sequence; //the sequence of RNA being analyzed
    MODS = calloc(strlen(mods)+1, sizeof(char));
    strcpy(MODS, mods); //the chemical modification data pertaining to the sequence. (i.e. must-pair nt, must not pair nt, etc.)
    WINDOW = window; //the span of the current window
    LENGTH = pseq->minLenOfHlix;
    TERMINAL_MISMATCHES = tmms;
    ASYMMETRY = asymm;
    MISMATCHES = pseq->maxNumMismatch;

    mylabs = labs;
    mylabsSize = labsSize;
    mylabsMax = labsMax;

    mylab = &mylabs[*mylabsSize];//malloc(sizeof(LabeledStructures));
//    initLabeledStructures(mylab);
    mylab->titlesize = 64;
    mylab->buffsize = 4096;
    mylab->title = (char*)calloc(mylab->titlesize, sizeof(char));
    mylab->structures = (char*)calloc(mylab->buffsize, sizeof(char));

    seq = pseq;

//printf("START %i, WINDOW %i, LENGTH %i, TMMS %i, ASYMM %i, MMS %i\n", START, WINDOW, LENGTH, TERMINAL_MISMATCHES, ASYMMETRY, MISMATCHES);
    set_args();
    char* locseq = calloc(strlen(subSeq)+1, sizeof(char));
//    fscanf(OPTIONS.infile, "%s", seq);
    strcpy(locseq, subSeq);
    int i;
    for(i=0; locseq[i]; i++)
        locseq[i] = toupper(locseq[i]);
 
    int *constraints;
    if (OPTIONS.constraints){
        char* consts = calloc(strlen(subSeq)+1, sizeof(char));
        strcpy(consts, subMod);
        constraints = interpreted_constraints(consts);
        free(consts);
//printf("I got subseq: %s\nsubmod: %s\n", seq, consts);
    } else {
        constraints = NULL;
    }
    //print( "%s\n", seq);
    start(locseq, constraints);
//#endif // _MPI

    //write out final labeled structures computed for this window
    if(mylab->title[0] != '\0') {
//      mylabs[(*mylabsSize)++] = mylab;
      if(*mylabsSize == *mylabsMax) {
        *mylabsMax += seq->strLen*100;
        mylabs = (LabeledStructures*)realloc(mylabs, *mylabsMax);
      }
      (*mylabsSize)++;
    } else {
//      freeLabeledStructures(&mylab);
      free(mylab->title);
      free(mylab->structures);
      mylab->title = NULL;
      mylab->structures = NULL;
      mylab->titlesize = 0;
      mylab->buffsize = 0;
    }

    if (OPTIONS.count){
        //print( "%d\n", OPTIONS.count);
    }

    if (OPTIONS.infile != stdin) fclose(OPTIONS.infile);
    if (OPTIONS.outfile != stdout) fclose(OPTIONS.outfile);

    free(MODS);
    free(SEQ);

    free(locseq);
    if (constraints) free(constraints);
    return 0;
}

void start(char *seq, int *constraints)
{
    int len = strlen(seq);
    state base;

    base.length = len;
    int struc[len];
    int i;

    for(i=0; i<len; i++)
        struc[i] = -1;

    base.structure = &struc[0];
    base.sequence = seq;
    base.intervals = NULL;
    base.constraints = constraints;
    make_interval(&base, 0, base.length, LOOP);
    refine_state(&base);
    unmake_interval(&base);
}


void refine_state(state *s){

    if(!allow_state(s))
        return;

    // This is an optimization step, changing the order
    // in which we walk the tree to improve pruning efficiency.
    // rearrange_intervals(s);

    OPTIONS.statecount ++;
    refine_state_locally(s);
}


//------------------------------ BACKTRACKING --------------------------------//

suboptinterval *make_interval(state *s, int ai, int aj, interval_flag aflag)
{
    suboptinterval *new = malloc(sizeof(suboptinterval));
    new->i = ai;
    new->j = aj;
    new->flag = aflag;
    new->next = s->intervals;
    s->intervals = new;
    return new;
}

void unmake_interval(state *s)
{
    suboptinterval *i = s->intervals->next;
    free(s->intervals);
    s->intervals = i;
}

void make_pair(state *s, int i, int j)
{
    s->structure[i] = j;
    s->structure[j] = i;
}

void unmake_pair(state *s, int i, int j)
{
    s->structure[i] = -1;
    s->structure[j] = -1;
}

void refine_state_locally(state *s)
{
    if (s->intervals == NULL){
        print_soln(s);
        return;
    }

    int i = s->intervals->i;
    int j = s->intervals->j;
    interval_flag flag = s->intervals->flag;
    suboptinterval *temp = s->intervals;
    s->intervals = s->intervals->next;

    // There are only two options
    // Either j pairs with something in the interval:

    int k;

    for (k = j-MIN_PAIR_DIST-2; k >= i; k--){

        if (allow_pair(s, k, j-1)){
            make_pair(s, k, j-1);

            make_interval(s, k+1, j-1, LOOP);
            make_interval(s, i, i, LOOP);

            refine_state(s);

            unmake_interval(s);
            unmake_interval(s);

            unmake_pair(s, k, j-1);
        }
    }

    // Or j does not pair at all:
      if (j - 1 >= i && j < s->length){
          make_interval(s, i, j-1, flag);
        
          refine_state(s);
        
          unmake_interval(s);

      } else {
          refine_state(s);
      }


    s->intervals = temp;
}

//----------------------------- OPTIMIZATION ---------------------------------//

/* void rearrange_intervals(state *s){ */
/*     // this is a mildly unintuitive optimization */
/*     // for no-lonely-pair filtering. NOT CURRENTLY IN USE */
/*     if (!s->intervals) return; */
/*     interval * curr = s->intervals; */
/*     if (curr->flag == SECUNDO){ */
/*         curr->flag = LOOP; */
/*         interval *temp; */
/*         temp = curr->next; */
/*         curr->next = temp->next; */
/*         temp->next = curr; */
/*         s->intervals = temp; */
/*     } */
/* } */

//---------------------------- FILTERING PAIRS -------------------------------//

int allow_pair(state *s, int i, int j)
// returns true if the pair is allowable.
{
    if (can_pair(s->sequence[i], s->sequence[j]) &&
        !constrained_to_not_pair(s, i, j) &&
        (!OPTIONS.constraints || chemical_modification_ok(s, i, j)) &&
        (OPTIONS.max_dist == 0 || j - i <= OPTIONS.max_dist)) {
        return TRUE;
    }
    return FALSE;
}


int can_pair(char a, char b)
{
    if (OPTIONS.allpair) return TRUE;

    if ((a == 'A' && b == 'U') ||
        (a == 'U' && b == 'A') ||

        (a == 'G' && b == 'C') ||
        (a == 'C' && b == 'G')) {
        return TRUE;
    }
    if (!OPTIONS.noGU &&
        ((a == 'G' && b == 'U') ||
         (a == 'U' && b == 'G'))) {
        return TRUE;
    }

    return FALSE;
}

int constrained_to_not_pair(state *s, int i, int j)
{
    if ((s->constraints[i] == UNPAIRED_BASE) || (s->constraints[j] == UNPAIRED_BASE)) return TRUE;
	
    return FALSE;
}

int death_triad(state *s, int i)
{
    // TODO Modify this to be a good little puppy about the following cases:
    // )(( ))( and GU pairs
    if(s->constraints[i] == MODIFIED_BASE &&
       s->structure[i] > i &&
       (i == s->length-1 || s->structure[i+1] > i) &&
       (i == 0 || s->structure[i-1] > i))
        return TRUE;
    return FALSE;
}

int chemical_modification_ok(state *s, int i, int j)
{
    make_pair(s, i, j);
    if (death_triad(s, i) ||
        (i != s->length-1 && death_triad(s, i+1)) ||
        (i != 0 && death_triad(s, i-1)) ||

        death_triad(s, j) ||
        (j != s->length-1 && death_triad(s, j+1)) ||
        (j != 0 && death_triad(s, j-1))) {

        unmake_pair(s, i, j);
        return FALSE;
    }

    unmake_pair(s, i, j);
    return TRUE;
}



//---------------------------- FILTERING STATES ------------------------------//

int allow_state(state *s)
// returns true if the state is allowable.
{
    if (!(OPTIONS.noLP && lonely_pair(s)) &&
        (!OPTIONS.helix_number || room_for_helices(s)) &&
        1){
        return TRUE;
    }
    return FALSE;
}


int lonely_pair(state *s)
{
    // This can totes be optimized, dudes.
    if (s->intervals == NULL) return FALSE;
    suboptinterval* curr = s->intervals;
    for(; curr; curr=curr->next){
        int i = s->intervals->i;
        int j = s->intervals->j;
        //print( "%d, %d, ", i, j);
        if (j-i < MIN_PAIR_DIST && i > 0){
            int a = i-1;
            int b = s->structure[i-1];
            if ((a == 0 ||
                 s->structure[a-1] < a ||
                 s->structure[a-1] == -1 ) &&
                (b == s->length-1 ||
                 s->structure[b+1] > b ||
                 s->structure[b+1] == -1) &&
                s->structure[a+1] == -1 &&
                s->structure[b-1] == -1) {
                return TRUE;
            }
        }
    }
    //print( "\n");
    return FALSE;
}




int room_for_helices(state *s) {
    int badnesses, sizei, sizej, number = helix_count(s);
    suboptinterval * inter;
    int i, j;
    //print( "\n\n");
    //print_soln(s);

    for(inter = s->intervals; inter; inter = inter->next){
        badnesses = 0;
        sizei = 0;
        sizej = 0;
        i = inter->i;
        j = inter->j+1;

        // First push out, 3' then 5'
        while(badnesses < OPTIONS.helix_badnesses &&
              sizej < OPTIONS.helix_size &&
              j < s->length){
            if(s->structure[j] == -1){
                badnesses++;
            } else if (s->structure[j] > j) {
                break; // no rooooOOOm
            } else {
                sizej++;
            }
            j++;
        }

        // and again on the other side, to taste.
        while(badnesses < OPTIONS.helix_badnesses &&
              sizei < sizej &&
              i < s->length){
            if(s->structure[i] == -1){
                badnesses++;
            } else if (s->structure[i] < i) {
                break;
            } else {
                sizei++;
            }
            i--;
        }

        sizei = sizei < sizej ? sizei : sizej;
        if (sizei >= OPTIONS.helix_size &&
            badnesses <= OPTIONS.helix_badnesses){
            number++;
        }

        int rest = inter->j - inter->i - MIN_PAIR_DIST;
        rest -= OPTIONS.helix_size - sizei;
        if (rest >= 0) {
            number++;
            number += rest % ((OPTIONS.helix_size*2) + MIN_PAIR_DIST);
        }
    }

    // print( "   %d possible helices of len 3", number);

    if (number >= OPTIONS.helix_number){
        //  print( " so I keep it!\n\n\n");
        return TRUE;
    }
    return FALSE;
}


int helix_count(state *s){
    int i = 0, nexti = 0, size, prej, number = 0, badnesses = 0;

    while(nexti < s->length){
        i = nexti;
        do {
            // Move forward to the next concentric possible-helix.
            for(; i<s->length && s->structure[i]<i; i++);
            if (i >= s->length) break;

            prej = s->structure[i];
            badnesses = 0;
            size = 1;
            nexti = i+1;

            // look at the current helix
            while(size < OPTIONS.helix_size && badnesses <= OPTIONS.helix_badnesses){
                i++;

                if (s->structure[i] > i){ // there's a pair
                    size++;
                    badnesses += prej - s->structure[i] - 1;
                    prej = s->structure[i];
                    nexti = prej;

                } else if (s->structure[i] == -1) { // no pair
                    badnesses++;
                }
            }

            // Increase count when a helix was found
            if (size >= OPTIONS.helix_size && badnesses <=  OPTIONS.helix_badnesses) {
                number++;
            }

        } while(i < prej && i<s->length);

        if (i >= s->length) break;
    }
    return number;
}



//-------------------------------- PRINTING ----------------------------------//

void print_state(state *s){
    // Prints a partial solution for debugging purposes.
    int i;
    print( "[");
    for(i=0; i<(s->length); i++){
        print( "%d, ", s->structure[i]);
    }
    print( "], ");
    print( "%d, ", s->length);

    char a[s->length+1];
    a[s->length] = (char)0;
    for(i=0; i<(s->length); i++){
        if(s->structure[i] == -1){
            a[i] = '.';
        } else if (s->structure[i] > i) {
            a[i] = '(';
        } else if (s->structure[i] < i) {
            a[i] = ')';
        }
    }

    print( "%s, ", a);
    suboptinterval *b = s->intervals;
    while(b){
        print( "(%d, %d, #%d), ", b->i, b->j, b->flag);
        b = b->next;
    }
    print( "\n");
}

void print_soln(state *s){

    if(OPTIONS.count){
        OPTIONS.count++;
        if (OPTIONS.count%10000 == 0) print( "%d \n", OPTIONS.count);
        return;
    }

    int i;
    char a[s->length+1];
    a[s->length] = (char)0;
    for(i=0; i<(s->length); i++){
        if(s->structure[i] == -1){
            a[i] = '.';
        } else if (s->structure[i] > i) {
            a[i] = '(';
        } else if (s->structure[i] < i) {
            a[i] = ')';
        }
    }
    //    float en = energy_of_struct(s->sequence, a);
//    print( "%s\n", a);
	label_struct(a, START);

}


//----------------------------- HANDLING INPUT -------------------------------//

void subopt_print_usage()
{
    print( "\n\n    usage: \n"
           "    subopt [-noLP] [-noGU] [-max-dist #] [-noML]\n"
           "           [-count] [-i] [-o] [-help | -h]\n\n"

           "    -noLP       : Restricts to only solutions with no lonely pairs. This\n"
           "                  greatly reduces the number of solutions.\n"
           "    -noGU       : Disallows G-U pairing.\n"
           "    -noML       : Restricts to only solutions without multiloops.\n"
           "    -count      : Does not print solutions --- only the number of total\n"
           "                  solutions found.\n\n"

           "    -max-dist # : Restricts to only solutions where all pairings\n"
           "                  are between bases # or fewer bases apart.\n"
           "    -e #        : Restricts to only solutions withing #kcal/mol of \n"
           "                  the MFE.\n\n"
           "    -i FILE     : Defines an input file, containing a sequence.\n\n"
           "    -o FILE     : Defines an output file, to overwrite with solutions.\n\n"

           "    -help | -h  : Prints this message and exits.\n"
           "\n\n");
    exit(0);
}


void set_args()
{

	OPTIONS.infile = stdin;
	OPTIONS.outfile = stdout;
	OPTIONS.noLP = 0;
	OPTIONS.noGU = 0;
	OPTIONS.max_dist = 0;
	OPTIONS.count = 0;
	OPTIONS.constraints = 1;
	OPTIONS.noML = 0;
	OPTIONS.allpair = 0;
	OPTIONS.helix_badnesses = 0;
	OPTIONS.helix_number = 0;
	OPTIONS.helix_size = 0;
	OPTIONS.statecount = 0;

}

int * interpreted_constraints(char *constraints)
{
    int i;
    int len = strlen(constraints);
    int *ret = calloc(sizeof(int), len);

    for (i=0; i<len; i++){
        if(constraints[i] == '.'){
            ret[i] = -1;
        } else if (constraints[i] == 'X') {
            ret[i] = UNPAIRED_BASE;
        } else if (constraints[i] == 'M') {
            ret[i] = MODIFIED_BASE;
        }
    }
    int start, end;
    while(1){
        end = 0;
        while (constraints[end] != ')' && end < len){
            end++;
        }
        if(end == len) break;
        start = end;
        while (constraints[start] != '('){
            start--;
        }
        ret[start] = end;
        ret[end] = start;
        constraints[start] = '.';
        constraints[end] = '.';
    }

    return ret;
}

//----------------------------- SERIALIZATION -------------------------------//

// l, struct, n, intervals
int * pack_state(state * S)
{
    int l, num_intervals;
    int * out = calloc((sizeof *out), (2+(4*S->length)));
    suboptinterval * i;

    out[0] = S->length;

    int j;
    for (j=0; j<S->length; j++){
        out[j+2] = S->structure[j];
    }

    l = S->length + 2;

    for ( num_intervals = 0, i = S->intervals; i; i = i->next, num_intervals ++ ) {
        out[l+(num_intervals*3)+0] = i->i;
        out[l+(num_intervals*3)+1] = i->j;
        out[l+(num_intervals*3)+2] = i->flag;
    }

    out[1] = num_intervals;

    return out;
}

state *interpret_message(int *ints)
{
    state * S = malloc(sizeof(state));
    S->intervals = NULL;


    S->length = ints[0];


    S->structure = malloc((sizeof(S->structure)*S->length));
    int i;
    for (i=0; i<S->length; i++){
        S->structure[i] = ints[i+2];
    }

    int num_ints = ints[1];

    for (i = S->length+2+(3*(num_ints-1)); i >= S->length+2; i -= 3) {
        make_interval(S, ints[i+0], ints[i+1], ints[i+2]);
    }

    return S;
}

//#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#//
//            Here marks the end of the code dedicated to the Subopt process for analyzing a window in the sequence.			   //
//            The rest of the code is dedicated to the labeling process for choosing which generated structures are			   //
//            adequate for reentry to Swellix.												   //
//#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#//


// Most of the functions written below are a quick, rough way to port the labeling code into C. It was originally written in Python
// and so there were many mechanics within Python used to facilitate the calculation that C doesn't support natively. I wrote these functions
// to imitate some of the functionality that was used in the Python script.

int min(int a, int b) {
	return a < b ? a : b;
}

bool arrIntCmp(int len, int arr1[len], int arr2[len]) {
	int i = 0;
	while(i < len) {
		if(arr1[i] != arr2[i]) return FALSE;
		i++;
	}
	return TRUE;
}

void arrIntCpy(int len, int (*dest)[len], int source[len]) {
	int i = 0;
	int* local = *dest;
	while(i < len) {
		local[i] = source[i];
		i++;
	}
}

void initArrInt(int len, int (*array)[len]) {
	int i = 0;
	int* local = *array;

	while(i < len) {
		local[i] = -1;
		i++;
	}
}

void intrev(int len, int (*array)[len], int size) {
	int* local = *array;
	int temp, i = 0, j = size-1;

	while(i < j) {
		temp = local[i];
		local[i] = local[j];
		local[j] = temp;
	
		i++;
		j--;
	}
}

char* slice(char* str, int start) {
	int len = strlen(str);
	if(start < 0) start += len;

	char* sliced = str + start;

	return sliced;
}

char* strrev(char* str) {
      char *p1, *p2;

      if (! str || ! *str)
            return str;

      int len = strlen(str);

      for (p1 = str, p2 = str + len - 1; p2 > p1; ++p1, --p2)
      {
            *p1 ^= *p2;
            *p2 ^= *p1;
            *p1 ^= *p2;
      }

      for (p1 = str; p1 < str + len; p1++)
      {
	    if(*p1 == '(') *p1 = ')';
	    else if(*p1 == ')') *p1 = '(';
      }
	
      return str;
}

int parse_conf(char* configfile) {
	FILE* config = fopen(configfile, "r");
	char* line = NULL;
	size_t len = 0;
	
	if(config != NULL) {
		char* tok;
		char* reg = " \"=\n";

		getline(&line, &len, config);
		tok = strtok(line, reg);
		tok = strtok(NULL, reg);
		strcpy(SEQ, tok);

		getline(&line, &len, config);
		tok = strtok(line, reg);
		tok = strtok(NULL, reg);
		strcpy(MODS, tok);

		getline(&line, &len, config);
		tok = strtok(line, reg);
		tok = strtok(NULL, reg);
		WINDOW = atoi(tok);

		getline(&line, &len, config);
		tok = strtok(line, reg);
		tok = strtok(NULL, reg);
		LENGTH = atoi(tok);

		getline(&line, &len, config);
		tok = strtok(line, reg);
		tok = strtok(NULL, reg);
		TERMINAL_MISMATCHES = atoi(tok);

		getline(&line, &len, config);
		tok = strtok(line, reg);
		tok = strtok(NULL, reg);
		ASYMMETRY = atoi(tok);

		getline(&line, &len, config);
		tok = strtok(line, reg);
		tok = strtok(NULL, reg);
		MISMATCHES = atoi(tok);
	}


	fclose(config);
	free(line);

	return 0;
}

int parse_structfile(FILE* infile, char*** linesinfile, int lineReadCap) {
	int i;

	if(infile == NULL) return -1;

	*linesinfile = (char**)calloc(lineReadCap, sizeof(char*));
	for(i = 0; i < lineReadCap; i++) (*linesinfile)[i] = (char*)calloc(WINDOW+2, sizeof(char));

	i = 0;
	while(i < lineReadCap && fgets((*linesinfile)[i], WINDOW+2, infile) != NULL) {
		(*linesinfile)[i][strlen((*linesinfile)[i])-1] = '\0';
		i++;
	}

	return i;
}

bool pairable(char c1, char c2) {
	if((c1 == 'G' && c2 == 'C') ||
	   (c1 == 'C' && c2 == 'G') ||

	   (c1 == 'A' && c2 == 'U') ||
	   (c1 == 'U' && c2 == 'A') ||

	   (c1 == 'G' && c2 == 'U') ||
	   (c1 == 'U' && c2 == 'G'))
		return TRUE;
	return FALSE;
}

void saveLabeledStructure(int end, int width, char* structure, int helices, char* helix_info, int wcpairs, int gupairs, int asymmetry, int inner_mismatches, int terminal_mismatches, int chemmod, int thermo) {
	char* titleString = calloc(256, sizeof(char));
	char saveString[256];
	sprintf(titleString, "%dx%d.lab", end, width);
//printf("new titleString: %s\n", titleString);
//printf("strcmp(%s, %s) = %d\n", mylab->title, titleString, strcmp(mylab->title, titleString));
        if(mylab->title && strcmp(mylab->title, titleString) != 0) {
//printf("making bundles for %s\n%s", mylab->title, mylab->structures);
          if(mylab->title[0] != '\0') {
//            mylabs[(*mylabsSize)++] = mylab;
            if(*mylabsSize == *mylabsMax) {
               *mylabsMax += seq->strLen*100;
              mylabs = (LabeledStructures*)realloc(mylabs, *mylabsMax);
            }
            (*mylabsSize)++;
          } else {
            //freeLabeledStructures(&mylab);
            free(mylab->title);
            free(mylab->structures);
          }
//printf("resetting mylab %s\n%s\n", mylab->title, mylab->structures);
          //mylab = malloc(sizeof(LabeledStructures));
          mylab = &mylabs[*mylabsSize];
          initLabeledStructures(mylab);
          resetLabeledStructures(mylab, titleString);
        }
        
	sprintf(saveString, "%d,%d,%s,%d,%s,%d,%d,%d,%d,%d,%d,%d\n", end, width, structure, helices, helix_info, wcpairs, gupairs, asymmetry, inner_mismatches, terminal_mismatches, chemmod, thermo);
//printf("saving to %s\n%s\n", mylab->title, saveString);
        if(strlen(mylab->structures)+256 > mylab->buffsize) {
          mylab->buffsize = mylab->buffsize << 1;
          mylab->structures = (char*)realloc(mylab->structures, mylab->buffsize);
        }
	strcat(mylab->structures, saveString);
        free(titleString);
}

int label_struct(char* structure, int start) {
  int prev_openOut[WINDOW], prev_openIn[WINDOW], prev_closeOut[WINDOW], prev_closeIn[WINDOW];
  int terminal = 0, width = 0;
  int tmmi, tmmi2, tmmi3, k, m, n;
  int structlen = strlen(structure);

  initArrInt(WINDOW, &prev_openOut);
  initArrInt(WINDOW, &prev_openIn);
  initArrInt(WINDOW, &prev_closeOut);
  initArrInt(WINDOW, &prev_closeIn);
	
  for(tmmi = 0; tmmi <= TERMINAL_MISMATCHES; tmmi++) {
    for(tmmi2 = 0; tmmi2 <= TERMINAL_MISMATCHES; tmmi2++) {
      for(tmmi3 = 0; tmmi3 <= TERMINAL_MISMATCHES; tmmi3++) {
        bool exit_loop = FALSE;
        for(k = 0; k <= MISMATCHES; k++) {
          for(m = 0; m <= MISMATCHES; m++) {
            for(n = 0; n <= MISMATCHES; n++) {
              if(tmmi > k || tmmi2 > m) continue;
					
              bool loop_continue = FALSE;
              int out_tmm = 0;
              int j = 0;
              while(j < structlen && structure[j] != ')') j++;
              if(j == structlen) return 0;

              int i = j;
              while(structure[i] != '(') i--;
					
              int loop = j - i - 1;
		
              int internal_mismatches = 0;
              int length = 0;
              int asymmetry = 0;
              int CM_score = 0;
              int mm_counts[3] = {k, m, n};
              int mms[WINDOW];
              mms[0] = 0;
              int mmsSize = 1;
	
              //label.py line 56
              if(((loop - 3)/2 - TERMINAL_MISMATCHES) >= 0 && TERMINAL_MISMATCHES <= mm_counts[0]) {
                length = tmmi;
                internal_mismatches = tmmi;
                mms[0] = tmmi;
              } else {
                if(((loop - 3)/2 - tmmi) >= 0) {
                  if(tmmi <= mm_counts[0]) {
                    length = tmmi;
                    internal_mismatches = tmmi;
                    mms[0] = tmmi;
                  } else {
                    length = mm_counts[0];
                    internal_mismatches = mm_counts[0];
                    mms[0] = mm_counts[0];
                  }
                } else if((loop - 3)/2 > 0) {
                  if(mm_counts[0] <= (loop - 3)/2) {
                    length = mm_counts[0];
                    internal_mismatches = mm_counts[0];
                    mms[0] = mm_counts[0];
                  } else {
                    length = (loop - 3)/2;
                    internal_mismatches = (loop - 3)/2;
                    mms[0] = (loop - 3)/2;
                  }
                }
              }

              int gu_pairs = 0;
              int wc_pairs = 0;

              int helices = 0;
              int openOut[WINDOW], openIn[WINDOW], closeIn[WINDOW], closeOut[WINDOW];
              int initIndx;
              for(initIndx = 0; initIndx < WINDOW; initIndx++) {
                openOut[initIndx] = -1;
                openIn[initIndx] = -1;
                closeIn[initIndx] = -1;
                closeOut[initIndx] = -1;
              }

              openOut[0] = 0;
              openIn[0] = i + internal_mismatches;
              closeIn[0] = j - internal_mismatches;
              closeOut[0] = 0;
              int openOutSize = 1, openInSize = 1, closeInSize = 1, closeOutSize = 1;

              bool reset = FALSE;
              bool fail = FALSE;
              while(j < structlen) {
                int modified_i = i;
                if(i < 0) modified_i += structlen;
                if(structure[modified_i] == '(' && structure[j] == ')') {
                  if(reset) {
                    length = 0;
                    asymmetry = 0;
                    internal_mismatches = 0;
                    reset = FALSE;
                    mms[mmsSize++] = 0;
                    int ins = 0;
                    if((openOut[helices] - i - 1) >= (j - closeOut[helices] - 1)) ins = j - closeOut[helices] - i - 1;
                    else ins = openOut[helices] - i - 1;
                    int next_tmmi = 0;
                    if(helices == 0) next_tmmi = tmmi2;
                    else next_tmmi = tmmi3;
                    //label.py line 111
                    if((ins - next_tmmi) >= 0) {
                      if(next_tmmi <= mm_counts[helices+1]) {
                        length = next_tmmi;
                        internal_mismatches = next_tmmi;
                        mms[helices+1] = next_tmmi;
                      } else {
                        length = mm_counts[helices+1];
                        internal_mismatches = mm_counts[helices+1];
                        mms[helices+1] = mm_counts[helices+1];
                      }
                    } else if(ins > 0) {
                      if(mm_counts[helices+1] <= ins) {
                        length = mm_counts[helices+1];
                        internal_mismatches = mm_counts[helices+1];
                        mms[helices+1] = mm_counts[helices+1];
                      } else {
                        length = ins;
                        internal_mismatches = ins;
                        mms[helices+1] = ins;
                      }
                    }
                    helices++;
                    if(helices > 2) {
                      fail = TRUE;
                      break;
                    }
                    openIn[openInSize++] = i + internal_mismatches;
                    closeIn[closeInSize++] = j - internal_mismatches;
                    openOut[openOutSize++] = i;
                    closeOut[closeOutSize++] = j;
                  }

                  openOut[helices] = i;
                  closeOut[helices] = j;

                  if(j == (structlen - 1) || i == 0) {
                    int mm = mm_counts[helices];
                    if((TERMINAL_MISMATCHES - internal_mismatches) >= 0 && (TERMINAL_MISMATCHES - internal_mismatches) <= mm) {
                      length += TERMINAL_MISMATCHES - internal_mismatches;
                      closeOut[helices] = j + TERMINAL_MISMATCHES - internal_mismatches;
                      openOut[helices] = i - TERMINAL_MISMATCHES + internal_mismatches;
                      out_tmm = TERMINAL_MISMATCHES - internal_mismatches;
                      mms[helices] += out_tmm;
                    } else {
                      length += mm - internal_mismatches;
                      closeOut[helices] = j + mm - internal_mismatches;
                      openOut[helices] = i - mm + internal_mismatches;
                      out_tmm = mm - internal_mismatches;
                      mms[helices] += out_tmm;
                    }
                  }

                  i--;
                  j++;
                  length++;
                } else if(structure[modified_i] == '(' && structure[j] == '.') {
                  asymmetry++;
                  j++;
                  if(asymmetry > ASYMMETRY) {
                    reset = TRUE;
                    if(length < LENGTH) {
                      fail = TRUE;
                      break;
                    }
                  }
                // label.py line 172
                } else if(structure[modified_i] == '.' && structure[j] == ')') {
                  asymmetry++;
                  i--;
                  if(asymmetry > ASYMMETRY) {
                    reset = TRUE;
                    if(length < LENGTH) {
                      fail = TRUE;
                      break;
                    }
                  }
                } else if(structure[modified_i] == '.' && structure[j] == '.') {
                  internal_mismatches++;
                  if(pairable(structure[modified_i], structure[j])) CM_score += 100;
                  if(internal_mismatches > mm_counts[helices]) {
                    reset = TRUE;
                    if(length < LENGTH) {
                      fail = TRUE;
                      break;
                    }
                  } else {
                    openOut[helices] = i;
                    closeOut[helices] = j;
                    mms[helices]++;
                    length++;
                  }
                  i--;
                  j++;
                }
              }
              exit_loop = TRUE;
              //label.py line 202
              if(fail) {
                if(k == MISMATCHES && m == MISMATCHES && n == MISMATCHES) loop_continue = TRUE;
                exit_loop = FALSE;
              } else {
                helices++;

                if(loop_continue) continue;

                width = j - i - 1;

                if((length < LENGTH)||(internal_mismatches > MISMATCHES)||(asymmetry > ASYMMETRY)||(CM_score > 1000)) continue;

                if(helices > 1) {
                  int linkedmms = 0;
                  int var;
                  for(var = 0; var < helices; var++){
                    if((openOut[var] - openIn[var+1]) == 1 && (closeIn[var+1] - closeOut[var]) == 1)
                    linkedmms += mms[var] + mms[var+1];
                  }
                  if(linkedmms > MISMATCHES) continue;
                }

                int extra_pairs_needed = LENGTH - length;
                int terminal_mismatches = min((int)((loop-3)/2), extra_pairs_needed);
                length += min(loop - 3, extra_pairs_needed);
                // label.py line 251
                if(length < LENGTH) {
                  terminal = extra_pairs_needed;
                  terminal_mismatches += terminal;
                }

                if(mylab->title) {
                  char* helix_info = (char*)calloc(256, sizeof(char));
                  char helperString[256];
                  char outputStruct[WINDOW+terminal+1];
                  char terminalDots[terminal+1];

                  int z = 0;
                  while(z < terminal) terminalDots[z] = '.';
                  terminalDots[terminal] = '\0';
                  while(openOut[helices-1] > -out_tmm) {
                    int index;
                    for(index = 0; index < helices; index++) {
                      openOut[index]--;
                      openIn[index]--;
                      closeOut[index]--;
                      closeIn[index]--;
                    }
                  }

                  intrev(WINDOW, &openOut, openOutSize);
                  intrev(WINDOW, &openIn, openInSize);
                  intrev(WINDOW, &closeOut, closeOutSize);
                  intrev(WINDOW, &closeIn, closeInSize);

                  if(arrIntCmp(WINDOW, prev_openOut, openOut) && arrIntCmp(WINDOW, prev_openIn, openIn) && arrIntCmp(WINDOW, prev_closeOut, closeOut) && arrIntCmp(WINDOW, prev_closeIn, closeIn))
                    continue;

                  arrIntCpy(WINDOW, &prev_openOut, openOut);
                  arrIntCpy(WINDOW, &prev_openIn, openIn);
                  arrIntCpy(WINDOW, &prev_closeOut, closeOut);
                  arrIntCpy(WINDOW, &prev_closeIn, closeIn);

                  int index;
                  for(index = 0; index < helices; index++) {
                    sprintf(helperString, "%d/%d|%d/%d", openOut[index], openIn[index], closeIn[index], closeOut[index]);
                    strcat(helix_info, helperString);
                    if(index != helices - 1)
                      strcat(helix_info, ",");
                  }
                  strcpy(outputStruct, slice(structure, i - terminal + 1));
                  strcat(outputStruct, terminalDots);
                  saveLabeledStructure(start+WINDOW+terminal, width+(terminal<<2), outputStruct, helices, helix_info, wc_pairs, gu_pairs, asymmetry, internal_mismatches, terminal_mismatches, CM_score, 0);
                  free(helix_info);
                }
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

int label(char* structure, char* outDir, int i) {

	label_struct(structure, i);

	return 0;
}
