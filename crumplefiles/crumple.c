#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <time.h>
#include "crumple.h"

#if defined _MPI
#include "parallel_crumple.c"
#endif

//----------------------------- CORE PROCESS ------------------------------//

/* These three functions define the overall shape of the crumpling
 process. The algo-rhythm is defined recursively: 

    I is a list of ranges [i, j); each range is a segment of the
    sequence that has not yet been examined.  S is a paritially formed
    secondary structure --- represented here as a set of paired bases.
    (Note: Here, + is append)
    
    crumple(I, S):
        if |I| == 0:
            print S
        else:
            pop (i, j) from I
            if i == j:
                crumple(I, S)
            for each possible cannonical pairing (k, j) where i<k<j
                crumple(([i, k-1), [k+1, j-1)) + I, S U {(k, j)})
            crumple(([i, j-1)) + I, S)
     
   crumple([0, |S|), o)

 However, for the sake of speed and parallelizability, I've implemented
 this algorithm iteratively. With each I,S pair is grouped information
 about the bases under examination currently -- that structure is a
 state. Each new state is pushed onto a stack, and dealt with in
 last-in-first-out order. 

 Because the depth of the stack is limited by the length of the
 sequence, the whole stack can be allocated from the get-go, resulting
 in a very fast calculation. */


#ifndef _TEST
// Initialization of options; allocation of the state-stack.

int main(int argc, char *argv[])
{
    START_TIME = clock();
#if defined _MPI

    MPI_Init(&argc, &argv);
    get_mpi_globals();
    get_args(argc, argv);
    char *seq = get_sequence();
    state_stack *stack = new_stack(seq);

    fold_parallel(stack);

    print("NODE FINAL COUNT:%ld\n", OPTION.count);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

#else

    get_args(argc, argv);
    char *seq = get_sequence();
    state_stack *stack = new_stack(seq);
    new_state(stack);
    print("%s\n", seq);

    fold_all(stack);

//    print("FINAL COUNT:%ld\n", OPTION.count);

#endif // _MPI

    // Tidying up.
    free_stack(stack);
    free(seq);
    if(OPTION.infile != stdin)
        fclose(OPTION.infile);
    if(OPTION.outfile != stdout)
        fclose(OPTION.outfile);

    return 0;
}
#endif // _TEST


// This loops through folding,
// It halts when there are no more states in the stack.
void fold_all(state_stack *stack)
{
    while(stack->current >= 0){        
        fold1(stack);
    }
}


// Performs a single step of folding.
void fold1(state_stack *stack)
{
    state *s = &(stack->states[stack->current]);    

    // If the state is no good, dump it from the stack.
    if(!allow_state(s, stack)){
        stack->current--;
    }

    // If there are no more intervals, then this state is a solution. Print and toss.
    if (s->current_interval == -1){
        if(allow_solution(stack)){
            OPTION.count++;
            if (OPTION.print){
//                if (OPTION.limit_output && ftell(OPTION.outfile) > 65535){
//                    rewind(OPTION.outfile);
//                }
                print_solution(stack, stack->current);
            }
            if  (OPTION.breadcrumb != 0 && 
                 ((clock() - START_TIME) / CLOCKS_PER_SEC) > OPTION.breadcrumb){
                START_TIME = clock();
                print_solution(stack, stack->current);
            }
        }
        stack->current--;
        return;
    }

    // Take a look at the current interval.
    int start = s->intervals[s->current_interval].start;
    int end = s->intervals[s->current_interval].end;
    int pairer = s->current_base;

    // If the interval is totally empty, then toss the interval.
    if(start >= end){
        s->current_interval--;
        if(s->current_interval >= 0){
            s->current_base = s->intervals[s->current_interval].start + MIN_PAIR_DIST + 1;
            s->max_base = s->intervals[s->current_interval].end;
        }
        return;
    }

    // First, we'll try every base in the interval with which the
    // 'start' base can pair.
    for(;pairer < s->max_base; pairer++){
        if (allow_pair(stack, start, pairer)) {
            // When a suitable pair is found, mark our place and 
            // make a copy of the current state, with the new pair.
            s->current_base = pairer+1; 
            
            s = duplicate_state(s, stack);
            s->unpair = 1;

            make_pair(s, start, pairer);

            // Replace the current interval with two new intervals:
            // one inside the pair, and one in the remaining space 3'-ward of the pair.
            s->current_interval--;

            if(OPTION.noML && end < stack->length){
                // (If we don't want  multiloops. make the second interval 0 length)
                new_interval(s, end, end);
                int i; for(i=pairer+1; i<end; i++) s->structure[i]=UNPAIRED_BASE; 
            } else {
                new_interval(s, pairer+1, end);
            }

            new_interval(s, start+1, pairer);

            // Now return; the old state will sit on the stack until we
            // come back and continue from here.
            return;
        }
    }

    // If there aren't any more pairs, then look at solutions where 'start' doesn't pair: 
    if(s->unpair && allow_unpair(stack, start)){
        s->current_interval--;
        s->structure[start] = UNPAIRED_BASE;
        new_interval(s, start+1, end);
    } else {
        stack->current--;
    }
}



//----------------------------- STATE MANIPULATION --------------------------//

void new_interval(state *s, int start, int end)
{
    s->current_interval++;
    s->intervals[s->current_interval].start = start;
    s->intervals[s->current_interval].end = end;
    s->current_base = start + MIN_PAIR_DIST + 1;
    s->max_base = end;
    s->unpair = 1;
}



void make_pair(state *s, int i, int j)
{
    s->structure[i] = j;
    s->structure[j] = i;
}



void unmake_pair(state *s, int i, int j)
{
    s->structure[i] = UNEXAMINED_BASE;
    s->structure[j] = UNEXAMINED_BASE;
}



state *duplicate_state(state *s, state_stack *stack)
{
    // Copies a state right onto the end of the stack.

    stack->current++;
    state *new = &(stack->states[stack->current]);

    new->current_interval = s->current_interval;
    new->current_base = s->current_base;
    new->max_base = s->max_base;
    new->unpair = s->unpair;

    memcpy(new->structure, s->structure, stack->length*sizeof(int));

    memcpy(new->intervals, s->intervals, (s->current_interval+1)*sizeof(interval));

    return new;
}



//------------------------- WHOLE-STACK MANIPULATION -------------------------//

state_stack * new_stack(char *seq)
{
    int length = strlen(seq);
    state_stack *s = malloc(sizeof(state_stack));
    state *states = calloc(length+3, sizeof(state));
    int i;
    for(i=0; i<length; i++){
        states[i].intervals = calloc(length, sizeof(interval));
        states[i].structure = calloc(length, sizeof(int));
    }
    s->states = states;
    s->length = length;
    s->current = -1;
    s->sequence = seq;
    if (OPTION.constraints){
        s->constraints = get_constraints(s->length);
    }
    return s;
}



void free_stack(state_stack *s){
    int i;
    for(i=0; i<s->length+3; i++){
        free(s->states[i].intervals);
        free(s->states[i].structure);
    }
    free(s->states);
    if (OPTION.constraints){
        free(s->constraints);
    }
    free(s);
}



void new_state(state_stack *s)
{
    s->current++;
    s->states[s->current].current_interval = -1;
    new_interval(&(s->states[s->current]), 0, s->length);
    int i;
    for(i=0; i<s->length; i++){
        s->states[s->current].structure[i] = -2;
    }
    s->states[s->current].unpair = 1;
    s->states[s->current].max_base = s->length;
    s->states[s->current].current_base = MIN_PAIR_DIST + 1;
}


void shift_stack_up(state_stack *s)
{
    //    print("SHIFTIN!");
    int i;
    for(i=1; i<=s->current; i++){
        s->states[i-1].current_interval = s->states[i].current_interval;
        s->states[i-1].current_base = s->states[i].current_base;
        s->states[i-1].max_base = s->states[i].max_base;
        s->states[i-1].unpair = s->states[i].unpair;
               
        memcpy(s->states[i-1].structure, s->states[i].structure, 
               s->length*sizeof(int));
        
        memcpy(s->states[i-1].intervals, s->states[i].intervals, 
               (s->states[i].current_interval+1)*sizeof(interval));
    }
    s->current--;  
}
//---------------------------- FILTERING PAIRS -------------------------------//

int allow_unpair(state_stack * s, int i)
{
    // returns true if it is acceptable for the base i to remain unpaired.
    return (!OPTION.constraints
            || s->constraints[s->length * s->length + i] != MUST_PAIR);
}



int allow_pair(state_stack * s, int i, int j)
{
    // returns true if the pair i, j is acceptable and violates no constraints.
    return (can_pair(s->sequence[i], s->sequence[j]) &&
            (!OPTION.constraints || !constrained_to_not_pair(s, i, j)) &&
            (!OPTION.constraints || chemical_modification_ok(s, i, j)) &&
            (OPTION.max_dist == 0 || j - i <= OPTION.max_dist));
}



int can_pair(char a, char b)
{
    // Canonical & GU pairing constraints.
    if (OPTION.allpair)
        return TRUE;

    if ((a == 'A' && b == 'U') ||
        (a == 'U' && b == 'A') ||
        (a == 'G' && b == 'C') || (a == 'C' && b == 'G')) {
        return TRUE;
    }
    if (!OPTION.noGU && ((a == 'G' && b == 'U') || (a == 'U' && b == 'G'))) {
        return TRUE;
    }

    return FALSE;
}



int gu_pair(state_stack * s, int i, int j)
{
    return ((s->sequence[i] == 'G' && s->sequence[j] == 'U') ||
            (s->sequence[i] == 'U' && s->sequence[j] == 'G'));
}



int death_triad(state_stack * s, int i)
{
    int *structure = s->states[s->current].structure;
    // Returns true if base i would never be solvent-accessible.

    // If the triad contains a GU pair, it's safe.
    if ((structure[i] >= 0 &&
         gu_pair(s, i, structure[i])) ||
        (i + 1 < s->length &&
         structure[i + 1] >= 0 &&
         gu_pair(s, i + 1, structure[i + 1])) ||
        (i - 1 >= 0 &&
         structure[i - 1] >= 0 &&
         gu_pair(s, i - 1, structure[i - 1]))) 
        return FALSE;

    // If any of the triad is unpaired, it's safe.
    if (s->constraints[s->length * s->length + i] == MODIFIED_BASE &&
        structure[i] > i &&
        (i == s->length - 1 || structure[i + 1] > i) &&
        (i == 0 || structure[i - 1] > i))
        return TRUE;

    return FALSE;
}



int chemical_modification_ok(state_stack * s, int i, int j)
{
    // Returns true if pair i, j causes no CM conflicts.
    state *sta = &(s->states[s->current]);
    make_pair(sta, i, j);
    if (death_triad(s, i) ||
        (i != s->length - 1 && death_triad(s, i + 1)) ||
        (i != 0 && death_triad(s, i - 1)) ||
        death_triad(s, j) ||
        (j != s->length - 1 && death_triad(s, j + 1)) ||
        (j != 0 && death_triad(s, j - 1))) {

        unmake_pair(sta, i, j);
        return FALSE;
    }

    unmake_pair(sta, i, j);
    return TRUE;
}



int constrained_to_not_pair(state_stack * s, int i, int j)
{
    // Returns true if i or j must not pair, or if i,j is forbidden.
    return (s->constraints[s->length * i + j] == RESTRICTED_PAIR ||
            s->constraints[s->length * s->length + i] == UNPAIRABLE_BASE ||
            s->constraints[s->length * s->length + j] == UNPAIRABLE_BASE);
}



//---------------------------- FILTERING STATES ------------------------------//

int allow_state(state *s, state_stack * stack)
{
    // returns true if the state is acceptable and violates no constraints.
    if (!(OPTION.noLP && lonely_pair(s, stack))
        && 1) {
        return TRUE;
    }
    return FALSE;
}



int lonely_pair(state *s, state_stack *stack)
{
    // THIS DEPENDS ON INTERVAL ORDER NOW.
    // When the first interval contains no bases, it is finally possible
    // to check the pair immediately outside it: .([]......).
    // If that pair is lonely now, it always will be. 
    if (s->current_interval == -1)
        return FALSE;

    //print_state(stack, stack->current);
    interval in = s->intervals[s->current_interval];
    if (in.end - in.start < MIN_PAIR_DIST &&
        in.end < stack->length) {
        int a = s->structure[in.end];
        int b = in.end;
        if ((a == 0 || // a i the first base or
             s->structure[a - 1] < a || // the base before a points the wrong way )(
             s->structure[a - 1] < 0) &&  // or the base before a is unpaired.
            (b == stack->length - 1 ||
             s->structure[b + 1] > b ||
             s->structure[b + 1] < 0) &&
            (s->structure[a + 1] == -1 && s->structure[b - 1] == -1)) {
            return TRUE;
        }
    }
    return FALSE;
}



// ---------------- FILTERING SOLUTIONS --------------------//

int allow_solution(state_stack * stack)
{
    // Returns true if the state is acceptable and violates no constraints.
    if ((!OPTION.helix_length || helix_filter(stack))
        && 1) {
        return TRUE;
    }
    return FALSE;
}



int helix_filter(state_stack * stack)
{
    // Counts helices, based on the helix constraints in OPTION
    int helices = 0;
    int current;
    int bulge = 0, loop = 0, pairs = 0, cloop = 0, cbulge = 0;
    state *s = &stack->states[stack->current];

    for (current = 0; current < stack->length; current++) {
        // Zoom forward until a ( is found.
        if (s->structure[current] > current) {
            // A helix has started!
            int jcurrent = s->structure[current];
            pairs = 0;
            bulge = 0;
            loop = 0;
            cbulge = 0;
            cloop = 0;
            while (1) {
                if (current >= jcurrent) {
                    // We must be done!
                    if (current == jcurrent)
                        cloop++;        // if there's one last .
                    if (bulge <= OPTION.helix_max_bulge &&
                        pairs + loop/2 >= (OPTION.helix_length -
                                           (OPTION.helix_max_mismatch - loop/2))){
                        helices++;
                    }
                    break;
                }
                if (s->structure[current] > current) {
                    if (s->structure[jcurrent] == current) {
                        // Found a pair.
                        pairs++;
                        bulge += cbulge;
                        loop += cloop;
                        cloop = cbulge = 0;

                    } else if (s->structure[jcurrent] == UNPAIRED_BASE) {
                        // Found a bulge.
                        current--;
                        cbulge++;

                    } else if (s->structure[jcurrent] != current) {
                        // Found a multiloop.
                        if (bulge <= OPTION.helix_max_bulge &&
                            pairs + loop/2 >= (OPTION.helix_length -
                                               (OPTION.helix_max_mismatch -
                                               loop/2))) {
                            helices++;
                        }
                        current--;
                        break;
                    }
                } else if (s->structure[current] == UNPAIRED_BASE) {
                    if (s->structure[jcurrent] == UNPAIRED_BASE) {
                        // open loop.
                        cloop += 2;
                    } else if (s->structure[jcurrent] > current) {
                        // open bulge.
                        cbulge++;
                        jcurrent++;
                    }
                }
                current++;
                jcurrent--;
            }
        }
    }

    return helices >= OPTION.helix_min_count;
}



// ---------------------------- Printing ----------------------------//

void print_a_state(state *s, int length)
{
    // Prints a partial solution for debugging purposes.
    int i;
    char a[length + 1];
    a[length] = (char) 0;
    for (i = 0; i < (length); i++) {
        if (s->structure[i] == UNPAIRED_BASE) {
            a[i] = '.';
        } else if (s->structure[i] == UNEXAMINED_BASE) {
            a[i] = '|';
        } else if (s->structure[i] > i) {
            a[i] = '(';
        } else if (s->structure[i] < i) {
            a[i] = ')';
        }
    }

    print("%s, %d intervals, current base:%d, max_base:%d, unpair: %d. ", 
          a,
          s->current_interval+1, 
          s->current_base,
          s->max_base, 
          s->unpair);

    interval *b = s->intervals;
    for(i=0; i<=s->current_interval; i++) {
        print("(%d-%d], ", b[i].start, b[i].end);
    }
    print("\n");
}



void print_state(state_stack *stack, int index){
    print_a_state(&stack->states[index], stack->length);
}



void print_stack(state_stack *stack){
    int i;
    for(i=0; i<=stack->current; i++){
        print_state(stack, i);
    }
    print("\n\n");
}



void print_solution(state_stack *stack, int index)
{
    // Prints a finished solution.
    if(OPTION.breadcrumb)
        print("%ld ", OPTION.count);
    state s = stack->states[index];
    int i;
    char a[stack->length + 1];
    a[stack->length] = (char) 0;
    for (i = 0; i < (stack->length); i++) {
        if (s.structure[i] == UNPAIRED_BASE) {
            a[i] = '.';
        } else if (s.structure[i] == UNEXAMINED_BASE) {
            a[i] = '|';
        } else if (s.structure[i] > i) {
            a[i] = '(';
        } else if (s.structure[i] < i) {
            a[i] = ')';
        }
    }

    print("%s\n", a);
}



//-------------------------------- User Input --------------------------------//

void print_usage()
{
    print("\n\n    usage: \n"
          "    crumple [-noLP] [-noGU] [-max-dist #] [-noML]\n"
          "           [-count [-b #]] [-i] [-o] [-help | -h]\n\n"
          "    -b #        : If -count is used, b determines the 'breadcrumb' rate.\n"
          "                  If -b 1000, solutions will print once ever 1000 seconds.\n"
          "    -C          : Uses constraints from input: M for solvent-\n"
          "                  accessible bases; X or x for modified, unpairing\n"
          "                  bases; ( and ) for covarying, must-pair bases; . for\n"
          "                  unconstrained bases.\n"
          "    -help | -h  : Prints this message and exits.\n"
          "    -i FILE     : Defines an input file, containing a sequence.\n\n"
          "    -max-dist # : Restricts to only solutions where all pairings\n"
          "                  are between bases # or fewer bases apart.\n"
          "    -noLP       : Restricts to only solutions with no lonely pairs.\n"
          "    -noGU       : Disallows G-U pairing.\n"
          "    -noML       : Restricts to only solutions without multiloops.\n"
          "    -o FILE     : Defines an output file, to overwrite with solutions.\n\n"
          "    --count      : Does not print solutions --- only the number of total\n"
          "                  solutions found.\n"
          "    --helix-count #      : Print only solutions with # helices.\n"
          "    --helix-length #     : Print only solutions with helices # bases long.\n"
          "    --helix-bulge #      : Print only solutions with <= # bulges per helix.\n"
          "    --helix-mismatches # : Print only solutions with <= # mismatches per helix.\n\n"
          "    --polling-interval # : MPI only; how many states to process before \n"
          "                           sending work to needy nodes.\n\n"
          "\n\n");
    exit(1);
}



void get_args(int argc, char *argv[])
{
    OPTION.noLP = 0;
    OPTION.noGU = 0;
    OPTION.max_dist = 0;
    OPTION.print = 1;
    OPTION.count = 0;
    OPTION.constraints = 0;
    OPTION.noML = 0;
    OPTION.allpair = 0;
    OPTION.infile = stdin;
    OPTION.outfile = stdout;
    OPTION.breadcrumb = 0;
    OPTION.helix_length = 0;
    OPTION.helix_min_count = 0;
    OPTION.helix_max_bulge = 0;
    OPTION.helix_max_mismatch = 0;
    OPTION.polling_interval = 64;
    OPTION.limit_output = FALSE;

    int arg;
    for (arg = 1; arg < argc; arg++) {
        if (!strcmp("-noLP", argv[arg]))
            OPTION.noLP = TRUE;

        else if (!strcmp("-C", argv[arg]))
            OPTION.constraints = TRUE;

        else if (!strcmp("-noGU", argv[arg]))
            OPTION.noGU = TRUE;

        else if (!strcmp("--max-dist", argv[arg]))
            OPTION.max_dist = atoi(argv[++arg]);

        else if (!strcmp("--allpair", argv[arg]))
            OPTION.allpair = 1;

        else if (!strcmp("--helix-count", argv[arg]))
            OPTION.helix_min_count = atoi(argv[++arg]);

        else if (!strcmp("--helix-mismatches", argv[arg]))
            OPTION.helix_max_mismatch = atoi(argv[++arg]);

        else if (!strcmp("--helix-bulge", argv[arg]))
            OPTION.helix_max_bulge = atoi(argv[++arg]);

        else if (!strcmp("--helix-length", argv[arg]))
            OPTION.helix_length = atoi(argv[++arg]);

        else if (!strcmp("--count", argv[arg]))
            OPTION.print = 0;

        else if (!strcmp("-b", argv[arg]))
            OPTION.breadcrumb = atoi(argv[++arg]);

        else if (!strcmp("--polling-interval", argv[arg]))
            OPTION.polling_interval = atoi(argv[++arg]);

        else if (!strcmp("-i", argv[arg])) {
            errno = 0;
            OPTION.infile = fopen(argv[++arg], "r");
            if (errno) {
              fprintf(stderr,
                      "\n\nnames flicker\nbut\none is missing\n\nSorry: %s", argv[arg]);
              exit(1);
            }
        } else if (!strcmp("-o", argv[arg])) {
            arg++;
            char *ofile = calloc(1, strlen(argv[arg]) + 15);
            memcpy(ofile, argv[arg], strlen(argv[arg]));
#if defined(_MPI)
            sprintf(ofile + strlen(argv[arg]), "-%d", mpi_rank);
#endif
            OPTION.outfile = fopen(ofile, "w");
            free(ofile);
        }

        else if (!strcmp("-noML", argv[arg]))
            OPTION.noML = TRUE;

        else if (!strcmp("-help", argv[arg]) || !strcmp("-h", argv[arg]))
            print_usage();

        else if (!strcmp("-limit-output", argv[arg]))
            OPTION.limit_output = TRUE;

        else {
            print("\n    Unknown option %s", argv[arg]);
            print_usage();
        }
    }
}



char *get_sequence()
{
    char seq[40000];

    fscanf(OPTION.infile, "%s\n", seq);

    char *seq2 = malloc(sizeof(char) * (strlen(seq) + 1));
    seq2 = memcpy(seq2, seq, sizeof(char) * (strlen(seq) + 1));
    int i;
    for (i = 0; seq[i]; i++)
        seq2[i] = toupper(seq2[i]);
    return seq2;
}



int *get_constraints(int len)
{
    // Toss into get_constraints; add error checking.
    int *constraints;
    if (OPTION.constraints) {
        char consts[len];
        fscanf(OPTION.infile, "%s", consts);
        constraints = interpret_constraints(consts);
    } else {
        constraints = NULL;
    }
    return constraints;
}



int *interpret_constraints(char *constraints)
{
    // TODO: make this a 2D array, please!
    int i, j;
    int len = strlen(constraints);
    int *ret = calloc(sizeof(int), len * (len + 1));

    for (i = 0; i < len; i++) {
        if (constraints[i] == '.') {
            ret[len * len + i] = UNCONSTRAINED_BASE;
        } else if (constraints[i] == 'X' || constraints[i] == 'x') {
            ret[len * len + i] = UNPAIRABLE_BASE;
        } else if (constraints[i] == 'M') {
            ret[len * len + i] = MODIFIED_BASE;
        } else if (constraints[i] == 'V') {
            ret[len * len + i] = MUST_PAIR;
        }
    }

    int start, end;

    while (1) {
        end = 0;
        while (constraints[end] != ')' && end < len) {
            end++;
        }
        if (end == len)
            break;
        start = end;
        while (constraints[start] != '(') {
            start--;
        }
        ret[len * len + start] = MUST_PAIR;
        ret[len * len + end] = MUST_PAIR;
        for (j = start; j <= end; j++) {
            for (i = 0; i <= start; i++) {
                ret[i * len + j] = RESTRICTED_PAIR;
                ret[j * len + i] = RESTRICTED_PAIR;
            }
            for (i = end; i < len; i++) {
                ret[j * len + i] = RESTRICTED_PAIR;
                ret[i * len + j] = RESTRICTED_PAIR;
            }
        }
        ret[start * len + end] = 0;
        ret[end * len + start] = 0;

        constraints[start] = '.';
        constraints[end] = '.';
    }

    // Must never pair addition
    while (1) {
        end = 0;
        while (constraints[end] != ']' && end < len) {
            end++;
        }
        if (end == len)
            break;
        start = end;
        while (constraints[start] != '[') {
            start--;
        }
        ret[start * len + end] = RESTRICTED_PAIR;
        ret[end * len + start] = RESTRICTED_PAIR;

        constraints[start] = '.';
        constraints[end] = '.';
    }

    return ret;
}
