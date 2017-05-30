#ifndef CRUMPLEH
#define CRUMPLEH

#define print(...) (fprintf(OPTION.outfile, __VA_ARGS__))
#define TRUE	  1
#define FALSE	  0
#define MIN_PAIR_DIST 3
// Flags and constants. 
enum {
    UNEXAMINED_BASE = -2,
    UNPAIRED_BASE = -1,
    UNCONSTRAINED_BASE,
    UNPAIRABLE_BASE,
    MODIFIED_BASE,
    MUST_PAIR,
    RESTRICTED_PAIR
};

 // options, filters, etc. See print_usage() and set_args().
struct options_struct {
    int noLP;
    int noGU;
    int max_dist;
    long count;
    int print;
    int constraints;
    int noML;
    int allpair;
    int breadcrumb;
    int helix_length;
    int helix_max_mismatch;
    int helix_max_bulge;
    int helix_min_count;
    int polling_interval;
    int limit_output;
    FILE *infile;
    FILE *outfile;
} OPTION;

typedef struct {
    int start;
    int end;
} interval;

enum {
    NEW_INTERVAL,
    PAIRING,
    NOT_PAIRING,
    COMPLETE
};

typedef struct {
    interval *intervals;
    int *structure;
    int current_interval;
    int current_base;
    int unpair;
    int max_base;
} state;

typedef struct{
    state *states;
    int length;
    int current;
    char *sequence;
    int *constraints;
} state_stack;

long START_TIME;

int main(int argc, char *argv[]);
void fold_all(state_stack *stack);
void fold1(state_stack *stack);

void new_interval(state *s, int start, int end);
void make_pair(state *s, int i, int j);
void unmake_pair(state *s, int i, int j);


state_stack * new_stack(char *seq);
void free_stack(state_stack *s);
void new_state(state_stack *s);
state *duplicate_state(state *s, state_stack *stack);
void new_interval(state *s, int start, int end);
void shift_stack_up(state_stack *s);

int allow_unpair(state_stack * s, int i);
int allow_pair(state_stack * s, int i, int j);
int can_pair(char a, char b);
int gu_pair(state_stack * s, int i, int j);
int death_triad(state_stack * s, int i);
int chemical_modification_ok(state_stack * s, int i, int j);
int constrained_to_not_pair(state_stack * s, int i, int j);

int allow_state(state *s, state_stack * stack);
int lonely_pair(state *s, state_stack * stack);

int allow_solution(state_stack * stack);
int helix_filter(state_stack * stack);

void print_a_state(state *s, int length);
void print_state(state_stack *stack, int index);
void print_stack(state_stack *stack);
void print_solution(state_stack *stack, int index);
void print_usage();
void get_args(int argc, char *argv[]);
char *get_sequence();
int *get_constraints(int len);
int *interpret_constraints(char *constraints);

#endif
