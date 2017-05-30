#ifndef __SWELLIX_STATS_H
#define __SWELLIX_STATS_H


float get_mfe_structure(char* sequence, char* structure);
void init_vrna(char* sequence);
void update_stats(config* seq);
void update_max_distance(config* seq);
void update_min_energy(config* seq);
void update_motif_count(config* seq);
void cleanup_stats();
void print_stats(config* seq);

#endif // __SWELLIX_STATS_H
