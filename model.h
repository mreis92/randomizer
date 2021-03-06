#ifndef MODEL_H
#define MODEL_H

#include "types.h"
#include "util.h"
#include "dataset.h"
#include "fasta.h"

/* If you want to allow N errors in your sequences, define NUM_ERRORS as N+1 */
#define NUM_ERRORS 1
#define PSEUDOCOUNTS 0.001

/* A statistical model that contain the Markov chains of order 0 and 1 for a given sequence */
typedef struct stat_model_t {
	long double p[NUCLEOTIDES];
	long double A[NUCLEOTIDES][NUCLEOTIDES];
} stat_model_t;

/* A set of statistical models, with an id associated with each model */
typedef struct stat_models_t {
	ushort modeln;
	stat_model_t **models;
	char **ids;
} stat_models_t;

/* The relevant information about a mirna, namely its sequence, id, and its background
	probability + probability of appearing on a given transcript */
typedef struct mirna_info_t {
	char *sequence;
	char *id;
	uint length;
	long double *background_prob;
	long double **gene_prob;

} mirna_info_t;

/* Information about a set of mirnas */
typedef struct mirnas_info_t {
	mirna_info_t **mirna;
	ushort seqn;
} mirnas_info_t;

stat_models_t *create_gene_models(dataset_t * dataset);
void destroy_gene_models(stat_models_t * models);
stat_model_t *create_background_model(dataset_t * background_ds);
void update_model(stat_model_t * model, char *seq);
void normalize_model(stat_model_t * model);
void pseudocounts(stat_model_t * model);
void print_model(stat_model_t * model);

mirna_info_t *create_mirna_info(ushort gene_seqn);
mirnas_info_t *create_mirnas_info(ushort seqn, ushort gene_seqn);
void destroy_mirna_info(mirna_info_t * mirna);
void destroy_mirnas_info(mirnas_info_t * mirnas);

#endif
