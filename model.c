#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "model.h"

/* Creates a set of models from a dataset. Each sequence
    of the dataset creates a new statistical model */
stat_models_t *create_gene_models(dataset_t * dataset)
{
	int i;
	stat_models_t *models = (stat_models_t *)safe_malloc(sizeof(stat_models_t));

	models->modeln = dataset->seqn;
	models->models =
	    (stat_model_t **) safe_malloc(sizeof(stat_model_t *) * models->modeln);
	models->ids = dataset->ids;

	for (i = 0; i < models->modeln; i++) {
		models->models[i] =
		    (stat_model_t *)safe_calloc(1, sizeof(stat_model_t));
		update_model(models->models[i], dataset->sequences[i]);
		normalize_model(models->models[i]);
	}

	return models;
}

/* Destroys the set of models given by "models" */
void destroy_gene_models(stat_models_t * models)
{
	int i;

	for (i = 0; i < models->modeln; i++)
		safe_free(models->models[i]);

	safe_free(models->models);
	safe_free(models);
}

/* Creates a model of the background, from the set of sequences given
    by "background_ds". Unlike a gene model, each of the sequences in
    the dataset will contribute to the same model               */
stat_model_t *create_background_model(dataset_t * background_ds)
{

	stat_model_t *model = (stat_model_t *) safe_calloc(1, sizeof(stat_model_t));
	int i;
	for (i = 0; i < background_ds->seqn; i++) {
		update_model(model, background_ds->sequences[i]);
	}
	normalize_model(model);

	return model;
}

/* Given a model and a sequence, it updates the model with the
    new sequence. This ways, there is no transition between the
    last character of a sequence and the first character of the
    next sequence                                               */
void update_model(stat_model_t * model, char *seq)
{
	int i;
	int line_length = strlen(seq);
	int previous = int_code(seq[0]);
	model->p[previous]++;
	model->p[_N_]++;

	for (i = 1; i < line_length; i++) {
		int pos = int_code(seq[i]);
		if (pos) {
			model->p[pos]++;
			model->p[_N_]++;
		}
		if (pos && previous) {
			model->A[previous][pos]++;
			model->A[previous][_N_]++;
		}
		previous = pos;
	}
}

/* Aplies pseudocounts to the model and changes the
    representation from frequencies to probabilities.   */
void normalize_model(stat_model_t * model)
{
	int i, j;

	pseudocounts(model);

	for (i = 1; i < NUCLEOTIDES; i++) {
		model->p[i] = model->p[i] / model->p[_N_];

		for (j = 1; j < NUCLEOTIDES; j++) {
			model->A[i][j] = model->A[i][j] / model->A[i][_N_];
		}
	}
}

/* Applies pseudocounts to the Markov chains of a statistical model */
void pseudocounts(stat_model_t * model)
{
	int i, j;

	for (i = 1; i < NUCLEOTIDES; i++) {
		model->p[i] += PSEUDOCOUNTS;
		model->p[0] += PSEUDOCOUNTS;

		for (j = 1; j < NUCLEOTIDES; j++) {
			model->A[i][0] += PSEUDOCOUNTS;
			model->A[i][j] += PSEUDOCOUNTS;
		}
	}
}

/* Prints relevant information about a statistical model */
void print_model(stat_model_t * model)
{
	int j, k;
	long double sum = 0;

	printf("\n---------MODEL---------\n");
	for (j = 1; j < 5; j++) {
		printf("p(%c) = %Lg\n", char_code(j), model->p[j]);
		sum += model->p[j];
	}
	printf("SUM: %Lg\n", sum);
	sum = 0;

	for (j = 1; j < 5; j++) {
		for (k = 1; k < 5; k++) {
			printf("p(%c)(%c) = %Lg\n", char_code(j), char_code(k), model->A[j][k]);
			sum += model->A[j][k];
		}
	}
	printf("SUM: %Lg\n", sum);
	printf("\n");
}

/* Initializes an entity that contains information about a mirna in this 
	context. Each mirna, will have a score for each gene and background, with a certain number
	of allowed errors */
mirna_info_t *create_mirna_info(ushort gene_seqn)
{
	int e;

	mirna_info_t *mirna = (mirna_info_t *)safe_malloc(sizeof(mirna_info_t));
	mirna->gene_prob =
	    (long double **)safe_malloc(sizeof(long double *) * NUM_ERRORS);
	mirna->background_prob =
	    (long double *)safe_malloc(sizeof(long double) * NUM_ERRORS);

	for (e = 0; e < NUM_ERRORS; e++)
		mirna->gene_prob[e] =
		    (long double *)safe_malloc(sizeof(long double) * gene_seqn);

	return mirna;
}

mirnas_info_t *create_mirnas_info(ushort seqn, ushort gene_seqn)
{
	int i;
	mirnas_info_t *mirnas = (mirnas_info_t *)safe_malloc(sizeof(mirnas_info_t));
	mirnas->mirna = (mirna_info_t **)safe_malloc(sizeof(mirna_info_t *) * seqn);
	mirnas->seqn = seqn;

	for (i = 0; i < seqn; i++) {
		mirnas->mirna[i] = create_mirna_info(gene_seqn);
	}

	return mirnas;
}

void destroy_mirna_info(mirna_info_t * mirna)
{
	int e;
	for (e = 0; e < NUM_ERRORS; e++)
		safe_free(mirna->gene_prob[e]);

	safe_free(mirna->gene_prob);
	safe_free(mirna->background_prob);
	safe_free(mirna);
}

void destroy_mirnas_info(mirnas_info_t * mirnas)
{
	int i;

	for (i = 0; i < mirnas->seqn; i++) {
		destroy_mirna_info(mirnas->mirna[i]);
	}

	safe_free(mirnas->mirna);
	safe_free(mirnas);
}
