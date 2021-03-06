#include <stdlib.h>
#include <stdio.h>

#include "dataset.h"

/* Creates a dataset with multiple sequences, each one with an associated id */
dataset_t *create_dataset(ushort seqn, char **sequences, char **ids)
{
	dataset_t *ds = (dataset_t *) safe_malloc(sizeof(dataset_t));
	ds->seqn = seqn;
	ds->ids = ids;
	ds->sequences = sequences;

	return ds;
}

/* Frees the allocated memory of the data structure that corresponds to a dataset */
void destroy_dataset(dataset_t * ds)
{
	int i = 0;

	if (ds == NULL)
		return;

	for (i = 0; i < ds->seqn; i++) {
		safe_free(ds->sequences[i]);
		safe_free(ds->ids[i]);
	}

	safe_free(ds->sequences);
	safe_free(ds->ids);
	safe_free(ds);
}
