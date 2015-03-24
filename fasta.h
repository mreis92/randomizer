#ifndef FASTA_H
#define FASTA_H

#include <stdio.h>

#include "dataset.h"
#include "util.h"

#define BUFFERLEN   100

dataset_t *parse_fasta(FILE *);

#endif
