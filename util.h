#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>
#include <stdio.h>
#include "types.h"

#define _N_ 0
#define _C_ 1
#define _G_ 2
#define _A_ 3
#define _T_ 4
#define NUCLEOTIDES 5

void *safe_malloc(size_t);
void *safe_calloc(size_t, size_t);
void *safe_realloc(void *, size_t);
void safe_free(void *);
FILE *safe_fopen(const char *, const char *);
uint int_code(char c);
char char_code(uint c);
char *reverse_complement(char *);
char *get_string(FILE *, uint);
char *str_concat(char *, char *);
char *prepend_char(char *, char);
float fmax3(float, float, float);
int count_occurrences(char *string, char c);

#endif
