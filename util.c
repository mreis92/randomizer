#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "util.h"

/* Safely allocates memory of space given by size */
void *safe_malloc(size_t size)
{
	void *mem = NULL;

	if (size == 0)
		printf("Allocating 0 bytes!");

	if ((mem = malloc(size)) != NULL)
		return mem;
	else {
		printf("Memory allocation error");
		return NULL;
	}
}

/*Safely allocates memory of space given by nmemb*size */
void *safe_calloc(size_t nmemb, size_t size)
{
	void *mem = NULL;

	if (size == 0 || nmemb == 0)
		printf("Allocating 0 bytes!");

	if ((mem = calloc(nmemb, size)) != NULL)
		return mem;
	else {
		printf("Memory allocation error");
		return NULL;
	}

}

/* Safely reallocates memory by the argument size, at the location pointed by
   ptr */
void *safe_realloc(void *ptr, size_t size)
{
	void *mem = NULL;

	if ((mem = realloc(ptr, size)) != NULL)
		return mem;
	else {
		printf("Memory re-allocation error");
		return NULL;
	}

}

/* Safely frees memory */
void safe_free(void *mem)
{
	if (mem != NULL)
		free(mem);
}

/* Safely opens a file */
FILE *safe_fopen(const char *path, const char *mode)
{
	FILE *file = NULL;
	file = fopen(path, mode);

	if (file == NULL) {
		printf("Could not open file: %s\n", path);
		return NULL;
	}

	return file;
}

/* Returns the code of a nucleotide char */
uint int_code(char c)
{
	switch (c) {
	case 'C':
	case 'c':
		return _C_;
	case 'G':
	case 'g':
		return _G_;
	case 'A':
	case 'a':
		return _A_;
	case 'T':
	case 't':
	case 'U':
	case 'u':
		return _T_;
	default:
		return 0;
	}
}

/* Returns the nucleotide char for the int code */
char char_code(uint c)
{
	switch (c) {
	case _C_:
		return 'C';
	case _G_:
		return 'G';
	case _A_:
		return 'A';
	case _T_:
		return 'T';
	case _N_: 
		return 'N';
	default:
		return '?';
	}
}

/* Returns a new sequence, but reversed */
char *reverse_complement(char *seq)
{
	int i;
	char *res = (char *)safe_malloc(sizeof(char) * (strlen(seq) + 1));
	int last_pos = strlen(seq) - 1;

	for (i = 0; i < strlen(seq); i++) {
		switch (seq[i]) {
		case 'C':
		case 'c':
			res[last_pos - i] = 'G';
			break;
		case 'G':
		case 'g':
			res[last_pos - i] = 'C';
			break;
		case 'A':
		case 'a':
			res[last_pos - i] = 'T';
			break;
		case 'T':
		case 't':
		case 'U':
		case 'u':
			res[last_pos - i] = 'A';
			break;
		default:
			break;
		}
	}
	res[i] = '\0';
	return res;
}

/* Gets a string from a file, that may be of variable size. The parameter
   size just gives a default size for the allocated memory */
char *get_string(FILE * fp, uint size)
{
	char c;
	int len = 0;
	char *str = safe_realloc(NULL, sizeof(char) * size);

	while (EOF != (c = fgetc(fp)) && c != '\n') {

		str[len++] = c;
		if (len == size) {
			str = safe_realloc(str, sizeof(char) * (size <<= 1));
		}
	}

	str[len++] = '\0';
	return safe_realloc(str, sizeof(char) * len);
}

/* Concatenates str2 to str1 */
char *str_concat(char *str1, char *str2)
{
	int len = strlen(str1) + strlen(str2) + 1;
	char *res = (char *)safe_malloc(sizeof(char) * len);
	snprintf(res, len, "%s%s", str1, str2);
	return res;
}

/* Appends a char to the beginning of str */
char *prepend_char(char *str, char c){
	int len = strlen(str);
	char *res = (char *)safe_malloc(sizeof(char) * (len+2));
	char *temp = res;
	*res++ = c;

	while((*res++ = *str++));

	return temp;
}

/* Determines the max between the three floats passed as arguments */
float fmax3(float f1, float f2, float f3){
	float max = f1;

	if(f2 > max)
		max = f2;
	if(f3 > max)
		max = f3;

	return max;
}

/* Counts the number of occurrences of c in string */
int count_occurrences(char *string, char c){
	int count;

	for (count = 0; string[count]; string[count] == c ? count++ : *string++);

	return count;
}

