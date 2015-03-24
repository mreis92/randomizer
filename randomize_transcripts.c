#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "types.h"
#include "model.h"
#include "fasta.h"
#include "util.h"

#define NUM_TRANSCRIPTS 1
#define FASTA_LEN 80

int round_int( double r ) {
    return (r > 0.0) ? (r + 0.5) : (r - 0.5); 
}

int calculate_mean_length(dataset_t* d){
	int i;	
	int sum = 0;

	for (i = 0; i < d->seqn; i++) {
		sum += strlen(d->sequences[i]);		
	}

	return round_int(sum/d->seqn);
}

/* String randomization using Fisher-Yates Shuffle */ 
char *randomize_string_1(char* str){
	int i;
	int len = strlen(str);
	char *res = (char*)safe_malloc(sizeof(char)*(len+1));
	strcpy(res, str);

	srand(time(NULL));

	for(i = 0; i < len; i++){
		int r = (rand() % (len-i)) + i;

		int temp = res[r];
		res[r] = res[i];
		res[i] = temp;
	}	

	return res;
}

/* String generation from a Markov model, preserving statistical properties */
char *randomize_string(int len, stat_model_t* model){
	int i, j, prev, r;
	char *res = (char*)safe_malloc(sizeof(char)*(len+1));
	long double prob, weight = 0;

	r = rand() % 100; /*percentages between 0 and 99 */
	prev = 0;

	for(i = 1; i < 5; i++){
		prob = model->p[i] * 100;

		if(weight + prob > r){
			res[0] = char_code(i);
			prev = i;
			break;
		}
		else 
			weight += prob;
	}
	


	for(i = 1; i < len; i++) {
		weight = 0;
		r = rand() % 100;

		for(j = 1; j < 5; j++){
			prob = model->A[prev][j] * 100;

			if(weight + prob > r){
				res[i] = char_code(j);
				prev = j;
				break;
			}
			else 
				weight += prob;
		}
	}

	res[i] = '\0';
	return res;	
}

int main(int argc, char **argv)
{
	int i, j, k;
	FILE *transcripts = safe_fopen(argv[1], "r");
	char *folder = (argv[2] != NULL ? argv[2] : "");

	dataset_t *tds = parse_fasta(transcripts);
	int num_chars = calculate_mean_length(tds);
	stat_model_t *model = create_background_model(tds);
	srand(time(NULL));

	for(i = 0; i < NUM_TRANSCRIPTS; i++){
		FILE *f;
		char buffer[50];
		
		sprintf(buffer, "%stranscript_file_%d.fa", folder, i);
		f = safe_fopen(buffer, "w");

		for (j = 0; j < 100000; j++) {
			char *random = randomize_string(num_chars, model);
			int num_lines = num_chars/FASTA_LEN + (num_chars%FASTA_LEN ? 1 : 0);

			fprintf(f, ">random_%d\n", j);

			for(k = 0; k < num_lines; k++){
				fprintf(f, "%.*s\n", FASTA_LEN, random + (k*FASTA_LEN));
			}

			safe_free(random);

		}

		fclose(f);
	}

	safe_free(model);	
	destroy_dataset(tds);
	return 0;
}
