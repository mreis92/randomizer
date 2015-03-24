CC = gcc
ADDITIONAL = -g -O3
OBJECTS = util.o dataset.o model.o fasta.o randomize_transcripts.o
CFLAGS = -Wall $(ADDITIONAL)
LIBS = -lm 
ALIBS = $(LIBS)

EXE = rand

rand: $(OBJECTS) 
	$(CC) $(CFLAGS) -o $(EXE) $(OBJECTS) -L. $(ALIBS)

clean:
	rm -f $(OBJECTS) $(EXE) *~

indent:
	indent -linux *.c *.h
	rm *~

test: $(EXE)
	./$(EXE)
