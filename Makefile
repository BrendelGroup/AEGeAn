CC=gcc
CFLAGS=-Wall -m64 -DWITHOUT_CAIRO
INCS=-I /usr/local/include/genometools/ -I inc/

obj:				obj/VangRelation.o obj/VangSchemaEntry.o obj/VangDegreeConstraint.o
				

obj/VangRelation.o:		inc/VangRelation.h src/VangRelation.c
				@- mkdir -p obj
				$(CC) $(CFLAGS) $(INCS) -c -o obj/VangRelation.o src/VangRelation.c

obj/VangSchemaEntry.o:		inc/VangSchemaEntry.h src/VangSchemaEntry.c
				@- mkdir -p obj
				$(CC) $(CFLAGS) $(INCS) -c -o obj/VangSchemaEntry.o src/VangSchemaEntry.c

obj/VangDegreeConstraint.o:	inc/VangDegreeConstraint.h src/VangDegreeConstraint.c
				@- mkdir -p obj
				$(CC) $(CFLAGS) $(INCS) -c -o obj/VangDegreeConstraint.o src/VangDegreeConstraint.c

clean:	
	rm -rf obj
