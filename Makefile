CC=		gcc

dyncut:
	gcc -o $@ dyncut.c thread_pool.c number.c fastq.c -pthread -lz

clean:
	rm -fr *.o *~ dyncut *.dSYM
