CC=		gcc

TN5dyncut:
	gcc -o dyncut dyncut.c thread_pool.c number.c -pthread -lz

clean:
	rm -fr *.o *~ dyncut *.dSYM
