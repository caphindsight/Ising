ising: Ising.c
	gcc -lm -g -O0 -o ising Ising.c

release: Ising.c
	gcc -lm -O2 -o ising Ising.c

clean:
	rm -f ising

.PHONY: build clean
