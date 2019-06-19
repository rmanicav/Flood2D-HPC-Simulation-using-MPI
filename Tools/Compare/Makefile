CC := $(if $(shell which CC),CC,mpic++)

EXECS= compare


all: $(EXECS)

compare:
	$(CC) -o $@ compare.cpp -g -I../../include -I../../domain-nd/include --std=c++11 -lz


.PHONY: clean
clean:
	rm -f *.so *.o *~ compare