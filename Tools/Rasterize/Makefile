CC := $(if $(shell which CC),CC,mpic++)

EXECS= rasterize


all: $(EXECS)

rasterize:
	$(CC) -o $@ rasterize.cpp -g -I../../include -I../../domain-nd/include -std=c++11 -IFreeImage/Dist ./FreeImage/Dist/libfreeimage.a


.PHONY: clean
clean:
	rm -f *.so *.o *~ rasterize