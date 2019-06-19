To use Freeimage first unzip it and then build the library
	unzip FreeImage.zip
	cd FreeImage
	make clean && make

To build:
	make clean && make

To make graphical image from the text outputs use any of the script
For example:
	./process_taum.sh

In titan you have to interactively login and use:
	qsub ./process_taum.sh
