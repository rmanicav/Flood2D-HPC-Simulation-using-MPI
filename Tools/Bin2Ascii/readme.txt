To build:
	make clean && make

To convert all binary files in a directory use any of the script
For example:
	./convert_taum.sh

In titan you have to interactively login and use:
	aprun ./convert_taum.sh

bin2ascii takes 5 argument
	./bin2ascii $INPUT_DIR $OUTPUT_DIR $Row $Col $Padding
	
$INPUT_DIR = Input folder directory of Binary files
$OUTPUT_DIR = Output folder directory of Ascii files
$Row = Number of row
$Col = Number of col
$Padding = Padding used as a ghost cell

For example (Taum Sauk):
	./bin2ascii "../../output/flood2d/bin" "../..output/flood2d/txt" 1136 624 1