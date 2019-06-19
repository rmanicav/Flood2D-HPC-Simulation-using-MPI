To build:
	make clean && make

To convert all output files between two directory use any of the script
For example:
	./compare.sh $Row $Col $Number $INPUT_DIR_1 $INPUT_DIR_2 $OUTPUT_DIR

In titan you have to interactively login and use:
	qsub ./compare_taum.sh $Row $Col $Number $INPUT_DIR_1 $INPUT_DIR_2 $OUTPUT_DIR
	
$INPUT_DIR_1 = Input folder directory of txt output files
$INPUT_DIR_1 = Input folder directory of another txt output files
$OUTPUT_DIR = Output folder directory of compared value
$Row = Number of row
$Col = Number of col
$Number = No of files to compare

For example (Taum Sauk):
	./compare_taum.sh 1136 624 60 ../../output/flood2d-output-pkg-1/txt ../../output/flood2d-output-pkg-2/txt ../../output/flood2d-compare