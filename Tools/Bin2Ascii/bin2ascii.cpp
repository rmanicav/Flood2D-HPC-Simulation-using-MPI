#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <dirent.h>
#include <cstring>
#include <sys/stat.h>

using namespace std;

void binToAscii(string str, string input_dir, string output_dir, int row, int col, int ghost_cell){
    string name1 = str, name2 = str.substr(0, str.length()-4) + ".txt";
	string input_file = input_dir + name1;
	string output_file = output_dir + name2;

	row = row + ghost_cell * 2;
	col = col + ghost_cell * 2;
	
    ifstream input(input_file.c_str(), ios::binary);

    float *arr = new float [row*col];

    input.read( (char*) arr, sizeof(float) * row * col );
    input.close();

    ofstream output(output_file.c_str());
    output << std::setprecision(4) << std::fixed;

    for(int i=ghost_cell; i<row-ghost_cell; i++){
        for(int j=ghost_cell; j<col-ghost_cell; j++){
            output << arr[i*col+j];
            if(j<col-1-ghost_cell){
                output << " ";
            }
        }
        output << endl;
    }
    output.close();

    delete[] arr;
}

int main(int argc, char* argv[])
{			
	string input_dir(argv[1]);
	string output_dir(argv[2]);
	int row = atoi(argv[3]);
	int col = atoi(argv[4]);
	int ghost_cell = atoi(argv[5]);
	
	DIR* dir_ = opendir(output_dir.empty() ? "." : output_dir.c_str());
	if( !dir_) {
		mkdir(output_dir.c_str(), S_IRWXU);
	}
	else closedir(dir_);
	

    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (input_dir.c_str())) != NULL) {
        while ((ent = readdir (dir)) != NULL) {
            if( strcmp(ent->d_name, ".") != 0 && strcmp(ent->d_name, "..") != 0 ){
                string str (ent->d_name);
                cout << "Reading " << str << endl;
                binToAscii(str, input_dir, output_dir, row, col, ghost_cell);
            }
        }
        closedir (dir);
    } else {
        perror ("");
        return 1;
    }

    return 0;
}
