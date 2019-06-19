//============================================================================
// Name        : flood-output-rasterize.cpp
// Author      : R. Marshall
// Version     :
// Copyright   : 2015
// Description : Generate animated gif from Flood2D output files
//============================================================================
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <dirent.h>
#include <getopt.h>
#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include "DEMFile.h"
#include "RasterMatrix.h"

#define BGTYPE ".bmp"
#define BASE_LAYER_NAME "_baselayer" BGTYPE

#define CONTAINS(X, Y) (std::find(X.begin(), X.end(), Y) != X.end())
#define HAS_FILE(X, Y) (X.find(Y) != std::string::npos)

enum Format {
  PLAIN_ASCII   = 0, // space-delimited ASCII grid 
  RC_ASCII      = 1, // plain ASCII, with first row: R C 
  ARCINFO_ASCII = 2  // ArcInfo/ASCII Grid
} ;



typedef float value_t;
typedef int dim;
typedef std::vector<std::string>::iterator file_list_itr;
typedef std::vector<std::string> file_list_t;

struct opts_t{
  int nrows, ncols, fgcolors, bgcolors;

  std::string datadir, bgdata, prefix;
  Format bgformat;

  void print() {
    fprintf(stderr, "Options:\n\t");
    fprintf(stderr, "nrows=%d\n\tncols=%d\n\tfgcolors=%d\n\tbgcolors=%d\n\t",
      nrows, ncols, fgcolors, bgcolors
    );
    
    fprintf(stderr, "datadir=%s\n\tbgfile=%s\n\tprefix=%s\n\tbgformat=%s\n",
      datadir.c_str(), bgdata.c_str(), prefix.c_str(), "ARCINFO_ASCII"
    );
  }
};

using namespace std;

file_list_t get_file_list(std::string dir_path, std::string prefix)
{
  file_list_t file_list;
  DIR *dir;
  struct dirent *ent;

  if ((dir = opendir (dir_path.c_str())) != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      bool match = false;

      for(unsigned i=0; i< prefix.size(); ++i)
      {
        match = ((ent->d_name)[i] == prefix[i]);

        if ( ! match) break;
      }

      if (match)
        file_list.push_back(ent->d_name);
    }
    closedir (dir);
  } else {
    perror ("could not open directory");
  }

  return file_list;
}


////////////////

void print_usage(opts_t opts) {
    fprintf(stderr, "Usage: f2d_rasterize ");
    fprintf(stderr, " --nrows num --ncols num");
    fprintf(stderr, "[fgcolors=%d bgcolors=%d datadir=%s bgfile=%s prefix=%s bgformat=%s]\n",
      opts.fgcolors, opts.bgcolors, 
      (opts.datadir).c_str(), (opts.bgdata).c_str(), "H", "ARCINFO_ASCII"
    );
}



opts_t
parseopts (int argc, char **argv)
{
  opts_t opts;
  int opt = 0, idx = 0;

  // default values
  opts.nrows = -1;
  opts.ncols = -1;
  opts.fgcolors = -1;
  opts.bgcolors = -1;

  opts.datadir = "../output";
  opts.bgdata = "../input/dem/ts_dam_break.dem";
  opts.prefix = "H";  
  opts.bgformat = ARCINFO_ASCII;

  static struct option long_options[] = {
    {"nrows",     required_argument, 0, 'R' },
    {"ncols",     required_argument, 0, 'C' },
    {"fgcolors",  optional_argument, 0, 'f' },
    {"bgcolors",  optional_argument, 0, 'c' },
    {"datadir",   optional_argument, 0, 'd' },
    {"bgdata",    optional_argument, 0, 'B' },
    {"prefix",    optional_argument, 0, 'P' },    
    {"bgformat",  optional_argument, 0, 't' },
    {0,           0,                 0,  0  }
  };

  while ( 
    ((opt = getopt_long(argc, argv,"R:C:fcdBt", long_options, &idx )) != -1)
  ) {
    switch (opt) {
      case 'R' : opts.nrows      = atoi(optarg);  break;
      case 'C' : opts.ncols      = atoi(optarg);  break;
      case 'f' : opts.fgcolors   = atoi(optarg);  break;
      case 'c' : opts.bgcolors   = atoi(optarg);  break;
      case 'd' : opts.datadir    = optarg;        break;
      case 'B' : opts.bgdata     = optarg;        break;
      case 'P' : opts.prefix     = optarg;        break;      
      case 't' : 
        switch (atoi(optarg)) {
          case ARCINFO_ASCII : opts.bgformat = ARCINFO_ASCII; break;
          case PLAIN_ASCII   : opts.bgformat = PLAIN_ASCII;   break;
          case RC_ASCII      : opts.bgformat = RC_ASCII;      break;
          default            : opts.bgformat = ARCINFO_ASCII;
        }
      default: print_usage(opts); 
      exit(EXIT_FAILURE);
    }
  }
  
  for (int index = optind; index < argc; index++)
    fprintf (stderr, "Non-option argument %s\n", argv[index]);

  return opts;
}




void gen_bg_img(opts_t opts )
{
  std::string empty("");
  DEMFile baselayer(opts.bgdata);
  RasterMatrix<value_t> M(opts.nrows, opts.ncols);

  for (int i = 0; i < opts.nrows; ++i)
    for (int j = 0; j < opts.ncols; ++j) {
      M.set_value(i, j, baselayer.matrix_.get_val(i,j));
    }

  gradient_t<value_t> grdt = M.get_gradient(255, true);
  std::string baselayer_path( opts.datadir + "/" + opts.prefix + BASE_LAYER_NAME );
      std::cerr << "baselayer_path:"<< baselayer_path << std::endl;
  M.rasterize(baselayer_path , grdt, empty, opts.bgcolors, opts.fgcolors);
}


///////////////

int main(int argc, char* argv[]) 
{
  opts_t opts = parseopts(argc, argv);

  std::string bgimg_name = opts.prefix + BASE_LAYER_NAME;
  file_list_t files = get_file_list(opts.datadir, opts.prefix);

  std::sort(files.begin(), files.end());

  bool has_bmp = (CONTAINS(files, bgimg_name)) ;

  if ( ! has_bmp)
  {
    gen_bg_img(opts);
    std::sort(files.begin(), files.end());
    files = get_file_list(opts.datadir, opts.prefix);
  }
  
  has_bmp = (CONTAINS(files, bgimg_name)) ;

  if(has_bmp)
  {
    std::cerr << "Rasterizing..." << std::endl;

    for(file_list_itr it = files.begin(); it != files.end(); it++)
    {
      if(*it != bgimg_name) {
        std::cerr << *it << std::endl;
        std::string fullpath( opts.datadir + "/" + (*it) );
        std::string baselayer_path( opts.datadir +"/" + bgimg_name );
        
        RasterMatrix<value_t> M;
        M.load_from_file(opts.nrows, opts.ncols, fullpath);
        gradient_t<value_t> grdt;
        grdt.begin = -10.0;
        grdt.end = 10.0;
        grdt.step_size = 255;
        //M.get_gradient(255,false);

        fullpath.replace((fullpath.size())-4, 4, ".gif");

        M.rasterize(fullpath , M.get_gradient(255,false), baselayer_path, 
           opts.bgcolors, opts.fgcolors, M.get_palette_vec(opts.fgcolors, VIOLET_RED, PALE_GREEN_4));         
          //opts.bgcolors, opts.fgcolors, M.get_palette_vec(opts.fgcolors, NAVY_BLUE, CORNFLOWER_BLUE));
      }
    }
  }
  else 
  {
    perror ("could not find image or data file for background");
  }

  return 0;
}
