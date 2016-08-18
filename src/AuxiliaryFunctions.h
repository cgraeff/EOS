//
//  AuxiliaryFunctions.h
//  hadrons EOS
//
//  Created by Clebson Graeff on 2016-08-16.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#ifndef AuxiliaryFunctions_h
#define AuxiliaryFunctions_h

#define FILENAME_MAX_SIZE 256
#define PATH_MAX_SIZE 1024
#define FILEPATH_MAX_SIZE PATH_MAX_SIZE + FILENAME_MAX_SIZE

#include <gsl/gsl_vector.h>

int WriteVectorsToFile(const char * filename, const char * header, int vectors_count, ...);
int WriteVectorsToFileUpToIndex(const char * filename, const char * header, int vector_index, int vectors_count, ...);
int WriteIndexedVectorsToFile(const char * filename, const char * header, int vectors_count, ...);

gsl_vector * VectorNewVectorFromDivisionElementByElement(gsl_vector * numerator, gsl_vector * denominator);

gsl_vector * VectorNewVectorFromDivisionElementByElement(gsl_vector * numerator, gsl_vector * denominator);

// Simple functions to open files and creating directories if necessary.
// The path should be defined by SetFilePath(). To unset, just call
// SetFilePath(NULL). The total size of path and filename combined
// must be smaller than FILEPATH_MAX_SIZE characters. Additionally,
// each path element must be smaller than PATH_ELEMENT_MAX_SIZE
// characters.
//
// Neither expansion of '~', nor of '.', '..' are supported.
// Examples of tested cases are
//
//      a/path/to/dir
//      a/path/to/dir/
//
// After using a path prefix, set it to NULL so any following code that
// expects an 'empty' path can work as expected.
void SetFilePath(const char path[]);
FILE * OpenFile(const char filename[]);

#endif /* AuxiliaryFunctions_h */
