//
//  AuxiliaryFunctions.h
//  hadrons EOS
//
//  Created by Clebson Graeff on 2016-08-16.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#ifndef AuxiliaryFunctions_h
#define AuxiliaryFunctions_h

#include <gsl/gsl_vector.h>

int WriteVectorsToFile(const char * filename, const char * header, int vectors_count, ...);
int WriteVectorsToFileUpToIndex(const char * filename, const char * header, int vector_index, int vectors_count, ...);
int WriteIndexedVectorsToFile(const char * filename, const char * header, int vectors_count, ...);

gsl_vector * VectorNewVectorFromDivisionElementByElement(gsl_vector * numerator, gsl_vector * denominator);

// Simple function to open files and creating directories if necessary
// it is tested only with paths with format:
//      "this/is/a/path/to/my_file.txt"
// and with no path at all:
//      "my_file.txt"
FILE * OpenFile(const char filename[]);

#endif /* AuxiliaryFunctions_h */
