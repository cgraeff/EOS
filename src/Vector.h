//
//  Vector.h
//  EOS
//
//  Created by Clebson Graeff on 2/16/16.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#ifndef Vector_h
#define Vector_h

#include <gsl/gsl_vector.h>

int WriteVectorsToFile(const char * filename, const char * header, int vectors_count, ...);
int WriteIndexedVectorsToFile(const char * filename, const char * header, int vectors_count, ...);
gsl_vector * VectorNewVectorFromDivisionElementByElement(gsl_vector * numerator, gsl_vector * denominator);


#endif /* Vector_h */
