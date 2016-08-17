//
//  AuxiliaryFunctions.c
//  hadrons EOS
//
//  Created by Clebson Graeff on 2016-08-16.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <sys/stat.h>

#include "AuxiliaryFunctions.h"

int WriteVectorsToFileUpToIndex(const char * filename, const char * header, int vector_index, int vectors_count, ...)
{
    FILE * output = OpenFile(filename);
    
    fprintf(output, "%s", header);
    
    va_list arg_list;
    va_start(arg_list, vectors_count);
    
    gsl_vector * vectors[vectors_count];
    
    for (int i = 0; i < vectors_count; i++){
        gsl_vector * v = va_arg(arg_list, gsl_vector *);
        vectors[i] = v;
    }
    
    va_end(arg_list);
    
    if (vectors_count > 1)
        for (int i = 0; i < vectors_count - 1; i++)
            if ((vectors[i])->size != (vectors[i + 1])->size) {
                
                printf("ERROR: Vectors have different sizes.\n");
                exit(EXIT_FAILURE);
            }
    
    for (int i = 0; i < vector_index; i++) {
        for (int j = 0; j < vectors_count; j++) {
            
            double x = gsl_vector_get(vectors[j], i);
            
            fprintf(output, "%20.15E", x);
            
            if (j != vectors_count - 1)
                fprintf(output, "\t");
        }
        
        fprintf(output, "\n");
    }
    
    fclose(output);
    
    return 0;
}

int WriteVectorsToFile(const char * filename, const char * header, int vectors_count, ...)
{
    FILE * output = OpenFile(filename);
	
	fprintf(output, "%s", header);
	
	va_list arg_list;
	va_start(arg_list, vectors_count);
	
	gsl_vector * vectors[vectors_count];
	
	for (int i = 0; i < vectors_count; i++){
		gsl_vector * v = va_arg(arg_list, gsl_vector *);
		vectors[i] = v;
	}
	
	va_end(arg_list);
	
	if (vectors_count > 1)
		for (int i = 0; i < vectors_count - 1; i++)
			if ((vectors[i])->size != (vectors[i + 1])->size) {
				
				printf("ERROR: Vectors have different sizes.\n");
				exit(EXIT_FAILURE);
			}
	
	for (int i = 0; i < (vectors[0])->size; i++) {
		for (int j = 0; j < vectors_count; j++) {
			
			double x = gsl_vector_get(vectors[j], i);
			
			fprintf(output, "%20.15E", x);
			
			if (j != vectors_count - 1)
				fprintf(output, "\t");
		}
		
		fprintf(output, "\n");
	}
	
	fclose(output);
	
	return 0;
}

int WriteIndexedVectorsToFile(const char * filename, const char * header, int vectors_count, ...)
{
    FILE * output = OpenFile(filename);
	
	fprintf(output, "%s", header);
	
	va_list arg_list;
	va_start(arg_list, vectors_count);
	
	gsl_vector * vectors[vectors_count];
	
	for (int i = 0; i < vectors_count; i++) {
		vectors[i] = va_arg(arg_list, gsl_vector *);
	}
	
	va_end(arg_list);
	
	if (vectors_count > 1)
		for (int i = 0; i < vectors_count - 1; i++)
			if ((vectors[i])->size != (vectors[i + 1])->size) {
				
				printf("ERROR: Vectors have different sizes.\n");
				exit(EXIT_FAILURE);
			}
	
	for (int i = 0; i < (vectors[0])->size; i++) {
		
		fprintf(output, "%d\t", i);
		
		for (int j = 0; j < vectors_count; j++) {
			
			fprintf(output, "%20.15E", gsl_vector_get(vectors[j], i));
			
			if (j != vectors_count - 1)
				fprintf(output, "\t");
		}
		
		fprintf(output, "\n");
	}
	
	fclose(output);
	
	return 0;
}

gsl_vector * VectorNewVectorFromDivisionElementByElement(gsl_vector * numerator, gsl_vector * denominator)
{
	if (numerator->size != denominator->size) {
		printf("The vectors must have the same size!\n");
		exit(EXIT_FAILURE);
	}
	
	gsl_vector * v = gsl_vector_alloc(numerator->size);
	
	for (int i = 0; i < numerator->size; i++){
		double value = gsl_vector_get(numerator, i) / gsl_vector_get(denominator, i);
		gsl_vector_set(v, i, value);
	}

	return v;
}

FILE * OpenFile(const char filename[])
{
    // If there is no slash, it means that
    // the file must me in the current dir
    char * last_slash = strrchr(filename, '/');
    
    if (last_slash == NULL){
        FILE * file = fopen(filename, "w");
        
        if (NULL == file) {
            printf("Could not open %s for writting.\n", filename);
            perror("Reason");
            exit(EXIT_FAILURE);
        }
        
        return file;
    }
    
    // Recursivelly create dirs in path
    char tmp[256];
    char *p = NULL;
    
    snprintf(tmp, sizeof(tmp),"%s", filename);

    for(p = tmp + 1; *p; p++)
        if(*p == '/') {
            *p = 0;
            struct stat st = {0};
            if (stat(tmp, &st) == -1)
                mkdir(tmp, S_IRWXU);
            *p = '/';
        }
    
    // Finally, open the file
    FILE * file = fopen(filename, "w");
    
    if (NULL == file) {
        printf("Could not open %s for writting.\n", filename);
        perror("Reason");
        exit(EXIT_FAILURE);
    }
    
    return file;
}
