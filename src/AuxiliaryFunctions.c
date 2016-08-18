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

static const char * open_file_prefix_path = NULL;

FILE * OpenFile(const char filename[])
{
    /*
     * Check filename for '/' or overflow
     */
    if (NULL != strrchr(filename, '/')){
        printf("OpenFile: Filename contains a path element. This is not supported.\n");
        abort();
    }

    if (strlen(filename) > FILENAME_MAX_SIZE){
        printf("OpenFile: Filename is too long.\n");
        abort();
    }

    // If there is not prefix path,
    // just open the file
    if (open_file_prefix_path == NULL){
        FILE * file = fopen(filename, "w");

        if (NULL == file) {
            printf("Could not open %s for writting.\n", filename);
            perror("Reason");
            exit(EXIT_FAILURE);
        }
        
        return file;
    }
    
    /*
     * Recursivelly create dirs in path
     */
    char tmp[PATH_MAX_SIZE];

    int length = snprintf(tmp, sizeof(tmp),"%s", open_file_prefix_path);

    // Check for final slash in path
    // add if it isn't there
    if (tmp[length - 1] != '/'){

        // Check if bounds will be exceded
        // (the + 1 account for the /)
        if (length + 1 >= PATH_MAX_SIZE){
            printf("OpenFile: Path length is too long.\n");
            abort();
        }
        tmp[length] = '/';
        tmp[length + 1] = 0; // End of string must be null
    }

    // Create dirs:
    //  Reads the string sequentially replacing occurences of '/' with null ('\0')
    //  (this will effectivelly limit the string to characters before the '/' being replaced)
    //  checking if the dir pointed by de resulting string exists, creating if it doesn't exist,
    //  restoring the modified '/' (and moving to the next occurence of '/'. The loop will
    //  continue ultil the pointer poinst to the final null character, whose value can be
    //  obtained by deference and is zero, ending the loop.
    char *pointer = NULL;
    for(pointer = tmp; *pointer; pointer++)
        if(*pointer == '/') {
            *pointer = 0;
            struct stat st = {0};
            if (stat(tmp, &st) == -1)
                mkdir(tmp, S_IRWXU);
            *pointer = '/';
        }
    
    // Append filename to path.
    // Use tmp as it will have the final '/'
    char filepath[FILEPATH_MAX_SIZE];
    strcpy(filepath, tmp);

    unsigned long path_length = strlen(filepath);
    unsigned long filename_length = strlen(filename);
    unsigned long total_length = path_length + filename_length;

    // Check if bounds will be exceded
    if (total_length >= FILEPATH_MAX_SIZE){
        printf("OpenFile: Path + filename length is too long.\n");
        abort();
    }

    strcat(filepath, filename);

    // Finally, open the file
    FILE * file = fopen(filepath, "w");
    
    if (NULL == file) {
        printf("OpenFile: Could not open %s for writting.\n", filename);
        perror("Reason");
        exit(EXIT_FAILURE);
    }
    
    return file;
}

void SetFilePath(const char path[])
{
    if (NULL == path){
      open_file_prefix_path = NULL;
      return;
    }

    switch (path[0]){
    case '~':
        printf("SetFilePath: ~ expansion is not supported.\n");
        abort();
        break;
    case '.':
        printf("SetFilePath: . is unnecessary, .. is unsuported.\n");
    }

    open_file_prefix_path = path;
}
