//
//  Parameters.c
//  hadrons EOS
//
//  Created by Clebson Graeff on 2016-02-16.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_integration.h>

#include "Constants.h"
#include "Parameters.h"
#include "CommandlineOptions.h"

// Chosen parameter set (globally accessible)
Parameters parameters;
static Parameters parameters_sets_list[256] = {0};
static int parameters_sets_list_count = 0;

Parameters NewCopyOfParametersSetFromTemplate();
void AppendParametersSetToList(Parameters a_set);

void ParametersSetup(void)
{
    // See header file for units and other relevant
    // information about the parameters below.

  	// START DECLARATION OF PARAMETERS SETS:
	//
  	// 	For declaration of parameters sets:
	// 		- copy from the default one,
	// 		- overwrite values as desired,
	// 		- append to the list of parameterizations
    //
    // The first declared (and appended) parameters set will be treated as the standard
    // case. To run others, use the "-p IDENTIFIER" commandline option, without
    // quotation marks) where IDENTIFIER is the content of the variable
    // parameters_set_identifier in the definitions bellow.

	Parameters p;

  	p = NewCopyOfParametersSetFromTemplate();

    p.parameters_set_identifier = "eNJL1";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.parameters.G_S = 4.855;         // (fm)^2
    p.parameters.G_V = 4.65;          // (fm)^2
    p.parameters.G_SV = -6.583;       // (fm)^8
    p.parameters.G_RHO = 0.5876;      // (fm)^2
    p.parameters.G_VRHO = 0.0;        // (fm)^8
    p.parameters.G_SRHO = 0.0;        // (fm)^8
    p.parameters.CUTOFF = 388.189;    // (MeV)
    p.parameters.bare_mass = 0.0;     // (MeV)

  	AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	
    p.parameters_set_identifier = "eNJL1OmegaRho1";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.parameters.G_S = 4.855;
    p.parameters.G_V = 4.65;
    p.parameters.G_SV = -6.583;
    p.parameters.G_RHO = 0.5976;
    p.parameters.G_VRHO = -1.0;
    p.parameters.G_SRHO = 0.0;
    p.parameters.CUTOFF = 388.189;
    p.parameters.bare_mass = 0.0;
    
    AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	    
    p.parameters_set_identifier = "eNJL1OmegaRho2";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.parameters.G_S = 4.855;
    p.parameters.G_V = 4.65;
    p.parameters.G_SV = -6.583;
    p.parameters.G_RHO = 0.6476;
    p.parameters.G_VRHO = -6.0;
    p.parameters.G_SRHO = 0.0;
    p.parameters.CUTOFF = 388.189;
    p.parameters.bare_mass = 0.0;
    
    AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	
    p.parameters_set_identifier = "eNJL2";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.parameters.G_S = 3.8;
    p.parameters.G_V = 3.8;
    p.parameters.G_SV = -4.228;
    p.parameters.G_RHO = 0.6313;
    p.parameters.G_VRHO = 0.0;
    p.parameters.G_SRHO = 0.0;
    p.parameters.CUTOFF = 422.384;
    p.parameters.bare_mass = 0.0;
    
    AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	            
    p.parameters_set_identifier = "eNJL2OmegaRho1";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.parameters.G_S = 3.8;
    p.parameters.G_V = 3.8;
    p.parameters.G_SV = -4.288;
    p.parameters.G_RHO = 0.6413;
    p.parameters.G_VRHO = -1.0;
    p.parameters.G_SRHO = 0.0;
    p.parameters.CUTOFF = 422.384;
    p.parameters.bare_mass = 0.0;
    
    AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	            
    p.parameters_set_identifier = "eNJL3";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.parameters.G_S = 1.93;
    p.parameters.G_V = 3.0;
    p.parameters.G_SV = -1.8;
    p.parameters.G_RHO = 0.65;
    p.parameters.G_VRHO = 0.0;
    p.parameters.G_SRHO = 0.0;
    p.parameters.CUTOFF = 534.815;
    p.parameters.bare_mass = 0.0;
     
    AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	           
    p.parameters_set_identifier = "eNJL3SigmaRho1";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.parameters.G_S = 1.93;
    p.parameters.G_V = 3;
    p.parameters.G_SV = -1.8;
    p.parameters.G_RHO = 0.0269;
    p.parameters.G_VRHO = 0.0;
    p.parameters.G_SRHO = 0.5;
    p.parameters.CUTOFF = 534.815;
    p.parameters.bare_mass = 0.0;
     
    AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	           
    p.parameters_set_identifier = "eNJL1m";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.parameters.G_S = 1.3833;
    p.parameters.G_V = 1.781;
    p.parameters.G_SV = -2.943;
    p.parameters.G_RHO = 0.7;
    p.parameters.G_VRHO = 0.0;
    p.parameters.G_SRHO = 0.0;
    p.parameters.CUTOFF = 478.248;
    p.parameters.bare_mass = 450.0;
    
    AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	            
    p.parameters_set_identifier = "eNJL1mSigmaRho1";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.parameters.G_S = 1.3833;
    p.parameters.G_V = 1.781;
    p.parameters.G_SV = -2.943;
    p.parameters.G_RHO = 0.0739;
    p.parameters.G_VRHO = 0.0;
    p.parameters.G_SRHO = 1.0;
    p.parameters.CUTOFF = 478.248;
    p.parameters.bare_mass = 450.0;
    
    AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	            
    p.parameters_set_identifier = "eNJL2m";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.parameters.G_S = 1.078;
    p.parameters.G_V = 1.955;
    p.parameters.G_SV = -2.74;
    p.parameters.G_RHO = 0.75;
    p.parameters.G_VRHO = 0.0;
    p.parameters.G_SRHO = 0.0;
    p.parameters.CUTOFF = 502.466;
    p.parameters.bare_mass = 500.0;
    
    AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	            
    p.parameters_set_identifier = "eNJL2mSigmaRho1";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.parameters.G_S = 1.078;
    p.parameters.G_V = 1.955;
    p.parameters.G_SV = -2.74;
    p.parameters.G_RHO = 0.1114;
    p.parameters.G_VRHO = 0.0;
    p.parameters.G_SRHO = 1.0;
    p.parameters.CUTOFF = 502.466;
    p.parameters.bare_mass = 500.0;
    
    AppendParametersSetToList(p);
    
  	// END DECLARATION OF PARAMETERS SETS
    
    // Verify that the set identifiers are unique
    for (int i = 0; i < parameters_sets_list_count; i++){
        for (int j = 0; j < parameters_sets_list_count; j++){
            
            if (i == j)
                continue;
            
            if(!strcasecmp(parameters_sets_list[i].parameters_set_identifier,
                           parameters_sets_list[j].parameters_set_identifier)){
                printf("Two parameters sets share the \"%s\" identifier, but it should be unique.\n",
                       parameters_sets_list[i].parameters_set_identifier);
                exit(EXIT_FAILURE);
            }
        }
    }

  	// If asked to, print parameters sets and exit
  	if (options.list_available_parameterizations){
		for (int i = 0; i < parameters_sets_list_count; i++){
		  	printf("Parameters set %s\n"
				   "\tOrigin: %s\n",
				   parameters_sets_list[i].parameters_set_identifier,
				   parameters_sets_list[i].parameters_set_origin);
		}
		exit(EXIT_SUCCESS);
	}
    
    // If there isn't at least one set, exit
    if (parameters_sets_list_count == 0){
        printf("There are not parameters set declared. Declare at least one set.\n");
        exit(EXIT_FAILURE);
    }
}

void SetParametersSet(char * parameters_set_identifier)
{
    // If the identifier is null, return default case
    if (parameters_set_identifier == NULL){
        parameters = parameters_sets_list[0];
        return;
    }
        
    for (int i = 0; i < parameters_sets_list_count; i++){
          
        Parameters p = parameters_sets_list[i];
                          
        if (!strcasecmp(p.parameters_set_identifier, parameters_set_identifier)){
			parameters = p;
            return;
		}
    }
    
    printf("Parameters set %s unrecognized.\n"
           "Use -l to list available parameterizations.\n",
           parameters_set_identifier);
    exit(EXIT_FAILURE);
}

Parameters NewCopyOfParametersSetFromTemplate()
{
  	Parameters p;

  	p.parameters_set_origin = "A standard parameters set. With theory parameters from eNJL1"
    						  "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
    p.parameters_set_identifier = "Template";

    parameters.nucleon_mass = 939.0; // (MeV)

	parameters.minimum_density = 0.01; // fm^-3
	parameters.maximum_density = 0.35; // fm^-3

	parameters.points_number = 1000;

	parameters.proton_fraction = 0.5;
	
	parameters.number_of_colors = 1;
	parameters.number_of_flavors = 1;
	
	parameters.degeneracy = (double)parameters.number_of_colors * (double)parameters.number_of_flavors / pow(M_PI, 2.0);
	
	parameters.lower_bound_gap_eq_solve = 1.0E-3; // Low, but not zero. In zero f(M) = 0;
	parameters.upper_bound_gap_eq_solve = 3000.0;	// MeV (about 150% of nucleon mass)
	
	parameters.G_S = 4.855;         // (fm)^2
    parameters.G_V = 4.65;          // (fm)^2
    parameters.G_SV = -6.583;       // (fm)^8
    parameters.G_RHO = 0.5876;      // (fm)^2
    parameters.G_VRHO = 0.0;        // (fm)^8
    parameters.G_SRHO = 0.0;        // (fm)^8
    parameters.CUTOFF = 388.189;    // (MeV)
    parameters.bare_mass = 0.0;     // (MeV)
    
  	return p;
}

void AppendParametersSetToList(Parameters a_set)
{
    // Append to list of parameterizations
    parameters_sets_list[parameters_sets_list_count] = a_set;
    parameters_sets_list_count++;
    
    return;
}

void PrintParametersToFile(FILE * file)
{
    fprintf(file, "\n--- PARAMETERS ---\n");
    fprintf(file, "parameters_set_origin = %s\n", parameters.parameters_set_origin);
    fprintf(file, "parameters_set_identifier = %s\n\n", parameters.parameters_set_identifier);
    
    fprintf(file, "nucleon_mass = %f\n", parameters.nucleon_mass);

	fprintf(file, "minimum_density = %f\n", parameters.minimum_density);
	fprintf(file, "maximum_density = %f\n", parameters.maximum_density);

	fprintf(file, "points_number = %d\n", parameters.points_number);

	fprintf(file, "proton_fraction = %f\n", parameters.proton_fraction);
	
	fprintf(file, "number_of_colors = %d\n", parameters.number_of_colors);
	fprintf(file, "number_of_flavors = %d\n", parameters.number_of_flavors);
	
	fprintf(file, "lower_bound_gap_eq_solve = %f\n", parameters.lower_bound_gap_eq_solve);
	fprintf(file, "upper_bound_gap_eq_solve = %f\n", upper_bound_gap_eq_solve);
	
	fprintf(file, "G_S = %f\n", parameters.G_S);
    fprintf(file, "G_V = %f\n", parameters.G_V);
    fprintf(file, "G_SV = %f\n", parameters.G_SV);
    fprintf(file, "G_RHO = %f\n", parameters.G_RHO);
    fprintf(file, "G_VRHO = %f\n", parameters.G_VRHO);
    fprintf(file, "G_SRHO = %f\n", parameters.G_SRHO);
    fprintf(file, "CUTOFF = %f\n", parameters.CUTOFF);
    fprintf(file, "bare_mass = %f\n", parameters.bare_mass);
    
    return;
}

