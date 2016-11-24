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
      	
    p.theory.G_S = 4.855;         // (fm)^2
    p.theory.G_V = 4.65;          // (fm)^2
    p.theory.G_SV = -6.583;       // (fm)^8
    p.theory.G_RHO = 0.5876;      // (fm)^2
    p.theory.G_VRHO = 0.0;        // (fm)^8
    p.theory.G_SRHO = 0.0;        // (fm)^8
    p.theory.cutoff = 388.189;    // (MeV)
    p.theory.bare_mass = 0.0;     // (MeV)

  	AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	
    p.parameters_set_identifier = "eNJL1OmegaRho1";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.theory.G_S = 4.855;
    p.theory.G_V = 4.65;
    p.theory.G_SV = -6.583;
    p.theory.G_RHO = 0.5976;
    p.theory.G_VRHO = -1.0;
    p.theory.G_SRHO = 0.0;
    p.theory.cutoff = 388.189;
    p.theory.bare_mass = 0.0;
    
    AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	    
    p.parameters_set_identifier = "eNJL1OmegaRho2";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      
    p.theory.G_S = 4.855;
    p.theory.G_V = 4.65;
    p.theory.G_SV = -6.583;
    p.theory.G_RHO = 0.6476;
    p.theory.G_VRHO = -6.0;
    p.theory.G_SRHO = 0.0;
    p.theory.cutoff = 388.189;
    p.theory.bare_mass = 0.0;
    
    AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	
    p.parameters_set_identifier = "eNJL2";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.theory.G_S = 3.8;
    p.theory.G_V = 3.8;
    p.theory.G_SV = -4.228;
    p.theory.G_RHO = 0.6313;
    p.theory.G_VRHO = 0.0;
    p.theory.G_SRHO = 0.0;
    p.theory.cutoff = 422.384;
    p.theory.bare_mass = 0.0;
    
    AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	            
    p.parameters_set_identifier = "eNJL2OmegaRho1";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.theory.G_S = 3.8;
    p.theory.G_V = 3.8;
    p.theory.G_SV = -4.288;
    p.theory.G_RHO = 0.6413;
    p.theory.G_VRHO = -1.0;
    p.theory.G_SRHO = 0.0;
    p.theory.cutoff = 422.384;
    p.theory.bare_mass = 0.0;
    
    AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	            
    p.parameters_set_identifier = "eNJL3";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.theory.G_S = 1.93;
    p.theory.G_V = 3.0;
    p.theory.G_SV = -1.8;
    p.theory.G_RHO = 0.65;
    p.theory.G_VRHO = 0.0;
    p.theory.G_SRHO = 0.0;
    p.theory.cutoff = 534.815;
    p.theory.bare_mass = 0.0;
     
    AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	           
    p.parameters_set_identifier = "eNJL3SigmaRho1";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.theory.G_S = 1.93;
    p.theory.G_V = 3.0;
    p.theory.G_SV = -1.8;
    p.theory.G_RHO = 0.0269;
    p.theory.G_VRHO = 0.0;
    p.theory.G_SRHO = 0.5;
    p.theory.cutoff = 534.815;
    p.theory.bare_mass = 0.0;
     
    AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	           
    p.parameters_set_identifier = "eNJL1m";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.theory.G_S = 1.3833;
    p.theory.G_V = 1.781;
    p.theory.G_SV = -2.943;
    p.theory.G_RHO = 0.7;
    p.theory.G_VRHO = 0.0;
    p.theory.G_SRHO = 0.0;
    p.theory.cutoff = 478.248;
    p.theory.bare_mass = 450.0;
    
    AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	            
    p.parameters_set_identifier = "eNJL1mSigmaRho1";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.theory.G_S = 1.3833;
    p.theory.G_V = 1.781;
    p.theory.G_SV = -2.943;
    p.theory.G_RHO = 0.0739;
    p.theory.G_VRHO = 0.0;
    p.theory.G_SRHO = 1.0;
    p.theory.cutoff = 478.248;
    p.theory.bare_mass = 450.0;
    
    AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	            
    p.parameters_set_identifier = "eNJL2m";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.theory.G_S = 1.078;
    p.theory.G_V = 1.955;
    p.theory.G_SV = -2.74;
    p.theory.G_RHO = 0.75;
    p.theory.G_VRHO = 0.0;
    p.theory.G_SRHO = 0.0;
    p.theory.cutoff = 502.466;
    p.theory.bare_mass = 500.0;
    
    AppendParametersSetToList(p);
    
  	////

  	p = NewCopyOfParametersSetFromTemplate();
  	            
    p.parameters_set_identifier = "eNJL2mSigmaRho1";
    p.parameters_set_origin = "Helena Pais, Débora P. Menezes, and Constança Providência, Phys. Rev. C 93, 065805 – Published 8 June 2016";
      	
    p.theory.G_S = 1.078;
    p.theory.G_V = 1.955;
    p.theory.G_SV = -2.74;
    p.theory.G_RHO = -0.1114;
    p.theory.G_VRHO = 0.0;
    p.theory.G_SRHO = 1.0;
    p.theory.cutoff = 502.466;
    p.theory.bare_mass = 500.0;
    
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

  	p.parameters_set_origin = "A standard parameters set. With no theory parameters.";
    p.parameters_set_identifier = "Template";

	p.minimum_density = 1.0E-2; // fm^-3
	p.maximum_density = 0.74; // fm^-3

	p.points_number = 3000;
	
    p.temperature = 0.0;		// (MeV)
	p.proton_fraction = 0.5;	// (no dimension)

	p.gap_equation_solver.max_iterations = 1000;
	p.gap_equation_solver.lower_bound = 1.0E-5; // Low, but not zero. In zero f(M) = 0;
	p.gap_equation_solver.upper_bound = 3000.0;	// MeV (much bigger than expected maximum mass)
	p.gap_equation_solver.abs_error = 1E-10;
	p.gap_equation_solver.rel_error = 1E-10;

	p.theory.G_S = 0;              // (fm)^2
    p.theory.G_V = 0;              // (fm)^2
    p.theory.G_SV = 0;             // (fm)^8
    p.theory.G_RHO = 0;            // (fm)^2
    p.theory.G_VRHO = 0.0;         // (fm)^8
    p.theory.G_SRHO = 0.0;         // (fm)^8
    p.theory.cutoff = 0.0;         // (MeV)
    p.theory.nucleon_mass = 939.0; // (MeV) // FIXME: this should be the bare_mass!
    p.theory.electron_mass = 0.511;// (MeV)
    p.theory.bare_mass = 0.0;      // (MeV)

    p.multiroot.max_iterations = 3000;
    p.multiroot.proton_fraction_mapping_scale = 100;
    p.multiroot.use_last_solution_as_guess = true;

    p.multiroot.abs_error = 1E-7;
    p.multiroot.rel_error = 1E-4;
    
    p.multiroot.guesses.mass = p.theory.nucleon_mass;
    p.multiroot.guesses.proton_fraction = 0.02;

    p.multiroot.special_case.mass_tolerance = 0.1;
    p.multiroot.special_case.lower_bound = 0.1; // The mapping with tanh(x) prevents
    p.multiroot.special_case.upper_bound = 0.9; // the use of extremes of the proton fraction
    
    p.finite_temperature.use_last_solution_as_guess = true;
    p.finite_temperature.guesses.mass = 1000.0; // MeV
    p.finite_temperature.guesses.proton_renormalized_chemical_potential = 1000.0; // MeV
    p.finite_temperature.guesses.neutron_renormalized_chemical_potential = 1000.0; // MeV
    p.finite_temperature.special_case.mass_tolerance = 0.1;
    
    p.fermi_dirac_integrals.max_interval_num = 1000;
    p.fermi_dirac_integrals.abs_error = 1E-5;
    p.fermi_dirac_integrals.rel_error = 1E-5;
    p.fermi_dirac_integrals.max_sub_interval = 1000;
    p.fermi_dirac_integrals.integration_key = GSL_INTEG_GAUSS61;
    
    p.entropy_integrals.interval_num = 1000;
    p.entropy_integrals.abserr = 1E-5;
    p.entropy_integrals.lower_limit = 0.0;
    p.entropy_integrals.upper_limit = parameters.theory.cutoff;
    p.entropy_integrals.abs_error = 1.0E-3;
    p.entropy_integrals.rel_error = 1.0E-3;
    p.entropy_integrals.max_sub_interval = 1000;
    p.entropy_integrals.integration_key = GSL_INTEG_GAUSS61;
    
  	return p;
}

void AppendParametersSetToList(Parameters a_set)
{
    // Append to list of parameterizations
    parameters_sets_list[parameters_sets_list_count] = a_set;
    parameters_sets_list_count++;
    
    return;
}

void SetParametersFromCommandline()
{
    if (options.temperature_override){
        if (options.temperature < 0){
            printf("Temperature must be a non-negative value (%f was provided).\n",
                   options.temperature);
            exit(EXIT_FAILURE);
        }

        parameters.temperature = options.temperature;
    }

    if (options.proton_fraction_override){
        if (options.proton_fraction < 0.0 || options.proton_fraction > 1.0){
            printf("The proton fraction must lie in the [0,1] range (%f was provided).\n",
                   options.proton_fraction);
            exit(EXIT_FAILURE);
        }

        parameters.proton_fraction = options.proton_fraction;
    }
}

void PrintParametersToFile(FILE * file)
{
    fprintf(file, "\n--- PARAMETERS ---\n");
    fprintf(file, "parameters_set_origin = %s\n", parameters.parameters_set_origin);
    fprintf(file, "parameters_set_identifier = %s\n\n", parameters.parameters_set_identifier);
    
    fprintf(file, "nucleon_mass = %f\n", parameters.theory.nucleon_mass);

	fprintf(file, "minimum_density = %f\n", parameters.minimum_density);
	fprintf(file, "maximum_density = %f\n", parameters.maximum_density);

	fprintf(file, "points_number = %d\n", parameters.points_number);

	fprintf(file, "proton_fraction = %f\n", parameters.proton_fraction);
	
	fprintf(file, "G_S = %f\n", parameters.theory.G_S);
    fprintf(file, "G_V = %f\n", parameters.theory.G_V);
    fprintf(file, "G_SV = %f\n", parameters.theory.G_SV);
    fprintf(file, "G_RHO = %f\n", parameters.theory.G_RHO);
    fprintf(file, "G_VRHO = %f\n", parameters.theory.G_VRHO);
    fprintf(file, "G_SRHO = %f\n", parameters.theory.G_SRHO);
    fprintf(file, "CUTOFF = %f\n", parameters.theory.cutoff);
    fprintf(file, "bare_mass = %f\n", parameters.theory.bare_mass);
    fprintf(file, "temperature = %fzn", parameters.temperature);
    
    return;
}

