//
//  Parameters.h
//  hadrons EOS
//
//  Created by Clebson Graeff on 2016-02-16.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#ifndef Parameters_h
#define Parameters_h

#include <stdbool.h>

typedef struct _parameters
{
    char * parameters_set_identifier;
    char * parameters_set_origin;

	double minimum_density;		// (fm^-3)
	double maximum_density;		// (fm^-3)

    double temperature;			// MeV
	double proton_fraction;		// (no dimension)

	int points_number;

    struct _theory{
	    double G_S;				// scalar-isoscalar coupling (fm^2)
	    double G_V;				// vector-isoscalar coupling (fm^2)
	    double G_RHO;			// vector-isovector coupling (fm^2)
	    double G_SV;			// scalar-isovector coupling (fm^8)
	    double G_VRHO;			// (fm^8)
        double G_SRHO;			// (fm^8)

	    double cutoff;			// (MeV)
	    double nucleon_mass;	// (MeV)
        double electron_mass;	// (MeV)
        double bare_mass;		// (MeV)
    } theory;

	// Solution for zero temperature with
	// nucleons only
    struct _gap_equation_solver{
        int max_iterations;
        double lower_bound;
        double upper_bound;
        double abs_error;
        double rel_error;
    } gap_equation_solver;

	// Solution for zero temperature with
	// beta equilibrium and charge neutrality
	// (ignores proton_fraction as it is a
	// variable in this case)
    struct _multiroot{
        int max_iterations;
        double proton_fraction_mapping_scale;
        bool use_last_solution_as_guess;

        double abs_error;
        double rel_error;

        struct _guesses{
            double mass;
            double proton_fraction;
        } guesses;

        // Zero mass special case:
        //  To masses less than mass_tolerance, the value of mass
        //  will be assumed to be zero and the proton fraction will
        //  be searched between the given bounds
        struct _special_case{
            double mass_tolerance;
            double lower_bound;
            double upper_bound;
        } special_case;
    } multiroot;
    
} Parameters;

extern Parameters parameters;

void ParametersSetup(void);
void SetParametersSet(char * parameters_set_identifier);
void SetParametersFromCommandline();

void PrintParametersToFile(FILE * file);

#endif /* Parameters_h */
