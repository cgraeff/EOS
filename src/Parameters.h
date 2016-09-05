//
//  Parameters.h
//  hadrons EOS
//
//  Created by Clebson Graeff on 2016-02-16.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#ifndef Parameters_h
#define Parameters_h

//typedef struct _parameters_
//{
//    double mass;
//    double proton_fraction;
//    
//} RootGuesses;

typedef struct _parameters
{
    char * parameters_set_identifier;
    char * parameters_set_origin;
    
	double G_S;		// scalar-isoscalar coupling (fm^2)
	double G_V;		// vector-isoscalar coupling (fm^2)
	double G_RHO; 	// vector-isovector vector_coupling (fm^2)
	double G_SV;	// (fm^8)
	double G_VRHO;	// (fm^8)
    double G_SRHO;  // (fm^8)

	double CUTOFF;	// \Lambda (MeV)
	double nucleon_mass; // (MeV)
    double electron_mass; // (MeV)
    double bare_mass;

	double minimum_density; // (fm^-3)
	double maximum_density; // (fm^-3)

	int points_number; // Number of points in which the above range will be divide into

	double proton_fraction; // (no dimension)
	
	double lower_bound_gap_eq_solve;
	double upper_bound_gap_eq_solve;
    
    double temperature;
	
    struct _multiroot{
        int max_iterations;
        
        struct _guesses{
            double mass;
            double proton_fraction;
        } guesses;
        
        double abs_error;
        double rel_error;
    } multiroot;
    
} Parameters;

extern Parameters parameters;

void ParametersSetup(void);
void SetParametersSet(char * parameters_set_identifier);

void PrintParametersToFile(FILE * file);

#endif /* Parameters_h */
