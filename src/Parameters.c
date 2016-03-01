//
//  Parameters.c
//  EOS
//
//  Created by Clebson Graeff on 2/16/16.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <math.h>

#include "Parameters.h"

Parameters parameters;

void ParametersSetup(void)
{
	// Parameters from Helena Pais, Débora P. Menezes, Constança Providência

	// eNJL1
 	parameters.G_S = 4.855;
  	parameters.G_V = 4.65;
	parameters.G_SV = -6.583;
	parameters.G_RHO = 0.587601034;
	parameters.G_VRHO = 0.0;

	parameters.CUTOFF = 388.189;
	parameters.nucleon_mass = 939.0;

	parameters.minimum_density = 0.01; // fm^-3
	parameters.maximum_density = 0.35; // fm^-3

	parameters.points_number = 1000;

	parameters.proton_fraction = 0.5;
	
	parameters.number_of_colors = 1;
	parameters.number_of_flavors = 1;
	
	parameters.degeneracy = (double)parameters.number_of_colors * (double)parameters.number_of_flavors / pow(M_PI, 2.0);
	
	parameters.lower_bound_gap_eq_solve = 1.0E-3; // Low, but not zero. In zero f(M) = 0;
	parameters.upper_bound_gap_eq_solve = 3000.0;	// MeV (about 150% of nucleon mass)
}
