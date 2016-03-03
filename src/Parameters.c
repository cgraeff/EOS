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
    
    typedef enum{
        eNJL1,
        eNJL1OmegaRho1,
        eNJL1OmegaRho2,
        eNJL2,
        eNJL2OmegaRho1,
        eNJL3,
        eNJL3SigmaRho1,
        eNJL1m,
        eNJL1mSigmaRho1,
        eNJL2m,
        eNJL2mSigmaRho1
    }Parameters_set;
    
    Parameters_set par = eNJL1;
    
    // Parameters from Helena Pais, Débora P. Menezes, Constança Providência
    switch (par) {
        case eNJL1:
            parameters.G_S = 4.855;
            parameters.G_V = 4.65;
            parameters.G_SV = -6.583;
            parameters.G_RHO = 0.5876;
            parameters.G_VRHO = 0.0;
            parameters.G_SRHO = 0.0;
            parameters.CUTOFF = 388.189;
            parameters.bare_mass = 0.0;
            break;

        case eNJL1OmegaRho1:
            parameters.G_S = 4.855;
            parameters.G_V = 4.65;
            parameters.G_SV = -6.583;
            parameters.G_RHO = 0.5976;
            parameters.G_VRHO = -1.0;
            parameters.G_SRHO = 0.0;
            parameters.CUTOFF = 388.189;
            parameters.bare_mass = 0.0;
            break;
            
        case eNJL1OmegaRho2:
            parameters.G_S = 4.855;
            parameters.G_V = 4.65;
            parameters.G_SV = -6.583;
            parameters.G_RHO = 0.6476;
            parameters.G_VRHO = -6.0;
            parameters.G_SRHO = 0.0;
            parameters.CUTOFF = 388.189;
            parameters.bare_mass = 0.0;
            break;

        case eNJL2:
            parameters.G_S = 3.8;
            parameters.G_V = 3.8;
            parameters.G_SV = -4.228;
            parameters.G_RHO = 0.6313;
            parameters.G_VRHO = 0.0;
            parameters.G_SRHO = 0.0;
            parameters.CUTOFF = 422.384;
            parameters.bare_mass = 0.0;
            break;
            
        case eNJL2OmegaRho1:
            parameters.G_S = 3.8;
            parameters.G_V = 3.8;
            parameters.G_SV = -4.288;
            parameters.G_RHO = 0.6413;
            parameters.G_VRHO = -1.0;
            parameters.G_SRHO = 0.0;
            parameters.CUTOFF = 422.384;
            parameters.bare_mass = 0.0;
            break;
            
        case eNJL3:
            parameters.G_S = 1.93;
            parameters.G_V = 3.0;
            parameters.G_SV = -1.8;
            parameters.G_RHO = 0.65;
            parameters.G_VRHO = 0.0;
            parameters.G_SRHO = 0.0;
            parameters.CUTOFF = 534.815;
            parameters.bare_mass = 0.0;
            break;
            
        case eNJL3SigmaRho1:
            parameters.G_S = 1.93;
            parameters.G_V = 3;
            parameters.G_SV = -1.8;
            parameters.G_RHO = 0.0269;
            parameters.G_VRHO = 0.0;
            parameters.G_SRHO = 0.5;
            parameters.CUTOFF = 534.815;
            parameters.bare_mass = 0.0;
            break;
            
        case eNJL1m:
            parameters.G_S = 1.3833;
            parameters.G_V = 1.781;
            parameters.G_SV = -2.943;
            parameters.G_RHO = 0.7;
            parameters.G_VRHO = 0.0;
            parameters.G_SRHO = 0.0;
            parameters.CUTOFF = 478.248;
            parameters.bare_mass = 450.0;
            break;
            
        case eNJL1mSigmaRho1:
            parameters.G_S = 1.3833;
            parameters.G_V = 1.781;
            parameters.G_SV = -2.943;
            parameters.G_RHO = 0.0739;
            parameters.G_VRHO = 0.0;
            parameters.G_SRHO = 1.0;
            parameters.CUTOFF = 478.248;
            parameters.bare_mass = 450.0;
            break;
            
        case eNJL2m:
            parameters.G_S = 1.078;
            parameters.G_V = 1.955;
            parameters.G_SV = -2.74;
            parameters.G_RHO = 0.75;
            parameters.G_VRHO = 0.0;
            parameters.G_SRHO = 0.0;
            parameters.CUTOFF = 502.466;
            parameters.bare_mass = 500.0;
            break;
            
        case eNJL2mSigmaRho1:
            parameters.G_S = 1.078;
            parameters.G_V = 1.955;
            parameters.G_SV = -2.74;
            parameters.G_RHO = 0.1114;
            parameters.G_VRHO = 0.0;
            parameters.G_SRHO = 1.0;
            parameters.CUTOFF = 502.466;
            parameters.bare_mass = 500.0;
            break;
        default:
            break;
    }
}
