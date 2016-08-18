//
//  Tests.c
//  hadrons EOS
//
//  Created by Clebson Graeff on 2016-08-16.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <stdbool.h>
#include <strings.h>
#include <stdio.h>
#include <math.h>

#include "CommandlineOptions.h"
#include "Tests.h"
#include "AuxiliaryFunctions.h"
#include "ZeroTemperatureEOS.h"
#include "Constants.h"
#include "Parameters.h"

void RunTests()
{
    printf("Running tests ...\n");
           
#pragma mark Zeroed Gap Equation
    
    // Write the zeroed gap equation for a set of proton and neutron
    // densities, for the eNJL1 parameters set. This serves for
    // verification purposes, as we should get a characteristic
    // curve. This tests functions used in the zero temperature case.
	if (true)
	{
        printf("\tZeroed Gap Equation\n");
        SetFilePath("tests/zeroed-gap-equation/data/");

        SetParametersSet("eNJL1");
        
        double neutron_density[9] = {0.0125, 0.025, 0.0375, 0.05, 0.0625, 0.075, 0.0875, 0.1, 0.1125};
        double proton_density[9] = {0.0125, 0.025, 0.0375, 0.05, 0.0625, 0.075, 0.0875, 0.1, 0.1125};
        
        double neutron_fermi_momentum[9];
        double proton_fermi_momentum[9];
        
        for (int i = 0; i < 9; i++){
            neutron_fermi_momentum[i] = pow(3.0 * pow(CONST_HBAR_C, 3.0) * pow(M_PI, 2.0) * neutron_density[i], 1.0 / 3.0);
            proton_fermi_momentum[i] = pow(3.0 * pow(CONST_HBAR_C, 3.0) * pow(M_PI, 2.0) * neutron_density[i], 1.0 / 3.0);
        }
		
        for (int i = 0; i < 9; i++){
            for (int j = 0; j < 9; j++){
            
                double m = 0;
                
                char filename[256];
                sprintf(filename, "gap_dens_%d_%d.dat", i, j);
                
                FILE * f = OpenFile(filename);
                
                gap_equation_input input;
                input.proton_density = proton_density[i];
                input.neutron_density = neutron_density[j];
                input.proton_fermi_momentum = proton_fermi_momentum[i];
                input.neutron_fermi_momentum = neutron_fermi_momentum[j];
                
                while (m < 1000.0) {
                    fprintf(f, "%20.15E\t%20.15E\n", m, GapEquation(m, &input));
                    m += 0.5;
                }
                
                fclose(f);
            }
        }
        
        SetFilePath("tests/zeroed-gap-equation/");
      	FILE * log_file = OpenFile("run.log");

        fprintf(log_file, "Write the zeroed gap equation for a set of proton and neutron\n"
                          "densities, for the eNJL1 parameters set. This serves for\n"
                          "verification purposes, as we should get a characteristic\n"
                          "curve. This tests functions used in the zero temperature case.\n\n");
        fprintf(log_file, "The following values of densities (and corresponding Fermi momenta)\n"
                          "were used:\n");
        fprintf(log_file, "neutron_density = {0.0125, 0.025, 0.0375, 0.05, 0.0625, 0.075, 0.0875, 0.1, 0.1125};\n"
                          "proton_density = {0.0125, 0.025, 0.0375, 0.05, 0.0625, 0.075, 0.0875, 0.1, 0.1125};\n");
        
        PrintParametersToFile(log_file);
        
		fclose(log_file);
	}
           
    printf("done.\n");
}
