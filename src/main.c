//
//  main.c
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-02-16.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_vector.h>

#include "CommandlineOptions.h"
#include "Tests.h"
#include "Parameters.h"
#include "Constants.h"
#include "ZeroTemperatureEOS.h"
#include "AuxiliaryFunctions.h"

int SolveZeroTemperatureEOS();
int SolveFiniteTemperatureEOS();

int main(int argc, char * argv[])
{
    CommandlineOptionsParse(argc, argv);
    ParametersSetup();

    if (options.tests){
  		RunTests();
        exit(EXIT_SUCCESS);
    }

  	// If option -p is used, set parameters set accordingly,
	// otherwise, use default set
  	SetParametersSet(options.parameterization);
    
    // If the temperature was chosen using
    // commandline options, use it
    // (-1.0 is a place holder value)
    if (options.temp != -1.0)
        parameters.temperature = options.temp;
    
    if (parameters.temperature == 0){
        SolveZeroTemperatureEOS();
    }
    else if (parameters.temperature > 0){
        //SolveFiniteTemperatureEOS();
        printf("Finite temperature is not yet implemented.\n");
        abort();
    }
    else{
        printf("Values of temperature must be non-negative.\n");
        printf("(%f was provided).\n", parameters.temperature);
        exit(EXIT_FAILURE);
    }

    return 0;
}

/* 
 * The EOS for zero temperature are solved in the function bellow. 
 * In this particular implementation, the barionic density was chosen
 * as the running parameter.
 */
int SolveZeroTemperatureEOS(){
    
    // Print name of parametrization
    if (options.verbose)
        printf("Calculation performed with %s parameters set.\n", parameters.parameters_set_identifier);
    
    gsl_vector * barionic_density_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * scalar_density_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * mass_vector = gsl_vector_alloc(parameters.points_number);
	
	gsl_vector * proton_chemical_potential_vector = gsl_vector_alloc(parameters.points_number);
	gsl_vector * neutron_chemical_potential_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * kinectic_energy_density_vector = gsl_vector_alloc(parameters.points_number);
	gsl_vector * termodynamic_potential_vector = gsl_vector_alloc(parameters.points_number);
	
	gsl_vector * pressure_vector = gsl_vector_alloc(parameters.points_number);
	gsl_vector * energy_density_vector = gsl_vector_alloc(parameters.points_number);
	
    /*
     * Main loop (on barionic density)
     */
    
    if (options.verbose)
        printf("Solving gap equation and equations of state ...\n");
    
    double barionic_density = parameters.minimum_density;
    
    // Define the density step. We subtract 1 from the number of points to
    // make sure that the last point corresponds to maximum_density
    double density_step = (parameters.maximum_density - parameters.minimum_density)
                          / (double)(parameters.points_number - 1);
    
    for (int i = 0; i < parameters.points_number; i++){
        
        gsl_vector_set(barionic_density_vector, i, barionic_density);
        
        if (options.verbose){
            printf("\r\tBarionic density: %f", barionic_density);
            fflush(stdout);
        }
        
        double proton_density = parameters.proton_fraction * barionic_density;
        double neutron_density = (1.0 - parameters.proton_fraction) * barionic_density;
        
        gsl_vector_set(barionic_density_vector, i, barionic_density);
		
		double proton_fermi_momentum = pow(3.0 * pow(CONST_HBAR_C, 3.0) * pow(M_PI, 2.0) * proton_density,
                                           1.0 / 3.0);
		double neutron_fermi_momentum = pow(3.0 * pow(CONST_HBAR_C, 3.0) * pow(M_PI, 2.0) * neutron_density,
                                            1.0 / 3.0);
		
		double mass = SolveGapEquation(proton_density,
                                       neutron_density,
                                       proton_fermi_momentum,
                                       neutron_fermi_momentum);

		double total_scalar_density = ScalarDensity(mass, proton_fermi_momentum, parameters.CUTOFF)
									  + ScalarDensity(mass, neutron_fermi_momentum, parameters.CUTOFF);

		
        gsl_vector_set(scalar_density_vector, i , total_scalar_density);
        gsl_vector_set(mass_vector, i, mass);
        
        double proton_chemical_potential =	ProtonChemicalPotential(proton_fermi_momentum,
																	total_scalar_density,
																	mass,
																	barionic_density,
																	proton_density,
																	neutron_density);
		
		gsl_vector_set(proton_chemical_potential_vector, i, proton_chemical_potential);
        
        double neutron_chemical_potential =	NeutronChemicalPotential(neutron_fermi_momentum,
																	 total_scalar_density,
																	 mass,
																	 barionic_density,
																	 proton_density,
																	 neutron_density);
		
		gsl_vector_set(neutron_chemical_potential_vector, i, neutron_chemical_potential);
        
		double kinectic_energy_density = KinecticEnergyDensity(mass, proton_fermi_momentum, neutron_fermi_momentum);
        gsl_vector_set(kinectic_energy_density_vector, i, kinectic_energy_density);

		
		double termodynamic_potential = TermodynamicPotential(total_scalar_density,
															  barionic_density,
															  proton_density,
															  neutron_density,
															  proton_chemical_potential,
															  neutron_chemical_potential,
															  kinectic_energy_density);
		
		gsl_vector_set(termodynamic_potential_vector, i, termodynamic_potential);
		
		double pressure = Pressure(termodynamic_potential);
		
		gsl_vector_set(pressure_vector, i, pressure);
		
		double energy_density = EnergyDensity(pressure,
											  proton_chemical_potential,
											  neutron_chemical_potential,
											  proton_density,
											  neutron_density);
		
		gsl_vector_set(energy_density_vector, i, energy_density);
        
        // Prepare next iteration
        barionic_density += density_step;
    }
    if (options.verbose)
        printf("\n"); // As print inside the loop doesn't use new line, we need one now
    
    // Calculate energy per particle
    gsl_vector * energy_by_nucleon_vector = VectorNewVectorFromDivisionElementByElement(energy_density_vector,
                                                                                        barionic_density_vector);
    
    for (int i = 0; i < energy_by_nucleon_vector->size; i++)
        gsl_vector_set(energy_by_nucleon_vector, i, gsl_vector_get(energy_by_nucleon_vector, i) - parameters.nucleon_mass);
    
    
    /*
     * Save results
     */
	
    if (options.verbose)
        printf("Saving results ...\n");
    
    if (options.dirs)
        SetFilePath("output/IR/data/");

    WriteVectorsToFile("mass.dat",
					   "# barionic density, mass\n",
					   2,
					   barionic_density_vector,
					   mass_vector);

	WriteVectorsToFile("scalar_density.dat",
					   "# barionic density, scalar density\n",
					   2,
					   barionic_density_vector,
					   scalar_density_vector);

	WriteVectorsToFile("thermodynamic_potential.dat",
					   "# density, grand canonical potential per unit volume\n",
					   2,
					   barionic_density_vector, termodynamic_potential_vector);

	WriteVectorsToFile("proton_chemical_potential.dat",
					   "# barionic dentisy, proton chemical potential\n",
					   2,
					   barionic_density_vector,
					   proton_chemical_potential_vector);

    WriteVectorsToFile("neutron_chemical_potential.dat",
                       "# barionic dentisy, neutron chemical potential\n",
                       2,
                       barionic_density_vector,
                       neutron_chemical_potential_vector);
    
    WriteVectorsToFile("kinectic_energy_density.dat",
                      "# barionic density, kinectic energy density\n",
                      2,
                      barionic_density_vector,
                      kinectic_energy_density_vector);

	if (options.dirs)
        SetFilePath("output/EOS/data/");
        	   
	WriteVectorsToFile("pressure.dat",
					   "# barionic density, pressure\n",
					   2,
					   barionic_density_vector,
					   pressure_vector);

	WriteVectorsToFile("energy_density.dat",
					   "# barionic density, energy per unit volume\n",
					   2,
					   barionic_density_vector,
					   energy_density_vector);

	WriteVectorsToFile("energy_density_per_particle.dat",
					   "# energy / barionic density = energy by nucleon \n"
					   "# barionic density, energy / barionic density\n",
					   2,
					   barionic_density_vector, energy_by_nucleon_vector);

    /*
     * Clean up
     */
    
    // Free vectors
    gsl_vector_free(barionic_density_vector);

    gsl_vector_free(scalar_density_vector);
    gsl_vector_free(mass_vector);
    
    gsl_vector_free(proton_chemical_potential_vector);
    gsl_vector_free(neutron_chemical_potential_vector);
    gsl_vector_free(termodynamic_potential_vector);

    gsl_vector_free(pressure_vector);
    gsl_vector_free(energy_density_vector);
    gsl_vector_free(energy_by_nucleon_vector);
    
    if (options.verbose)
        printf("Done!\n");

    return 0;
}



