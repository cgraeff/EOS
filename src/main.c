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
int SolveZeroTemperatureStarEOS();

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
    
    // Set parameters overrides from commandline
    SetParametersFromCommandline();

    if (parameters.temperature == 0){
        if (options.stars){
            SolveZeroTemperatureStarEOS();
            return 0;
        }
        
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
	gsl_vector * thermodynamic_potential_vector = gsl_vector_alloc(parameters.points_number);
	
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

		double total_scalar_density = ScalarDensity(mass, proton_fermi_momentum, parameters.theory.cutoff)
									  + ScalarDensity(mass, neutron_fermi_momentum, parameters.theory.cutoff);


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

		
		double thermodynamic_potential = ThermodynamicPotential(total_scalar_density,
                                                                barionic_density,
                                                                proton_density,
                                                                neutron_density,
                                                                proton_chemical_potential,
                                                                neutron_chemical_potential,
                                                                kinectic_energy_density);
		
		gsl_vector_set(thermodynamic_potential_vector, i, thermodynamic_potential);
		
		double pressure = Pressure(thermodynamic_potential);
		
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
        gsl_vector_set(energy_by_nucleon_vector, i, gsl_vector_get(energy_by_nucleon_vector, i) - parameters.theory.nucleon_mass);
    
    
    /*
     * Save results
     */
	
    if (options.verbose)
        printf("Saving results ...\n");
    
    if (options.dirs)
        SetFilePath("output/IR/data/");

    WriteVectorsToFile("mass.dat",
					   "# barionic density (fm^{-3}), mass (MeV)\n",
					   2,
					   barionic_density_vector,
					   mass_vector);

	WriteVectorsToFile("scalar_density.dat",
					   "# barionic density (fm^{-3}), scalar density (fm^{-3})\n",
					   2,
					   barionic_density_vector,
					   scalar_density_vector);

	WriteVectorsToFile("thermodynamic_potential.dat",
					   "# density (fm^{-3}), grand canonical potential per unit volume (MeV/fm^{3})\n",
					   2,
					   barionic_density_vector, thermodynamic_potential_vector);

	WriteVectorsToFile("proton_chemical_potential.dat",
					   "# barionic density (fm^{-3}), proton chemical potential (MeV)\n",
					   2,
					   barionic_density_vector,
					   proton_chemical_potential_vector);

    WriteVectorsToFile("neutron_chemical_potential.dat",
                       "# barionic density (fm^{-3}), neutron chemical potential (MeV)\n",
                       2,
                       barionic_density_vector,
                       neutron_chemical_potential_vector);
    
    WriteVectorsToFile("kinectic_energy_density.dat",
                      "# barionic density (fm^{-3}), kinectic energy density (MeV)\n",
                      2,
                      barionic_density_vector,
                      kinectic_energy_density_vector);

    if (options.dirs)
        SetFilePath("output/Rep_Tsue/data/");
    
    WriteVectorsToFile("pressure_chem_pot.dat",
                       "# chemical potential (MeV), pressure (MeV/fm^{3})\n",
                       2,
                       proton_chemical_potential_vector,
                       pressure_vector);
    
	if (options.dirs)
        SetFilePath("output/EOS/data/");
        	   
	WriteVectorsToFile("pressure.dat",
					   "# barionic density (fm^{-3}), pressure (MeV/fm^{3})\n",
					   2,
					   barionic_density_vector,
					   pressure_vector);

	WriteVectorsToFile("energy_density.dat",
					   "# barionic density (fm^{-3}), energy per unit volume (MeV/fm^{3})\n",
					   2,
					   barionic_density_vector,
					   energy_density_vector);

	WriteVectorsToFile("energy_density_per_particle.dat",
					   "# energy density / barionic density = energy by nucleon \n"
					   "# barionic density (fm^{-3}), energy / barionic density  (MeV)\n",
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
    gsl_vector_free(thermodynamic_potential_vector);

    gsl_vector_free(pressure_vector);
    gsl_vector_free(energy_density_vector);
    gsl_vector_free(energy_by_nucleon_vector);
    
    if (options.verbose)
        printf("Done!\n");

    return 0;
}


int SolveZeroTemperatureStarEOS(){

    // Print name of parametrization
    if (options.verbose)
        printf("Calculation performed with %s parameters set.\n", parameters.parameters_set_identifier);

    gsl_vector * barionic_density_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * scalar_density_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * mass_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * proton_fraction_vector = gsl_vector_alloc(parameters.points_number);

	gsl_vector * proton_chemical_potential_vector = gsl_vector_alloc(parameters.points_number);
	gsl_vector * neutron_chemical_potential_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * kinectic_energy_density_vector = gsl_vector_alloc(parameters.points_number);
	gsl_vector * thermodynamic_potential_vector = gsl_vector_alloc(parameters.points_number);

	gsl_vector * pressure_vector = gsl_vector_alloc(parameters.points_number);
	gsl_vector * energy_density_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * electron_energy_density_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * electron_pressure_vector = gsl_vector_alloc(parameters.points_number);

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

        double mass;
        double proton_fraction;
        SolveMultiRoots(barionic_density, &mass, &proton_fraction);
        
        double proton_density = barionic_density * proton_fraction;
        double neutron_density = barionic_density * (1.0 - proton_fraction);
        
        // Ensure charge neutrality
        double electron_density = proton_density;
        
        double proton_fermi_momentum = FermiMomentum(proton_density);
        double neutron_fermi_momentum = FermiMomentum(barionic_density * (1.0 - proton_fraction));
        double electron_fermi_momentum = FermiMomentum(electron_density);

		double total_scalar_density = ScalarDensity(mass, proton_fermi_momentum, parameters.theory.cutoff)
									  + ScalarDensity(mass, neutron_fermi_momentum, parameters.theory.cutoff);


        gsl_vector_set(scalar_density_vector, i , total_scalar_density);
        gsl_vector_set(mass_vector, i, mass);
        gsl_vector_set(proton_fraction_vector, i, proton_fraction);

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


        // BEWARE: The thermodynamic potential don't have the
        // electron contribution
		double thermodynamic_potential = ThermodynamicPotential(total_scalar_density,
                                                                barionic_density,
                                                                proton_density,
                                                                neutron_density,
                                                                proton_chemical_potential,
                                                                neutron_chemical_potential,
                                                                kinectic_energy_density);

		gsl_vector_set(thermodynamic_potential_vector, i, thermodynamic_potential);

        // Electron contributions to energy and pressure
        double electron_pressure = ElectronPressure(electron_fermi_momentum);
        double electron_energy_density = ElectronEnergyDensity(electron_fermi_momentum);
        
        gsl_vector_set(electron_pressure_vector, i, electron_pressure);
        gsl_vector_set(electron_energy_density_vector, i, electron_energy_density);
        
		double pressure = electron_pressure + Pressure(thermodynamic_potential);

		gsl_vector_set(pressure_vector, i, pressure);

		double energy_density = electron_energy_density
                                + EnergyDensity(pressure,
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
        gsl_vector_set(energy_by_nucleon_vector, i, gsl_vector_get(energy_by_nucleon_vector, i) - parameters.theory.nucleon_mass);


    /*
     * Save results
     */

    if (options.verbose)
        printf("Saving results ...\n");

    if (options.dirs)
        SetFilePath("output/IR/data/");

    WriteVectorsToFile("mass.dat",
					   "# barionic density (fm^{-3}), mass (MeV)\n",
					   2,
					   barionic_density_vector,
					   mass_vector);

	WriteVectorsToFile("scalar_density.dat",
					   "# barionic density (fm^{-3}), scalar density (fm^{-3})\n",
					   2,
					   barionic_density_vector,
					   scalar_density_vector);

	WriteVectorsToFile("thermodynamic_potential.dat",
					   "# barionic density (fm^{-3}), grand canonical potential per unit volume (MeV/fm^{3})\n",
					   2,
					   barionic_density_vector,
                       thermodynamic_potential_vector);

	WriteVectorsToFile("proton_chemical_potential.dat",
					   "# barionic density (fm^{-3}), proton chemical potential (MeV)\n",
					   2,
					   barionic_density_vector,
					   proton_chemical_potential_vector);

    WriteVectorsToFile("neutron_chemical_potential.dat",
                       "# barionic density (fm^{-3}), neutron chemical potential (MeV)\n",
                       2,
                       barionic_density_vector,
                       neutron_chemical_potential_vector);

    WriteVectorsToFile("kinectic_energy_density.dat",
                      "# barionic density (fm^{-3}), kinectic energy density (MeV/fm^{3})\n",
                      2,
                      barionic_density_vector,
                      kinectic_energy_density_vector);
    
    WriteVectorsToFile("proton_fraction.dat",
                       "# barionic density (fm^{-3}), proton fraction\n",
                       2,
                       barionic_density_vector,
                       proton_fraction_vector);
    
    WriteVectorsToFile("electron_energy_density.dat",
                       "# barionic density (fm^{-3}), electron energy density contribution (MeV/fm^{3})\n",
                       2,
                       barionic_density_vector,
                       electron_energy_density_vector);
    
    WriteVectorsToFile("electron_pressure.dat",
                       "# barionic density (fm^{-3}), electron pressure contribution (MeV/fm^{3})\n",
                       2,
                       barionic_density_vector,
                       electron_pressure_vector);

    if (options.dirs)
        SetFilePath("output/Rep_Tsue/data/");

    WriteVectorsToFile("pressure_chem_pot.dat",
                       "# chemical potential (MeV), pressure (MeV/fm^{3})\n",
                       2,
                       proton_chemical_potential_vector,
                       pressure_vector);

	if (options.dirs)
        SetFilePath("output/EOS/data/");

	WriteVectorsToFile("pressure.dat",
					   "# barionic density (fm^{-3}), pressure (MeV/fm^{3})\n",
					   2,
					   barionic_density_vector,
					   pressure_vector);

	WriteVectorsToFile("energy_density.dat",
					   "# barionic density (fm^{-3}), energy per unit volume (MeV/fm^{3})\n",
					   2,
					   barionic_density_vector,
					   energy_density_vector);

	WriteVectorsToFile("energy_density_per_particle.dat",
					   "# energy density / barionic density = energy by nucleon \n"
					   "# barionic density (fm^{-3}), energy / barionic density (MeV)\n",
					   2,
					   barionic_density_vector, energy_by_nucleon_vector);
    
    WriteVectorsToFile("rho_energy_pressure.dat",
                       "# barionic density (fm^{-3}), energy density (MeV/fm^{3}), pressure (MeV/fm^{3})\n",
                       3,
                       barionic_density_vector,
                       energy_density_vector,
                       pressure_vector);

    /*
     * Clean up
     */

    // Free vectors
    gsl_vector_free(barionic_density_vector);

    gsl_vector_free(scalar_density_vector);
    gsl_vector_free(mass_vector);

    gsl_vector_free(proton_chemical_potential_vector);
    gsl_vector_free(neutron_chemical_potential_vector);
    gsl_vector_free(thermodynamic_potential_vector);

    gsl_vector_free(pressure_vector);
    gsl_vector_free(energy_density_vector);
    gsl_vector_free(energy_by_nucleon_vector);

    if (options.verbose)
        printf("Done!\n");

    return 0;
}
