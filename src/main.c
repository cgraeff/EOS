//
//  main.c
//  EOS
//
//  Created by Clebson Graeff on 2/16/16.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>

#include "Parameters.h"
#include "Constants.h"
#include "EOS.h"
#include "Vector.h"

int main(int argc, char * argv[])
{
    ParametersSetup();
    
    double barionic_density = parameters.minimum_density;
	
	// Define the density step. We subtract 1 from the number of points to
	// make sure that the last point corresponds to maximum_density
    double density_step = (parameters.maximum_density - parameters.minimum_density)
    / (parameters.points_number - 1);
    
    gsl_vector * barionic_density_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * proton_density_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * neutron_density_vector = gsl_vector_alloc(parameters.points_number);
	
	gsl_vector * proton_fermi_momentum_vector = gsl_vector_alloc(parameters.points_number);
	gsl_vector * neutron_fermi_momentum_vector = gsl_vector_alloc(parameters.points_number);
	
    gsl_vector * scalar_density_vector = gsl_vector_alloc(parameters.points_number);
    gsl_vector * mass_vector = gsl_vector_alloc(parameters.points_number);
	
	gsl_vector * proton_chemical_potential_vector = gsl_vector_alloc(parameters.points_number);
	gsl_vector * neutron_chemical_potential_vector = gsl_vector_alloc(parameters.points_number);
	gsl_vector * kinectic_energy_density_vector = gsl_vector_alloc(parameters.points_number);
	gsl_vector * termodynamic_potential_vector = gsl_vector_alloc(parameters.points_number);
	
	gsl_vector * pressure_vector = gsl_vector_alloc(parameters.points_number);
	gsl_vector * energy_density_vector = gsl_vector_alloc(parameters.points_number);
	
    for (int i = 0; i < parameters.points_number; i++){
        
        barionic_density += density_step;
        
        double proton_density = parameters.proton_fraction * barionic_density;
        double neutron_density = (1.0 - parameters.proton_fraction) * barionic_density;
        
        gsl_vector_set(barionic_density_vector, i, barionic_density);
        gsl_vector_set(proton_density_vector, i , proton_density);
        gsl_vector_set(neutron_density_vector, i, neutron_density);
		
		double proton_fermi_momentum = pow(3.0 * pow(CONST_HBAR_C, 3.0) * pow(M_PI, 2.0) * proton_density, 1.0 / 3.0);
		double neutron_fermi_momentum = pow(3.0 * pow(CONST_HBAR_C, 3.0) * pow(M_PI, 2.0) * neutron_density, 1.0 / 3.0);
		
		gsl_vector_set(proton_fermi_momentum_vector, i, proton_fermi_momentum);
		gsl_vector_set(neutron_fermi_momentum_vector, i, neutron_fermi_momentum);
		
		
		double m = 0;
		
		char filename[256];
		sprintf(filename, "Gap/Gap_dens_%f.dat", barionic_density);
		FILE * f = fopen(filename, "w");
		gap_equation_input input;
		input.barionic_density = barionic_density;
		input.neutron_fermi_momentum = neutron_fermi_momentum;
		input.proton_fermi_momentum = proton_fermi_momentum;
		
		while (m < 2000) {
			fprintf(f, "%20.15E\t%20.15E\n", m, GapEquation(m, &input));
			m += 0.5;
		}
		
		double mass = SolveGapEquation(barionic_density, proton_fermi_momentum, neutron_fermi_momentum);

		double total_scalar_density = scalar_density(mass, proton_fermi_momentum, parameters.CUTOFF)
									  + scalar_density(mass, neutron_fermi_momentum, parameters.CUTOFF);

		
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
 
    }
	
    // Write pressure, energy, chemical_potential, density, ...
	
	WriteVectorsToFile("densities.dat",
					   "# barionic, proton, and neutron densities\n",
					   3,
					   barionic_density_vector,
					   proton_density_vector,
					   neutron_density_vector);
	
	WriteVectorsToFile("fermi_momentum.dat",
					   "# barionic density, proton fermi momentum, neutron fermi momentum\n",
					   3,
					   barionic_density_vector,
					   proton_fermi_momentum_vector,
					   neutron_fermi_momentum_vector);
	
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
	
	WriteVectorsToFile("kinectic_energy.dat",
					   "# density, kinectic energy per unit volume\n",
					   2,
					   barionic_density_vector,
					   kinectic_energy_density_vector);
	
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
	
	WriteVectorsToFile("chemical_potentials.dat",
					   "# barionic dentisy, proton chemical potential, neutron chemical potential",
					   3,
					   barionic_density_vector,
					   proton_chemical_potential_vector,
					   neutron_chemical_potential_vector);
	
	WriteVectorsToFile("proton_chemical_potential-mass.dat",
					   "# proton chemical potential, mass\n",
					   2,
					   proton_chemical_potential_vector,
					   mass_vector);
	
	gsl_vector * energy_by_nucleon_vector = VectorNewVectorFromDivisionElementByElement(energy_density_vector,
																						barionic_density_vector);
	
	WriteVectorsToFile("energy_by_nucleon.dat",
					   "# energy / barionic density = energy by nucleon \n"
					   "# barionic density, energy / barionic density\n",
					   2,
					   barionic_density_vector, energy_by_nucleon_vector);
    return 0;
}



