//
//  EOS.c
//  EOS
//
//  Created by Clebson Graeff on 2016-02-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include "Parameters.h"
#include "Constants.h"
#include "ZeroTemperatureEOS.h"

double SolveGapEquation(double proton_density,
                        double neutron_density,
                        double proton_fermi_momentum,
                        double neutron_fermi_momentum)
{
	const gsl_root_fsolver_type * T	= gsl_root_fsolver_bisection; // Maybe this would be better: gsl_root_fsolver_brent

	gsl_root_fsolver * s = gsl_root_fsolver_alloc(T);

	gap_equation_input input;
	
	input.proton_density = proton_density;
    input.neutron_density = neutron_density;
	input.proton_fermi_momentum = proton_fermi_momentum;
	input.neutron_fermi_momentum = neutron_fermi_momentum;
	
	gsl_function F;
	F.function = &GapEquation;
	F.params = &input;
    
    // Test if the limits straddle the root
    // If they don't, we will assume that
    // the mass = 0 root was found
    if (GSL_SIGN(GSL_FN_EVAL(&F, parameters.lower_bound_gap_eq_solve))
        == GSL_SIGN(GSL_FN_EVAL(&F, parameters.upper_bound_gap_eq_solve)))
        return 0;
	
	gsl_root_fsolver_set(s, &F, parameters.lower_bound_gap_eq_solve, parameters.upper_bound_gap_eq_solve);
	
	int i = 0;
	double x_lower;
	double x_upper;
	do{
		i++;
		
		int status = gsl_root_fsolver_iterate(s);
		
		if (status != GSL_SUCCESS){
			printf("ERROR: No solution to the gap equation was found!\n");
			exit(EXIT_FAILURE);
		}
		
		x_lower = gsl_root_fsolver_x_lower(s);
		x_upper = gsl_root_fsolver_x_upper(s);
	} while(GSL_CONTINUE == gsl_root_test_interval(x_lower, x_upper, CONST_ABS_ERROR_GAP_EQ_SOLVING, CONST_REL_ERROR_GAP_EQ_SOLVING)
			&& i <= CONST_MAX_ITERATIONS_GAP_EQ_SOLVING);

	double result = gsl_root_fsolver_root(s);

	void gsl_root_fsolver_free(gsl_root_fsolver * S);
		
	return result;
}

double GapEquation(double mass, void * input)
{
	gap_equation_input * param = (gap_equation_input *)input;
    
    double barionic_density = param->proton_density + param->neutron_density;
    double rho_3 = param->proton_density - param->neutron_density;
	
	double scalar_density = scalar_density_function(mass, param->proton_fermi_momentum, parameters.CUTOFF)
                            + scalar_density_function(mass, param->neutron_fermi_momentum, parameters.CUTOFF);
	
	double gap_1st_term = 2.0  * CONST_HBAR_C * parameters.G_S * scalar_density;
	double gap_2nd_term = - 2.0 * CONST_HBAR_C * parameters.G_SV * scalar_density * pow(barionic_density, 2.0);
    double gap_3rd_term = - 2.0  * CONST_HBAR_C * parameters.G_SRHO * scalar_density * pow(rho_3, 2.0);

	return mass + gap_1st_term + gap_2nd_term + gap_3rd_term - parameters.bare_mass;
}

double scalar_density_function(double mass, double fermi_momentum, double cutoff)
{
	return pow(CONST_HBAR_C, -3.0) * (mass / pow(M_PI, 2.0)) * (F0(mass, fermi_momentum) - F0(mass, cutoff));
}

double vacuum_scalar_density()
{
	return 2.0 * pow(CONST_HBAR_C, -3.0) * (parameters.nucleon_mass / pow(M_PI, 2.0))
			* (F0(parameters.nucleon_mass, 0.0) - F0(parameters.nucleon_mass, parameters.CUTOFF));
}

double ProtonChemicalPotential(double proton_fermi_momentum,
							   double scalar_density,
							   double mass,
							   double barionic_density,
							   double proton_density,
							   double neutron_density)
{
	double chemical_potential = sqrt(pow(mass, 2.0) + pow(proton_fermi_momentum, 2.0));
	
	double rho_3 = (proton_density - neutron_density);
	
	chemical_potential += 2.0 * CONST_HBAR_C * (parameters.G_V * barionic_density
												+ parameters.G_SV * barionic_density * pow(scalar_density, 2.0)
												+ parameters.G_RHO * rho_3
												+ parameters.G_VRHO * pow(rho_3, 2.0) * barionic_density
												+ parameters.G_VRHO * pow(barionic_density, 2.0) * rho_3
                                                + parameters.G_SRHO * pow(scalar_density, 2.0) * rho_3);
	
	return chemical_potential;

}

double NeutronChemicalPotential(double neutron_fermi_momentum,
								double scalar_density,
								double mass,
								double barionic_density,
								double proton_density,
								double neutron_density)
{
	double chemical_potential = sqrt(pow(mass, 2.0) + pow(neutron_fermi_momentum, 2.0));
	
	double rho_3 = (proton_density - neutron_density);
	
	chemical_potential += 2.0 * CONST_HBAR_C * (parameters.G_V * barionic_density
                                                + parameters.G_SV * barionic_density * pow(scalar_density, 2.0)
                                                - parameters.G_RHO * rho_3
                                                + parameters.G_VRHO * pow(rho_3, 2.0) * barionic_density
                                                - parameters.G_VRHO * pow(barionic_density, 2.0) * rho_3
                                                - parameters.G_SRHO * pow(scalar_density, 2.0) * rho_3);
	
	return chemical_potential;
}


double KinecticEnergyDensity(double mass, double proton_fermi_momentum, double neutron_fermi_momentum)
{
	double proton_kinectic_energy = parameters.degeneracy * (F2(mass, proton_fermi_momentum) - F2(mass, parameters.CUTOFF));
	double neutron_kinectic_energy = parameters.degeneracy * (F2(mass, neutron_fermi_momentum) - F2(mass, parameters.CUTOFF));
	
	return (proton_kinectic_energy + neutron_kinectic_energy) / (pow(CONST_HBAR_C, 3.0));
}

double VacuumKinecticEnergyDensity()
{
	return 2.0 * parameters.degeneracy * (pow(CONST_HBAR_C, -3.0))
		   * (F2(parameters.nucleon_mass, 0) - F2(parameters.nucleon_mass, parameters.CUTOFF));
}

double VacuumEnergyDensity()
{
    double scalar_density_0 = vacuum_scalar_density();

	return VacuumKinecticEnergyDensity() + parameters.bare_mass * scalar_density_0 - parameters.G_S * pow(scalar_density_0, 2.0) * CONST_HBAR_C;
}

double EnergyDensity(double pressure,
					 double proton_chemical_potential,
					 double neutron_chemical_potential,
					 double proton_density,
					 double neutron_density)
{
	return - pressure + proton_chemical_potential * proton_density + neutron_chemical_potential * neutron_density;
}

double Pressure(double termodynamic_potential)
{
	double vacuum_energy_density = VacuumEnergyDensity();
	
	return - termodynamic_potential + vacuum_energy_density;
}

double TermodynamicPotential(double scalar_density,
							 double barionic_density,
							 double proton_density,
							 double neutron_density,
							 double proton_chemical_potential,
							 double neutron_chemical_potential,
							 double kinectic_energy_density)
{
	double rho_3 = proton_density - neutron_density;
	
	double omega = kinectic_energy_density + parameters.bare_mass * scalar_density;
	
	omega += - proton_chemical_potential * proton_density
			 - neutron_chemical_potential * neutron_density;
	
	omega += (- parameters.G_S * pow(scalar_density, 2.0)
			  + parameters.G_V * pow(barionic_density, 2.0)
			  + parameters.G_SV * pow(scalar_density * barionic_density, 2.0)
			  + parameters.G_RHO * pow(rho_3, 2.0)
			  + parameters.G_VRHO * pow(barionic_density * rho_3, 2.0)
              + parameters.G_SRHO * pow(scalar_density * rho_3, 2.0)) * CONST_HBAR_C;
	
	return omega;
}

double F0(double mass, double momentum)
{
	double E = sqrt(pow(mass, 2.0) + pow(momentum, 2.0));
	
	return (1.0 / 2.0) * (momentum * E - pow(mass, 2.0) * log((momentum + E) / mass));
}

double F2(double mass, double momentum)
{
	double E = sqrt(pow(mass, 2.0) + pow(momentum, 2.0));
				   
	return (1.0 / 8.0) * (-3.0 * pow(mass, 2.0) * momentum + 2.0 * pow(momentum, 3.0)) * E
           + (3.0 / 8.0) * pow(mass, 4.0) * log((momentum + E) / mass);
}
