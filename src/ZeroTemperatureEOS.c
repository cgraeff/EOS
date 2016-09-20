//
//  ZeroTemperatureEOS.c
//  hadrons EOS
//
//  Created by Clebson Graeff on 2016-02-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>

#include "Parameters.h"
#include "Constants.h"
#include "ZeroTemperatureEOS.h"
#include "RootFinding.h"

typedef struct _multi_dim_root_params{
    double barionic_density;
} multi_dim_root_params;

double F_E(double mass, double momentum);

int MultiDimensionalRootFinderHelperFunction(const gsl_vector *x,
                                             void             *p,
                                             gsl_vector       *return_values);

double ZeroMassSpecialCaseHelperFunction(double x,
                                         void  *par);

    // Mappings:
    //      The variables for the root finding are assumed to cover the range
    //      (-\infty, +\infty), but that is not the case for the variables
    //      that we are trying to solve. Here the proton fraction $y_p$ is
    //      such that
    //          $y_p \in [0, 1]$
    //      and the mass $m$ is such that
    //          $m \in [0, +\infty)$.
    //      To solve that, we use the mappings:
    //          $m = x^2$
    //      and
    //          $y_p = \frac{\tanh(x / a) + 1}{2}$.
    //      The initial guesses must be transformed by inverting the relations
    //      above
    // Handle special case: Zero mass case
    //  The zero mass case is important as most calculations will be
    //  performed at this particular case, which due to characteristics
    //  of the multidimensional root-finding algorithm, may be problematic to
    //  solve (it works most of the time, but sometimes calculations result in NaNs).
    //  This is due to problems in the calculation of derivatives of
    //  the function with respect to mass which arise from low variability
    //  of the function near zero mass.
    //
    //  This case, however is not the one reached at the start of calculations.
    //  In addition to that, once it is reached, all subsequent calculations are
    //  performed at approximatelly zero mass.
    //
    //  We take these characteristics into account and do the zero mass case
    //  with a special path, where we just assume mass = 0 (this effectivelly
    //  reduces the dimension of the system). This will avoid
    //  any calculation of potentially problematic derivatives. The special
    //  path is triggered by the condition
    //      parameters.multiroot.guesses.mass < parameters.multiroot.mass_tolerance
    //  which must be true. The tolerance should be adusted in the setup of
    //  parameters.
void SolveMultiRoots(double barionic_density,
double * return_mass,
double * return_proton_fraction)
{
    // Set dimension (number of equations|variables to solve|find)
    const int dimension = 2;

    // Set up parameters to be passed to helper function
    multi_dim_root_params p;
    p.barionic_density = barionic_density;

    // Check for zero mass special case. As mass != 0 is the
    // case that appears first, it is implemented first.
    if (parameters.multiroot.guesses.mass > parameters.multiroot.special_case.mass_tolerance){

        gsl_multiroot_function f;
        f.f = &MultiDimensionalRootFinderHelperFunction;
        f.n = dimension;
        f.params = (void *)&p;

        gsl_vector * initial_guess = gsl_vector_alloc(dimension);
        gsl_vector * return_results = gsl_vector_alloc(dimension);

        gsl_vector_set(initial_guess, 0, sqrt(parameters.multiroot.guesses.mass));
        gsl_vector_set(initial_guess, 1, parameters.multiroot.proton_fraction_mapping_scale
                                         * atanh(2.0 * parameters.multiroot.guesses.proton_fraction - 1.0));

        int status = MultidimensionalRootFinder(dimension,
        										&f,
                                                initial_guess,
                                                parameters.multiroot.abs_error,
                                                parameters.multiroot.rel_error,
                                                parameters.multiroot.max_iterations,
                                                return_results);

        if (status != 0){
            printf("Something is wrong with the rootfinding.\n");
            abort();
        }

        // Save results in return variables,
        // taking care of the mappinps
        *return_mass = pow(gsl_vector_get(return_results, 0), 2.0);
        *return_proton_fraction = (tanh(gsl_vector_get(return_results, 1)
                                        / parameters.multiroot.proton_fraction_mapping_scale) + 1.0) / 2.0;

        // Free vectors
        gsl_vector_free(initial_guess);
        gsl_vector_free(return_results);

        // Save solution as guess for next iteration
        if (parameters.multiroot.use_last_solution_as_guess == true){
            parameters.multiroot.guesses.mass = *return_mass;
            parameters.multiroot.guesses.proton_fraction = *return_proton_fraction;
        }

        return;

    } 
    else{ // Handle special case: Zero mass case

        gsl_function F;
        F.function = &ZeroMassSpecialCaseHelperFunction;
        F.params = &p;

        // Set root bounds observing the mappings
        double lower_bound = parameters.multiroot.proton_fraction_mapping_scale
                             * atanh(2.0 * parameters.multiroot.special_case.lower_bound - 1.0);
        double upper_bound = parameters.multiroot.proton_fraction_mapping_scale
                             * atanh(2.0 * parameters.multiroot.special_case.upper_bound - 1.0);

        // As we are left with just one variable and one equation to solve,
        // now an one-dimensional algorithm may be employed. Otherwise,
        // the dimension ought to be decreased by one an the multidimensional
        // employed again.
        double return_result;

        int status = UnidimensionalRootFinder(&F,
                                              lower_bound,
                                              upper_bound,
                                              parameters.multiroot.abs_error,
                                              parameters.multiroot.rel_error,
                                              parameters.multiroot.max_iterations,
                                              &return_result);

        if (status != 0){
            printf("Something is wrong with the rootfinding.\n");
            abort();
        }

        // Save results in return variables,
        // taking care of the mappings
        *return_mass = 0.0;
        *return_proton_fraction = (tanh(return_result / parameters.multiroot.proton_fraction_mapping_scale) + 1.0)
                                  / 2.0;
        return;
    }
}

double ZeroMassSpecialCaseHelperFunction(double x, void * par)
{
    const int dimension = 2;

    gsl_vector * input_values = gsl_vector_alloc(dimension);
    gsl_vector * return_values = gsl_vector_alloc(dimension);

    // Set mass = 0, which is our special case
    gsl_vector_set(input_values, 0, 0);

    // Pass value selected by the root finding routine
    gsl_vector_set(input_values, 1, x);

    MultiDimensionalRootFinderHelperFunction(input_values, par, return_values);

    double return_value =gsl_vector_get(return_values, 1);

    gsl_vector_free(input_values);
    gsl_vector_free(return_values);

    return return_value;
}

int MultiDimensionalRootFinderHelperFunction(const gsl_vector * x,
                                             void * p,
                                             gsl_vector * return_values)
{
    multi_dim_root_params * params = (multi_dim_root_params *)p;
    double barionic_density = params->barionic_density;

    // Mappings:
    //      The variables for the root finding are assumed to cover the range
    //      (-\infty, +\infty), but that is not the case for the variables
    //      that we are trying to solve. Here the proton fraction $y_p$ is
    //      such that
    //          $y_p \in [0, 1]$
    //      and the mass $m$ is such that
    //          $m \in [0, +\infty)$.
    //      To solve that, we use the mappings:
    //          $m = x^2$
    //      and
    //          $y_p = \frac{\tanh(x / a) + 1}{2}$.
    //      The initial guesses must be transformed by inverting the relations
    //      above

   	const double mass = pow(gsl_vector_get(x, 0), 2.0);
    const double proton_fraction = (tanh(gsl_vector_get(x, 1)
                                         / parameters.multiroot.proton_fraction_mapping_scale) + 1) / 2.0;

    double proton_density = proton_fraction * barionic_density;
    double neutron_density = (1.0 - proton_fraction) * barionic_density;

    // Ensure charge neutrality
    double electron_density = proton_density;

    double proton_fermi_momentum = FermiMomentum(proton_density);
	double neutron_fermi_momentum = FermiMomentum(neutron_density);
    double electron_fermi_momentum = FermiMomentum(electron_density);

    double total_scalar_density = ScalarDensity(mass, proton_fermi_momentum, parameters.theory.cutoff)
                                  + ScalarDensity(mass, neutron_fermi_momentum, parameters.theory.cutoff);
    
    double proton_chemical_potential = ProtonChemicalPotential(proton_fermi_momentum,
                                                               total_scalar_density,
                                                               mass,
                                                               barionic_density,
                                                               proton_density,
                                                               neutron_density);

    double neutron_chemical_potential = NeutronChemicalPotential(neutron_fermi_momentum,
                                                                 total_scalar_density,
                                                                 mass,
                                                                 barionic_density,
                                                                 proton_density,
                                                                 neutron_density);

    double electron_chemical_potential = sqrt(pow(electron_fermi_momentum, 2.0)
                                              + pow(parameters.theory.electron_mass, 2.0));
    
    // Set up parameters for gap equation
	gap_equation_input gap_params;
    gap_params.proton_fermi_momentum = proton_fermi_momentum;
    gap_params.neutron_fermi_momentum = neutron_fermi_momentum;
    gap_params.proton_density = proton_density;
    gap_params.neutron_density = neutron_density;

    double zeroed_gap_eq = GapEquation(mass, &gap_params);
    double zeroed_beta_equilibrium_relation = neutron_chemical_potential
                                         - proton_chemical_potential
                                         - electron_chemical_potential;
    
    // Prepare return vector
   	gsl_vector_set(return_values, 0, zeroed_gap_eq);
   	gsl_vector_set(return_values, 1, zeroed_beta_equilibrium_relation);

    return GSL_SUCCESS;
}

double FermiMomentum(double density)
{
    return CONST_HBAR_C * pow(3.0 * pow(M_PI, 2.0) * density, 1.0 / 3.0);
}

double SolveGapEquation(double proton_density,
                        double neutron_density,
                        double proton_fermi_momentum,
                        double neutron_fermi_momentum)
{
	gap_equation_input input;

	input.proton_density = proton_density;
    input.neutron_density = neutron_density;
	input.proton_fermi_momentum = proton_fermi_momentum;
	input.neutron_fermi_momentum = neutron_fermi_momentum;
	
	gsl_function F;
	F.function = &GapEquation;
	F.params = &input;

    double return_result;
    int status = UnidimensionalRootFinder(&F,
                                          parameters.gap_equation_solver.lower_bound,
                                          parameters.gap_equation_solver.upper_bound,
                                          parameters.gap_equation_solver.abs_error,
                                          parameters.gap_equation_solver.rel_error,
                                          parameters.gap_equation_solver.max_iterations,
                                          &return_result);

    // If the bounds don't straddle the root (-1 return status),
    // we assume that the only solution is mass = 0
    if (status == -1){
        return 0;
    }

	return return_result;
}

double GapEquation(double mass, void * input)
{
	gap_equation_input * param = (gap_equation_input *)input;

    double barionic_density = param->proton_density + param->neutron_density;
    double rho_3 = param->proton_density - param->neutron_density;

	double scalar_density = ScalarDensity(mass, param->proton_fermi_momentum, parameters.theory.cutoff)
                            + ScalarDensity(mass, param->neutron_fermi_momentum, parameters.theory.cutoff);

	double gap_1st_term = 2.0  * CONST_HBAR_C * parameters.theory.G_S * scalar_density;
	double gap_2nd_term = - 2.0 * CONST_HBAR_C * parameters.theory.G_SV * scalar_density * pow(barionic_density, 2.0);
    double gap_3rd_term = - 2.0  * CONST_HBAR_C * parameters.theory.G_SRHO * scalar_density * pow(rho_3, 2.0);

	return mass + gap_1st_term + gap_2nd_term + gap_3rd_term - parameters.theory.bare_mass;
}

double ScalarDensity(double mass, double fermi_momentum, double cutoff)
{
    if (mass == 0.0){
        return 0.0;
    }

	return pow(CONST_HBAR_C, -3.0) * (mass / pow(M_PI, 2.0)) * (F0(mass, fermi_momentum) - F0(mass, cutoff));
}

double VacuumScalarDensity()
{
	return 2.0 * pow(CONST_HBAR_C, -3.0) * (parameters.theory.nucleon_mass / pow(M_PI, 2.0))
			* (F0(parameters.theory.nucleon_mass, 0.0) - F0(parameters.theory.nucleon_mass, parameters.theory.cutoff));
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

	chemical_potential += 2.0 * CONST_HBAR_C * (parameters.theory.G_V * barionic_density
												+ parameters.theory.G_SV * barionic_density * pow(scalar_density, 2.0)
												+ parameters.theory.G_RHO * rho_3
												+ parameters.theory.G_VRHO * pow(rho_3, 2.0) * barionic_density
												+ parameters.theory.G_VRHO * pow(barionic_density, 2.0) * rho_3
                                                + parameters.theory.G_SRHO * pow(scalar_density, 2.0) * rho_3);

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
	
	chemical_potential += 2.0 * CONST_HBAR_C * (parameters.theory.G_V * barionic_density
                                                + parameters.theory.G_SV * barionic_density * pow(scalar_density, 2.0)
                                                - parameters.theory.G_RHO * rho_3
                                                + parameters.theory.G_VRHO * pow(rho_3, 2.0) * barionic_density
                                                - parameters.theory.G_VRHO * pow(barionic_density, 2.0) * rho_3
                                                - parameters.theory.G_SRHO * pow(scalar_density, 2.0) * rho_3);

	return chemical_potential;
}


double KinecticEnergyDensity(double mass, double proton_fermi_momentum, double neutron_fermi_momentum)
{
	double proton_kinectic_energy = (CTE_NUM_COLORS  / pow(M_PI, 2.0))
                                    * (F2(mass, proton_fermi_momentum) - F2(mass, parameters.theory.cutoff));
	double neutron_kinectic_energy = (CTE_NUM_COLORS / pow(M_PI, 2.0))
                                     * (F2(mass, neutron_fermi_momentum) - F2(mass, parameters.theory.cutoff));

	return (proton_kinectic_energy + neutron_kinectic_energy) / (pow(CONST_HBAR_C, 3.0));
}

double VacuumKinecticEnergyDensity()
{
	return 2.0 * (CTE_NUM_COLORS / pow(M_PI, 2.0)) * (pow(CONST_HBAR_C, -3.0))
		   * (F2(parameters.theory.nucleon_mass, 0) - F2(parameters.theory.nucleon_mass, parameters.theory.cutoff));
}

double VacuumEnergyDensity()
{
    double scalar_density_0 = VacuumScalarDensity();

	return VacuumKinecticEnergyDensity()
           + parameters.theory.bare_mass * scalar_density_0
           - parameters.theory.G_S * pow(scalar_density_0, 2.0) * CONST_HBAR_C;
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

	double omega = kinectic_energy_density + parameters.theory.bare_mass * scalar_density;
	
	omega += - proton_chemical_potential * proton_density
			 - neutron_chemical_potential * neutron_density;
	
	omega += (- parameters.theory.G_S * pow(scalar_density, 2.0)
			  + parameters.theory.G_V * pow(barionic_density, 2.0)
			  + parameters.theory.G_SV * pow(scalar_density * barionic_density, 2.0)
			  + parameters.theory.G_RHO * pow(rho_3, 2.0)
			  + parameters.theory.G_VRHO * pow(barionic_density * rho_3, 2.0)
              + parameters.theory.G_SRHO * pow(scalar_density * rho_3, 2.0)) * CONST_HBAR_C;
	
	return omega;
}

double ElectronPressure(double electron_fermi_momentum)
{
    return F2(parameters.theory.electron_mass, electron_fermi_momentum) / (3.0 * pow(M_PI, 2.0) * pow(CONST_HBAR_C, 3.0));
}

double ElectronEnergyDensity(double electron_fermi_momentum)
{
    return F_E(parameters.theory.electron_mass, electron_fermi_momentum) / (pow(M_PI, 2.0) * pow(CONST_HBAR_C, 3.0));
}

double F0(double mass, double momentum)
{
	double E = sqrt(pow(mass, 2.0) + pow(momentum, 2.0));
	
    double first_term = momentum * E;
    double second_term = - pow(mass, 2.0) * log(momentum + E);
    double third_term = pow(mass, 2.0) * log(mass);

    // In the limit of mass -> 0, the logarithm may give a NaN,
    // so when that happens, just make sure the third term
    // have the correct value
    if (third_term != third_term)
        third_term = 0;
    
    return (1.0 / 2.0) * (first_term +second_term + third_term);
}

double F2(double mass, double momentum)
{
    if (mass == 0.0){
        return pow(momentum, 4.0) / 4.0;
    }
    
	double E = sqrt(pow(mass, 2.0) + pow(momentum, 2.0));
				   
	return (1.0 / 8.0) * (-3.0 * pow(mass, 2.0) * momentum + 2.0 * pow(momentum, 3.0)) * E
           + (3.0 / 8.0) * pow(mass, 4.0) * log((momentum + E) / mass);
}

double F_E(double mass, double momentum)
{
    if (mass == 0)
        return pow(momentum, 4.0) / 4.0;
    
    double E = sqrt(pow(mass, 2.0) + pow(momentum, 2.0));
    
    return (momentum * pow(E, 3.0)
            - 0.5 * pow(mass, 2.0) * momentum * E
            - 0.5 * pow(mass, 4.0) * log ((momentum + E) / mass))
    / 4.0;
}
