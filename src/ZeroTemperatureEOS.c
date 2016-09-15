//
//  EOS.c
//  EOS
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

typedef struct _multi_dim_root_params{
    double barionic_density;
} multi_dim_root_params;

double F_E(double mass, double momentum);

int MultiDimensionalRootFinderHelperFunction(const gsl_vector * x,
                                             void * p,
                                             gsl_vector * return_values);

void SolveMultiRoots(double barionic_density, double * return_mass, double * return_proton_fraction)
{
    
    // Set up parameters to be passed to helper function
    multi_dim_root_params p;
    p.barionic_density = barionic_density;
    
    // Handle special case: mass near zero
    //  When running using last solution as guess, we may
    //  avoid the calculation of a two dimensional problem
    //  when mass is (near) zero
    if (parameters.multiroot.use_last_solution_as_guess == true){
        if (parameters.multiroot.guesses.mass < parameters.multiroot.mass_tolerance){

            const gsl_root_fsolver_type * T	= gsl_root_fsolver_bisection;
            gsl_root_fsolver * s = gsl_root_fsolver_alloc(T);
            
            gsl_function F;
            F.function = &GapEquation;
            F.params = &p;
            
            double lower_bound = parameters.multiroot.proton_fraction_mapping_scale
                                * atanh(2.0 * parameters.multiroot.special_case.lower_bound - 1.0);
            double upper_bound = parameters.multiroot.proton_fraction_mapping_scale
                                * atanh(2.0 * parameters.multiroot.special_case.upper_bound - 1.0);
            
            gsl_root_fsolver_set(s, &F, lower_bound, upper_bound);
            
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
            } while(GSL_CONTINUE == gsl_root_test_interval(x_lower,
                                                           x_upper,
                                                           CONST_ABS_ERROR_GAP_EQ_SOLVING, // change this as we are reusing
                                                           CONST_REL_ERROR_GAP_EQ_SOLVING) // error variables
                    && i <= CONST_MAX_ITERATIONS_GAP_EQ_SOLVING); // this too
            
            double result = gsl_root_fsolver_root(s);
            
            void gsl_root_fsolver_free(gsl_root_fsolver * S);
            
            *return_mass = 0.0;
            *return_proton_fraction = (tanh(result / parameters.multiroot.proton_fraction_mapping_scale) + 1.0) / 2.0;
            return;
        }
    }
    
    // Set dimension (number of equations|variables to solve|find)
    const int dimension = 2;

    gsl_multiroot_function f;
    f.f = &MultiDimensionalRootFinderHelperFunction;
    f.n = dimension;
    f.params = (void *)&p;

    gsl_vector * initial_guess = gsl_vector_alloc(dimension);
    
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
    gsl_vector_set(initial_guess, 0, sqrt(parameters.multiroot.guesses.mass));
    gsl_vector_set(initial_guess, 1, parameters.multiroot.proton_fraction_mapping_scale
                                     * atanh(2.0 * parameters.multiroot.guesses.proton_fraction - 1.0));

    int status;
    size_t iter = 0;
    
    const gsl_multiroot_fsolver_type * solver_type = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver * solver = gsl_multiroot_fsolver_alloc(solver_type,
                                                                 dimension);
    
    gsl_multiroot_fsolver_set(solver, &f, initial_guess);
    
    do {
        iter++;
        status = gsl_multiroot_fsolver_iterate(solver);
        
        if (status == GSL_EBADFUNC){
            printf("TwodimensionalRootFinder: Error: Infinity or division by zero.\n");
            abort();
        }
        else if (status == GSL_ENOPROG){
            printf("TwodimensionalRootFinder: Error: Solver is stuck. Try a different initial guess.\n");
            abort();
        }
        
        // Check if the root is good enough:
        // tests for the convergence of the sequence by comparing the last step dx with the
        // absolute error epsabs and relative error epsrel to the current position x. The test
        // returns GSL_SUCCESS if the following condition is achieved,
        //
        // |dx_i| < epsabs + epsrel |x_i|
        
        gsl_vector * x = gsl_multiroot_fsolver_root(solver); // current root
        gsl_vector * dx = gsl_multiroot_fsolver_dx(solver); // last step
        
        status = gsl_multiroot_test_delta(dx,
                                          x,
                                          parameters.multiroot.abs_error,
                                          parameters.multiroot.rel_error);
        
    } while (status == GSL_CONTINUE
             && iter < parameters.multiroot.max_iterations);
    
    // Save results in return variables
    *return_mass = pow(gsl_vector_get(gsl_multiroot_fsolver_root(solver), 0), 2.0);;
    *return_proton_fraction = (tanh(gsl_vector_get(gsl_multiroot_fsolver_root(solver), 1)
                                    / parameters.multiroot.proton_fraction_mapping_scale) + 1.0) / 2.0;
    
    // Free vectors
    gsl_multiroot_fsolver_free(solver);
    gsl_vector_free(initial_guess);
    
    // Save solution as guess for next iteration
    if (parameters.multiroot.use_last_solution_as_guess == true){
        parameters.multiroot.guesses.mass = *return_mass;
        parameters.multiroot.guesses.proton_fraction = *return_proton_fraction;
    }

    return;
}

double SpecialCaseHelperFunction(double x, void * par)
{
    gsl_vector * input_values = gsl_vector_alloc(2);
    gsl_vector * return_values = gsl_vector_alloc(2);
    
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
    
    double total_scalar_density = ScalarDensity(mass, proton_fermi_momentum, parameters.CUTOFF)
                                  + ScalarDensity(mass, neutron_fermi_momentum, parameters.CUTOFF);
    
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
    
    double electron_chemical_potential = sqrt(pow(electron_fermi_momentum, 2.0) + pow(parameters.electron_mass, 2.0));
    
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
    // Maybe this would be better: gsl_root_fsolver_brent
	const gsl_root_fsolver_type * T	= gsl_root_fsolver_bisection;

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
	} while(GSL_CONTINUE == gsl_root_test_interval(x_lower,
                                                   x_upper,
                                                   CONST_ABS_ERROR_GAP_EQ_SOLVING,
                                                   CONST_REL_ERROR_GAP_EQ_SOLVING)
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
	
	double scalar_density = ScalarDensity(mass, param->proton_fermi_momentum, parameters.CUTOFF)
                            + ScalarDensity(mass, param->neutron_fermi_momentum, parameters.CUTOFF);
	
	double gap_1st_term = 2.0  * CONST_HBAR_C * parameters.G_S * scalar_density;
	double gap_2nd_term = - 2.0 * CONST_HBAR_C * parameters.G_SV * scalar_density * pow(barionic_density, 2.0);
    double gap_3rd_term = - 2.0  * CONST_HBAR_C * parameters.G_SRHO * scalar_density * pow(rho_3, 2.0);

	return mass + gap_1st_term + gap_2nd_term + gap_3rd_term - parameters.bare_mass;
}

double ScalarDensity(double mass, double fermi_momentum, double cutoff)
{
    if (mass == 0.0)
        return 0.0;
    
	return pow(CONST_HBAR_C, -3.0) * (mass / pow(M_PI, 2.0)) * (F0(mass, fermi_momentum) - F0(mass, cutoff));
}

double VacuumScalarDensity()
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
	double proton_kinectic_energy = (CTE_NUM_COLORS  / pow(M_PI, 2.0))
                                    * (F2(mass, proton_fermi_momentum) - F2(mass, parameters.CUTOFF));
	double neutron_kinectic_energy = (CTE_NUM_COLORS / pow(M_PI, 2.0))
                                     * (F2(mass, neutron_fermi_momentum) - F2(mass, parameters.CUTOFF));
	
	return (proton_kinectic_energy + neutron_kinectic_energy) / (pow(CONST_HBAR_C, 3.0));
}

double VacuumKinecticEnergyDensity()
{
	return 2.0 * (CTE_NUM_COLORS / pow(M_PI, 2.0)) * (pow(CONST_HBAR_C, -3.0))
		   * (F2(parameters.nucleon_mass, 0) - F2(parameters.nucleon_mass, parameters.CUTOFF));
}

double VacuumEnergyDensity()
{
    double scalar_density_0 = VacuumScalarDensity();

	return VacuumKinecticEnergyDensity()
           + parameters.bare_mass * scalar_density_0
           - parameters.G_S * pow(scalar_density_0, 2.0) * CONST_HBAR_C;
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

double ElectronPressure(double electron_fermi_momentum)
{
    return F2(parameters.electron_mass, electron_fermi_momentum) / (3.0 * pow(M_PI, 2.0) * pow(CONST_HBAR_C, 3.0));
}

double ElectronEnergyDensity(double electron_fermi_momentum)
{
    return F_E(parameters.electron_mass, electron_fermi_momentum) / (pow(M_PI, 2.0) * pow(CONST_HBAR_C, 3.0));
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
    if (mass == 0.0)
        return pow(momentum, 4.0) / 4.0;
    
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