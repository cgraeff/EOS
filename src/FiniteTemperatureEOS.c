//
//  FiniteTemperatureEOS.c
//  hadrons EOS
//
//  Created by Clebson Graeff on 2016-11-08.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <gsl/gsl_integration.h>
#include "FiniteTemperatureEOS.h"
#include "Constants.h"
#include "RootFinding.h"
#include "Parameters.h"

typedef struct _ft_multi_dim_root_params{
    double proton_barionic_density;
    double neutron_barionic_density;
} ft_multi_dim_root_params;

double OnedimensionalIntegrator(gsl_function * F, double lower_limit, double upper_limit);

double BarionicDensityAtFiniteTemperature(double mass,
                                          double chemical_potential,
                                          double upper_limit);

int FiniteTemperatureMultiDimensionalRootFinderFunction(const gsl_vector   *x,
                                                        void               *p,
                                                        gsl_vector         *return_values);
int FiniteTemperatureZeroMassSpecialCaseHelperFunction(const gsl_vector   *x,
                                                       void               *p,
                                                       gsl_vector         *return_values);

double FermiDiracDistributionForParticles(double energy,
                                          double chemical_potential,
                                          double temperature);
double FermiDiracDistributionForAntiparticles(double energy,
                                              double chemical_potential,
                                              double temperature);

int SolveMultiRootsFiniteTemp(double  proton_density,
                              double  neutron_density,
                              double *return_mass,
                              double *return_proton_renormalized_chemical_potential,
                              double *return_neutron_renormalized_chemical_potential)
{
    // Set up parameters to be passed to helper function
    ft_multi_dim_root_params p;
    p.proton_barionic_density = proton_density;
    p.neutron_barionic_density = neutron_density;
    
    // Check for zero mass special case. As mass != 0 is the
    // case that appears first, it is implemented first.
    if (parameters.finite_temperature.guesses.mass > parameters.multiroot.special_case.mass_tolerance){
        // FIXME: this is really wrong, the test
        // should not test the initial guess, but the current value of the guess. This will work only if we are
        // using the last solution as guess. The simplest way to solve this is always use the last solution as guess
        
        // Set dimension (number of equations|variables to solve|find)
        const int dimension = 3;
        
        gsl_multiroot_function f;
        f.f = &FiniteTemperatureMultiDimensionalRootFinderFunction;
        f.n = dimension;
        f.params = (void *)&p;
        
        gsl_vector * initial_guess = gsl_vector_alloc(dimension);
        gsl_vector * return_results = gsl_vector_alloc(dimension);

        gsl_vector_set(initial_guess,
                       0,
                       sqrt(parameters.finite_temperature.guesses.mass));
        gsl_vector_set(initial_guess,
                       1,
                       sqrt(parameters.finite_temperature.guesses.proton_renormalized_chemical_potential));
        gsl_vector_set(initial_guess,
                       2,
                       sqrt(parameters.finite_temperature.guesses.proton_renormalized_chemical_potential));
        
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
        // taking care of the mappings
        *return_mass = pow(gsl_vector_get(return_results, 0), 2.0);
        *return_proton_renormalized_chemical_potential = pow(gsl_vector_get(return_results, 1), 2.0);
        *return_neutron_renormalized_chemical_potential = pow(gsl_vector_get(return_results, 2), 2.0);

        
        // Free vectors
        gsl_vector_free(initial_guess);
        gsl_vector_free(return_results);
        
        // Save solution as guess for next iteration
        if (parameters.finite_temperature.use_last_solution_as_guess == true){
            parameters.finite_temperature.guesses.mass = *return_mass;
            parameters.finite_temperature.guesses.proton_renormalized_chemical_potential =
                *return_proton_renormalized_chemical_potential;
            parameters.finite_temperature.guesses.neutron_renormalized_chemical_potential =
                *return_neutron_renormalized_chemical_potential;
        }
        
        return 0;
        
    }
    else{ // Handle special case: Zero mass case
        
        // Set dimension (number of equations|variables to solve|find)
        const int dimension = 2;

        gsl_multiroot_function f;
        f.f = &FiniteTemperatureZeroMassSpecialCaseHelperFunction;
        f.n = dimension;
        f.params = (void *)&p;
        
        gsl_vector * initial_guess = gsl_vector_alloc(dimension);
        gsl_vector * return_results = gsl_vector_alloc(dimension);
        
        gsl_vector_set(initial_guess,
                       0,
                       sqrt(parameters.finite_temperature.guesses.mass));
        gsl_vector_set(initial_guess,
                       1,
                       sqrt(parameters.finite_temperature.guesses.proton_renormalized_chemical_potential));
        gsl_vector_set(initial_guess,
                       2,
                       sqrt(parameters.finite_temperature.guesses.proton_renormalized_chemical_potential));
        
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
        
        
        // Free vectors
        gsl_vector_free(initial_guess);
        gsl_vector_free(return_results);
        
        // Save solution as guess for next iteration
        if (parameters.finite_temperature.use_last_solution_as_guess == true){
            parameters.finite_temperature.guesses.mass = *return_mass;
            parameters.finite_temperature.guesses.proton_renormalized_chemical_potential =
                *return_proton_renormalized_chemical_potential;
            parameters.finite_temperature.guesses.neutron_renormalized_chemical_potential =
                *return_neutron_renormalized_chemical_potential;
        }
        

    }

    return 0;
}

int FiniteTemperatureZeroMassSpecialCaseHelperFunction(const gsl_vector   *x,
                                                       void               *p,
                                                       gsl_vector         *return_values)
{
    const int dimension = 3;
    
    gsl_vector * input_values = gsl_vector_alloc(dimension);
    gsl_vector * results = gsl_vector_alloc(dimension);
    
    // Set mass = 0, which is our special case
    gsl_vector_set(input_values, 0, 0);
    
    // Pass value selected by the root finding routine
    gsl_vector_set(input_values, 1, gsl_vector_get(x, 0));
    gsl_vector_set(input_values, 2, gsl_vector_get(x, 1));
    
    FiniteTemperatureMultiDimensionalRootFinderFunction(input_values, p, results);
    
    gsl_vector_set(return_values, 0, gsl_vector_get(results, 1));
    gsl_vector_set(return_values, 1, gsl_vector_get(results, 2));
    
    gsl_vector_free(input_values);
    gsl_vector_free(results);
    
    return 0;
}

int FiniteTemperatureMultiDimensionalRootFinderFunction(const gsl_vector   *x,
                                                        void               *p,
                                                        gsl_vector         *return_values)
{
    // receive constant values
    ft_multi_dim_root_params * params = (ft_multi_dim_root_params *)p;
    double proton_barionic_density = params->proton_barionic_density;
    double neutron_barionic_density = params->neutron_barionic_density;
    double total_barionic_density = proton_barionic_density + neutron_barionic_density;
    double rho_3 = proton_barionic_density - neutron_barionic_density;
    
    // take care of mappings
   	const double mass = pow(gsl_vector_get(x, 0), 2.0);
    const double proton_renormalized_chemical_potential = pow(gsl_vector_get(x, 1), 2.0);
    const double neutron_renormalized_chemical_potential = pow(gsl_vector_get(x, 2), 2.0);
    
    // calculate things to be zeroed
    double proton_scalar_density = ScalarDensityAtFiniteTemperature(mass,
                                                                    proton_renormalized_chemical_potential,
                                                                    parameters.theory.cutoff);
    double neutron_scalar_density = ScalarDensityAtFiniteTemperature(mass,
                                                                     neutron_renormalized_chemical_potential,
                                                                     parameters.theory.cutoff);
    double total_scalar_density = proton_scalar_density + neutron_scalar_density;
    
    // the zeroed gap eq is probably the same as in the zero temp case
    // try to reuse by using a function
    double zeroed_gap = mass
                        - parameters.theory.bare_mass
                        - 2.0 * parameters.theory.G_S * total_scalar_density
                        + 2.0 * parameters.theory.G_SV * pow(total_scalar_density, 2.0) * total_barionic_density
                        + 2.0 * parameters.theory.G_SRHO * pow(rho_3, 2.0) * total_scalar_density;
    
    double proton_barionic_density_integral = BarionicDensityAtFiniteTemperature(mass,
                                                                                 proton_renormalized_chemical_potential,
                                                                                 parameters.theory.cutoff);
    double neutron_barionic_density_integral = BarionicDensityAtFiniteTemperature(mass,
                                                                                  neutron_renormalized_chemical_potential,
                                                                                  parameters.theory.cutoff);
    
    double zeroed_p_dens_gap = proton_barionic_density - proton_barionic_density_integral;
    double zeroed_n_dens_gap = neutron_barionic_density - neutron_barionic_density_integral;
    
    // Prepare return vector
   	gsl_vector_set(return_values, 0, zeroed_gap);
   	gsl_vector_set(return_values, 1, zeroed_p_dens_gap);
    gsl_vector_set(return_values, 2, zeroed_n_dens_gap);
    
    return GSL_SUCCESS;
}

typedef struct _fermi_dirac_distrib_integrand{
    double mass;
    double chemical_potential;
    double temperature;
}fermi_dirac_distrib_integrand;

double ScalarDensityAtFiniteTemperatureIntegrand(double momentum, void * params);

double ScalarDensityAtFiniteTemperature(double mass,
                                        double renormalized_chemical_potential,
                                        double upper_limit)
{
    fermi_dirac_distrib_integrand p;
    p.mass = mass;
    p.chemical_potential = renormalized_chemical_potential;
    p.temperature = parameters.temperature;
    
    gsl_function F;
    F.function = &ScalarDensityAtFiniteTemperatureIntegrand;
    F.params = &p;
    
    double integral = OnedimensionalIntegrator(&F, 0.0, upper_limit);
    
    return CTE_NUM_COLORS * CTE_NUM_FLAVORS * mass * integral / (pow(M_PI, 2.0) * pow(CONST_HBAR_C, 3.0));
}

// FIXME: move this
double FermionNumberDistributionPlus(double energy, double chemical_potential, double temperature);
double FermionNumberDistributionMinus(double energy, double chemical_potential, double temperature);


double ScalarDensityAtFiniteTemperatureIntegrand(double momentum, void * params)
{
    fermi_dirac_distrib_integrand * p = (fermi_dirac_distrib_integrand *) params;
    
    double E = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));
    
    double plus_term = FermionNumberDistributionPlus(E,
                                                     p->chemical_potential,
                                                     p->temperature);
    double minus_term = FermionNumberDistributionMinus(E,
                                                       p->chemical_potential,
                                                       p->temperature);
    return (plus_term - minus_term) * pow(momentum, 2.0) / E;
}

double BarionicDensityAtFiniteTemperatureIntegrand(double momentum, void * params);

double BarionicDensityAtFiniteTemperature(double mass,
                                          double chemical_potential,
                                          double upper_limit)
{
    fermi_dirac_distrib_integrand p;
    p.mass = mass;
    p.chemical_potential = chemical_potential;
    p.temperature = parameters.temperature;
    
    gsl_function F;
    F.function = &BarionicDensityAtFiniteTemperatureIntegrand;
    F.params = &p;
    
    double integral = OnedimensionalIntegrator(&F, 0.0, upper_limit);
    
    return CTE_NUM_COLORS * CTE_NUM_FLAVORS * mass * integral / (pow(M_PI, 2.0) * pow(CONST_HBAR_C, 3.0));
}

double BarionicDensityAtFiniteTemperatureIntegrand(double momentum, void * params)
{
    fermi_dirac_distrib_integrand * p = (fermi_dirac_distrib_integrand *) params;
    
    double E = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));
    
    double plus_term = FermionNumberDistributionPlus(E,
                                                     p->chemical_potential,
                                                     p->temperature);
    double minus_term = FermionNumberDistributionMinus(E,
                                                       p->chemical_potential,
                                                       p->temperature);
    return (plus_term + minus_term - 1.0) * pow(momentum, 2.0) / E;
}

// FIXME: Rewrite functions that depend on these distributions so that they use
// the others named FermiDiracDistributionFor{Particles|Antiparticles}
double FermionNumberDistributionPlus(double energy, double chemical_potential, double temperature)
{
    return 1.0 / (1.0 + exp((energy - chemical_potential) / temperature));
}

double FermionNumberDistributionMinus(double energy, double chemical_potential, double temperature)
{
    return 1.0 / (1.0 + exp(-(energy + chemical_potential) / temperature));
}

typedef struct _entropy_integrand_parameters{
    double mass;
    double temperature;
    double renormalized_chemical_potential;
} entropy_integrand_parameters;

double g(double t, double e, double c)
{
    double a = - (e - c)/t;
    return log1p(exp(a)) - a * exp(a) / (1.0 + exp(a));
}

double EntropyIntegrand(double momentum, void * parameters)
{
    entropy_integrand_parameters * p = (entropy_integrand_parameters *)parameters;
    
    double energy = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));
    
    double first_term = g(p->temperature, energy, p->renormalized_chemical_potential);
    double second_term = g(p->temperature, energy, -p->renormalized_chemical_potential);
    
    return pow(momentum, 2.0) * (first_term + second_term);
}

double Entropy(double mass, double temperature, double renormalized_chemical_potential)
{
    int interval_num = 1000;
    
    entropy_integrand_parameters p;
    p.mass = mass;
    p.temperature = temperature;
    p.renormalized_chemical_potential = renormalized_chemical_potential;
    
    gsl_function F;
    F.function = &EntropyIntegrand;
    F.params = &p;
    
    gsl_integration_workspace * workspace =
    gsl_integration_workspace_alloc(interval_num);
    
    // FIXME: move parameters to Parameters.[h,c]
    double integral = 0;
    double abserr;
    double lower_limit = 0.0;
    double upper_limit = parameters.theory.cutoff;
    double abs_error = 1.0E-3;
    double rel_error = 1.0E-3;
    int max_sub_interval = interval_num;
    int integration_key = GSL_INTEG_GAUSS61;
    
    gsl_integration_qag(&F,
                        lower_limit,
                        upper_limit,
                        abs_error,
                        rel_error,
                        max_sub_interval,
                        integration_key,
                        workspace,
                        &integral,
                        &abserr);
    
    gsl_integration_workspace_free(workspace);
    
    return CTE_NUM_COLORS * CTE_NUM_FLAVORS * pow(CONST_HBAR_C, -3.0) * integral / pow(M_PI, 2.0);
}

// FIXME: Remove once I'm sure the other version is working ok
double EntropyIntegrandFromDerivative(double momentum, void * parameters)
{
    // This is the expression I obtained from $- \partial \omega / \partial T$
    
    entropy_integrand_parameters * p = (entropy_integrand_parameters *)parameters;
    
    double energy = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));
    
    double fd_dist_part = FermiDiracDistributionForParticles(energy,
                                                             p->renormalized_chemical_potential,
                                                             p->temperature);
    
    double fd_dist_antipart = FermiDiracDistributionForAntiparticles(energy,
                                                                     p->renormalized_chemical_potential,
                                                                     p->temperature);
    
    double particles_term = - log1p(-fd_dist_part)
    + (energy - p->renormalized_chemical_potential) * fd_dist_part / p->temperature;
    
    double antiparticles_term = - log1p(-fd_dist_antipart)
    + (energy + p->renormalized_chemical_potential) * fd_dist_antipart / p->temperature;
    
    return pow(momentum, 2.0) * (particles_term + antiparticles_term);
}

double EntropyIntegrandArt(double momentum, void * parameters)
{
    entropy_integrand_parameters * p = (entropy_integrand_parameters *)parameters;
    
    double energy = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));
    
    double fd_dist_part = FermiDiracDistributionForParticles(energy,
                                                             p->renormalized_chemical_potential,
                                                             p->temperature);
    
    double fd_dist_antipart = FermiDiracDistributionForAntiparticles(energy,
                                                                     p->renormalized_chemical_potential,
                                                                     p->temperature);
    
    double particles_term = fd_dist_part * log(fd_dist_part)
    + (1.0 - fd_dist_part) * log1p(-fd_dist_part);
    
    double antiparticles_term = fd_dist_antipart * log(fd_dist_antipart)
    + (1.0 - fd_dist_antipart) * log1p(-fd_dist_antipart);
    
    return -pow(momentum, 2.0) * (particles_term + antiparticles_term);
}

double FermiDiracDistributionForParticles(double energy,
                                          double chemical_potential,
                                          double temperature)
{
    return 1.0 / (1.0 + exp((energy - chemical_potential)/temperature));
}

double FermiDiracDistributionForAntiparticles(double energy,
                                              double chemical_potential,
                                              double temperature)
{
    return 1.0 / (1.0 + exp((energy + chemical_potential)/temperature));
}

double OnedimensionalIntegrator(gsl_function * F, double lower_limit, double upper_limit)
{
    gsl_integration_workspace * workspace =
    gsl_integration_workspace_alloc(parameters.fermi_dirac_integrals_max_interval_num);
    
    double integral = 0;
    double abserr;
    gsl_integration_qag(F,
                        lower_limit,
                        upper_limit,
                        parameters.fermi_dirac_integrals_abs_error,
                        parameters.fermi_dirac_integrals_rel_error,
                        parameters.fermi_dirac_integrals_max_sub_interval,
                        parameters.fermi_dirac_integrals_integration_key,
                        workspace,
                        &integral,
                        &abserr);

    gsl_integration_workspace_free(workspace);
    
    return integral;
}

double FiniteTemperatureProtonChemicalPotential(double renormalized_proton_chemical_potential,
                                                double mass,
                                                double total_scalar_density,
                                                double proton_barionic_density,
                                                double neutron_barionic_density)
{
    double total_barionic_density = proton_barionic_density + neutron_barionic_density;
    double rho_3 = proton_barionic_density - neutron_barionic_density;

    return renormalized_neutron_chemical_potential
           - 2.0 * parameters.theory.G_V * total_barionic_density
           - 2.0 * parameters.theory.G_SV * total_barionic_density * pow(total_scalar_density, 2.0)
           - 2.0 * parameters.theory.G_RHO * rho_3;
           - 2.0 * parameters.theory.G_VRHO * total_barionic_density * pow(rho_3, 2.0)
           - 2.0 * parameters.theory.G_VRHO * rho_3 * pow(total_barionic_density, 2.0)
           - 2.0 * parameters.theory.G_SRHO * rho_3 * pow(total_scalar_density, 2.0);
}

double FiniteTemperatureNeutronChemicalPotential(double renormalized_neutron_chemical_potential,
                                                 double mass,
                                                 double total_scalar_density,
                                                 double proton_barionic_density,
                                                 double neutron_barionic_density)
{
    double total_barionic_density = proton_barionic_density + neutron_barionic_density;
    double rho_3 = proton_barionic_density - neutron_barionic_density;

    return renormalized_neutron_chemical_potential
           - 2.0 * parameters.theory.G_V * total_barionic_density
           - 2.0 * parameters.theory.G_SV * total_barionic_density * pow(total_scalar_density, 2.0)
           + 2.0 * parameters.theory.G_RHO * rho_3;
           - 2.0 * parameters.theory.G_VRHO * total_barionic_density * pow(rho_3, 2.0)
           + 2.0 * parameters.theory.G_VRHO * rho_3 * pow(total_barionic_density, 2.0)
           + 2.0 * parameters.theory.G_SRHO * rho_3 * pow(total_scalar_density, 2.0);
}

// FIXME: just call this FermiDiracDistributionParameters would
// be much better.
typedef struct _ft_temp_kinectic_energy_integrand_params{
    double mass;
    double chemical_potential;
    double temperature;
} ft_temp_kinectic_energy_integrand_params;

double FiniteTemperatureKinecticEnergyDensity(double mass,
                                              double proton_renormalized_chemical_potential,
                                              double neutron_renormalized_chemical_potential)
{
    gsl_function F;
    F.function = &FiniteTemperatureKinecticEnergyDensityIntegrand;

    ft_temp_kinectic_energy_integrand_params *params;
    params.mass = mass;
    params.chemical_potential = proton_renormalized_chemical_potential;
    params.temperature = &(parameters.temperature);

    F.params = params;

    double integral = OneDimensionalIntegrator(&F, 0.0, parameters.theory.cutoff);

    double proton_kinectic_energy = CTE_NUM_COLORS * integral / pow(M_PI, 2.0);

    params.chemical_potential = proton_renormalized_chemical_potential;

    double integral = OneDimensionalIntegrator(&F, 0.0, parameters.theory.cutoff);

    double neutron_kinectic_energy = CTE_NUM_COLORS * integral / pow(M_PI, 2.0);

    return proton_kinectic_energy + neutron_kinectic_energy;
}

double FiniteTemperatureKinecticEnergyDensityIntegrand(double momentum, void *par)
{
    ft_temp_kinectic_energy_integrand_params *params = (ft_temp_kinectic_energy_integrand_params *)par;

    double energy = sqrt(pow(momentum, 2.0) + pow(params->mass, 2.0));

    double dist_plus = FermionNumberDistributionPlus(energy,
                                                     params->chemical_potential,
                                                     params->temperature);

    double dist_minus = FermionNumberDistributionMinus(energy,
                                                       params->chemical_potential,
                                                       params->temperature);

    return pow(momentum, 2.0) * (dist_plus - dist_minus) / energy;
}

double FiniteTemperatureThermodynamicPotential()
{
    
    return 0;
}
