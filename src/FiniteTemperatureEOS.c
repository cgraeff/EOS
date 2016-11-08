//
//  FiniteTemperatureEOS.c
//  hadrons EOS
//
//  Created by Clebson Graeff on 2016-11-08.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

int SolveMultiRootsFiniteTemp(double  proton_density,
                              double  neutron_density,
                              double *return_mass,
                              double *return_proton_renormalized_chemical_potential,
                              double *return_neutron_renormalized_chemical_potential)
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
        
        return 0;
    }

    return 0;
}

double FiniteTemperatureZeroMassSpecialCaseHelperFunction(double  x,
                                                          void   *par)
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

int FTMultiDimensionalRootFinderHelperFunction(const gsl_vector   *x,
                                               void               *p,
                                               gsl_vector         *return_values)
{
    // receive constant values
    multi_dim_root_params * params = (multi_dim_root_params *)p;
    double proton_barionic_density = params->proton_barionic_density;
    double neutron_barionic_density = params->neutron_barionic_density;
    double rho_3 = proton_barionic_density - neutron_barionic_density;
    
    // take care of mappings
   	const double mass = pow(gsl_vector_get(x, 0), 2.0);
    const double proton_renormalized_chemical_potential = pow(gsl_vector_get(x, 1), 2.0);
    const double neutron_renormalized_chemical_potential = pow(gsl_vector_get(x, 2), 2.0);
    
    // calculate things to be zeroed
    double proton_scalar_density = ScalarDensityAtFiniteTemperature(mass, proton_renormalized_chemical_potential, par);
    double neutron_scalar_density = ScalarDensityAtFiniteTemperature(mass, neutron_renormalized_chemical_potential, par);
    double total_scalar_density = proton_scalar_density + neutron_scalar_density;
    
    // the zeroed gap eq is probably the same as in the zero temp case
    // try to reuse by using a function
    double zeroed_gap = mass
                        - parameters.bare_mass
                        - 2.0 * parameters.G_S * total_scalar_density
                        + 2.0 * parameters.G_SV * pow(total_scalar_density, 2.0) * total_barionic_density
                        + 2.0 * parameters.G_SRHO * pow(rho_3, 2.0) * total_scalar_density;
    
    double proton_barionic_density_integral = BarionicDensityAtFiniteTemperature(proton_renormalized_chemical_potential,
                                                                                 temperature);
    double neutron_barionic_density_integral = BarionicDensityAtFiniteTemperature(neutron_renormalized_chemical_potential,
                                                                                  temperature);
    
    double zeroed_p_dens_gap = proton_barionic_density - proton_barionic_density_integral;
    double zeroed_n_dens_gap = neutron_barionic_density - neutron_barionic_density_integral;
    
    // Prepare return vector
   	gsl_vector_set(return_values, 0, zeroed_gap_eq);
   	gsl_vector_set(return_values, 1, zeroed_p_dens_gap);
    gsl_vector_set(return_values, 2, zeroed_n_dens_gap);
    
    return GSL_SUCCESS;
}

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
    
    return NUM_COLORS * NUM_FLAVORS * mass * integral / (pow(M_PI, 2.0) * pow(CTE_HC, 3.0));
}

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

double BarionicDensityAtFiniteTemperature(double mass,
                                          double chemical_potential,
                                          double upper_limit)
{
    fermi_dirac_distrib_integrand p;
    p.mass = mass;
    p.chemical_potential = renormalized_chemical_potential;
    p.temperature = parameters.temperature;
    
    gsl_function F;
    F.function = &BarionicDensityAtFiniteTemperatureIntegrand;
    F.params = &p;
    
    double integral = OnedimensionalIntegrator(&F, 0.0, upper_limit);
    
    return NUM_COLORS * NUM_FLAVORS * mass * integral / (pow(M_PI, 2.0) * pow(CTE_HC, 3.0));
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

double FiniteTemperatureThermodynamicPotential()
{
    
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
    double upper_limit = parameters.cutoff;
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
    
    return NUM_COLORS * NUM_FLAVORS * pow(CONST_HBAR_C, -3.0) * integral / pow(M_PI, 2.0);
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
