//
//  FiniteTemperatureEOS.h
//  hadrons EOS
//
//  Created by Clebson Graeff on 2016-11-08.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

double ScalarDensityAtFiniteTemperature(double mass,
                                        double renormalized_chemical_potential,
                                        double upper_limit);

int SolveMultiRootsFiniteTemp(double  proton_density,
                              double  neutron_density,
                              double *return_mass,
                              double *return_proton_renormalized_chemical_potential,
                              double *return_neutron_renormalized_chemical_potential);

double FiniteTemperatureProtonChemicalPotential(double renormalized_proton_chemical_potential,
                                                double mass,
                                                double total_scalar_density,
                                                double proton_barionic_density,
                                                double neutron_barionic_density);

double FiniteTemperatureNeutronChemicalPotential(double renormalized_neutron_chemical_potential,
                                                 double mass,
                                                 double total_scalar_density,
                                                 double proton_barionic_density,
                                                 double neutron_barionic_density);
double FiniteTemperatureKinecticEnergyDensity();
