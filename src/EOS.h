//
//  EOS.h
//  EOS
//
//  Created by Clebson Graeff on 2/17/16.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#ifndef EOS_h
#define EOS_h

typedef struct _gap_equation_input{
	
	double proton_fermi_momentum;
	double neutron_fermi_momentum;
	double barionic_density;
	
} gap_equation_input;

double SolveGapEquation(double barionic_density, double proton_fermi_momentum, double neutron_fermi_momentum);

double scalar_density(double mass, double fermi_momentum, double cutoff);

double ProtonChemicalPotential(double proton_fermi_momentum,
                               double scalar_density,
                               double mass,
                               double barionic_density,
                               double proton_density,
                               double neutron_density);

double NeutronChemicalPotential(double neutron_fermi_momentum,
                                double scalar_density,
                                double mass,
                                double barionic_density,
                                double proton_density,
                                double neutron_density);

double EnergyDensity(double pressure,
                     double proton_chemical_potential,
                     double neutron_chemical_potential,
                     double proton_density,
                     double neutron_density);

double Pressure(double termodynamic_potential);

double TermodynamicPotential(double scalar_density,
                             double barionic_density,
                             double proton_density,
                             double neutron_density,
                             double proton_chemical_potential,
                             double neutron_chemical_potential,
                             double kinectic_energy_density);

double KinecticEnergyDensity(double mass, double proton_fermi_momentum, double neutron_fermi_momentum);

double F0(double mass, double momentum);
double F2(double mass, double momentum);
double GapEquation(double mass, void * input);

#endif /* EOS_h */
