# Hadrons EOS

This program calculates the Equations of State for quarks according the SU(2)
version of the extended NJL model (eNJL). As references we use (Pais 2016).

1. **(Pais 2016)** Helena Pais, Débora P. Menezes, and Constança Providência,
   *Neutron stars: From the inner crust to the core with the (extended)
   Nambu–Jona-Lasinio model*, Phys. Rev. C 93, 065805 – Published 8 June 2016,
   [DOI](http://dx.doi.org/10.1103/PhysRevC.93.065805)

## Obtaining this program
The latest version of this code is available at [github.com/cgraeff/quarks_EOS](https://github.com/cgraeff/quarks_EOS).
A copy can be easily obtained if your system have `git` installed, just issue the following
command in a terminal:
 ```
   git clone https://github.com/cgraeff/quarks_EOS.git
 ```

## Requisites

To build and run this code `make`, a C compiler (default is `gcc`) and the
GSL (GNU Scientific Library) must be installed on the system. For plotting
graphics, `gnuplot` is also required.

On Linux the installation varies from distribution to distribution, but generally
there are two packages, one for regular use and one for developing.
The development version is the one needed for compilation of this code.

On OSX, the location of GSL headers and libraries is assumed to be
`/usr/local/include` and `/usr/local/lib`. This is the default if GSL is
installed from Homebrew. To change, edit the variables at `src/Makefile`.

## Build and run instructions

Basic build (generate `heos`executable):
* Build with `make` in the top dir;

Running:
* When executed, `heos` will calculate the equations of state with the default
  parameterization;
* The following options are available:
 * `-p par`: uses `par` parameters set;
 * `-t val`: uses `val` for temperature value. Must be a `double` value;
 * `-y val`: uses `val` for proton fraction. Must be a `double` value.
 * `-l`: list available parameters set;
 * `-q`: quiet (supress information written to standard out);
 * `-d`: write results using a dir structure;
 * `-a`: run tests. Automatically sets `-d`;
 * `-s`: run with charge neutrality and beta equilibrium;
 * `-u`: prints usage;
 * `-h`: also prints usage;

For easier running and testing, the following commands are provided:
* Run with `make run` (implies `-d`);
* Remove product files with `make clean`;
* Arguments may be passed with `make ARGS="-p Set" run`,
  where `-p Set` stand for example arguments;
* Run for many parameters sets at once using `make multirun` or `make ARGS="-t 10" multirun`,
where `-t 10` stands for example arguments. The sets must be listed in the `MULTIRUN_SETS`
variable in the `Makefile` at the root dir;
* Run tests with `make tests` (it is a shortcut to `make ARGS="-a" run` with
  the default parameterization);
* As a shorthand to `make ARGS="-s" run` or `make ARGS="-s" multirun`,
  `make srun` and `make smultirun` may be used.

When running on the default tree (that is, on the cloned or downloade dir), the
results can be plotted with
* Plot results with `make graph`
* Plot results for multirun with `make mgraph`. This also plots the `multioutput/gnuplot.gpi`
which may access files for multiple parameters sets and is useful for making plots comparing
results for different sets;
* Plot tests with `make tgraph`;

Two files for each plot will be created, a `png` file and a `tex` file. The second
is a version created using the `tiks` terminal for `gnuplot` and can be used
directly into `(pdf)latex`. The advantage is a much cleaner plot and proper
`latex` equations and text in labels.

## Code structure

### Parameters

Parameters sets must be declared in `ParametersSetup()` in `Parameters.c`.
Each parameters set is cloned from a reference set with
`NewCopyOfParametersSetFromTemplate()` then the variables which differ from the
template are set. Each parameters set must have a unique identifier. Finally,
the set must be appendend to the list using `AppendParametersSetToList()`.

To choose a parameters set, the function `SetParametersSet()` must be used.
This is used before the calculation of the main results and if no parameters set
is explicitly requested in the command line, the first set declared is used.
Each test in `Tests()` (in `Tests.c`) should declare a set.

Once set using `SetParametersSet()`, use `parameters.a_parameter` to access the
parameters values, where `a_parameter` is one of the parameters declared in the
struct in `Parameters.h`.

### Commandline options

The commandline options set in `CommandlineOptions.*` are globally available
using `options.an_option`.

### Tests

Tests should be declared in the function `RunTests()` in `Tests.c`. This
function is executed when the executable is called with `-a` option.

### Files and paths

Paths can be easily set with `SetFilePath()`. After that, all files created with
`OpenFile()` will be created at that path. After a section of code that uses
`SetFilePath()`, the path should be reset to the working dir with
`SetFilePath(NULL)`. This is important as other sections of code may expect to
write in the working dir while using `OpenFile()` and this is accomplished with
a `NULL` path.
