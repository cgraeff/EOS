//
//  RootFinding.c
//  hadrons EOS
//
//  Created by Clebson Graeff on 2016-09-20.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include "RootFinding.h"

int MultidimensionalRootFinder(const int                dimension,
                               gsl_multiroot_function  *f,
                               gsl_vector              *initial_guess,
                               double                   abs_error,
                               double                   rel_error,
                               int                      max_iterations,
                               gsl_vector              *results)
{
    int status;
    size_t iter = 0;

    const gsl_multiroot_fsolver_type * solver_type = gsl_multiroot_fsolver_broyden;
    gsl_multiroot_fsolver * solver = gsl_multiroot_fsolver_alloc(solver_type,
                                                                 dimension);

    gsl_multiroot_fsolver_set(solver, f, initial_guess);

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
                                          abs_error,
                                          rel_error);

    } while (status == GSL_CONTINUE
             && iter < max_iterations);

    // Save results in return variables
    gsl_vector_memcpy(results, gsl_multiroot_fsolver_root(solver));

    // Free vectors
    gsl_multiroot_fsolver_free(solver);

    return 0;
}


int UnidimensionalRootFinder(gsl_function	*F,
                             double          lower_bound,
                             double          upper_bound,
                             double          abs_error,
                             double          rel_error,
                             int             max_iterations,
                             double         *return_result)
{

    const gsl_root_fsolver_type * T	= gsl_root_fsolver_bisection;
    gsl_root_fsolver * s = gsl_root_fsolver_alloc(T);

    // Test if the limits straddle the root,
    // if they don't, we will return -1.
    if (GSL_SIGN(GSL_FN_EVAL(F, lower_bound)) == GSL_SIGN(GSL_FN_EVAL(F, upper_bound)))
        return -1;

    gsl_root_fsolver_set(s, F, lower_bound, upper_bound);

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
                                                   abs_error,
                                                   rel_error)
            && i <= max_iterations);

    double result = gsl_root_fsolver_root(s);

    void gsl_root_fsolver_free(gsl_root_fsolver * S);

    *return_result = result;

    return 0;
}
