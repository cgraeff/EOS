//
//  Constants.h
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-02-23.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#ifndef Constants_h
#define Constants_h


#endif /* Constants_h */

#define CONST_HBAR_C 197.326
#define CTE_NUM_COLORS 1
#define CTE_NUM_FLAVORS 1

// Parameters for root finding in the solution of the Gap Equation
//		CONST_ABS_ERROR_GAP_EQ_SOLVING: An absolute error in MeV for the solution of the Gap Equation
//		CONST_REL_ERROR_GAP_EQ_SOLVING: A fraction of the value
// The root is accepted when
//		|x_lower - x_upper| < CONST_ABS_ERROR_GAP_EQ_SOLVING + CONST_REL_ERROR_GAP_EQ_SOLVING * MIN(|x_upper|, |x_lower|)
//
#define CONST_ABS_ERROR_GAP_EQ_SOLVING 0.05		
#define CONST_REL_ERROR_GAP_EQ_SOLVING 5.0E-4
#define CONST_MAX_ITERATIONS_GAP_EQ_SOLVING 1000
