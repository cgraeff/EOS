//
//  CommandlineOptions.h
//  hadrons EOS
//
//  Created by Clebson Graeff on 2016-08-16.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#ifndef CommandlineOptions_h
#define CommandlineOptions_h

#include <stdbool.h>

typedef struct _options{
	// List options and flags that will be acessible
	// during the execution
	bool verbose;
    bool dirs;
  	bool tests;
    bool stars;
  	bool list_available_parameterizations;
	char * parameterization;
    double temp;
} Options;

extern Options options;

int CommandlineOptionsParse(int argc, char * argv[]);
void CommandlineOptionsPrintUsage();

#endif /* CommandlineOptions.h */
