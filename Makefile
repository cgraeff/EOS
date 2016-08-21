SHELL := /bin/bash # Use bash as shell
TARGET = heos

.PHONY: all run graph tests tgraph clean

all:
	cd src; make
run:
	./$(TARGET) -d $(ARGS)
graph:
	cd output; \
	for dir in `echo */`; do \
		cd "$$dir"; \
		gnuplot gnuplot.gpi; \
		cd ..; \
	done
multirun:
	for key in eNJL1 eNJL2 eNJL3 eNJL1m eNJL2m; do \
		./$(TARGET) -d -p "$$key"; \
		if [ -d multioutput/"$$key" ]; then rm -r multioutput/"$$key"; fi; \
		cp -r output multioutput/"$$key"; \
	done
mgraph:
	# in the following we exploit the fact that every parameterization starts with an 'e'
	# and run their graphics routines, then we run the general one.
	# Aditional sets may be listed in the 'echo' bellow with the form
	# multioutput/a_set multioutput/another multioutput/prefix*/
	# (Note that when using the wildcard *, the slash '/' at the end is
	# necessary 
	for dir in `echo multioutput/e*/`; do \
		cd "$$dir"; \
		for subdir in `echo */`; do \
			cd "$$subdir"; \
			gnuplot gnuplot.gpi; \
			cd ..; \
		done; \
		cd ../..; \
	done; \
	cd multioutput; gnuplot gnuplot.gpi
tests:
	./$(TARGET) -a $(ARGS)
tgraph:
	for dir in `echo tests/*/`; do cd "$$dir" && gnuplot gnuplot.gpi && cd ../..; done
clean:
	-rm -f $(TARGET)
	cd src; make clean
	find . -name "*.dat" -type f -delete
	find . -name "*.log" -type f -delete
	find . -name "*.png" -type f -delete
	find . -name "*.tex" -type f -delete
	-rm -rf multioutput/e*/
