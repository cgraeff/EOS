SHELL := /bin/bash # Use bash as shell
TARGET = heos

# List set for multirun
MULTIRUN_SETS = eNJL1 eNJL2 eNJL3 eNJL1m eNJL2m

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
	for key in $(MULTIRUN_SETS); do \
		./$(TARGET) -d -p "$$key"; \
		if [ -d multioutput/"$$key" ]; then rm -r multioutput/"$$key"; fi; \
		cp -r output multioutput/"$$key"; \
	done
mgraph:
	for dir in $(MULTIRUN_SETS); do \
		cd "multioutput/$$dir"; \
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
	cd multioutput; rm -rf $(MULTIRUN_SETS)
