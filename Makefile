SHELL := /bin/bash # Use bash as shell
TARGET = heos

# List set for multirun
MULTIRUN_SETS = eNJL1 eNJL2 eNJL3 eNJL1m eNJL2m eNJL1OmegaRho1 eNJL2OmegaRho1 eNJL3SigmaRho1

.PHONY: all run graph tests tgraph clean

all:
	@echo "[Building ...]"
	@cd src; make
	@echo "[done.]"
run:
	@./$(TARGET) -d $(ARGS)
srun:
	@./$(TARGET) -s -d $(ARGS)
graph:
	@echo "[Plotting ...]"
	@cd output; \
	for dir in `echo */`; do \
		cd "$$dir"; \
		gnuplot gnuplot.gpi; \
		cd ..; \
	done
	@echo "[done.]"
multirun:
	@echo "[Running for multiple parameterizations ...]"
	@for key in $(MULTIRUN_SETS); do \
		./$(TARGET) -d -p "$$key" $(ARGS); \
		if [ -d multioutput/"$$key" ]; then rm -r multioutput/"$$key"; fi; \
		cp -r output multioutput/"$$key"; \
	done
	@echo "[done.]"
smultirun:
	@echo "[Running for multiple parameterizations ...]"
	@for key in $(MULTIRUN_SETS); do \
		./$(TARGET) -s -d -p "$$key" $(ARGS); \
		if [ -d multioutput/"$$key" ]; then rm -r multioutput/"$$key"; fi; \
		cp -r output multioutput/"$$key"; \
	done
	@echo "[done.]"
mgraph:
	@echo "[Plotting for multiple parameterizations ...]"
	@for dir in $(MULTIRUN_SETS); do \
		echo $$dir; \
		cd "multioutput/$$dir"; \
		for subdir in `echo */`; do \
			cd "$$subdir"; \
			gnuplot gnuplot.gpi; \
			cd ..; \
		done; \
		cd ../..; \
	done; \
	cd multioutput; gnuplot gnuplot.gpi
	@echo "[done.]"
tests:
	@echo "[Running tests ...]"
	@./$(TARGET) -a $(ARGS)
	@echo "[done.]"
tgraph:
	@echo "[Plotting tests ...]"
	@for dir in `echo tests/*/`; do cd "$$dir" && gnuplot gnuplot.gpi && cd ../..; done
	@echo "[done.]"
clean:
	@echo "[Cleaning ...]"
	@-rm -f $(TARGET)
	@cd src; make clean
	@find . -name "*.dat" -type f -delete
	@find . -name "*.log" -type f -delete
	@find . -name "*.png" -type f -delete
	@find . -name "*.tex" -type f -delete
	@cd multioutput; rm -rf $(MULTIRUN_SETS)
	@echo "[done.]"
