all:
	cd src; make
run:
	mkdir -p data/gap
	./eos
graph:
	cd data; make
clean:
	-rm -f eos
	cd src; make clean
	cd data; make clean
