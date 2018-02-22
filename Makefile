CFLAGS=-Wall -Wextra -Wshadow -Winline -std=c++11 -lm -O3
FASTFLAGS=-fopenmp -fwhole-program

SHELL=/bin/bash

.PHONY: all clean distclean

all: preprocessor processor

preprocessor: preprocessor.cpp oboparser.cpp gafparser.cpp termbuilder.cpp linkbuilder.cpp
	-@echo "  PROG   $@"
	@g++ $(CFLAGS) $^ -o $@

processor: processor.cpp
	-@echo "  PROG   $@"
	@g++ $(CFLAGS) $(FASTFLAGS) $^ -o $@

formatter_test: formatters.cpp
	@g++ $(CFLAGS) -DTEST=1 $^ -o $@

gaf.bin: preprocessor official.gaf
	-@echo "  GAF"
	-@time ./preprocessor gaf official.gaf gaf.bin
obo.bin: preprocessor gene_ontology_ext.obo
	-@echo "  OBO"
	-@time ./preprocessor obo gene_ontology_ext.obo obo.bin
test.out: processor obo.bin gaf.bin input_genes.txt
	-@echo "  PROCESSOR"
	-@time ./processor 0.2 2 obo.bin gaf.bin input_genes.txt > test.out

clean:
	-rm preprocessor processor formatter_test
distclean: clean
