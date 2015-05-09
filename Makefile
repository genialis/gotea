OBJS=gafparser oboparser preprocessor termbuilder linkbuilder
GOMPOBJS=processor

CFLAGS=-Wall -Wextra -Wshadow -Winline -std=c++11 -lm -g
FASTFLAGS=-fopenmp

SHELL=/bin/bash

.PHONY: all clean distclean

all: preprocessor processor

%.d: %.cpp
	-@echo "  DEP    $@"
	@g++ -std=c++11 -MM $< > $@
-include $(OBJS:%=%.d) $(GOMPOBJS:%=%.d)

$(OBJS:%=%.o): %.o: %.cpp
	-@echo "  C++    $@"
	@g++ $(CFLAGS) -c $< -o $@

$(GOMPOBJS:%=%.o): %.o: %.cpp
	-@echo "  MPC++  $@"
	@g++ $(CFLAGS) $(FASTFLAGS) -c $< -o $@

preprocessor: preprocessor.o oboparser.o gafparser.o termbuilder.o linkbuilder.o
	-@echo "  LD     $@"
	@g++ $(CFLAGS) $^ -o $@

processor: processor.o
	-@echo "  LD     $@"
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
	-rm $(OBJS:%=%.o) $(GOMPOBJS:%=%.o) preprocessor processor formatter_test
distclean: clean
	-rm $(OBJS:%=%.d) $(GOMPOBJS:%=%.d)
