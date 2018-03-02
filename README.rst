======================================
Gene Ontology Term Enrichment Analysis
======================================

Functionality and Output
========================

Gene ontology term enrichment analysis (GOTEA) is used to test weather a
GO term is statistically enriched (over- or under- represented) for the
given set of genes. The input set of genes can be functionally profiled
by determining which GO terms appear more or less frequently than
expected by chance when examining the set of input genes annotated to
the GO term.

P-value
-------

A p-value is calculated to determine weather the GO term over- or under-
representation is significant or not.

For example, out of total 20,000 genes in the human genome, 440 map to
the GO term `induction of apoptosis`. Therefore, 2.2% (440 divided by
20,000) of the genes are involved in the induction of apoptosis. If a
set of input genes contains 500 genes, 11 genes (500 multiplied by 2.2%)
would be expected to be involved in induction of apoptosis. If more than
11 genes are involved in induction of apoptosis than this GO term would
be over-represented in the given gene set.

A p-value is calculated with hypergeometric test. It uses
`hypergeometric distribution
<https://en.wikipedia.org/wiki/Hypergeometric_distribution>`__
to measure the statistical significance weather among a set of
input genes (*n* draws) a proportion of genes that are in intersection
with genes assigned to specific GO term (*k* successes) taken from a
population with size of all genes form chosen Species (*N* size of
population) containing a proportion of genes assigned to GO term (*K*
successes).

Scores
------

Scores are showing weather there is over-representation or
under-representation of the input set of genes relative to the reference
list. The score equals 1 when the proportion of input genes involved in
GO term is as expected (same to proportion of genes assigned to GO term
among total number of genes).

**Scores more than 1** for GO term means that more genes are observed in
the input gene set than expected, so there is an **over-representation**
of input genes involved in this term.

**Score between 0 and 1** for GO term means fewer genes are observed in
the input gene set than expected, there is un **under-representation**
of input genes involved in this term.

High score number means that more genes are over-represented, closer the
value is to 0 higher is under-representation of input genes for selected
enriched term.

Score is a ratio between number of genes from input gene set that are
assigned to specific term (*k* successes) divided by the number of genes
in input gene set (*n* draws) and reference list which is a number of
genes in term (*K* successes) relatively to whole population of genes in
Species (*N* size of population).
