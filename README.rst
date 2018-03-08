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

Technical Details
=================

The process is split into two parts: an offline slow pre-processing part
and a fast parallelized analyzer. The preprocessor marshalls the textual
data in GAF and OBO files into a binary form which can be quickly loaded
by the analyzer; the analyzer then, given a set of input genes, computes
hypergeometric probabilities for all terms in the ontology and outputs a
JSON structure detailing the per-term scores and gene associations.

OBO Preprocessor
----------------

- Input file: text file in standard OBO format

- Output file: binary OBO graph digest

- Usage::

    preprocessor obo <input_file.obo> <output_file>

Terms in the constructed graph are linked with identity and relationship
information as given in the input file. The attributes that determine
the relationships are:

- ``is_a: <term_id>``,

- ``relationship: part_of <term_id>``,

- ``relationship: regulates <term_id>``,

- ``relationship: positively_regulates <term_id>``,

- ``relationship: negatively_regulates <term_id>``.

For all of these, ``<term_id>`` denotes the *parent* of the term being
processed.

Relationships given with ``relationship: has_part <term_id>`` mean the
*reverse*: i.e. the term being processed is the parent and term
``<term_id>`` is its child.

For each term, a binary record is written to the output file. The record
contains the term's name, identifier, namespace identifier, lists of
parent and child terms, and scratch space for various metadata (this is
used exclusively by the processor and is here to alleviate any
allocation and setup overhead).

GAF Preprocessor
----------------

- Input file: text file in standard GAF format

- Output file: binary GAF link digest

- Usage::

    preprocessor gaf <input_file.gaf> <output_file>

GAF files simply associate genes with ontology terms. For each such
association, a link is added to the output file between the gene and the
ontology term, with the following information:

- the name of the gene is expected to be in the second column,

- the name of the term is expected to be in the fifth column,

- the link is only added if the fourth column is either:

  * empty, or

  * contains none of ``NOT``, ``colocalizes_with``, and
    ``contributes_to``.

Processor
---------

- Input files: binary OBO digest, binary GAF digest and a list of
  input genes

- Standard out: a result JSON with p-value-filtered terms and their
  scores and gene associations

- Usage::

    processor <p-val-thresh> <gene-count-thresh> <obo_file> <gaf_file> <input_file>

The OBO and GAF digests are memory-mapped and run through an index fixup
step to correct intra-file references and convert string indexes into
pointers.

The rest of the processing proceeds as follows:

1) OBO term ids are pushed into a caching array, which is used to speed
   up term lookup by id. Graph roots are found and referenced in a
   separate array.

2) Gene names from the input GAF are pushed into a caching array for
   faster lookup. Each gene is added to the gene list of the term it is
   linked to, and the process recurses up the graph: each of the term's
   parent nodes is also associated with the gene, up to (and including)
   one of the graph roots.

3) Input genes are 'intersected' with the gene lists of each term. That
   is: proceeding recursively from graph roots down through all
   child terms, intersections are computed between the set of input
   genes and the sets of genes associated through the GAF. Each input
   gene which has been associated with the term via the GAF, is added
   to an intersection list.

4) The following steps are done for each term. Because processing for
   each term is entirely self-contained and has no effect on anything
   outside the term's metadata (including the output buffer, storage
   for which was allocated beforehand), the loop can be embarrassingly
   parallelized.

   a) Hypergeometric probabilities are computed for every term in the
      ontology. Denote with ``k`` the size of gene set intersection,
      with ``n`` the total number of input genes, with ``K`` the
      number of GAF gene associations on this term (recursive!), and
      with ``N`` the number of all genes which had any associations at
      all with the ontology. The probability is then::

          e ^ (
              (log K + log (N-K) + log n + log (N-n)) -
              (log k + log (K-k) + log (n-k) + log (N-K-n+k) + log N)
          )

      Where ``e`` is Euler's number and ``log`` is the natural
      logarithm.  Computation is done with logarithms to avoid
      overflow problems. The equation in its original form::

          [K (N-K) n (N-n)] / [k (K-k) (n-k) (N-K-n+k) N]

   b) For each term which has at least one associated gene, a further
      'score' is computed. Using the same notation as above::

          (k / n) / (K / N)

   c) For each term, as much of the output textual representation is
      generated as possible. After this tep, the serial output code will
      only need to construct the surrounding JSON.

5) For each root term in the ontology graph, a tree of descendant terms
   is rendered and printed on standard out. Obsolete roots are ignored;
   children are filtered according to the p-value and gene count
   thresholds specified on the command line. Eligible terms (that is,
   terms whose p-value is *at most* the threshold value and the gene
   intersection size is *at least* the threshold value) are printed;
   ineligible terms *and the entire ontology subtree below them* are
   skipped.

6) After the main tree structure, for every dumped term its gene
   associations are also printed.

Output structure (random scores, fictitious OBO graph)::

    {
        "total_genes": 9,   (* GAF genes with any term associations *)
        "tree": {
            "BP": [   (* root; array are all immediate descendants *)
                {
                    "matched": 1,   (* intersection size *)
                                    (* size of gene_ids below *)
                    "pval": 0.025958803815505518,
                    "score": 1.515797517532959,
                    "total": 3,  (* all GAF associated genes *)
                    "term_id": "GMM:30.2.22",
                    "term_name": "signalling.receptor kinases.proline extensin like",
                    "gene_ids": [
                        "PGSC0003DMG400002675"
                        (* intersected genes: subset of input list *)
                    ],
                    "children": [
                        (* eligible descendants of this term *)
                        (* same structure as above, array of dicts *)
                    ]
                }
            ],
            "CC": [
                (* cellular component; same structure as above *)
            ],
            "MF": [
                (* molecular function; same structure as above *)
            ]
        },
        "gene_associations": {
            "GMM:30.2.22": [
                "PGSC0003DMG400000082",
                "MICRO.7728.C3"
            ],
            (* so on; every term appearing in "tree" above with *)
            (* all GAF-specified gene associations *)
        }
    }

Any of the root terms may be omitted from the output if it does not
appear in the OBO. The OBO is expected to have at most three
non-obsolete root terms, at most one each from the "biological_process",
"cellular_function" and "molecular_function" namespaces.

Binary Blob Format
------------------

The aim is to be as close to memory layout as possible; no marshalling
is done, except for the bare minimum required to re-link pointers and
get an actual memory layout. The processor does some further copying in
order to associate strings from the OBO blog to strings in the GAF blob.

A side effect of this is that the blob is very platform-dependent; blobs
produced on a given platform can only be used on platforms with the same
integer and pointer sizes and the same endianness.

GAF Blob
^^^^^^^^

Header:

- ``tsize [int]``: the size, in bytes, of the first section of the file:
  header + string nest. The header comprises three ``int`` fields; the
  strings in the nest are serialized as described above. Each string
  item is padded so that ``tsize`` up to the end of it is a multiple of
  ``sizeof(size_t)``. At the end, ``tsize`` is further rounded up to the
  next multiple of 4096, so that the next section of the file begins on
  a page boundary.

- ``gene_count [int]``: the number of gene names in the string nest.

- ``link_count [int]``: the number of gene<>term associations.

Following the header is the string nest which contains uniqued strings,
(all gene names, followed by term ids). Each string is laid out as
follows:

- ``len [size_t]``: the length, in bytes, of the string (excluding this
  metadata).

- ``hash [size_t]``: the precomputed hash of the string, used for quick
  equality checks.

- ``string [char[]]``: the actual string. There is extra padding at the
  end so that the next string begins at a multiple of ``size_t`` from
  the beginning of the file.

The link list begins at the next page (4K) boundary from the end of the
string nest. It is an array of ``int_link_t`` structs, with string
pointers changed into file offsets.

Immediately following the link list is an array of precomputed log
values (used for p-value computation in the processor), from 0 to the
gene count (inclusive), of type ``float_type`` (defined in
``processor.h``, usually ``double``).

OBO Blob
^^^^^^^^

Header:

- ``str_link_size [int]`` the length, in bytes, of the first section of
  the file, comprising the header (two ``int`` fields), the string nest
  and the link array. As in the GAF blob, each string is padded up to a
  multiple of ``size_t``, the nest as a whole is rounded up to page
  size.

- ``term_count [int]``: the number of terms in the ontology.

Following the header is the string nest, containing all term-related
strings (id and name), laid out the same as in the GAF blob.

Following the string nest and beginning on a page boundary from the
beginning of the file is a topology description of the ontology graph.
For each term, there are two arrays: the list of parent terms and the
list of descendant terms. The arrays have no sizing information, as that
is provided in the term structures that follow. Each term association is
a ``ptrdiff_t`` pointing to the offset, from the beginning of the file,
to the ``int_term_t`` structure of the term being associated to the
current one by the edge.

Following the link list, beginning on a page boundary, is the array of
term structures (``int_term_t``). The structures have some scratch space
that has nothing to do with the ontology, but is there to provide
quickly-accessible work memory for the processor (which does no extra
allocations for the terms, other than for the output buffer).

Following this is a buffer segment. This could be allocated by the
processor, but it is quicker to "precompute" it and just map it into
memory. There is a buffer for each term, containing enough space for 200
bytes (enough for the static parts of the term's JSON representation)
and the id and name strings.
