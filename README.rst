PathCORE
--------
Python 3 implementation of methods described in
Chen et al.'s 2017 PathCORE paper for identifying pathway-pathway
interactions using features constructed from transcriptomic data.

This code has been tested on Python 3.5.

Installation
----------------
To install the current PyPI version (recommended), run::

    pip install PathCORE

For the latest GitHub version, run::

    pip install git+https://github.com/greenelab/PathCORE.git#egg=PathCORE

Package contents
----------------

=====================================
feature_pathway_overrepresentation.py
=====================================
The methods in this module are used to identify the pathways
overrepresented in features extracted from a transcriptomic dataset
of genes-by-samples. Features must preserve the genes in the dataset
and assign weights to these genes based on some distribution.

===========
network.py
===========
Contains the data structure ``CoNetwork`` that stores information
about the pathway co-occurrence network. The output from
a pathway enrichment analysis in ``feature_pathway_overrepresentation.py``
serves as input into the ``CoNetwork`` constructor.

============================
network_permutation_test.py
============================
The methods in this module are used to filter the constructed
co-occurence network. We implement a permutation test that evaluates
and removes edges (pathway-pathway relationships) in the network
that cannot be distinguished from a null model of random associations.
The null model is created by generating _N_ permutations of the network.

Acknowledgements
----------------
This work was supported by the Penn Institute for Bioinformatics
