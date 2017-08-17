PathCORE
--------
Python 3 implementation of methods described in
`Chen et al.'s 2017 PathCORE paper <https://doi.org/10.1101/147645>`_ for
identifying pathway-pathway interactions using features constructed from
transcriptomic data.

This code has been tested on Python 3.5.
The documentation for the modules in the package can be
`accessed here <http://pathcore-demo.herokuapp.com/static/data/docs_pathcore/index.html>`_.

Installation
----------------
To install the current PyPI version (recommended), run::

    pip install PathCORE

For the latest GitHub version, run::

    pip install git+https://github.com/greenelab/PathCORE.git#egg=PathCORE

Examples
---------
We recommend that users of the PathCORE software begin by reviewing the
examples in the `PathCORE-analysis <https://github.com/greenelab/PathCORE-analysis>`_
repository. The analysis repository contains shell scripts and wrapper
analysis scripts that demonstrate how to run the methods in this package
on features constructed from a broad compendium according to the 
`workflow we describe in our paper <https://github.com/greenelab/PathCORE-analysis#the-pathcore-analysis-workflow>`_.

Specifically, `this Jupyter notebook <https://github.com/greenelab/PathCORE-analysis/blob/master/jupyter-notebooks/Supplemental_PAO1_FastICA_example.ipynb>`_
is a simple example of the workflow and a great place to start.

Package contents
----------------

=====================================
feature_pathway_overrepresentation.py
=====================================
The methods in this module are used to identify the pathways
overrepresented in features extracted from a transcriptomic dataset
of genes-by-samples. Features must preserve the genes in the dataset
and assign weights to these genes based on some distribution.
[`feature_pathway_overrepresentation documentation. <http://pathcore-demo.herokuapp.com/static/data/docs_pathcore/source/pathcore.html#module-pathcore.feature_pathway_overrepresentation>`_]

===========
network.py
===========
Contains the data structure ``CoNetwork`` that stores information
about the pathway co-occurrence network. The output from
a pathway enrichment analysis in ``feature_pathway_overrepresentation.py``
serves as input into the ``CoNetwork`` constructor.
[`CoNetwork documentation. <http://pathcore-demo.herokuapp.com/static/data/docs_pathcore/source/pathcore.html#module-pathcore.network>`_]

============================
network_permutation_test.py
============================
The methods in this module are used to filter the constructed
co-occurence network. We implement a permutation test that evaluates
and removes edges (pathway-pathway relationships) in the network
that cannot be distinguished from a null model of random associations.
The null model is created by generating *N* permutations of the network.
[`network_permutation_test documentation. <http://pathcore-demo.herokuapp.com/static/data/docs_pathcore/source/pathcore.html#module-pathcore.network_permutation_test>`_]

Acknowledgements
----------------
This work was supported by the Penn Institute for Bioinformatics
