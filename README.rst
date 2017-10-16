PathCORE-T
----------
Python 3 implementation of methods described in
`Chen et al.'s 2017 PathCORE-T paper <https://doi.org/10.1101/147645>`_.

Note that this software was renamed from PathCORE to PathCORE-T in Oct 2017.
The T specifies that pathway co-occurrence relationships are identified using
features extracted from **transcriptomic** data. 
The module itself is still named `pathcore` to maintain backwards
compatibility for users of the original PathCORE software package. 

This code has been tested on Python 3.5.
The documentation for the modules in the package can be
`accessed here <http://pathcore-demo.herokuapp.com/static/data/docs_pathcore/index.html>`_.

Installation
----------------
To install the current PyPI version (recommended), run::

    pip install PathCORE-T

For the latest GitHub version, run::

    pip install git+https://github.com/greenelab/PathCORE-T.git#egg=PathCORE-T

Examples
---------
We recommend that users of the PathCORE-T software begin by reviewing the
examples in the `PathCORE-T-analysis <https://github.com/greenelab/PathCORE-T-analysis>`_
repository. The analysis repository contains shell scripts and wrapper
analysis scripts that demonstrate how to run the methods in this package
on features constructed from a broad compendium according to the 
`workflow we describe in our paper <https://github.com/greenelab/PathCORE-T-analysis#the-pathcore-analysis-workflow>`_.

Specifically, `this Jupyter notebook <https://github.com/greenelab/PathCORE-T-analysis/blob/master/jupyter-notebooks/Supplemental_PAO1_FastICA_example.ipynb>`_
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
