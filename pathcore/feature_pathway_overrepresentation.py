"""Determine the pathways overrepresented in a constructed feature.
"""
from crosstalk_correction import crosstalk_correction
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests


# handle the floating point errors thrown in the Fisher's exact test
# manually.
np.seterr(all="raise")


def pathway_enrichment_with_overlap_correction(feature_weight_vector,
                                               pathway_definitions,
                                               gene_signature_definition,
                                               alpha=0.05,
                                               correct_all_genes=True,
                                               metadata=False):
    """Identify pathways overrepresented in a constructed feature
    according to a user-specified criterion (see `gene_signature_definitions`
    in the parameters list below) for identifying the feature's gene
    signature. Donato et al.'s (2013) algorithm for pathway crosstalk
    correction (removes gene overlap between pathway definitions) is applied to
    the pathway definitions before overrepresentation analysis.
    Here, we refer to it as *overlap-correction*.

    Parameters
    -----------
    feature_weight_vector : pandas.Series(float), shape = n
      A vector containing gene weights
    pathway_definitions : dict(str -> set(str))
      A pathway (key) is defined by a set of genes (value).
    gene_signature_definition : functools.partial callable,
                                returns (set(), set())
      Accepts the `feature_weight_vector` as input. Provide a function to
      distinguish positive and/or negative gene signatures.
      Both a positive & a negative signature may be appropriate if the
      feature's gene weight distribution spans positive and negative values.
      If this is not the case, specify a single gene signature by returning
      one of the sides as an empty set.
    alpha : float (default=0.05)
      Significance level for pathway enrichment.
    correct_all_genes : bool (default=True)
      The overlap correction procedure is applied independently to both
      the gene signature (union of positive and negative signature when
      applicable) _and_ the genes outside of the signature (termed
      *remaining* genes).
      If not `correct_all_genes`, overlap correction is not applied to
      the set of remaining genes.
    metadata : bool (default=False)
      Gather information to store in a MongoDB-backed Flask web application.
      Users can interact with the PathCORE-T-produced network and analyze the
      genes underlying a pair of pathways linked in the network.

    Returns
    -----------
    tup([pandas.DataFrame|None], dict())
    tup[0] : pandas.DataFrame: dataframe of significant pathways
           | None if the gene signature does not contain any genes in the
                  pathway definitions
    tup[1] : if `metadata`:
             {"positive signature": <set() positive gene signature>,
              "negative signature": <set() negative gene signature>,
              "pathway definitions": <dict(str -> set())
                overlap-corrected definitions--only signature genes>}
             else: {}
    """
    genes_in_pathway_definitions = set.union(*pathway_definitions.values())
    positive_gene_signature, negative_gene_signature = \
        gene_signature_definition(feature_weight_vector)
    gene_signature = ((positive_gene_signature | negative_gene_signature) &
                      genes_in_pathway_definitions)
    if not gene_signature:
        return (None, {})

    additional_information = {}
    n_genes = len(feature_weight_vector)

    if metadata:
        additional_information["positive_signature"] = (
          positive_gene_signature & gene_signature)
        additional_information["negative_signature"] = (
          negative_gene_signature & gene_signature)
    corrected_pathway_definitions = crosstalk_correction(
        pathway_definitions,
        gene_set=gene_signature,
        all_genes=correct_all_genes)

    if metadata:
        collect_signature_pathway_definitions = {}
        for pathway, (signature_defn, remaining_defn) in \
                corrected_pathway_definitions.items():
            collect_signature_pathway_definitions[pathway] = signature_defn
        additional_information["pathway_definitions"] = \
            collect_signature_pathway_definitions

    pathway_positive_series = single_side_pathway_enrichment(
        corrected_pathway_definitions, positive_gene_signature, n_genes)
    pathway_negative_series = single_side_pathway_enrichment(
        corrected_pathway_definitions, negative_gene_signature, n_genes)

    pvalue_information = pathway_positive_series.append(
        pathway_negative_series)
    side_information = _pathway_side_information(
        pathway_positive_series, pathway_negative_series,
        pvalue_information.index)
    significant_pathways = _significant_pathways_dataframe(
        pvalue_information, side_information, alpha)
    return significant_pathways, additional_information


def _pathway_side_information(pathway_positive_series,
                              pathway_negative_series,
                              index):
    """Create the pandas.Series containing the side labels that correspond
    to each pathway, based on the user-specified gene signature definition.
    """
    positive_series_label = pd.Series(["pos"] * len(pathway_positive_series))
    negative_series_label = pd.Series(["neg"] * len(pathway_negative_series))
    side_information = positive_series_label.append(
        negative_series_label)
    side_information.index = index
    side_information.name = "side"
    return side_information


def _significant_pathways_dataframe(pvalue_information,
                                    side_information,
                                    alpha):
    """Create the significant pathways pandas.DataFrame.
    Given the p-values corresponding to each pathway in a feature,
    apply the FDR correction for multiple testing and remove those that
    do not have a q-value of less than `alpha`.
    """
    significant_pathways = pd.concat(
        [pvalue_information, side_information], axis=1)
    # fdr_bh: false discovery rate, Benjamini & Hochberg (1995, 2000)
    below_alpha, qvalues, _, _ = multipletests(
        significant_pathways["p-value"], alpha=alpha, method="fdr_bh")
    below_alpha = pd.Series(
        below_alpha, index=pvalue_information.index, name="pass")
    qvalues = pd.Series(
        qvalues, index=pvalue_information.index, name="q-value")
    significant_pathways = pd.concat(
        [significant_pathways, below_alpha, qvalues], axis=1)
    significant_pathways = significant_pathways[significant_pathways["pass"]]
    significant_pathways.drop("pass", axis=1, inplace=True)
    significant_pathways.loc[:, "pathway"] = significant_pathways.index
    return significant_pathways


def single_side_pathway_enrichment(pathway_definitions,
                                   gene_signature,
                                   n_genes):
    """Identify overrepresented pathways using the Fisher's exact test for
    significance on a given pathway definition and gene signature.
    (FDR correction for multiple testing is applied in
    `_significant_pathways_dataframe`).

    Parameters
    -----------
    pathway_definitions : dict(str -> set(str))
      Pathway definitions, *post*-overlap-correction if this function
      is called from `pathway_enrichment_with_overlap_correction`.
      A pathway (key) is defined by a set of genes (value).
    gene_signature : set(str)
      The set of genes we consider to be enriched in a feature.
    n_genes : int
      The total number of genes for which we have assigned weights in the
      features of an unsupervised model.

    Returns
    -----------
    pandas.Series, for each pathway, the p-value from applying the Fisher's
      exact test.
    """
    if not gene_signature:
        return pd.Series(name="p-value")
    pvalues_list = []
    for pathway, definition in pathway_definitions.items():
        if isinstance(definition, tuple):
            definition = set.union(*definition)

        both_definition_and_signature = len(definition & gene_signature)
        in_definition_not_signature = (len(definition) -
                                       both_definition_and_signature)
        in_signature_not_definition = (len(gene_signature) -
                                       both_definition_and_signature)
        neither_definition_nor_signature = (n_genes -
                                            both_definition_and_signature -
                                            in_definition_not_signature -
                                            in_signature_not_definition)
        contingency_table = np.array(
            [[both_definition_and_signature, in_signature_not_definition],
             [in_definition_not_signature, neither_definition_nor_signature]])
        try:
            _, pvalue = stats.fisher_exact(
                contingency_table, alternative="greater")
            pvalues_list.append(pvalue)
        # FPE can occur when `neither_definition_nor_signature` is very
        # large and `both_definition_and_signature` is very small (near zero)
        except FloatingPointError:
            pvalues_list.append(1.0)
    pvalues_series = pd.Series(
        pvalues_list, index=pathway_definitions.keys(), name="p-value")
    return pvalues_series


def pathway_enrichment_no_overlap_correction(feature_weight_vector,
                                             pathway_definitions,
                                             gene_signature_definition,
                                             alpha=0.05,
                                             metadata=False):
    """Identify pathways overrepresented in a constructed feature
    according to a user-specified criterion (see `gene_signature_definitions`
    in the parameters list below) for identifying the feature's gene
    signature. Overlap-correction is not applied to the pathway definitions.

    Parameters
    -----------
    feature_weight_vector : pandas.Series(float), shape = n
      A vector containing gene weights
    pathway_definitions : dict(str -> set(str))
      A pathway (key) is defined by a set of genes (value).
    gene_signature_definition : functools.partial callable,
                                returns (set(), set())
      Accepts the `feature_weight_vector` as input. Provide a function to
      distinguish positive and/or negative gene signatures.
      Both a positive & negative signature may be appropriate if the feature's
      gene weight distribution spans positive and negative values. If this
      is not the case, a user can just specify a single gene signature by
      returning one or the other as an empty set.
    alpha : float (default=0.05)
      Significance level for pathway enrichment.
    metadata : bool (default=False)
      Return information about the gene signature(s)

    Returns
    -----------
    tup([pandas.DataFrame|None], dict())
    tup[0] : pandas.DataFrame: dataframe of significant pathways
           | None if the gene signature does not contain any genes in the
                  pathway definitions
    tup[1] : if `metadata`:
             {"positive signature": <set() positive gene signature>,
              "negative signature": <set() negative gene signature>,
              "pathway definitions": <dict(str -> set()) the pathway genes
               that are in the gene signature(s)>}
             else: {}
    """
    genes_in_pathway_definitions = set.union(*pathway_definitions.values())

    positive_gene_signature, negative_gene_signature = \
        gene_signature_definition(feature_weight_vector)
    gene_signature = ((positive_gene_signature | negative_gene_signature) &
                      genes_in_pathway_definitions)

    if not gene_signature:
        return (None, {})

    additional_information = {}
    n_genes = len(feature_weight_vector)

    if metadata:
        additional_information["positive_signature"] = positive_gene_signature
        additional_information["negative_signature"] = negative_gene_signature

        collect_signature_pathway_definitions = {}
        for pathway, definition in pathway_definitions.items():
            signature_definition = gene_signature & definition
            if signature_definition:
                collect_signature_pathway_definitions[pathway] = (
                    signature_definition)
        additional_information["pathway_definitions"] = (
            collect_signature_pathway_definitions)

    pathway_positive_series = single_side_pathway_enrichment(
        pathway_definitions, positive_gene_signature, n_genes)
    pathway_negative_series = single_side_pathway_enrichment(
        pathway_definitions, negative_gene_signature, n_genes)
    pvalue_information = pathway_positive_series.append(
        pathway_negative_series)
    side_information = _pathway_side_information(
        pathway_positive_series, pathway_negative_series,
        pvalue_information.index)
    significant_pathways = _significant_pathways_dataframe(
        pvalue_information, side_information, alpha)
    return significant_pathways, additional_information
