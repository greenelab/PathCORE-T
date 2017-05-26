"""Methods for permuting the pathways across constructed features of one
or more unsupervised models.
A permutation maintains the number of unique pathways that correspond to a
feature. If applicable, it also preserves the side (positive or negative)
in which pathway was overrepresented as well.
"""
import csv

from statsmodels.sandbox.stats.multicomp import multipletests


def aggregate_permuted_network(observed_networks):
    """This method handles the case where multiple observed networks are
    generated (e.g. from models produced by different random seed
    initializations of an unsupervised feature construction algorithm).
    We handle analysis of multiple networks by aggregating them; likewise,
    we require that the permutation test generates N aggregate permuted
    networks.

    Parameters
    -----------
    observed_networks : list(CoNetwork)
      the list of observed networks, generated from models produced by
      different random seed initializations of the same unsupervised feature
      construction algorithm

    Returns
    -----------
    CoNetwork, the aggregate permuted network created by generating
      a permutation for each individual observed network and then aggregating
      them into a single network
    """
    aggregate_permuted_network = None
    for pathcore_network in observed_networks:
        permuted_network = pathcore_network.permute_pathways_across_features()
        if not aggregate_permuted_network:
            aggregate_permuted_network = permuted_network
        else:
            aggregate_permuted_network.aggregate(permuted_network)
    return aggregate_permuted_network


def network_edges_permutation_test(observed_network,
                                   permuted_networks,
                                   alpha,
                                   n_networks=1,
                                   output_edges_to_file=None,
                                   output_network_to_file=None):
    """Given the observed network and N permutations of that network,
    determine the significance of each observed, weighted edge. Edges that
    are not distinguishable from the null model of random associations will
    be removed from the network file output. Additionally, weights of the
    significant edges are updated to an odds ratio that quantifies the
    "unexpectedness" of an observed pathway co-occurrence relationship
    based on the generated null model.

    Parameters
    -----------
    observed_network : CoNetwork
      the original observed network (this can be an aggregate of
      multiple networks, see CoNetwork's `aggregate` function in
      `network.py`)
    permuted_networks : list(CoNetwork)
      the list of N permuted networks generated from the original observed
      CoNetwork by calling its `permute_pathways_across_features` function
      N times and collecting the output in a list.
    alpha : float
      specify the threshold for significance testing
    n_networks : int (default=1)
      in the case where `observed_network` is an aggregate of multiple
      co-occurrence networks, `n_networks` can be greater than 1.
    output_edges_to_file : str|None (default=None)
      specify the filepath to write the edges that had q-values
      below `alpha`
    output_network_to_file : str|None (default=None)
      specify the filepath to write the final, filtered network

    Returns
    -----------
    CoNetwork, the observed network after the permutation test (filtered
    to only report significant edges, weighted by their odds ratios)
    """
    n_permutations = len(permuted_networks)
    n_features = observed_network.n_features

    edge_weight_distributions = _get_edge_weight_distributions(
        observed_network, permuted_networks, n_networks * n_features)

    pvalues = []
    edge_expected_weight = []
    for edge_id, edge_obj in observed_network.edges.items():
        observed_weight = edge_obj.weight
        weight_distribution = edge_weight_distributions[edge_id]
        eq_or_above_observed = sum(weight_distribution[observed_weight:])
        pvalue = eq_or_above_observed / float(n_permutations)
        pvalues.append(pvalue)
        edge_expected_weight.append((edge_id, _edge_expected_weight(
            weight_distribution[1:])))
    # fdr_bh: false discovery rate, Benjamini & Hochberg (1995, 2000)
    below_alpha, qvalues, _, _ = multipletests(
        pvalues, alpha=alpha, method="fdr_bh")

    significant_edge_indices = [index for index, passed in
                                enumerate(below_alpha) if passed]
    whitelist = set([edge_expected_weight[index][0] for index in
                     significant_edge_indices])

    observed_network.weight_by_edge_odds_ratios(
        edge_expected_weight, whitelist)

    n_significant_edges = len(significant_edge_indices)
    print("{0} edges are significant under the null distribution, generated "
          "from {1} permutations, for alpha = {2}.".format(n_significant_edges,
                                                           n_permutations,
                                                           alpha))
    if output_edges_to_file:
        _output_significant_edges_to_file(
            output_edges_to_file, observed_network, edge_expected_weight,
            pvalues, qvalues, significant_edge_indices)

    if output_network_to_file:
        observed_network_df = observed_network.to_dataframe(
            whitelist=whitelist)
        observed_network_df.to_csv(
            output_network_to_file, sep="\t", index=False)

    return observed_network


def _get_edge_weight_distributions(observed_network,
                                   permuted_networks,
                                   max_possible_edge_weight):
    """Helper function that returns the null distribution of each
    observed, weighted edge.
    """
    edge_weight_distributions = {}
    for edge in observed_network.edges.keys():
        edge_weight_distributions[edge] = [0] * max_possible_edge_weight

    for permuted_network in permuted_networks:
        vertex_conversion = observed_network.convert_pathway_mapping(
            permuted_network.pathways)
        for permuted_edge_id, edge in permuted_network.edges.items():
            edge_id = observed_network.remapped_edge(
                vertex_conversion, permuted_edge_id)
            if edge_id in edge_weight_distributions:
                edge_weight_distributions[edge_id][edge.weight] += 1
    return edge_weight_distributions


def _edge_expected_weight(edge_weight_distribution):
    """Helper function that computes the expected weight of an edge from
    the given distribution
    """
    edge_exists_in_n_permutations = sum(edge_weight_distribution)
    if not edge_exists_in_n_permutations:
        return 1
    sum_weight_counts = 0
    for weight_index, n_permutations in enumerate(edge_weight_distribution):
        sum_weight_counts += (weight_index + 1) * n_permutations

    return sum_weight_counts / float(edge_exists_in_n_permutations)


def _output_significant_edges_to_file(filepath,
                                      observed_network,
                                      edge_expected_weight,
                                      pvalues,
                                      qvalues,
                                      significant_edge_indices):
    """Helper function to report the results of the permutation test. Writes
    the table of significant edges and their p-value & q-values to the
    specified filepath.
    """
    with open(filepath, "w") as fp:
        writer = csv.writer(fp, delimiter="\t")
        writer.writerow(["pw0", "pw1",
                         "pvalue", "qvalue", "odds_ratio"])
        for index in significant_edge_indices:
            edge_id, _ = edge_expected_weight[index]
            pathway0, pathway1 = observed_network.get_edge_pathways(
                edge_id)
            odds_ratio = observed_network.edges[edge_id].weight
            output = (pathway0, pathway1,
                      pvalues[index], qvalues[index], odds_ratio)
            writer.writerow(output)
