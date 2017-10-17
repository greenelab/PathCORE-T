"""CoNetwork is the data structure that supports the pathway co-occurrence
network produced by the PathCORE-T software.

The classes Vertex and Edge are used in CoNetwork.
"""
import itertools
import random

import pandas as pd


class Vertex:
    """A pathway is represented as a Vertex object"""

    def __init__(self, vertex_id):
        self._id = vertex_id
        self.edges = {}

    def degree(self):
        return len(self.edges)

    def get_adjacent_vertex_ids(self):
        return self.edges.keys()


class Edge:
    """A co-occurrence relationship between two pathways is
    represented by an Edge object"""

    def __init__(self, vertex0_id, vertex1_id, which_features=[]):
        """
        vertex0_id : int
        vertex1_id : int
        which_features : list (default=[])
          The features in which this pathway-pathway relationship appears.
        """
        self.edge = (vertex0_id, vertex1_id)
        # If `which_features` is initialized as the empty list, assume this is
        # an unweighted network and set the weight to 1.
        self.weight = 1 if not which_features else len(which_features)
        self.which_features = which_features

        # This flag is set to a bool value after a permutation test
        # has been applied to the network. See `network_permutation_test.py`.
        self.is_significant = None

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        else:
            return sorted(self.edge) == sorted(other.edge)

    def __ne__(self, other):
        return not self.__eq__(other)

    def features_to_string(self):
        return " ".join(map(str, self.which_features))

    def connected_to(self, vertex_id):
        """
        Parameters
        -----------
        vertex_id : int
          Get what `vertex_id` is connected to.

        Returns
        -----------
        int|None, the vertex id connected to the
          input `vertex_id` in this edge, as long as `vertex_id` is
          one of the vertices connected by this edge.
        """
        if vertex_id not in self.edge:
            return None
        connected_to = (self.edge[1] if self.edge[0] == vertex_id
                        else self.edge[0])
        return connected_to


class CoNetwork:

    def __init__(self,
                 model_n_features,
                 significant_pathways=None,
                 permutation_max_iters=30000):
        """
        model_n_features : int
          the number of features the network is built from.
        significant_pathways : str|pandas.DataFrame|None
        - str: the path to a significant pathways file, which is then
            loaded into a pandas.DataFrame. Must be tab-separated
            and contain columns "feature," "pathway," and "side."
        - pandas.DataFrame: the PathCORE-T analysis produces a significant
            pathways pandas.DataFrame that a user can directly pass
            to the CoNetwork constructor as well.
        - None: Use the method `read_network_file` to load a PathCORE-T
            network file as a CoNetwork object after a CoNetwork object
            has been initialized with
            CoNetwork(model_n_features, significant_pathways=None, ...).
        permutation_max_iters : int (default=30000)
          Sets an upper bound to the number of times we can randomly assign
          a pathway to a feature during the side-preserving model permutation.
        """
        self.permutation_max_iters = permutation_max_iters

        self.vertices = {}  # vertex id -> vertex object
        self.edges = {}  # (vertex 0 id, vertex 1 id) -> edge object

        self.n_features = model_n_features
        self.features = set()

        self.n_pathways = 0
        self.pathways = {}  # pathway -> vertex id

        if isinstance(significant_pathways, str):
            self.feature_pathway_df = _load_significant_pathways_file(
                significant_pathways)
            self._construct_from_dataframe()
        elif isinstance(significant_pathways, pd.DataFrame):
            self.feature_pathway_df = significant_pathways
            self._construct_from_dataframe()
        elif isinstance(significant_pathways, dict):
            self._construct_from_permutation(significant_pathways)

    def __getitem__(self, search_id):
        for pathway, vertex_id in self.pathways.items():
            if search_id == vertex_id:
                return pathway
        return None

    def read_network_file(self, path_to_network_file):
        """
        Parameters
        -----------
        path_to_network_file : str
          Expects a network file with columns "pw0" and "pw1." A "features"
          column that specifies the features where the (pw0, pw1) edge is
          present will assign a weight to the edge, though it is not required
          (edge will have weight 1 if no "features" column exists).
              "features" format: space-separated feature numbers,
              e.g. 0.0 1.0 2.0
        """
        network_df = pd.read_table(path_to_network_file)
        network_edges = {}
        for _, row in network_df.iterrows():
            vertex0_id = self.add_pathway(row["pw0"])
            vertex1_id = self.add_pathway(row["pw1"])
            edge_id = self.edge_tuple(vertex0_id, vertex1_id)
            if "features" in row:
                network_edges[edge_id] = \
                  [int(float(f)) for f in row["features"].split(" ")]
            else:
                network_edges[edge_id] = []
        self._augment_network(network_edges)

    def _construct_from_dataframe(self):
        """Build the network using the significant pathways dataframe
        by identifying pairwise direct (same-side) relationships.
        """
        direct_edges = {}
        feature_side_grouping = self.feature_pathway_df.groupby(
            ["feature", "side"])
        for (feature, _), group in feature_side_grouping:
            co_occurring_pathways = group["pathway"].tolist()
            pairings = list(itertools.combinations(co_occurring_pathways, 2))
            if not pairings:
                continue
            self.features.add(feature)
            for pathway0, pathway1 in pairings:
                vertex0_id = self.add_pathway(pathway0)
                vertex1_id = self.add_pathway(pathway1)
                new_edge = self.edge_tuple(vertex0_id, vertex1_id)
                if new_edge not in direct_edges:
                    direct_edges[new_edge] = []
                direct_edges[new_edge].append(feature)
        self._augment_network(direct_edges)

    def _construct_from_permutation(self, significant_pathways):
        """Build the network from a dictionary of (side -> tuple lists),
        where the side is specified as "pos" and/or "neg" (from the feature
        gene signature(s)) and mapped to a tuple list of [(pathway, feature)].
        Used during the PathCORE-T permutation test by applying the method
        `permute_pathways_across_features` to an existing CoNetwork.
        """
        for side, pathway_feature_tuples in significant_pathways.items():
            feature_pathway_dict = self._collect_pathways_by_feature(
                pathway_feature_tuples)
            self._edges_from_permutation(feature_pathway_dict)

    def permute_pathways_across_features(self):
        """Return a permutation of the network. Requires that the
        significant pathways file has been specified during CoNetwork
        initialization (or set `self.feature_pathway_df` to the
        pandas.DataFrame afterwards).
        """
        side_labels = self.feature_pathway_df.side.unique()
        pathway_features = {}
        for side in side_labels:
            pathway_features[side] = []

        feature_grouping = self.feature_pathway_df.groupby(
            ["side", "feature"])

        for (side, feature), group in feature_grouping:
            pathways = group["pathway"].tolist()
            for pathway in pathways:
                pathway_features[side].append((pathway, feature))

        permuted_network = {}
        for side in side_labels:
            pathway_side_permutation = _pathway_feature_permutation(
                pathway_features[side], self.permutation_max_iters)
            while pathway_side_permutation is None:
                pathway_side_permutation = _pathway_feature_permutation(
                    pathway_features[side], self.permutation_max_iters)
            assert _permutation_correctness(
                pathway_side_permutation, pathway_features[side]), \
                ("Permutation on side {0} did not preserve the "
                 "expected invariants").format(side)
            permuted_network[side] = pathway_side_permutation
        return CoNetwork(self.n_features,
                         significant_pathways=permuted_network)

    def weight_by_edge_odds_ratios(self,
                                   edges_expected_weight,
                                   flag_as_significant):
        """Applied during the permutation test. Update the edges in the
        network to be weighted by their odds ratios. The odds ratio measures
        how unexpected the observed edge weight is based on the expected
        weight.

        Parameters
        -----------
        edges_expected_weight : list(tup(int, int), float)
          A tuple list of (edge id, edge expected weight) generated from the
          permutation test step.
        flag_as_significant : [set|list](tup(int, int))
          A set or list of edge ids that are considered significant against
          the null model of random associations generated in the permutation
          test
        """
        for edge_id, expected_weight in edges_expected_weight:
            edge_obj = self.edges[edge_id]
            edge_obj.weight /= expected_weight
            if edge_id in flag_as_significant:
                edge_obj.significant = True
            else:
                edge_obj.significant = False

    def aggregate(self, merge):
        """Combine this network with another network. The aggregation step
        takes the union of the edges in the two networks, where we take the
        sum of weights for edges common to both networks.

        Parameters
        -----------
        merge : CoNetwork
          the CoNetwork object being merged into the current network.
        """
        self.features = set()
        self.n_features += merge.n_features

        vertex_id_conversion = self.convert_pathway_mapping(merge.pathways)
        for edge_id, edge in merge.edges.items():
            edge_key = self.remapped_edge(
                vertex_id_conversion, edge_id)
            if edge_key in self.edges:
                if self.edges[edge_key].which_features:
                    self.edges[edge_key].which_features = []
                self.edges[edge_key].weight += edge.weight
            else:
                vertex0_id, vertex1_id = edge_key
                new_edge_obj = Edge(vertex0_id, vertex1_id, [])
                new_edge_obj.weight = edge.weight
                self.edges[edge_key] = new_edge_obj
                self._add_edge_to_vertex(vertex0_id, new_edge_obj)
                self._add_edge_to_vertex(vertex1_id, new_edge_obj)

    def convert_pathway_mapping(self, other_pathway_mapping):
        """Used to convert the pathway-to-vertex id mapping in one CoNetwork
        to the one used in the current CoNetwork (`self`).
        The following tasks are carried out in the remapping:
        (1) If `self.pathways` contains the pathway to be merged, map the
            vertex id in `other_pathway_mapping` to the vertex id
            in `self.pathways`.
        (2) If not, create a vertex in `self` and then add the key-value pair
            to `self.pathways` accordingly.

        Parameters
        -----------
        other_pathway_mapping : dict(str -> int)
          the `pathways` field in the CoNetwork class. This is a
          (pathway -> vertex id) map.

        Returns
        -----------
        dict(int -> int), the (other vertex id -> `self` vertex id)
          conversion map
        """
        vertex_id_conversion = {}
        for pathway, vertex_id in other_pathway_mapping.items():
            if pathway in self.pathways:
                vertex_id_conversion[vertex_id] = self.pathways[pathway]
            else:
                self_vertex_id = self.add_pathway(pathway)
                self.vertices[self_vertex_id] = Vertex(self_vertex_id)
                vertex_id_conversion[vertex_id] = self_vertex_id
        return vertex_id_conversion

    def remapped_edge(self, vertex_id_conversion, other_edge_id):
        """Given an edge (`vertex0_id`, `vertex1_id`) from a
        different CoNetwork, return the corresponding edge id from this
        CoNetwork (`self`).

        Parameters
        -----------
        vertex_id_conversion : dict (int -> int)
          the dict produced by `convert_pathway_mapping,` which maps
          the other CoNetwork's vertex ids to this CoNetwork's vertex ids.
        other_edge_id : tup(int, int)
          the edge id for the edge in the other CoNetwork

        Returns
        -----------
        tup(int, int), the corresponding edge id
        """
        vertex0_id, vertex1_id = other_edge_id
        self_vertex0 = vertex_id_conversion[vertex0_id]
        self_vertex1 = vertex_id_conversion[vertex1_id]
        edge_id = self.edge_tuple(self_vertex0, self_vertex1)
        return edge_id

    def edge_tuple(self, vertex0_id, vertex1_id):
        """To avoid duplicate edges where the vertex ids are reversed,
        we maintain that the vertex ids are ordered so that the corresponding
        pathway names are alphabetical.

        Parameters
        -----------
        vertex0_id : int
          one vertex in the edge
        vertex1_id : int
          the other vertex in the edge

        Returns
        -----------
        tup(int, int)|None, the edge id or None if the vertices do not
        exist in the network or they map to the same pathway (there should not
        be any self-loops in the network)
        """
        pw0 = self.__getitem__(vertex0_id)
        pw1 = self.__getitem__(vertex1_id)

        if not pw0 or not pw1:
            return None
        if pw0 < pw1:
            return (vertex0_id, vertex1_id)
        elif pw0 > pw1:
            return (vertex1_id, vertex0_id)
        else:
            return None

    def add_pathway(self, pathway):
        """Updates `self.pathways` and `self.n_pathways.`

        Parameters
        -----------
        pathway : str
          the pathway to add to the network.
        """
        if pathway not in self.pathways:
            self.pathways[pathway] = self.n_pathways
            self.n_pathways += 1
        return self.pathways[pathway]

    def get_pathway_from_vertex_id(self, vertex_id):
        """Return the pathway string corresponding to a vertex id

        Parameters
        -----------
        vertex_id : int

        Returns
        -----------
        str|None, the pathway name if the vertex is in this network
        """
        return self.__getitem__(vertex_id)

    def get_edge_pathways(self, edge_id):
        """Get the pathways associated with an edge.

        Parameters
        -----------
        edge_id : tup(int, int)

        Returns
        -----------
        tup(str, str)|None, the edge as a pair of 2 pathways if the edge id
          is in this network
        """
        vertex0_id, vertex1_id = edge_id
        pw0 = self.get_pathway_from_vertex_id(vertex0_id)
        pw1 = self.get_pathway_from_vertex_id(vertex1_id)
        if not pw0 or not pw1:
            return None
        return (pw0, pw1)

    def get_vertex_obj(self, vertex_id):
        """Get the Vertex object that corresponds to a vertex id

        Parameters
        -----------
        vertex_id : int

        Returns
        -----------
        Vertex|None, the Vertex obj if the vertex id is in this network
        """
        if vertex_id in self.vertices:
            return self.vertices[vertex_id]
        return None

    def get_vertex_obj_from_pathway(self, pathway):
        """Get the vertex object that corresponds to a pathway name

        Parameters
        -----------
        pathway : str

        Returns
        -----------
        Vertex|None, the Vertex obj if the pathway is in this network
        """
        if pathway in self.pathways:
            vertex_id = self.pathways[pathway]
            return self.vertices[vertex_id]
        else:
            return None

    def get_adjacent_pathways(self, pathway):
        """Get the pathways adjacent to this pathway in the network

        Parameters
        -----------
        pathway : str

        Returns
        -----------
        list(str), a list of pathways adjacent to the input pathway
        """
        vertex_id = self.pathways[pathway]
        adjacent = self.vertices[vertex_id].get_adjacent_vertex_ids()
        adjacent_pathways = []
        for adjacent_id in adjacent:
            adjacent_pathways.append(self.get_pathway_from_vertex_id(
                adjacent_id))
        return adjacent_pathways

    def to_dataframe(self, drop_weights_below=0, whitelist=None):
        """ Conversion of the network to a pandas.DataFrame.

        Parameters
        -----------
        drop_weights_below : int (default=0)
          specify an edge weight threshold - remove all edges with weight
          below this value
        whitelist : [set|list](tup(int, int))|None (default=None)
          option to pass in a set/list of edge ids (tup(int, int)) that should
          be kept in the resulting dataframe

        Returns
        -----------
        pandas.DataFrame
          a pandas.DataFrame containing the network edge information.
          columns = [pw0, pw1, weight]. an additional "features" column is
          returned if this network is not an aggregate of multiple networks.
        """
        network_df_cols = ["pw0", "pw1", "weight"]
        if self.features:
            network_df_cols.append("features")
        network_df = pd.DataFrame(columns=network_df_cols)
        idx = 0
        edge_pathways = set()
        for (v0, v1), edge_obj in self.edges.items():
            if (edge_obj.weight > drop_weights_below and
                    (whitelist is None or (v0, v1) in whitelist)):
                row = [self.__getitem__(v0),
                       self.__getitem__(v1),
                       edge_obj.weight]
                edge_pathways.add(v0)
                edge_pathways.add(v1)
                if self.features:
                    features = edge_obj.features_to_string()
                    row.append(features)
                network_df.loc[idx] = row
                idx += 1  # faster to append by index.
        network_df = network_df.sort_values(by=["weight"],
                                            ascending=False)
        print("The pathway co-occurrence network "
              "contains {0} pathways.".format(
                len(edge_pathways)))
        return network_df

    def _add_edge_to_vertex(self, vertex_id, edge):
        """Adds the edge to the Vertex object's `edges` dictionary
        """
        connected_to = edge.connected_to(vertex_id)
        if vertex_id not in self.vertices:
            vertex_obj = Vertex(vertex_id)
            self.vertices[vertex_id] = vertex_obj
        self.vertices[vertex_id].edges[connected_to] = edge.weight

    def _augment_network(self, edge_dict):
        """Given a dictionary of edges (edge id -> feature list), add all of
        these to the CoNetwork object
        """
        for (vertex0, vertex1), feature_list in edge_dict.items():
            edge_obj = Edge(vertex0, vertex1, feature_list)
            self.edges[(vertex0, vertex1)] = edge_obj
            self._add_edge_to_vertex(vertex0, edge_obj)
            self._add_edge_to_vertex(vertex1, edge_obj)

    def _edges_from_permutation(self, feature_pathway_dict):
        """Given a dictionary mapping each feature to the pathways
        overrepresented in the feature, build a CoNetwork by
        creating edges for every pairwise combination of pathways in a feature.
        """
        network_edges = {}
        for feature, pathway_list in feature_pathway_dict.items():
            for i in range(len(pathway_list)):
                for j in range(i + 1, len(pathway_list)):
                    vertex_i = pathway_list[i]
                    vertex_j = pathway_list[j]
                    new_edge = self.edge_tuple(vertex_i, vertex_j)
                    if new_edge not in network_edges:
                        network_edges[new_edge] = []
                    network_edges[new_edge].append(feature)
        self._augment_network(network_edges)

    def _collect_pathways_by_feature(self, pathway_feature_tuples):
        """Given a tuple list [(pathway, feature)], create a dictionary
        mapping each feature to the pathways overrepresented in the feature.
        """
        pathways_in_feature = {}
        for (pathway, feature) in pathway_feature_tuples:
            vertex = self.add_pathway(pathway)
            if feature not in pathways_in_feature:
                pathways_in_feature[feature] = []
            pathways_in_feature[feature].append(vertex)
        return pathways_in_feature


def _load_significant_pathways_file(path_to_file):
    """Read in the significant pathways file as a
    pandas.DataFrame.
    """
    feature_pathway_df = pd.read_table(
        path_to_file, header=0,
        usecols=["feature", "side", "pathway"])
    feature_pathway_df = feature_pathway_df.sort_values(
        by=["feature", "side"])
    return feature_pathway_df

################################################################
# PATHWAY-FEATURE PERMUTATION HELPER FUNCTIONS.
################################################################


def _pathway_feature_permutation(pathway_feature_tuples,
                                 permutation_max_iters):
    """Permute the pathways across features for one side in the
    network. Used in `permute_pathways_across_features`

    Parameters
    -----------
    pathway_feature_tuples : list(tup(str, int))
      a tuple list [(pathway, feature)] where the pathway, feature pairing
      indicates that a pathway was overrepresented in that feature
    permutation_max_iters : int
      specify the maximum number of iterations, limit the number of attempts
      we have to generate a permutation

    Returns
    -----------
    list(tup(str, int)), the list of pathway, feature pairings after the
      permutation
    """

    pathways, features = [list(elements_at_position)
                          for elements_at_position in
                          zip(*pathway_feature_tuples)]

    original_pathways = pathways[:]

    random.shuffle(pathways)
    feature_block_locations = {}
    i = 0
    while i < len(pathways):
        starting_index = i
        current_feature = features[i]
        pathway_set = set()
        # input is grouped by feature, so we want to keep track of the start
        # and end of a given "block" of the same feature--this corresponds
        # to all the pathways overrepresented in that feature.
        while i < len(pathways) and features[i] == current_feature:
            # check the results of the permutation. if `pathway_set` does
            # not contain the current pathway, we are maintaining the
            # necessary invariants in our permutation thus far.
            if pathways[i] not in pathway_set:
                pathway_set.add(pathways[i])
            else:
                k = 0
                random_pathway = None
                while True:
                    # select another random pathway from the list
                    # and get the feature to which it belongs
                    j = random.choice(range(0, len(pathways)))
                    random_pathway = pathways[j]
                    random_feature = features[j]
                    if (random_pathway != pathways[i] and
                            random_pathway not in pathway_set):
                        # if this is a feature we have not already seen,
                        # we are done.
                        if random_feature not in feature_block_locations:
                            break
                        # otherwise, look at the indices that correspond
                        # to that feature's block of pathways
                        feature_block_start, feature_block_end = \
                            feature_block_locations[random_feature]
                        pathway_block = pathways[feature_block_start:
                                                 feature_block_end]
                        # make sure that the current pathway is not in
                        # that block--ensures that we maintain the invariant
                        # after the swap
                        if pathways[i] not in pathway_block:
                            break
                    k += 1
                    if k > permutation_max_iters:
                        print("Permutation step: reached the maximum "
                              "number of iterations {0}.".format(
                                permutation_max_iters))
                        return None
                pathway_set.add(random_pathway)
                pathways[j] = pathways[i]
                pathways[i] = random_pathway
            i += 1
        ending_index = i
        feature_block_locations[current_feature] = (
            starting_index, ending_index)

    if original_pathways == pathways:
        return None
    return list(zip(pathways, features))


def _permutation_correctness(pathway_feature_tuples, original):
    """Determine whether the generated permutation is a valid permutation.
    Used in `permute_pathways_across_features`.
    (Cannot be identical to the original significant pathways list,
    and a feature should map to a distinct set of pathways.)
    """
    if pathway_feature_tuples:
        if set(pathway_feature_tuples) == set(original):
            return False
        pathways_in_feature = {}
        for (pathway, feature) in pathway_feature_tuples:
            if feature not in pathways_in_feature:
                pathways_in_feature[feature] = set()
            if pathway in pathways_in_feature[feature]:
                return False
            else:
                pathways_in_feature[feature].add(pathway)
    return True
