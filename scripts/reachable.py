#!/usr/bin/env python3


def reachable(G, X, Z):
    """Returns set of nodes reachable from X given Z via active trails, for Bayes Net G.

    G: Graph of Bayes Net
    X: Source Random Variable
    Z: Observation Random Variables

class Digraph:
    """A simple directed graph."""
    def __init__(self):
        self.nodes = set()
        self.forward_edges = dict()
        self.reverse_edges = dict()
        self.edges = self.forward_edges

    def add(self, node, children):
        """Adds a node and its outward edges to the graph.

        node: Identifies the node.  Must be hashable.

        children: An iterable of nodes corresponding to the outward edges.  Must also be hashable."""

        new_nodes = set.union({node}, children) - self.nodes
        self.nodes |= new_nodes
        for n in new_nodes:
            self.forward_edges[n] = set()
            self.reverse_edges[n] = set()

        self.forward_edges[node].update(children)
        for child in children:
            self.reverse_edges[child].add(node)

    def parents(self, node):
        """Returns the parents of the node."""
        return self.reverse_edges[node]

    def children(self, node):
        """Returns the children of the node."""
        return self.forward_edges[node]

    """

    pass
