#!/usr/bin/env python3

from enum import Enum

class Direction(Enum):
    UP = 1
    DOWN = 2


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


def ancestors(G, X):
    """Returns the ancestors of node X in Digraph G.  Includes X itself."""

    nodes_to_visit = {X}
    ancestors = set()

    while nodes_to_visit:
        node = nodes_to_visit.pop()
        if node not in ancestors:
            nodes_to_visit |= G.parents(node)
            ancestors.add(node)

    return ancestors


def reachable(G, source, observations):
    """Returns the set of nodes reachable via active trails from a source given a list of observations, for Bayes Net G.

    G: Digraph of Bayes Net
    source: Source node in G
    observations: List of Observation nodes in G

    See page 75 of Koller and Friedman, Probabilistic Graphical Models: Principles and Techniques
    """

    source_ancestors = ancestors(G, source)

    to_visit = {(source, Direction.UP)}
    visited = set()
    reachable_nodes = set()

    while to_visit:
        node, the_dir = to_visit.pop()
        if (node, the_dir) not in visited:
            visited.add((node, the_dir))
            if node not in observations:
                reachable_nodes.add(node)
            if the_dir is Direction.UP and node not in observations:
                to_visit |= {(p, Direction.UP) for p in G.parents(node)}
                to_visit |= {(c, Direction.DOWN) for c in G.children(node)}
            elif the_dir is Direction.DOWN:
                if node not in observations:
                    to_visit |= {(c, Direction.DOWN) for c in G.children(node)}
                if node in source_ancestors:
                    to_visit |= {(p, Direction.UP) for p in G.parents(node)}
    return reachable_nodes
