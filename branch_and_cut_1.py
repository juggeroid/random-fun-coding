import argparse
import itertools
import math
import multiprocessing
import operator
import sys
import threading
from collections import defaultdict
from time import perf_counter

import docplex
import networkx as nx
from docplex.mp.model import Model
from docplex.mp.model import Model
from networkx.algorithms.approximation import ramsey


###################################################################################################################

class Utilities:
    @staticmethod
    def parse_dimacs(file):
        edges = []
        with open(file, 'r') as file:
            for line in file:
                if line.startswith('p'):
                    print(line)
                if line.startswith('e'):
                    dummy, v_1, v_2 = line.split()
                    edges.append((v_1, v_2))
                else:
                    continue
        g = nx.Graph(edges)
        g = nx.relabel.convert_node_labels_to_integers(g)
        return g

    class reversible_uniq_dict(dict):
        def reversed(self):
            reversed = {}
            for k, v in self.items():
                reversed.setdefault(v, []).append(k)
            return reversed

    @staticmethod
    def is_integer(x: float):
        result_bool = math.isclose(x, round(x), rel_tol = 1e-7)
        return result_bool

    @staticmethod
    def epsilon_rounding(x):
        if Utilities.is_integer(x):
            ret = round(x)
        else:
            ret = x
        return ret

class Algorithms:

    @staticmethod
    def is_clique(graph, clique):
        subgraph = graph.subgraph(clique)
        N = len(clique)
        if subgraph.size() == N * (N - 1) / 2:
            return None
        for (u, v) in itertools.combinations(clique, 2):
            if not graph.has_edge(u, v):
                return (u, v)

    @staticmethod
    def coloring_set_batch_generate(graph, iterations):
        colorings = []
        independent_sets = []
        for _ in range(0, iterations):
            coloring = nx.greedy_color(graph, strategy = "random_sequential")
            colorings.append(coloring)
            independent_sets.append(Utilities.reversible_uniq_dict(coloring).reversed())
        return colorings, independent_sets

    @staticmethod
    def heuristic_by_removal(graph_reference):
        graph = graph_reference.copy()
        largest_clique, independent = nx.algorithms.approximation.ramsey.ramsey_R2(graph)
        sets = [independent]
        while graph:
            graph.remove_nodes_from(largest_clique)
            largest_clique, independent = nx.algorithms.approximation.ramsey.ramsey_R2(graph)
            if independent:
                sets.append(independent)
        return sets

    @staticmethod
    def initial_clique_heuristic(graph):
        sets = Algorithms.heuristic_by_removal(nx.complement(graph))
        return max(sets, key = len)

###################################################################################################################

PRECOMPUTED_COLORINGS = 5000
TAILING_THRESHOLD = 5

###################################################################################################################


class Solver:

    cplex_model: Model

    def __init__(self, path):
        self.cplex_model = docplex.mp.model.Model()
        self.cplex_model.parameters.threads = multiprocessing.cpu_count()
        print("CPLEX Model has been successfully initialized.")

        t0 = perf_counter()
        self.graph = Utilities.parse_dimacs(path)
        t1 = perf_counter()
        print("Preprocessing took {} s.".format(t1 - t0))

        t2 = perf_counter()
        self.colorings, self.independent_sets = Algorithms.coloring_set_batch_generate(self.graph, PRECOMPUTED_COLORINGS)
        self.initial_clique = Algorithms.initial_clique_heuristic(self.graph)
        clique_size = len(self.initial_clique)
        self.best_found = clique_size
        t3 = perf_counter()
        print("Found the clique of size {}, initializing the variables.".format(clique_size))
        print("Initial heuristic and coloring took {} s.".format(t3 - t2))

        self.variables = {node: self.cplex_model.continuous_var(name='x_{}'.format(node)) for node in self.graph.nodes()}


    ###################################################################################################################

    def separation(self, xs, top_k = None):

        if len(self.colorings) == 0:
            coloring = nx.algorithms.greedy_color(self.graph, strategy = "random_sequential")
            independent_sets = Utilities.reversible_uniq_dict(coloring).reversed()
        else:
            coloring = self.colorings.pop()
            independent_sets = self.independent_sets.pop()

        weights = defaultdict(int)
        for k, v in coloring.items():
            weights[v] += xs[k]

        constraints = []
        sorted_batch = sorted(weights, key = weights.get, reverse = True)
        if top_k is not None:
            sorted_batch = sorted_batch[:top_k] # Set it manually for bigger graphs.

        for k in sorted_batch:
            if Utilities.epsilon_rounding(weights[k]) <= 1:
                break
            else:
                lhs = sum(self.variables[i] for i in independent_sets[k])
                constraint = self.cplex_model.linear_constraint(lhs, 1, "le")
                constraints.append(constraint)

        return constraints

    ###################################################################################################################

    def initialize_constraints(self):
        self.cplex_model.add_constraints([self.variables[i] <= 1 for i in self.graph.nodes()])
        self.cplex_model.maximize(self.cplex_model.sum(self.variables))

    ###################################################################################################################

    def manage_branch_ordering(self, order):
        for constraint in order:
            self.cplex_model.add_constraint(constraint)
            self.branch_and_cut()
            self.cplex_model.remove_constraint(constraint)

    def branching_decision(self, cplex_variable_to_branch, actual_value):
        l = self.cplex_model.linear_constraint(cplex_variable_to_branch, 0, "le")
        r = self.cplex_model.linear_constraint(cplex_variable_to_branch, 1, "ge")
        if actual_value >= 0.5:
            self.manage_branch_ordering([r, l])
        else:
            self.manage_branch_ordering([l, r])

    def get_branching_candidate(self, xs):
        candidates = [(i, value) for i, value in enumerate(xs) if not Utilities.is_integer(value)]
        try:
            found_float = max(candidates, key=operator.itemgetter(1))[0]
        except (ValueError, TypeError):
            found_float = None
        return found_float

    def has_improved(self, s1):
        if Utilities.is_integer(s1):
            s1 = round(s1)
        else:
            s1 = math.floor(s1)
        return s1 > self.best_found

    def branch_and_cut(self):

        current = self.cplex_model.solve()
        if current is None:
            return
        xs = current.get_all_values()
        objective = current.get_objective_value()
        if not self.has_improved(objective):
            return

        enhanced_objective = 0
        iteration = 0
        while iteration < TAILING_THRESHOLD:
            constraints = self.separation(xs)
            if not constraints:
                break
            self.cplex_model.add_constraints(constraints)
            current = self.cplex_model.solve()
            if current is None:
                return
            xs = current.get_all_values()
            objective = current.get_objective_value()
            if not self.has_improved(objective):
                return
            if objective - enhanced_objective <= 1e-2:
                iteration += 1
            else:
                enhanced_objective = objective
                iteration = 0

        branching_candidate = self.get_branching_candidate(xs)
        if branching_candidate is None:
            missing_edges = Algorithms.is_clique(self.graph, [i for i, value in enumerate(xs) if value == 1.0])
            if missing_edges:
                source, target = missing_edges
                self.cplex_model.add_constraint(self.variables[source] + self.variables[target] <= 1)
                self.branch_and_cut()
            else:
                self.best_found = Utilities.epsilon_rounding(objective)
                return
            return

        self.branching_decision(self.cplex_model.get_var_by_index(branching_candidate), xs[branching_candidate])

    ###################################################################################################################

def main(args):
    solver = Solver(args.graph)
    solver.initialize_constraints()
    t5 = perf_counter()
    solver.branch_and_cut()
    t6 = perf_counter()
    print("[{}]: Time: {}, found: {}".format(args.graph, t6 - t5, solver.best_found))

if __name__ == "__main__":
    sys.setrecursionlimit(100000)
    threading.stack_size(200000000)
    parser = argparse.ArgumentParser()
    parser.add_argument("graph")
    args = parser.parse_args()
    main(args)

# find ./ -maxdepth 1 -type f -name *.clq -exec python branch_and_cut.py {} \;