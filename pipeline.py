import sys
import re
import subprocess
import numpy as np
import time
from enum import Enum
from Condition import Condition

class Solver(Enum):
    PLINGELING = 0
    GLUCOSE_SYRUP = 1

def parse_input(path):
    with open(path, "r") as f:
        lines = np.asarray(f.readlines())
        start_indices = np.asarray([i for i, line in enumerate(lines) if "//" in line] + [len(lines) + 1])
        m = int(lines[start_indices[0] + 1][len("segsites: "):])
        n = start_indices[1] - start_indices[0] - 4
        mats = np.ndarray((len(start_indices) - 1, n, m))

        for i, j in enumerate(start_indices[:-1] + 3):
            mats[i] = np.asarray([list(line.rstrip()) for line in lines[j:(start_indices[i + 1] - 1)]])

    return mats


def gen_tree_node_conditions(n, m, total_nodes):
    # Each leaf must be included in the tree (I(l, k, l) = True)
    root_i_clause = [1]
    node_i_clause = list(range(2, n+m+2))

    i_condition = Condition([root_i_clause] + [node_i_clause], True, n*m, total_nodes)

    leaf_i_condition = Condition(list(), True, m, n*total_nodes)

    for l in range(n):
        for i, val in enumerate(range(n+m+2, total_nodes+1)):
            leaf_i_condition.add_clause([(-1 if i != l else 1) * (l * (total_nodes) + val)])

    final_i_val = n * m * total_nodes

    root_t_condition = Condition([[final_i_val + 1]], True, m, total_nodes)
    leaf_t_condition = Condition([[x + final_i_val + n + m + 2] for x in range(n)], True, m, total_nodes)
    # m commodities and n+m internal nodes
    internal_node_t_condition = Condition(list(), True, n+m, 1)

    for k in range(m):
        for l in range(n):
            internal_node_t_condition.add_clause([-1 * (2 + (k*n+l)*total_nodes), final_i_val + 2 + k*total_nodes])
    
        internal_node_t_condition.add_clause([x for x in range(2, n*total_nodes+2, total_nodes)] + [-1 * (final_i_val + 2 + k*total_nodes)])
    
    node_conditions = [i_condition, leaf_i_condition, root_t_condition, internal_node_t_condition, leaf_t_condition]

    max_val = final_i_val + m*total_nodes

    return node_conditions, max_val

def gen_tree_edge_conditions(input):
    edge_conditions = Condition()
    return edge_conditions

def gen_tree_conditions(input):
    # Have adjacency mat for edges?
    # characters = columns
    # taxa = rows
    n = input.shape[0] # number of rows (and leaves)
    m = input.shape[1] # number of columns
    num_internal_nodes = n + m
    total_nodes = 1 + num_internal_nodes + n # root + internal + leaves
    # one edge from root to each internal node, all internal nodes connected (w/o cycles, so same as undirected handshake problem), one node from each internal node to each leaf
    num_edges = num_internal_nodes + (num_internal_nodes * (num_internal_nodes - 1))//2 + num_internal_nodes * n 

    node_conditions, final_node_var = gen_tree_node_conditions(n, m, total_nodes)

    for l in range(n):
        for c in range(m):
            for i in range(1, num_internal_nodes + 1):
                continue
                # I(i, c, l) variables

    for l in range(n):
        for c in range(m):
            for i in range(1, num_internal_nodes):
                for j in range(i, num_internal_nodes + 1):
                    continue
                    # f(i,j,c,l) variables

    for j in range (num_internal_nodes):
        continue
        # R(j) variables (objective)
        # Should this go here or in its own function like in prototein problem?
    

    conditions = list()
    return conditions

def gen_subtree_conditions(input):
    conditions = list()
    return conditions

def minimize_sat(conditions, solver):
    bound = 0
    time = 0

    results = {"time":time, "solver":solver, "bound":bound, "runs_required":runs_required}
    return results

def minimize_ilp(file):
    sol_file = "./gurobi_output/" + file + ".sol"
    lp_file = "./input/" + file + ".lp"

    # Generate lp file

    start = time.time()
    subprocess.run(["gurobi_cl", "ResultFile=" + sol_file, lp_file], capture_outut=True)
    end = time.time()

    if (result.returncode == 1):
        print(str(result.stdout))
        return
    else:
        with open("./gurobi_output/" + file + ".sol") as f:
            lines = f.readlines()
            if "value =" in lines[0]:
                bound = re.search(r"\d+", lines[0]).group()
            else:
                bound = 0

    results = {"time":end-start, "bound":bound}
    return results

def main(argv):
    outdir = "./output"
    solver = Solver.GLUCOSE_SYRUP

    if len(argv) < 2:
        print("Error: usage\n\tpython3 pipeline.py -o {output directory} -s solver [input files]")
        return
    if "-o" not in argv and "-s" not in argv:
        input_files = argv[1:]
    elif "-o" not in argv or "-s" not in argv:
        input_files = argv[3:]
    else:
        input_files = argv[5:]
    if "-o" in argv:
        outdir = argv[argv.index("-o") + 1]
    if "-s" in argv:
        solver = argv[argv.index("-s") + 1]

    for in_file in input_files:
        input_path = "./input/" + in_file
        input_matrices = parse_input(input_path)

        for mat in input_matrices:
            tree_conditions = gen_tree_conditions(mat)
            subtree_conditions = gen_subtree_conditions(mat, tree_conditions)
            conditions = tree_conditions + subtree_conditions
            sat_results = minimize_sat(conditions, Solver.GLUCOSE_SYRUP)
            ilp_results = minimize_ilp()

        print_results(sat_results, ilp_results)
    
    return

main(["pipeline.py", "-o", "test_output", "-s", "glucose", "test"])