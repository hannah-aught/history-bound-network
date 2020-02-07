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


def gen_tree_node_conditions(n, m):
    # Each leaf must be included in the tree (I(l, k, l) = True)
    total_nodes = 2*n+m+1
    root_i_condition = Condition([[1]], True, n*m, total_nodes)

    i_condition = Condition([list(range(2, n+m+2))], True, n*m, total_nodes - n + 1)
    leaf_i_condition = Condition([n+m+2], True, m, n*(total_nodes - n + 1))

    final_i_val = n * m * (total_nodes - n + 1)

    # m commodities and n+m internal nodes
    # commodities can't be repeated over using the Condition (requires adding a different number to first elements than last), but each node can be
    t_condition = Condition([[-1 * (final_i_val + 1)]], True, 2*n+m, 1)

    for k in range(m):
        for l in range(n):
            t_condition.add_clause([-1 * (2 + (k*n+l)*(total_nodes - n + 1)), final_i_val + 2 + k*(2*n+m)])
    
        t_condition.add_clause([x for x in range(k*n*(total_nodes - n + 1)+2, 3 + (k*n+l)*(total_nodes - n + 1), total_nodes - n + 1)] + [-1 * (final_i_val + 2 + k*(2*n+m))])
    
    node_conditions = [root_i_condition, i_condition, leaf_i_condition, t_condition]

    max_val = final_i_val + m*(2*n+m)

    return node_conditions, max_val

def gen_tree_edge_conditions(n, m, total_nodes, total_edges, num_node_vars):
    f_condition_1 = Condition([list(range(num_node_vars+1, num_node_vars+n+m+1))], True, m*n, total_edges)
    f_condition_2 = Condition(list(), False)
    f_condition_4 = Condition(list(), False)

    f_vars = list(range(num_node_vars + 1, num_node_vars + total_edges - (m+n)*(n-1) + 1))

    for k in range(m):
        for l in range(n):
            for j in range(1, m+n+1):
                current_node_f_vars = [num_node_vars + (k*n+l)*(total_edges-(m+n)*(n-1)) + j]
                current_i_var = (n+m + 2)*(l+k*n) + j + 1
                start_i_var = current_i_var - j
                f_condition_2.add_clause([-1*current_node_f_vars[-1], current_i_var])
                f_condition_4.add_clause([-1*current_node_f_vars[-1], start_i_var])


                for i in range(j-1):
                    current_node_f_vars.append(current_node_f_vars[-1] + m + n - max(i,1))
                    f_condition_2.add_clause([-1*current_node_f_vars[-1], current_i_var])
                    f_condition_4.add_clause([-1*current_node_f_vars[-1], (n+m + 2)*l+k*n + i + 2])

                f_condition_2.add_clause(current_node_f_vars + [-1*current_i_var])

            # Condition saying that there must be flow from at least one internal node to the leaf we're currently concerned with
            # The condition for *only* one edge going to the leaf comes in F condition 3
            f_condition_2.add_clause([x+1 for x in current_node_f_vars[1:]] + [current_node_f_vars[-1] + 2])

            for i, x in enumerate(current_node_f_vars[1:]):
                f_condition_4.add_clause([-1*(x+1), (m+n*2)*(l+k*n)+i + 2])
            f_condition_4.add_clause([-1*(current_node_f_vars[-1] + 2), (m+n*2)*(l+k*n)+i+3])

    f_condition_3 = Condition(list(), True, m*n, total_edges-n+1)

    for j in range(2, m+n+2):
        start_f_var = num_node_vars + j

        for start_i in range(j - 1):
            next_f_var = start_f_var + m + n - max(start_i, 1)

            if start_i == 0 and j == m+n+1:
                # no edge from the root to any leaf, so continue
                start_f_var = start_f_var + m + n - max(start_i, 1)
                continue

            for i in range(start_i + 1, j):
                f_condition_3.add_clause([-1*start_f_var, -1*next_f_var])
                next_f_var = next_f_var + m + n - max(i, 1)

            start_f_var = start_f_var + m + n - max(start_i, 1)

    last_f_var = num_node_vars + m*n*(f_vars[-1] - num_node_vars)

    f_condition_5 = Condition(list(), True, n, num_edges*m)
    f_condition_6 = Condition(list(), True, m, num_edges)

    for f in f_vars:
        f_condition_5.add_clause([x * f for x in range(1, n+1)])

    x_condition = Condition([[-1 * (num_node_vars + 1), num_f_vars + 1], [num_node_vars + 1, -1 * (num_f_vars + 1)]], True, total_nodes*n*m, total_nodes - 1)

    step = 4
    clause_index = 4

    while step <= m+n+2:
        for i in range(1, len(f_condition_1.clauses[clause_index])-1):
            for j in range(j+1, len(f_condition_1.clauses[clause_index])):
                f_condition_2.add_clause([-1*f_condition_1.clauses[clause_index][i], -1*f_condition_1.clauses[clause_index][j]])
            x_condition.add_clause([-1*f_condition_1.clauses[clause_index][i], x_var_num])
        x_condition.add_clause([-1*f_condition_1.clauses[clause_index][i+1], x_var_num])
        x_condition.add_clause(f_condition_1.clauses[clause_index][1:] + [-1 * x_var_num])
        clause_index = clause_index + step
        x_var_num = x_var_num + 1
        step = step + 1

    for x in f_condition_1.clauses[-1]:
        x_condition.add_clause([-1*x, x_var_num])
    x_condition.add_clause(f_condition_1.clauses[-1][1:] + [x_var_num])

    edge_conditions = [f_condition_1, f_condition_2, x_condition]
    max_val = num_node_vars + num_f_vars + m*(x_var_num - num_f_vars)
    return edge_conditions, x_var_num

def gen_tree_conditions(n, m):
    # Have adjacency mat for edges?
    # characters = columns
    # taxa = rows
    num_internal_nodes = n + m
    total_nodes = 1 + num_internal_nodes + n # root + internal + leaves
    total_edges = (n+m)*(1 + (n+m-1)//2 + n)

    node_conditions, final_node_var = gen_tree_node_conditions(n, m)
    edge_conditions, final_edge_var = gen_tree_edge_conditions(n, m, total_nodes, total_edges, final_node_var)

    conditions = [node_conditions, edge_conditions]
    return conditions, final_node_var, final_edge_var

def gen_subtree_conditions(input, n, m, final_node_var, final_edge_var):
    conditions = list()
    total_nodes = 2 * n + m + 1
    i_offset = n*m*(n+m)
    rct_condition_1 = Condition([[i_offset + 1, -1*(final_edge_var + 2)]], True, m*(n+m), 1)
    rct_condition_2 = Condition(list(), True, m, n+m+1)

    # Root can't be the root of the subtree
    rct_condition_2.add_clause([-1*(final_edge_var + 1)])

    for i in range(2, n+m+1):
        for j in range(i+1, n+m+2):
            rct_condition_2.add_clause([-1*(final_edge_var + i), -1*(final_edge_var + j)])

    final_rct_var = final_edge_var + m*(n+m+1)

    #CT vars exist for leaves too
    ct_condition = Condition([[-1*(final_rct_var + 1)]], True, m, total_nodes)
    
    for i in range(1, n+m+2):
        for j in range(2, n+m+2):
            ct_condition.add_clause([final_edge_var + j])
    
    leaf_ct_condition = Condition(list(), False)

    last_ct_internal_var = final_rct_var + m+n+1
    non_zero_indices = input.nonzero()

    for k in range(m):
        true_ct_vars = non_zero_indices[0][np.where(non_zero_indices[1] == k)]
        
        for l in range(n):
            if l in true_ct_vars:
                leaf_ct_condition.add_clause([last_ct_internal_var + k*n + l + 1])
            else:
                leaf_ct_condition.add_clause([-1*last_ct_internal_var + k*n + l + 1])

    return [rct_condition_1, rct_condition_2, ct_condition, leaf_ct_condition]

def gen_reticulation_conditions():
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
            n = mat.shape[0]
            m = mat.shape[1]
            tree_conditions, final_node_var, final_edge_var = gen_tree_conditions(n, m)
            subtree_conditions = gen_subtree_conditions(mat, n, m, final_node_var, final_edge_var)
            conditions = tree_conditions + subtree_conditions
            sat_results = minimize_sat(conditions, Solver.GLUCOSE_SYRUP)
            ilp_results = minimize_ilp()

        print_results(sat_results, ilp_results)
    
    return

main(["pipeline.py", "-o", "test_output", "-s", "glucose", "test"])