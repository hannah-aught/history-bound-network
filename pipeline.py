import sys
import re
import subprocess
import numpy as np
import time
from enum import Enum
from Condition import Condition
from sympy.logic.boolalg import to_cnf, Equivalent, Implies
from sympy import symbols

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

def gen_i_conditions(n, m, total_nodes):
    # Each leaf must be included in the tree (I(l, k, l) = True)
    root_i_condition = Condition([[1]], True, n*m, total_nodes-(n-1))

    i_condition = Condition([list(range(2, n+m+2))], True, n*m, total_nodes - n + 1)
    leaf_i_condition = Condition([[n+m+2]], True, n*m, total_nodes - n + 1)

    final_i_val = n * m * (total_nodes - n + 1)

    return [root_i_condition, i_condition, leaf_i_condition], final_i_val

def gen_t_conditions(n, m, total_nodes, final_i_var):
    # m commodities and n+m internal nodes
    # commodities can't be repeated over using the Condition (requires adding a different number to first elements than last), but each node can be
    root_t_condition = Condition([[-1*(final_i_var + 1)]], True, m, total_nodes)
    t_condition = Condition(list(), True, m, total_nodes)
    leaf_t_condition = Condition([[x] for x in range(final_i_var + n + m + 2, final_i_var + n + m + n + 2)], True, m, total_nodes)

    for j in range(n + m):
        for l in range(n):
            t_condition.add_clause([-1 * (2 + (l)*(total_nodes - n + 1) + j), final_i_var + 2 + j])
    
        t_condition.add_clause([x for x in range((l-1)*(total_nodes - n + 1) + j + 2, 3 + j + (n-1)*(total_nodes - n + 1), total_nodes - n + 1)] + [-1 * (final_i_var + 2 + j)])

    final_t_var = final_i_var + m*total_nodes

    return [root_t_condition, t_condition, leaf_t_condition], final_t_var

def gen_f_conditions(n, m, total_edges, final_t_var):
    f_condition_1 = Condition([list(range(final_t_var+1, final_t_var+n+m+1))], True, m*n, total_edges - (m+n)*(n-1))
    f_condition_2 = Condition(list(), False)

    f_vars = list(range(final_t_var + 1, final_t_var + total_edges - (m+n)*(n-1) + 1))

    for k in range(m):
        for l in range(n):
            for j in range(1, m+n+1):
                current_node_f_vars = [final_t_var + (k*n+l)*(total_edges-(m+n)*(n-1)) + j]
                current_i_var = (n+m + 2)*(l+k*n) + j + 1
                start_i_var = current_i_var - j
                f_condition_2.add_clause([-1*current_node_f_vars[-1], current_i_var])

                for i in range(j-1):
                    current_node_f_vars.append(current_node_f_vars[-1] + m + n - max(i,1))
                    f_condition_2.add_clause([-1*current_node_f_vars[-1], current_i_var])

                f_condition_2.add_clause(current_node_f_vars + [-1*current_i_var])

            # Condition saying that there must be flow from at least one internal node to the leaf we're currently concerned with
            # The condition for *only* one edge going to the leaf comes in F condition 3
            f_condition_2.add_clause([x+1 for x in current_node_f_vars[1:]] + [current_node_f_vars[-1] + 2])

    f_condition_3 = Condition(list(), True, m*n, total_edges-(m+n)*(n-1))

    for j in range(2, m+n+2):
        start_f_var = final_t_var + j

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
    
    f_condition_4 = Condition(list(), True, m*n, total_edges-(m+n)*(n-1))
    current_f_vars = f_vars

    for i in range(m+n): # don't need to look at edges out of last internal node bc there's only one, going to leaf l
        for start_j in range(min(m+n, m+n+1-i)):
            first_f_var = current_f_vars[start_j]

            for j in range(start_j + 1, min(m+n, m+n+1-i)):
                second_f_var = current_f_vars[j]
                f_condition_4.add_clause([-1*first_f_var, -1*second_f_var])
        
        current_f_vars = current_f_vars[min(m+n, m+n+1-i):]


    f_condition_5 = Condition(list(), True, n*m, total_edges-(m+n)*(n-1))
    current_f_vars = f_vars

    for i_prime in range(m+n):
        start_f_var = f_vars[i_prime]
        current_f_i_prime_vars = f_vars[m+n:]
        current_f_vars = current_f_vars[min(m+n, m+n+1-i_prime):]
        f_i_prime_vars = [f_vars[i_prime]]


        for i in range(i_prime):
            f_i_prime_vars.append(current_f_i_prime_vars[i_prime-i-1])
            current_f_i_prime_vars = current_f_i_prime_vars[min(m+n, m+n+1-i):]

        for f_i in current_f_vars[:min(m+n, m+n-i_prime)]: # second arg to min is simplified from m+n+2-(i+1)               
            f_condition_5.add_clause([-1*(f_i)] + [x for x in f_i_prime_vars])


    last_f_var = final_t_var + m*n*(f_vars[-1] - final_t_var)

    return [f_condition_1, f_condition_2, f_condition_3, f_condition_4, f_condition_5], last_f_var, f_vars

def gen_x_conditions(n, m, total_edges, last_f_var, f_vars):
    x_condition = Condition(list(), True, m, total_edges)
    x_vars = list(range(last_f_var + 1, last_f_var + total_edges + 1))
    leaf_f_vars = [0 for x in range(m+n)]
    last_index = len(f_vars) - 1

    for f in range(m+n):
        leaf_f_vars[m+n-f-1] = f_vars[-1-f*(f+1)//2]

    last_x_var = last_f_var + m*(x_vars[-1] - last_f_var)
    x_i = 0

    for i, f_var in enumerate(f_vars):
        if f_var in leaf_f_vars:
            for x in range(n):
                x_condition.add_clause([x*(total_edges-(n+m)*(n-1)) + f_var, -1*x_vars[x_i + i +x]])
                x_condition.add_clause([-1*(x*(total_edges-(n+m)*(n-1)) + f_var), x_vars[x_i + i +x]])

            x_i += x
            continue

        for l in range(n):
            x_condition.add_clause([-1*(l*(total_edges-(n+m)*(n-1)) + f_var), x_vars[x_i+i]])

        x_condition.add_clause([x*(total_edges-(n+m)*(n-1)) + f_var for x in range(n)] + [-1*x_vars[x_i+i]])

    
    return x_condition, last_x_var, x_vars

def gen_d_conditions(n, m, total_edges, last_x_var, x_vars):
    d_condition = Condition(list(), True, m, total_edges)
    d_vars = list(range(last_x_var + 1, last_x_var + total_edges + 1))

    for i, x_var in enumerate(x_vars):
        for k in range(m):
            d_condition.add_clause([-1*(k*(total_edges-(n+m)*(n-1)) + x_var), d_vars[i]])
        d_condition.add_clause([x_vars[i] + x*(total_edges-(m+n)*(n-1)) for x in range(m)] + [-1*d_vars[i]])

    return d_condition, d_vars[-1]

def gen_tree_conditions(n, m):
    # Have adjacency mat for edges?
    # characters = columns
    # taxa = rows
    num_internal_nodes = n + m
    total_nodes = 1 + num_internal_nodes + n # root + internal + leaves
    total_edges = (n+m)*(1 + (n+m-1)//2 + n)

    i_conditions, final_i_var = gen_i_conditions(n, m, total_nodes) # Conditions for each node being included in a commodity tree with flow passing through them to a specific leaf
    t_conditions, final_t_var = gen_t_conditions(n, m, total_nodes, final_i_var) # Conditions for each node being included in a commodity tree
    f_conditions, final_f_var, f_vars = gen_f_conditions(n, m, total_edges, final_t_var) # Conditions for each edge being included in a commodity tree going to a specific leaf
    x_conditions, final_x_var, x_vars = gen_x_conditions(n, m, total_edges, final_f_var, f_vars) # Condition for edges being included in a commodity tree
    d_conditions, final_d_var = gen_d_conditions(n, m, total_edges, final_x_var, x_vars) # Condition for edges being included in the DAG

    conditions = i_conditions + t_conditions + f_conditions + [x_conditions] + [d_conditions]
    return conditions, final_t_var, final_d_var

def gen_subtree_conditions(input, n, m, num_edges, final_node_var, final_edge_var, final_r_var):
    conditions = list()
    total_nodes = 2 * n + m + 1
    i_offset = n*m*(n+m+2)
    final_rct_var = final_r_var + m*total_nodes
    final_ct_var = final_rct_var + m*total_nodes
    rct_vars = symbols([str(r) for r in range(final_r_var + 1, final_rct_var + 1)])
    ct_vars = symbols([str(c) for c in range(final_rct_var + 1, final_ct_var + 1)])
    x_vars = symbols([str(x) for x in range(final_edge_var - (m+1)*num_edges + 1, final_edge_var - num_edges + 1)])
    rct_condition_1 = Condition([[i_offset + 1, -1*(final_r_var + 1)]], True, m*total_nodes, 1)
    rct_condition_2 = Condition([[x for x in range(final_r_var + 2, final_r_var + m+2*n+2)]], True, m, total_nodes)
    rct_condition_3 = Condition(list(), True, m, n+m+1)

    # Root can't be the root of the subtree
    rct_condition_3.add_clause([-1*(final_r_var + 1)])

    for i in range(2, total_nodes):
        for j in range(i+1, total_nodes + 1):
            rct_condition_3.add_clause([-1*(final_r_var + i), -1*(final_r_var + j)])


    #CT vars exist for leaves too
    ct_condition = Condition([[-1*(final_rct_var + 1)]], True, m, total_nodes)

    offset = 0
    clauses = list()

    for j in range(1,total_nodes):
        dnf_clause = rct_vars[j]
        offset = j-1

        for i in range(min(j, m+n+1)):
            if i == 0 and j > n + m:
                offset += n+m-1
                continue
            
            dnf_clause = dnf_clause | (ct_vars[i] & x_vars[offset])

            if i == 0:
                offset += m+n-1
            else:
                offset += m+n-(i-1)
        
        sympy_clauses = to_cnf(Equivalent(ct_vars[j], dnf_clause))
        clauses = sympy_to_dimacs(str(sympy_clauses))

        for clause in clauses:
            ct_condition.add_clause(clause)
    
    leaf_ct_condition = Condition(list(), False)

    last_ct_internal_var = final_rct_var + m+n+1
    non_zero_indices = input.nonzero()

    for k in range(m):
        true_ct_vars = non_zero_indices[0][np.where(non_zero_indices[1] == k)]
        
        for l in range(n):
            if l in true_ct_vars:
                leaf_ct_condition.add_clause([last_ct_internal_var + k*n + l + 1])
            else:
                leaf_ct_condition.add_clause([-1*(last_ct_internal_var + k*n + l + 1)])

    return [rct_condition_1, rct_condition_2, ct_condition, leaf_ct_condition]

def sympy_to_dimacs(expr):
    clauses = expr.split('&')

    for i, clause in enumerate(clauses):
        clause = clause.strip().replace(" | ", " ").replace("(", "").replace(")", "").replace("~", "-")
        clauses[i] = [int(x) for x in clause.split(" ")]

    return clauses

def gen_reticulation_conditions(n, m, goal_count, num_edges, final_d_var):
    r_vars = symbols([str(x) for x in range(final_d_var +  1, final_d_var + 2*n + m + 2)])
    d_vars = symbols([str(x) for x in range(final_d_var - num_edges + 1, final_d_var + 1)])
    final_r_var = final_d_var + 2*n + m + 1

    c_vars = symbols([str(x) for x in range(final_r_var + 1, final_r_var + (2*n + m + 1)*(goal_count + 2) + 1)])
    r_condition = Condition(list(), False)
    c_condition = Condition(list(), False)

    r_condition.add_clause([-1*(final_d_var + 1)])

    j_node = 2
    offset = 1

    for r in r_vars[1:-1]:
        clauses = list()
        vars = list()

        for i in range(min(j_node, m+n+1)):
            if i == 0 and j_node > n + m:
                offset += n+m-1
                continue # no edge from 0 to leaves

            vars.append(d_vars[offset])

            if i == 0:
                offset += m+n-1
            else:
                offset +=  m+n-(i-1)
        
        dnf_clause = False

        for x, var in enumerate(vars[:-1]):
            for y in range(x + 1, len(vars)):
                dnf_clause = dnf_clause | (var & vars[y])
        

        sympy_clauses = to_cnf(Equivalent(r, dnf_clause))
        clauses = sympy_to_dimacs(str(sympy_clauses))
        j_node += 1
        offset = j_node-1


        for clause in clauses:
            r_condition.add_clause(clause)

    
    for k in range(goal_count + 1):
        for i in range(len(r_vars)-1):
            current_c_var_index = k*(m+2*n+1) + i
            sympy_clauses = to_cnf(Implies(c_vars[current_c_var_index], c_vars[current_c_var_index + 1]))
            sympy_clauses = sympy_clauses & to_cnf(Implies(c_vars[current_c_var_index] & r_vars[i], c_vars[(k+1)*(m+2*n+1) + i + 1]))
            clauses = sympy_to_dimacs(str(sympy_clauses))

            for clause in clauses:
                c_condition.add_clause(clause)


    final_c_var = final_d_var + len(r_vars) + len(c_vars)
    c_condition.add_clause([-1*final_c_var])
    conditions = [r_condition, c_condition]


    return conditions, final_c_var

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
            num_edges = (n+m)*(1 + (n+m-1)//2 + n)
            c = 0

            tree_conditions, final_node_var, final_edge_var = gen_tree_conditions(n, m) #, final_node_var, final_edge_var 
            reticulation_conditions, final_r_var = gen_reticulation_conditions(n, m, c, num_edges, final_edge_var)
            subtree_conditions = gen_subtree_conditions(mat, n, m, num_edges, final_node_var, final_edge_var, final_r_var)
            conditions = tree_conditions + reticulation_conditions + subtree_conditions

            with open('test', 'w+') as f:
                for condition in conditions:
                    condition.write_condition(f)

            sat_results = minimize_sat(conditions, Solver.GLUCOSE_SYRUP)
            ilp_results = minimize_ilp()

        print_results(sat_results, ilp_results)
    
    return

main(["pipeline.py", "-o", "test_output", "-s", "glucose", "test"])