import sys
import numpy as np

def parse_input(path):
    with open(path, "r") as f:
        lines = np.asarray(f.readlines())

        for i in range(len(lines)):
            if "//" in lines[i]:
                i += 1
            else:

    return mat

def gen_tree_conditions(input):
    conditions = list()
    return conditions

def gen_subtree_conditions(input):
    conditions = list()
    return conditions

def minimize_sat(conditions, solver):
    bound = 0
    time = 0

    results = {"time":time, "solver":solver, "bound":bound}
    return results

def minimize_ilp():
    time = 0
    bound = 0

    results = {"time":time, "bound":bound}
    return results

def main(argv):
    if len(argv) < 2:
        print("Error: usage\n\tpython3 pipeline.py -o {output directory} -s solver [input files]")
        return
    if "-o" not in argv:
        outdir = "./output"
        input_files = argv[1:]
    else:
        outdir = argv[-1]
        input_files = argv[1:len(argv) - 1]

    for input_path in input_files:
        with open("./input/" + input_path, "r") as f:
            input_sequences = parse_input(input_path)
            for sequence in input_sequences:
                tree_conditions = gen_tree_conditions(sequence)
                subtree_conditions = gen_subtree_conditions(sequences, tree_conditions)
                conditions = tree_conditions + subtree_conditions
                sat_results = minimize_sat(conditions, solver)
                ilp_results = minimize_ilp()

            print_results(sat_results, ilp_results)
    
    return

main(sys.argv)