import sys

def main(args):
    file_path = args[1]
    start_var = int(args[2])
    end_var = int(args[3])
    num_rows = int(args[4])
    num_cols = int(args[5])
    d_vars = list()
    var_strings = list()
    node_degrees = [0 for x in range(2*num_rows + num_cols)]

    with open(file_path) as f:
        lines = f.readlines()
    
        for line in lines:
            if line[0] != "v":
                continue
            line = [int(x) for x in line[2:].split(" ")]

            for x in line:
                if abs(x) >= start_var and abs(x) <= end_var:
                    d_vars.append(x)

    i = 0
    j = 1
    current_var = start_var

    for var in d_vars:
        while current_var < abs(var):
            var_strings.append("D(" + str(i) + "," + str(j) + ")")
            node_degrees[j-1] += 1

            j += 1

            if (i == 0 and j == num_rows + num_cols) or j == 2*num_rows + num_cols:
                i += 1
                j = i+1
                print(var_strings)
                var_strings = list()
            

            current_var += 1


        if var > 0:
            var_strings.append("D(" + str(i) + "," + str(j) + ")")
            node_degrees[j-1] += 1


        j += 1
        current_var += 1

        if (i == 0 and j == num_rows + num_cols + 1) or j == 2*num_rows + num_cols + 1:
            i += 1
            j = i+1
            print(var_strings)
            var_strings = list()



    print(", ".join(var_strings))
    print([x if x > 1 else "" for x in node_degrees])



if __name__ == "__main__":
    #main(["parseDag.py", "vars10x10.txt", "29611", "30020", "10", "10"])
    main(sys.argv)