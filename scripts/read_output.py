import os

# Opening file
# dir_path = 'data/output/'
dir_path = "data/output/"
for path in os.listdir(dir_path):
    # check if current path is a file
    if os.path.isfile(os.path.join(dir_path, path)):
        
        if "output_data" in path:
            print("== path:", path)

            # file_to_open = "output_data70_b750000_q4_mc100_r150.csv"
            file1 = open(dir_path+path, 'r')
            count = 0

            counting_dict_max_coverage = {}
            counting_dict_max_profit = {}
            counting_dict_percentage_cycles = {}
            for line in file1:
                count += 1
                # print("Line{}: {}".format(count, line.strip()))
                line_array = line.split(" | ")
                # print(line_array)
                name_algorithm = line_array[0].split(" B: ")[0] + "-" + (line_array[0].split(" B: ")[1]).split(" #iterazioni: ")[0]

                # profitto massimo
                if name_algorithm in counting_dict_max_profit:
                    counting_dict_max_profit[name_algorithm] += float(line_array[4].split(" p_max: ")[1])
                else:
                    counting_dict_max_profit[name_algorithm] = float(line_array[4].split(" p_max: ")[1])

                # maximum coverage
                if name_algorithm in counting_dict_max_coverage:
                    counting_dict_max_coverage[name_algorithm] += float(line_array[5].split(" - max: ")[1].split("%)")[0])
                else:
                    counting_dict_max_coverage[name_algorithm] = float(line_array[5].split(" - max: ")[1].split("%)")[0])

                # percentage cycles
                if name_algorithm in counting_dict_percentage_cycles:
                    counting_dict_percentage_cycles[name_algorithm] += float(line_array[7].split("unici: ")[1])/400
                else:
                    counting_dict_percentage_cycles[name_algorithm] = float(line_array[7].split("unici: ")[1])/400

            # Closing files
            file1.close()
    
            # print(counting_dict_max_coverage)

            for x in counting_dict_max_coverage.keys():
                print(x, "| max_coverage:", float(counting_dict_max_coverage[x])/10, "| max_profit:",float(counting_dict_max_profit[x])/10, "| percentage cycles:",float(counting_dict_percentage_cycles[x])/10 )
            print()

            # for x in counting_dict_max_coverage.keys():
            #     print(float(counting_dict_max_coverage[x])/10)
            
            print("-------------------------------------")
