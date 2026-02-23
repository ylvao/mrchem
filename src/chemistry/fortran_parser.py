import re

# def generate_header(input_fortran, output_header):
#     num_pattern = re.compile(r"[-+]?\d*\.\d+[Dd][-+]?\d+|[-+]?\d*\.\d+|[-+]?\d+")

#     with open(input_fortran, 'r') as f:
#         content = f.read()
#         matches = num_pattern.findall(content)
#         values = [m.replace('D', 'E').replace('d', 'e') for m in matches]

#     # Skip metadata (first two numbers) and group into 5s
#     entries = [values[i:i+5] for i in range(2, len(values), 5)]
#     valid_entries = [e for e in entries if len(e) == 5]

#     with open(output_header, 'w') as out:
#         out.write("#ifndef DISPERSION_PARAMS_H\n#define DISPERSION_PARAMS_H\n\n")
#         out.write("struct DispersionEntry {\n    double c6;\n    int z1;\n    int z2;\n    double cn1;\n    double cn2;\n};\n\n")
        
#         # ADD THIS LINE: Define the size for the C++ loop
#         out.write(f"static const int DISPERSION_DATA_SIZE = {len(valid_entries)};\n\n")
        
#         out.write("static const DispersionEntry DISPERSION_DATA[] = {\n")
#         for e in valid_entries:
#             z1 = int(float(e[1]))
#             z2 = int(float(e[2]))
#             out.write(f"    {{{e[0]}, {z1}, {z2}, {e[3]}, {e[4]}}},\n")
#         out.write("};\n\n")
#         out.write("#endif\n")

# generate_header('pars.f', 'dispersion_params.h')
# print("Finished")



def generate_header(input_fortran, output_header):
    # num_pattern = re.compile(r"pars\(\s+\d+:\s+\d+\)=\(\/\n(\s+\.[,\s](0\.\d+D[+-]\d+),(0\.\d+D[+-]\d+),(0\.\d+D[+-]\d+),(0\.\d+D[+-]\d+),(0\.\d+D[+-]\d+)\n)+")
    num_pattern = re.compile(r"(0\.\d+E[+-]\d+),(0\.\d+E[+-]\d+),(0\.\d+E[+-]\d+),(0\.\d+E[+-]\d+),(0\.\d+E[+-]\d+)")

    with open(input_fortran, 'r') as f:
        content = f.read()
        content = content.replace("D", "E")
        matches = re.findall(num_pattern, content)
        print(matches[0])
        

    # # Skip metadata (first two numbers) and group into 5s
    # entries = [values[i:i+5] for i in range(2, len(values), 5)]
    # valid_entries = [e for e in entries if len(e) == 5]

    with open(output_header, 'w') as out:
        out.write("#ifndef DISPERSION_PARAMS_H\n#define DISPERSION_PARAMS_H\n\n")
        out.write("struct DispersionEntry {\n    double c6;\n    int z1;\n    int z2;\n    double cn1;\n    double cn2;\n};\n\n")
        
        # ADD THIS LINE: Define the size for the C++ loop
        string = re.compile(r"nlines=\s+(\d+)")
        n_entries = string.search(content).group(1)
        out.write(f"static const int DISPERSION_DATA_SIZE = {n_entries};\n\n")
        
        out.write("static const DispersionEntry DISPERSION_DATA[] = {\n")
        for m in matches:
            c6 = m[0]
            z1 = int(float(m[1]))
            z2 = int(float(m[2]))
            cn1 = m[3]
            cn2 = m[4]
            
            #     z1 = int(float(e[1]))
            #     z2 = int(float(e[2]))
            out.write(f"    {{{c6}, {z1}, {z2}, {cn1}, {cn2}}},\n")
            #     out.write(f"    {{{e[0]}, {z1}, {z2}, {e[3]}, {e[4]}}},\n")
        out.write("};\n\n")
        out.write("#endif\n")

generate_header('pars.f', 'dispersion_params.h')
print("Finished")





