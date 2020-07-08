fin = open("../aux_func.C")
auxfunc = fin.readlines()
fin.close()

fout = open("aux_func_print.C", "w")

for line in auxfunc:
    if "all functions that make up the cross-section model used" in line:
        fout.write(line[:-4]+"\n")
        fout.write('"print" means its the version used to store the cosalpha distribution */\n')
    else:
        fout.write(line)
        if "double alpha_pow = " in line:
            fout.write("\n")
            fout.write('  // part of "print" version that stores cosalpha\n')
            fout.write('  ofstream fout;\n')
            fout.write('  fout.open("cosa_scan.txt", ofstream::app);\n')
            fout.write('  fout << par[1] << " " << par[2] << " " << par[3] << " " << cosa << " " << fx1 * fx2 * dsdt * jac << endl;\n')
            fout.write('  fout.close();\n')
fout.close()

fin = open("../fitModel.C")
auxfunc = fin.readlines()
fin.close()

fout = open("fitModel_print.C", "w")

for line in auxfunc:
    if "also contains the minuitFunction wrapper" in line:
        fout.write(line[:-4]+"\n")
        fout.write('"print" means its the version used to store the cosalpha distribution */\n')
    elif "import" in line:
        fout.write('#import "aux_func_print.C"\n')
    else:
        fout.write(line)
fout.close()
