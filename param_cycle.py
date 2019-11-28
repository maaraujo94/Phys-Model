fname = 'param_list.txt'

fin = open(fname)
fins = fin.readlines()
fin.close()
# cycle in lambda
for lam in range(-10, 11, 1):
    fout = open('l_%s.txt'%str(lam),'w')
    for line in fins:
        if " lambda " in line:
            line = line.replace('1.0', str(lam/10.))
        fout.write(line)
    fout.close()
