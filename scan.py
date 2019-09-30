fname = 'param_list.txt'

fin = open(fname)
fins = fin.readlines()
fin.close()
# cycle in b
for i in range(16, 21, 1):
    b = i/10.
    fout = open('b_%s.txt'%str(i),'w')
    for line in fins:
        if " b " in line:
            line = line.replace('1.8', str(b))
        fout.write(line)
    fout.close()

# cycle in c
for i in range(4, 9, 1):
    c = i
    fout = open('c_%s.txt'%str(i),'w')
    for line in fins:
        if " c " in line:
            line = line.replace('6', str(c))
        fout.write(line)
    fout.close()
