fname = "full_state_list.txt"

fin = open(fname)
fins = fin.readlines()
fin.close()

fout = open("state_list.txt", "w")
for line in fins:
    if not ("jpsi" in line or "psi2" in line):
        line = '#'+line;
    fout.write(line)
fout.close()
