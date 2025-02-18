s = set()
r = set()

for line in open("cells.out"):
    op = line[0]
    c = line[2:].strip()
    if op == 'r':
        r.add(c)
    if op == 's':
        s.add(c)

for c in s-r:
    print(c)

print(len(s))
print(len(r-s))
print(len(s-r))


