import sys

def parse(file):
    return [
        line.strip() for line in open(file).read().split('\n')
        if len(line.strip()) > 0
    ]



d1 = parse(sys.argv[1])
d2 = parse(sys.argv[2])

s1 = set()
s2 = set()

for i,d,s in zip(range(2),[d1,d2],[s1,s2]):
    for j,line in enumerate(d):
        if line.count(' ') != 2:
            print('bug',i,j,line)
        if line in s:
            print('dup',i,j,line)
        
        s.add(line)
        
print('expected')
for s in s1 - s2:
    print(f'\t{s}')
print('unexpected')
for s in s2 - s1:
    print(f'\t{s}')
