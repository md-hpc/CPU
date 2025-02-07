import sys

def parse(file):
    return [
        line.strip() for line in open(file).read().strip().split('\n')
        if len(line.strip()) > 0
    ]

d1 = set(parse(sys.argv[1]))
d2 = set(parse(sys.argv[2]))

print('expected')
print(d1 - d2)
print('unexpected')
print(d2 - d1)
