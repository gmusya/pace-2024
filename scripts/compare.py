#!/bin/python3

import sys
  
# Arguments passed
names = sys.argv[1:]
content = [open(name, "r").read().split('\n') for name in names]

print('\t\t'.join(names))
for i in range(100):
    scores_as_str = [c[i].split(' ')[2] for c in content]
    print('\t\t'.join(scores_as_str))
