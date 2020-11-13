#!/bin/python
import random

def get_U():
        return sum(2**-i * random.getrandbits(1) for i in range(1, 25))

C = 362436.0 / 16777216.0
CD = 7654321.0 / 16777216.0
CM = 16777213.0 / 16777216.0
IP = 97
JP = 33
with open('fort.8', 'w') as f:
        f.write(' '.join([str(get_U()) for i in range(97)]))
        f.write('\n')
        f.write('{} {} {} {} {}'.format(C, CD, CM, IP, JP))

