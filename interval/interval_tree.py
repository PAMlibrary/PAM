import time
import random
import sys
from intervaltree import Interval, IntervalTree

if len(sys.argv) != 2:
   print ("Please provide the number of intervals");
   sys.exit(0);

n=int(sys.argv[1]);
random.seed();
intervals=[ Interval(k, k + random.randint(1,n)) for k in range(n) ];

t0=time.time();
t=IntervalTree(intervals);
t1=time.time();

print ("Construction time for n=", n, ": ", t1 - t0);

