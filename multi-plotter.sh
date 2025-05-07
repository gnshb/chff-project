#!/bin/bash

# Different cases
# A     P    N   SR
# 100  1e17  20  yes
# 100  1e17  20  no
# 10   1e17  20  yes
# 100  1e6   20  yes
# 100  1e6   20  no

./main.sh -S -A 100 -P 1e17 -N 20;
echo "Saved Plot for ./main.sh -S -A 100 -P 1e17 -N 20";
./main.sh -A 100 -P 1e17 -N 20;
 echo "Saved Plot for ./main.sh -A 100 -P 1e17 -N 20";
./main.sh -S -A 10 -P 1e17 -N 20;
 echo "Saved Plot for ./main.sh -S -A 10 -P 1e17 -N 20";
./main.sh -S -A 100 -P 1e6 -N 20;
 echo "Saved Plot for ./main.sh -S -A 100 -P 1e6 -N 20";
./main.sh -A 100 -P 1e6 -N 20;
 echo "Saved Plot for ./main.sh -A 100 -P 1e6 -N 20";
