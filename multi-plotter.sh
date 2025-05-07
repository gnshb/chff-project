#!/bin/bash

# Different cases
# A     P    N   SR
# $area  $pressure  $gridpoints  yes
# $area  $pressure  $gridpoints  no
# 10   $pressure  $gridpoints  yes
# $area  1e6   $gridpoints  yes
# $area  1e6   $gridpoints  no

while getopts "A:P:N:" opt; do
  case ${opt} in
    A)
      area="$OPTARG"
      ;;
    P)
      pressure="$OPTARG"
      ;;
    N)
      gridpoints="$OPTARG"
      ;;
    *)
      echo "Usage: $0 -A area_of_inlet/area_of_outlet -P pressure_inlet"
      echo "Optional: -S for Special Relativistic Calculation"
      exit 1
      ;;
  esac
done

rm ex.txt;

./main.sh -S -A $area -P $pressure -N $gridpoints >> ex.txt;
./main.sh -A $area -P $pressure -N $gridpoints >> ex.txt;
python analytic.py --SR -A $area -P $pressure -N $gridpoints >> ex.txt;
python analytic.py -A $area -P $pressure -N $gridpoints >> ex.txt;

python compare.py;

rm ex.txt;
