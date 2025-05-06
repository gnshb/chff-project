#!/bin/bash

# Create a data file from the provided arrays
cat <<EOF > data.dat
# Index  MassFlow  u       p
0   0.4472  4.4721  -0.0000
1   0.4472  4.4721  -0.0000
2   0.4472  4.4721  -0.0000
3   0.4472  4.4721  -0.0000
4   0.4472  4.4721   0.0000
5   0.4472  4.4721   0.0000
6   0.4472  4.4721   0.0001
7   0.4472  4.4721   0.0001
8   0.4472  4.4721   0.0001
9   0.4472  4.4721   0.0001
10  0.4472  4.4721   0.0001
11  0.4472  4.4721   0.0001
12  0.4472  4.4721   0.0001
13  0.4472  4.4721   0.0001
14  0.4472  4.4721   0.0001
15  0.4472  4.4721   0.0001
16  0.4472  4.4721   0.0001
17  0.4472  4.4721   0.0000
18  0.4472  4.4721   0.0000
19  0.4472  4.4721   0.0000
EOF

# Gnuplot script
gnuplot -persist <<-EOFMarker
    set title "Flow Quantities vs Index"
    set xlabel "Index"
    set ylabel "Value"
    set key outside
    plot "data.dat" using 1:2 with lines title "Mass Flow Rate", \
         "data.dat" using 1:3 with lines title "Velocity u", \
         "data.dat" using 1:4 with lines title "Pressure p"
EOFMarker

