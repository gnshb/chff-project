#!/bin/bash

g++ -lm -fopenmp main_sr.cpp;
./a.out | tee ad.txt;

# Input file
INPUT_FILE="ad.txt"

# Normalize whitespace to space
u=$(grep -A1 "^u:" "$INPUT_FILE" | tail -n1 | tr -s '\t ' ' ')
p=$(grep -A1 "^p:" "$INPUT_FILE" | tail -n1 | tr -s '\t ' ' ')

# Convert to arrays
IFS=' ' read -r -a u_array <<< "$u"
IFS=' ' read -r -a p_array <<< "$p"

nu=${#u_array[@]}
np=${#p_array[@]}

# Create data file for u
{
  echo "# Index u"
  for ((i=0; i<nu; i++)); do
    echo "$i ${u_array[$i]}"
  done
} > u_data.dat

# Create data file for p (full length)
{
  echo "# Index p"
  for ((i=0; i<np; i++)); do
    echo "$i ${p_array[$i]}"
  done
} > p_data.dat

# Gnuplot: 2 plots in one window (side-by-side)
gnuplot -persist <<-EOF
    set key outside
    set style line 1 lt rgb "green" lw 2 pt 7 ps 1.2
    set style line 2 lt rgb "red" lw 2 pt 7 ps 1.2

    set multiplot layout 1,2 title "Velocity and Pressure Plots"

    # Velocity plot
    set title "Velocity (u)"
    set xlabel "Index"
    set ylabel "u"
    plot "u_data.dat" using 1:2 with linespoints linestyle 1 title "u"

    # Pressure plot
    set title "Pressure (p)"
    set xlabel "Index"
    set ylabel "p"
    plot "p_data.dat" using 1:2 with linespoints linestyle 2 title "p"

    unset multiplot
EOF

rm ad.txt u_data.dat p_data.dat
