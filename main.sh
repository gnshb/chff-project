#!/bin/bash

g++ -lm -fopenmp main.cpp;
./a.out | tee ad.txt;
#!/bin/bash

# Input file
INPUT_FILE="ad.txt"

# Normalize whitespace to space
massflow=$(grep -A1 "Mass Flow Rate" "$INPUT_FILE" | tail -n1 | tr -s '\t ' ' ')
u=$(grep -A1 "^u:" "$INPUT_FILE" | tail -n1 | tr -s '\t ' ' ')
p=$(grep -A1 "^p:" "$INPUT_FILE" | tail -n1 | tr -s '\t ' ' ')

# Convert to arrays
IFS=' ' read -r -a mf_array <<< "$massflow"
IFS=' ' read -r -a u_array <<< "$u"
IFS=' ' read -r -a p_array <<< "$p"

# Trim p_array to match mf_array length
n=${#mf_array[@]}
if [[ ${#u_array[@]} -lt $n ]]; then
    echo "Error: 'u' array is shorter than 'Mass Flow Rate'."
    exit 1
fi
if [[ ${#p_array[@]} -gt $n ]]; then
    p_array=("${p_array[@]:0:$n}")
fi

# Create data file
{
  echo "# Index MassFlow u p"
  for ((i=0; i<n; i++)); do
    echo "$i ${mf_array[$i]} ${u_array[$i]} ${p_array[$i]}"
  done
} > plot_data.dat

# Gnuplot plot
gnuplot -persist <<-EOF
    set title "Flow Quantities vs Index"
    set xlabel "Index"
    set ylabel "Value"
    set key outside
    plot "plot_data.dat" using 1:2 with lines title "Mass Flow Rate", \
         "plot_data.dat" using 1:3 with lines title "u", \
         "plot_data.dat" using 1:4 with lines title "p"
EOF

rm ad.txt;
