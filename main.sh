#!/bin/bash

srflag=""
while getopts ":J:SA:P:N:" opt; do
  case ${opt} in
    A)
      area="$OPTARG"
      ;;
    P)
      pressure="$OPTARG"
      ;;
    S)
      sr="sr"
      srflag="--SR"
      ;;
    N)
      gridpoints="$OPTARG"
      ;;
    J)
      plot="true"
      ;;
    *)
      echo "Usage: $0 -A area_of_inlet/area_of_outlet -P pressure_inlet"
      echo "Optional: -S for Special Relativistic Calculation, -J for plotting"
      exit 1
      ;;
  esac
done

if [[ -z "$area" || -z "$pressure" ]]; then
  echo "Error: Both -A and -P arguments are required."
  echo "Usage: $0 -A area_of_outlet/area_of_inlet -P pressure_inlet"
  echo "Optional: -S for Special Relativistic Calculation"
  exit 1
fi

sed -i "s/double Pi = .*;/double Pi = ${pressure};/g" SIMPLE-RF.cpp;
sed -i "s/double Ai = .*;/double Ai = ${area};/g" SIMPLE-RF.cpp;
sed -i "s/int N = .*;/int N = ${gridpoints};/g" SIMPLE-RF.cpp;

if [[ "$sr" == "sr" ]]; then
    sed -i "s/\/\/ #define SR/#define SR/g" SIMPLE-RF.cpp;
else
    :
fi

g++ -lm -fopenmp SIMPLE-RF.cpp;

if [[ "$plot" == "true" ]]; then
	./a.out > ad.txt;
else
	./a.out;
fi

INPUT_FILE="ad.txt"

u=$(grep -A1 "^u:" "$INPUT_FILE" | tail -n1 | tr -s '\t ' ' ')
p=$(grep -A1 "^p:" "$INPUT_FILE" | tail -n1 | tr -s '\t ' ' ')

IFS=' ' read -r -a u_array <<< "$u"
IFS=' ' read -r -a p_array <<< "$p"

nu=${#u_array[@]}
np=${#p_array[@]}

{
  echo "# Index u"
  for ((i=0; i<nu; i++)); do
    echo "$i ${u_array[$i]}"
  done
} > u_data.dat

{
  echo "# Index p"
  for ((i=0; i<np; i++)); do
    echo "$i ${p_array[$i]}"
  done
} > p_data.dat

gnuplot -persist <<-EOF
    set terminal pngcairo size 1200,600 enhanced font 'Verdana,10'
    set output "output_images/SIMPLE_${sr}_A${area}_P${pressure}.png"

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

if [[ "$sr" == "sr" ]]; then
    sed -i "s/#define SR/\/\/ #define SR/g" SIMPLE-RF.cpp;
else
    :
fi
