#!/bin/bash
#
# plot output using gnuplot
#
# first run jonswap and output to out.txt
#
start=$(awk '/Spectrum/{ print NR; exit }' out.txt)
end=$(awk '/Time History/{ print NR; exit }' out.txt)

# spectra
sed -n "$((4+${start})),$((-1+${end}))p" out.txt > tmp.txt
gnuplot -p -e "plot 'tmp.txt' using 1:3 w l, '' using 1:4 w l"

# time trace
sed -n "$((4+${end})),999999p" out.txt > tmp.txt
gnuplot -p -e "plot 'tmp.txt' using 1:2 w l"

