#!/bin/bash
#
# plot output using gnuplot
#
# first run jonswap and output to out.txt
#
# Note that to plot 2 line in same plot, need to save
# output to file. But for 1 plot only, can just pipe it.
#

start=$(awk '/Spectrum/{ print NR; exit }' out.txt)
end=$(awk '/Time History/{ print NR; exit }' out.txt)

# spectra
sed -n "$((2+${start})),$((-1+${end}))p" out.txt > tmp.txt
gnuplot -p -e "set terminal dumb size $COLUMNS $LINES; \
               plot 'tmp.txt' using \"T\":\"PM\" w l, '' using \"T\":\"JS\" w l"
rm tmp.txt

# time trace
sed -n "$((2+${end})),999999p" out.txt  | \
gnuplot -p -e "set terminal dumb size $COLUMNS $LINES; \
               plot '-' using \"Time\":\"Elevation\" w l"


