#
# If using jonswap with the -ps or -pt options, just
# pipe output to this script and pass which column to plot.
#

if [ -z "$1" ]; then
    j=2
else
    j=$1
fi

gnuplot -p -e "set terminal dumb; \
               plot '-' using 1:$j w l"
