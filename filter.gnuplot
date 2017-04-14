set terminal svg enhanced
set output 'filter.svg'
set datafile separator "\t"
set ylabel 'Gain (dB)'
set xlabel 'Frequency (Mhz)'
set grid
plot 'lpf-140MHz.txt' using 1:2 with lines lw 2 lt 3 title 'Filter response'
