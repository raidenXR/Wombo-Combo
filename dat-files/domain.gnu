set term png size 840,480    
set output '../images/domain.png'
set title 'NTUA Campus'
set xlabel 'latitude'
set ylabel 'longitude'

plot 'coordinates_dense_tabbed.dat' using 2:3 title 'domain' with lines lw 2 lc rgb 'black'
