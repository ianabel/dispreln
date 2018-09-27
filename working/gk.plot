
# Gnuplot file

set log
set xtics
set ytics
set yrange [.1:35]
set xrange [.1:100]
plot for [i=0:1] 'gk.out' index i using 1:2 w lines
pause mouse 
set yrange [.0001:10]
plot for [i=0:1] 'gk.out' index i using 1:(-$3) w lines
pause mouse


