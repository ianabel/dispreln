
set log
set xrange [.1:100]
set yrange [.1:35]
plot for [i=0:3] 'Howes.out' index i using 1:2 w lines
pause mouse
set yrange [.0001:10]
plot for [i=0:3] 'Howes.out' index i using 1:(-$3) w lines
pause mouse

set yrange [.1:35]
plot for [i=0:3] 'Howes2.out' index i using 1:2 w lines
pause mouse
set yrange [.0001:10]
plot for [i=0:3] 'Howes2.out' index i using 1:(-$3) w lines
pause mouse

