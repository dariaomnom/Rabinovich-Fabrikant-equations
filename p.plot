set xlabel "X"
set ylabel "Y"
set zlabel "Z"
set title "Trajectory"

splot 'data.txt' using 1:2:3 with lines title "Rabinovich–Fabrikant attractor"

pause -1