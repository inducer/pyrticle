set terminal postscript eps
set output "kv-plot.eps"
set title "Kapchinskij-Vladimirskij Beam Evolution, no space charge"
set ylabel "Beam Radius [m]"
set xlabel "s [m]"
plot "particle-r-t.dat" title "Beam Radius (simulation, 10000 particles)" with lines, \
  "particle-r-t-theory.dat" title "Beam Radius (theory)" with lines
#pause mouse "Click the mouse"
