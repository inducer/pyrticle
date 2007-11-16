set yrange [0.0012:0.0015]
set xrange [0:0.1]
set title "6 radial subdiv"
plot \
  "sweep2-r6-p5000-2007-11-15-1725/beam-rad-rms-theory-0.dat" title "theoretical (?)" with lines, \
  "sweep2-r6-p500-2007-11-15-1725/beam-rad-rms-sim.dat" title "6/500", \
  "sweep2-r6-p1000-2007-11-15-1725/beam-rad-rms-sim.dat" title "6/1000", \
  "sweep2-r6-p3000-2007-11-15-1725/beam-rad-rms-sim.dat" title "6/3000", \
  "sweep2-r6-p5000-2007-11-15-1725/beam-rad-rms-sim.dat" title "6/5000"
pause -1

set title "10 radial subdiv"
plot \
  "sweep2-r10-p5000-2007-11-15-1725/beam-rad-rms-theory-0.dat" title "theoretical (?)" with lines, \
  "sweep2-r10-p500-2007-11-15-1725/beam-rad-rms-sim.dat" title "10/500", \
  "sweep2-r10-p1000-2007-11-15-1725/beam-rad-rms-sim.dat" title "10/1000", \
  "sweep2-r10-p3000-2007-11-15-1725/beam-rad-rms-sim.dat" title "10/3000", \
  "sweep2-r10-p5000-2007-11-15-1725/beam-rad-rms-sim.dat" title "10/5000"
pause -1

set title "12 radial subdiv"
plot \
  "sweep2-r10-p5000-2007-11-15-1725/beam-rad-rms-theory-0.dat" title "theoretical (?)" with lines, \
  "sweep2-r12-p500-2007-11-15-1725/beam-rad-rms-sim.dat"  title "12/500", \
  "sweep2-r12-p1000-2007-11-15-1725/beam-rad-rms-sim.dat" title "12/1000", \
  "sweep2-r12-p3000-2007-11-15-1725/beam-rad-rms-sim.dat" title "12/3000", \
  "sweep2-r12-p5000-2007-11-15-1725/beam-rad-rms-sim.dat" title "12/5000"
pause -1

set title "500 particle roundup"
plot \
  "sweep2-r10-p500-2007-11-15-1725/beam-rad-rms-theory-0.dat" title "theoretical (?)" with lines, \
  "sweep2-r6-p500-2007-11-15-1725/beam-rad-rms-sim.dat"  title "6/500", \
  "sweep2-r10-p500-2007-11-15-1725/beam-rad-rms-sim.dat" title "10/500", \
  "sweep2-r12-p500-2007-11-15-1725/beam-rad-rms-sim.dat" title "12/500"
pause -1
set title "1k particle roundup"
plot \
  "sweep2-r10-p1000-2007-11-15-1725/beam-rad-rms-theory-0.dat" title "theoretical (?)" with lines, \
  "sweep2-r6-p1000-2007-11-15-1725/beam-rad-rms-sim.dat"  title "6/1000", \
  "sweep2-r10-p1000-2007-11-15-1725/beam-rad-rms-sim.dat" title "10/1000", \
  "sweep2-r12-p1000-2007-11-15-1725/beam-rad-rms-sim.dat" title "12/1000"
pause -1

set title "3k particle roundup"
plot \
  "sweep2-r10-p3000-2007-11-15-1725/beam-rad-rms-theory-0.dat" title "theoretical (?)" with lines, \
  "sweep2-r6-p3000-2007-11-15-1725/beam-rad-rms-sim.dat"  title "6/3000", \
  "sweep2-r10-p3000-2007-11-15-1725/beam-rad-rms-sim.dat" title "10/3000", \
  "sweep2-r12-p3000-2007-11-15-1725/beam-rad-rms-sim.dat" title "12/3000"
pause -1

set title "5k particle roundup"
plot \
  "sweep2-r10-p5000-2007-11-15-1725/beam-rad-rms-theory-0.dat" title "theoretical (?)" with lines, \
  "sweep2-r6-p5000-2007-11-15-1725/beam-rad-rms-sim.dat"  title "6/5000", \
  "sweep2-r10-p5000-2007-11-15-1725/beam-rad-rms-sim.dat" title "10/5000", \
  "sweep2-r12-p5000-2007-11-15-1725/beam-rad-rms-sim.dat" title "12/5000"
pause -1

! cat "elcounts.dat"
pause -1

! cat "dt.dat"
pause -1

! cat "runtimes.dat"
pause -1
