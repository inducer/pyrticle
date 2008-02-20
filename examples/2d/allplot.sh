#! /bin/sh


for dir in meth*1753*; do
  id=`echo $dir | cut -f2-3 -d-`
  cmdline="$cmdline $dir/2d.dat prefix $id: plot t_sim,$1"
done
logtool --legend-expr $cmdline
