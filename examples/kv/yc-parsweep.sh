#! /bin/bash

for rsubdiv in 6 10 12; do
  for nparticles in 500 1000 3000 5000; do
    start-computation sweep1-r$rsubdiv-p$nparticles -- csub python2.5 ../with-charge.py \
      --nparticles=$nparticles \
      --radial-subdiv=$rsubdiv
  done
done
