#! /bin/sh
# Instructions: Run this script in an existing screen session.
# Caution: The status line is not automatically updated,
# switch windows once to see all runs.

for rec in advective shape normshape; do
  for push in monomial average; do
    echo $rec-$push
    screen -t "$rec-$push" /bin/bash -l start-computation methodtest-$rec-$push -- \
      python2.5 ../with-charge.py $rec $push
    sleep 1
  done
done
