#! /usr/bin/env gnuplot

l=2
f(r,alpha)= (l-r**2/l)**alpha 
plot [r=0:l] f(r, 1), f(r,2), f(r,3)
pause -1 "Hit return to exit"
