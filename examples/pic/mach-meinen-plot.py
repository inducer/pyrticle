import numpy

from matplotlib.pyplot import title, plot, show
for filename in db.q("select filename from runs"):
    data = numpy.array(db.q("select $t_sim,$W_part+$W_el where filename = ?", filename).fetchall())
    title(filename)
    plot(data[:,0], data[:,1])
show()
