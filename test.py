import lhsmdu
import numpy as np
from bokeh.models import ColumnDataSource
from scipy.stats import qmc

sampler = qmc.LatinHypercube(d=2)
Lcube = sampler.random(10)

print(Lcube)
print(type(Lcube))

print(Lcube[0])
print(Lcube[0,0])

print(Lcube[5,0])

