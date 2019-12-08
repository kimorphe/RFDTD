%matplotlib inline
import matplotlib.pyplot as plt
from ipywidgets import interact
import numpy as np

def f(k):
    x= np.linspace(0,10,num=1000)
    y=np.sin(k*x)
    plt.plot(x,y)
    plt.show()


interact(f,k=(0.0,5.0,0.1))
