import numpy as np
#import matplotlib as mpl
#mpl.use('AGG')
import matplotlib.pyplot as plt
#from matplotlib import gridspec
#import matplotlib.path as mpltPath
#import matplotlib.patches as mpatches
#from mpl_toolkits.basemap import Basemap
from scipy.stats import norm
from scipy.stats import entropy as KL_divergence

plt.close('all')
x = np.linspace(-10,15,num=50)
y = norm.pdf(x, loc = 2, scale = 3)
fig = plt.figure()

for i in range(9):

    y2 = norm.pdf(x, loc = 2+np.random.rand(1)*10, scale = 3+np.random.rand(1)*5)+np.random.rand(1)*0.05*np.random.rand(50)
    ax = fig.add_subplot(3,3,i+1)
    score = KL_divergence(y, y2, base=None)
    ax.plot(x,y,'r-')
    ax.plot(x,y2,'b-')
    ax.set_title(str(score))

