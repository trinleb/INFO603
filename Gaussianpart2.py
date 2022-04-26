import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def gaus2d(x, y, mx, my, s):
    return  np.exp(-((x - mx)**2. / (2. * s**2.) + (y - my)**2. / (2. * s**2.)))
    
x = np.linspace(0, 1,200)
y = np.linspace(0, 1,200)
x, y = np.meshgrid(x, y) # get 2D variables instead of 1D

#Data frame for Layout
dataf_Layout = pd.read_csv('Layout.txt', delimiter = "\t" )
#Layout=dataf_Layout.sort_values(by=['sampleID'], ascending=True)
print(dataf_Layout)

#data frame for Expression
dataf_Expression = pd.read_csv('Expression.txt', delimiter = "\t" )#,  sep=',')
print(dataf_Expression)

#X and Y from Layout
xx= dataf_Layout['X']
yy=dataf_Layout['Y']

#merge two dataframes

z_t=np.zeros((200,200))
np.random.seed(164)


#amp=abs(5*np.random.randn(11))
###A from Expression, Genes alphabetically arranged###
amp=dataf_Expression['A']


for i in range(len(xx)):
  z = gaus2d(x, y,xx[i],1-yy[i],.05)
  z_t=z_t+amp[i]*z
fig, ax = plt.subplots()
ax.pcolormesh(x, y, z_t, shading='gouraud', cmap=plt.cm.jet)
#im = ax.imshow(z_t,cmap=plt.cm.jet)
plt.show()

#remove by pd.merge

