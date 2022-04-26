
#Import Packages 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Code for Gaussian equation
def gaus2d(x, y, mx, my, s):
    return  np.exp(-((x - mx)**2. / (2. * s**2.) + (y - my)**2. / (2. * s**2.)))

#Code for plot
x = np.linspace(0, 1,200)
y = np.linspace(0, 1,200)
x, y = np.meshgrid(x, y) # get 2D variables instead of 1D

#Code to generate random numpy seeds z_t
z_t=np.zeros((200,200))
np.random.seed(164)

#Data frame for Layout
dataf_Layout = pd.read_csv('Layout.txt', delimiter = "\t" )

np.shape(dataf_Layout)
for col in dataf_Layout.columns:
    print(col) 
print(dataf_Layout)

#data frame for TCGA
dataf_TCGA = pd.read_csv('TCGA.csv')
print(dataf_TCGA)

#merge two dataframes
Transposed_dataf_TCGA = dataf_TCGA.transpose()
Transposed_dataf_TCGA["sampleID"]=Transposed_dataf_TCGA.index
#Transposed_dataf_TCGA.columns()
Transposed_dataf_TCGA=Transposed_dataf_TCGA.iloc[1: , :]

for col in Transposed_dataf_TCGA.columns:
    print(col) 

print(Transposed_dataf_TCGA)

Merged_TCGA=pd.merge(dataf_Layout, Transposed_dataf_TCGA, how="inner",left_on="sampleID", right_on="sampleID")
finalTransposed_dataf_TCGA=Merged_TCGA.iloc[: , :]
print(finalTransposed_dataf_TCGA)

Layout1=finalTransposed_dataf_TCGA.iloc[:,0:3]
print(Layout1)
#Code to read expression from final transposed csv
amp=finalTransposed_dataf_TCGA.iloc[0:,4:393]
print(amp)
amp.columns=dataf_TCGA["sampleID"]
print (amp)


def GT_Create(name, amp):      
#Code to create heatmap
    #amp1=amp[name]
    xx_Layout=Layout1["X"]
    yy_Layout=Layout1["Y"]
    for i in range(0,len(xx_Layout)):
        z = gaus2d(x, y,xx_Layout[i],1-yy_Layout[i],.05)
        z_t=z_t+amp[i]*z
    fig, ax = plt.subplots()
    ax.pcolormesh(x, y, z_t, shading='gouraud', cmap=plt.cm.jet)
    plt.show()

for names in amp.columns:
    GT_Create(names, amp[names])
    
#print(amp)
#Create function for gene terrain then put function in loop
