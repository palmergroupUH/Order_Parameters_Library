import numpy as np 
import matplotlib
#matplotlib.use('Agg') 
import matplotlib.colors as colors
import matplotlib.cm as cmx
import colormaps as cmaps
import matplotlib.pyplot as plt 
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator) 
from matplotlib import rcParams 
from matplotlib.ticker import MaxNLocator
import os

fig_width = 3.0  # width in inches 
fig_height = fig_width/1.333   # height in inches 
fig_size =  [fig_width,fig_height] 
params = {'backend': 'Agg',
          'axes.labelsize': 8,
          'axes.titlesize': 8,
		  'font.size': 8, 
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'figure.figsize': fig_size,
          'savefig.dpi' : 600,
          'font.family': 'sans-serif',
          'axes.linewidth' : 0.5,
          'xtick.major.size' : 2,
          'ytick.major.size' : 2,
          'font.size' : 8,
          'svg.fonttype' : 'none',
          'pdf.fonttype' : 42
          }

rcParams.update(params) 

lwidth=0.8
msize=4

 # colormap
n_curves = 11
values = list(range(n_curves))

plt.register_cmap(name='magma', cmap=cmaps.magma)
jet = cm = plt.get_cmap('magma')
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

# Data
# column 1: q, column 2: P(q) 
q_3375k=np.loadtxt('q_normalized_3375.txt')
q_3400k=np.loadtxt('q_normalized_3400.txt')
q_3425k=np.loadtxt('q_normalized_3425.txt')
q_3450k=np.loadtxt('q_normalized_3450.txt')
q_3500k=np.loadtxt('q_normalized_3500.txt')
q_3550k=np.loadtxt('q_normalized_3550.txt')
q_3600k=np.loadtxt('q_normalized_3600.txt')
q_3700k=np.loadtxt('q_normalized_3700.txt')
q_3800k=np.loadtxt('q_normalized_3800.txt')

# Two vertical meadian q lines: 3350k and 3800k
#median_q = [0.646,0.603] 
median_q = [0.643786707,0.603 ] 

x=np.linspace(0,1,200) 
y=np.linspace(0,6.0,200)


#Colors
col=list(range(n_curves))
col[0] = scalarMap.to_rgba(values[0])
col[1] = scalarMap.to_rgba(values[1])
col[2] = scalarMap.to_rgba(values[2])
col[3] = scalarMap.to_rgba(values[3])
col[4] = scalarMap.to_rgba(values[4])
col[5] = scalarMap.to_rgba(values[5])
col[6] = scalarMap.to_rgba(values[6])
col[7] = scalarMap.to_rgba(values[7]) 
col[8] = scalarMap.to_rgba(values[8]) 
col[9] = scalarMap.to_rgba(values[9]) 

#Line styles
#lstyle = [':', '-.', '--', '-', '--', '-.',':']
lstyle = ['_', '-', '-', '-', '-', '-','-','-','-','-']

#Labels
plot_ID  = ['3350 K', '3375 K','3400 K','3425 K','3450 K','3500 K','3550 K','3600 K','3700 K','3800 K']  

# Plot P(q) vs q: 
#plt.plot(q_3350k[:,0],q_3350k[:,1],label=plot_ID[0],color=col[0],linewidth=lwidth, linestyle = lstyle[0]) # plot only q 0 to 1 
plt.plot(q_3375k[:,0],q_3375k[:,1],label=plot_ID[1],color=col[1],linewidth=lwidth, linestyle = lstyle[1])	
plt.plot(q_3400k[:,0],q_3400k[:,1],label=plot_ID[2],color=col[2],linewidth=lwidth, linestyle = lstyle[2]) 
plt.plot(q_3425k[:,0],q_3425k[:,1],label=plot_ID[3],color=col[3],linewidth=lwidth, linestyle = lstyle[3])
plt.plot(q_3450k[:,0],q_3450k[:,1],label=plot_ID[4],color=col[4],linewidth=lwidth, linestyle = lstyle[4])
plt.plot(q_3500k[:,0],q_3500k[:,1],label=plot_ID[5],color=col[5],linewidth=lwidth, linestyle = lstyle[5])
plt.plot(q_3550k[:,0],q_3550k[:,1],label=plot_ID[6],color=col[6],linewidth=lwidth, linestyle = lstyle[6])
plt.plot(q_3600k[:,0],q_3600k[:,1],label=plot_ID[7],color=col[7],linewidth=lwidth, linestyle = lstyle[7]) 
plt.plot(q_3700k[:,0],q_3700k[:,1],label=plot_ID[8],color=col[8],linewidth=lwidth, linestyle = lstyle[8]) 
plt.plot(q_3800k[:,0],q_3800k[:,1],label=plot_ID[9],color=col[9],linewidth=lwidth, linestyle = lstyle[9]) 



# Plot vertical lines for two medians q
x[:] = median_q[0] #3350 k  
plt.plot(x,y,label='_no_legend_',color=col[0],linewidth=lwidth,linestyle = '--') # plot median of 3350 k 
x[:] = median_q[1] #3800 k 
plt.plot(x,y,label='_no_legend_',color=col[9],linewidth=lwidth, linestyle = '--') # plot median of 3800 k 

# Set axis properties
ax = plt.subplot(111)
minorLocator =  MultipleLocator(0.25) 

majorLocator = MultipleLocator(0.1) 

ax.yaxis.set_minor_locator( minorLocator) 
ax.xaxis.set_minor_locator( majorLocator)

ax.tick_params(axis='x',direction='in', pad=2.0)
ax.tick_params(axis='y',direction='in', pad=2.0)
ax.tick_params(axis='y',which='minor',length=1)
ax.tick_params(axis='x',which='minor',length=1)

plt.xlim([0,1])
plt.ylim([0,2])
plt.xticks([0.1,0.3,0.5,0.7,0.9])
plt.yticks([0,1,2]) 

# Set legend properties
handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1],labels[::-1],loc=(0.02,0.35),fontsize=7,frameon=False, labelspacing=0.07,ncol=1)

#plt.legend(loc=(0.1,0.385),fontsize=7,frameon=False, labelspacing=0.15,ncol=1)

#plt.legend(loc=(0.035,0.35),fontsize=7,frameon=False, labelspacing=0.15,ncol=2)
plt.text(0.025,0.925,"(a)", transform=ax.transAxes)

plt.subplots_adjust(left=0.15, bottom=0.15, right=0.99, top=0.95, wspace=0.0, hspace=0.0)

#plt.subplots_adjust(left=0.15, bottom=0.15, right=0.99, top=0.99, wspace=0.0, hspace=0.0)

plt.xlabel(r'$q$',labelpad=-1)
plt.ylabel(r'$P(q)$',labelpad=0) 

plt.savefig('fig2a.pdf',transparent=True)
plt.savefig('fig2a.eps',transparent=True)
plt.savefig('fig2a.png',transparent=False)

plt.show()

