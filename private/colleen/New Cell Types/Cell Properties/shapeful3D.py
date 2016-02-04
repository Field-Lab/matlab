#!/usr/bin/python

#import importlib
#importlib.import_module('mpl_toolkits.mplot3d').__path__
#import matplotlib 
#print importlib.import_module('matplotlib').__version__

import csv
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
#from mpl_toolkits.mplot3D import Axes3D

import importlib
importlib.import_module('mpl_toolkits.mplot3d').Axes3D
# plt.style.use('bmh')

with open('ARVOtest', 'r') as csvfile:
  spamreader = csv.reader(csvfile, delimiter=',')
  shape = []
  color = []
  color_names = ['red', 'green', 'blue', 'black']

  x = []
  y = []
  z = [];
  for row in spamreader:
    shape.append(int(row[0]))
    color.append(int(row[1]))
    x.append(float(row[2]))
    y.append(float(row[3]))
    z.append(float(row[4]))


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

markers = [u'D', u's', u'|', u'x', u'_', u'^', u'd', u'h', u'+', u'*', u'D', u'o', u'.', u'1', u'p', u'3', u'2', u'4', u'H', u'v', u'8', u'<', u'>']
for m in Line2D.markers:
    try:
        if len(m) == 1 and m != ' ':
            markers.append(m)
    except TypeError:
        pass


legend = color_names[:]

plot_color = map(lambda i:color_names[i-1],color)
plot_shape = map(lambda i:markers[i-1], shape)
for i in range(len(x)):
  shape_i = plot_shape[i]
  h=ax.scatter(x[i], y[i], z[i], c=plot_color[i], edgecolors = 'face', marker =shape_i)
  if shape_i == 'o' or shape_i == '^' or shape_i == 'x':
    legend[color_names.index(plot_color[i])] = h

plt.legend(legend, ['ON type 1', 'ON type 2', 'OFF type 1', 'OFF type 2'], scatterpoints=1)

plt.show();
# plt.gca().set_xlim([1/lim,lim])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#ax.xaxis.set_ticks_position('bottom')
#ax.yaxis.set_ticks_position('left')
ax.grid(b=False)
ax.set_axis_bgcolor('white')
plt.xlabel('PC 1')
plt.ylabel('PC 2')
plt.title('They cluster!!')
# plt.subplots_adjust(left=0.13, bottom=0.18)
#plt.savefig('ARVOabstract.pdf')

