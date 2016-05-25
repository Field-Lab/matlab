#!/usr/bin/python

import csv
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# plt.style.use('bmh')

with open('ARVOfiguredata', 'r') as csvfile:
  spamreader = csv.reader(csvfile, delimiter=',')
  shape = []
  color = []
  color_names = ['red', 'green', 'purple', 'black']

  x = []
  y = []
  z = [];
  for row in spamreader:
    shape.append(int(row[0]))
    color.append(int(row[1]))
    x.append(float(row[2]))
    y.append(float(row[3]))
    z.append(float(row[4]))
fig, ax = plt.subplots()
markers = [u'D', u's', u'|', u'x', u'^', u'd', u'h', u'+', u'*', u'D', u'v',  u'>', u'p', u'<', u'8', u'2', u'H', u'o', u'4', u'3', u'1']

# for m in Line2D.markers:
#     try:
#         if len(m) == 1 and m != ' ':
#             markers.append(m)
#     except TypeError:
#         pass
#
# print markers

legend = color_names[:]

plot_color = map(lambda i:color_names[i-1],color)
plot_shape = map(lambda i:markers[i-1], shape)
for i in range(len(x)):
  shape_i = plot_shape[i]
  h=plt.scatter(x[i], y[i], c=plot_color[i], edgecolors = 'face', marker =shape_i, s=15)
  if shape_i == 'o':
    legend[color_names.index(plot_color[i])] = h

plt.legend(legend, ['ON type 1', 'ON type 2', 'OFF type 1', 'OFF type 2'], scatterpoints=1)

#plt.show();

plt.gca().set_ylim([-3, 2])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.grid(b=False)
ax.set_axis_bgcolor('white')
plt.xlabel('PC 1')
plt.ylabel('PC 2')
plt.title('They cluster!!')
# plt.subplots_adjust(left=0.13, bottom=0.18)
plt.savefig('ARVOabstract.pdf')


