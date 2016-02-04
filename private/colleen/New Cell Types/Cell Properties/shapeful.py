#!/usr/bin/python

import csv
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


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
fig, ax = plt.subplots()
#fig.set_size_inches(6,6)
markers = [u'D', u's', u'|', u'x', u'^', u'd', u'h', u'+', u'*', u'D', u'v',  u'>', u'p', u'<', u'8', u'2', u'H', u'o', u'4', u'3', u'1', u'D', u's', u'|', u'x', u'^', u'd', u'h', u'+', u'*', u'D', u'v',  u'>', u'p', u'<', u'8', u'2', u'H', u'o', u'4', u'3', u'1']

# for m in Line2D.markers:
#     try:
#         if len(m) == 1 and m != ' ':
#             markers.append(m)
#     except TypeError:
#         pass
#
# print markers
print color
print shape
legend = color_names[:]

plot_color = map(lambda i:color_names[i-1],color)

plot_shape = map(lambda i:markers[i-1], shape)
for i in range(len(x)):
  shape_i = plot_shape[i]
  h=plt.scatter(x[i], y[i], c=plot_color[i], edgecolors = 'face', marker =shape_i, s=15)
  if shape_i == 'D' or shape_i == 'd' or shape_i == '^' :
    legend[color_names.index(plot_color[i])] = h

plt.legend(legend, ['ON Type 1', 'ON Type 2', 'OFF Type 1', 'OFF Type 2', 'ON LBC', 'ON Type 3', 'ON Type 4'], scatterpoints=1, loc = "upper left", fontsize = 10)

plt.show()
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())

#plt.gca().set_aspect(6/5)
#plt.axis('square')

plt.gca().set_ylim([-2, 1.5])
plt.gca().set_xlim([-3, 1.5])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#ax.set_xticklabels([])

ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

ax.grid(b=False)
ax.set_axis_bgcolor('white')
plt.xlabel('First Principal Component (a.u.)', fontsize = 14)
plt.ylabel('Second Principal Component (a.u.)', fontsize = 14)
plt.legend(loc = 'lower left')
#plt.title('They cluster!!')
# plt.subplots_adjust(left=0.13, bottom=0.18)
#plt.savefig('ARVOfigure_new.png', dpi= 600)


