import numpy as np
import pylab as plt

### Pie chart where the slices are ordered and plotted counter-clockwise

sizes = [100./3, 100./3, 100./3]		# defining 3 slices for the pie chart
text_diagram_1 = r'$P_{rot} < 50$' + ' ' + r'$d$' + '\n' + r'$0.1\% < S_{ph} < 1\%$' + '\n' + r'$0.18 < S_{CaII \, H&K} < 0.3$' + '\n' + r'$0.7 < S_{H\alpha} < 0.75$' + '\n' + r'$0.55 < S_{CaII \, IR} < 0.62$' + '\n' + r'$-0.5 < \Delta m_{NUV} < 1$' 				# legend in the 1st slice
text_diagram_2 = r'$P_{rot} < 50$' + ' ' + r'$d$' + '\n' + r'$S_{ph} > 1\%$' + '\n' + r'$S_{CaII \, H&K} > 0.3$' + '\n' + r'$S_{H\alpha} > 0.75$' + '\n' + r'$S_{CaII \, IR} > 0.62$' + '\n' + r'$\Delta m_{NUV} < -0.5$' + '\nPossibly flares'						# legend in the 2nd slice
text_diagram_3 =  r'$P_{rot} > 50$' + ' ' + r'$d$' + '\n' + r'$0.1\% < S_{ph} < 1\%$' + '\n' + r'$0.18 < S_{CaII \, H&K} < 0.3$' + '\n' + r'$0.7 < S_{H\alpha} < 0.75$' + '\n' + r'$0.55 < S_{CaII \, IR} < 0.62$' + '\n' + r'$-0.5 < \Delta m_{NUV} < 1$'				# legend in the 3rd slice


## Plotting the pie chart

plt.figure(figsize=(7., 7.))
plt.pie(sizes, labeldistance=.5, startangle=90, colors=['b', 'g', 'r'], wedgeprops={"alpha": 0.4})
plt.annotate('Fast-rotating\nsingle RGs', (-0.62, 0.55), weight='bold', fontsize=12)
plt.annotate(text_diagram_1, (-0.75, 0.02), fontsize=12)
plt.annotate('Fast-rotating\nclose binary RGs', (0.1, 0.55), weight='bold', fontsize=12)
plt.annotate(text_diagram_2, (0.2, -0.05), fontsize=12)
plt.annotate('Slowly-rotating\nRGs: single (85%)\nor binary (15%)', (-0.3, -0.4), weight='bold', fontsize=12)
plt.annotate(text_diagram_3, (-0.29, -0.92), fontsize=12)
plt.title('Classification of active red giants', fontsize='x-large')
plt.axis('equal')
plt.savefig('./Pie_chart_classification.pdf', format='pdf')
plt.close()
