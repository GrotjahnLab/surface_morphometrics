from matplotlib import pyplot as plt
import numpy as np
purple = [.5, 0, .5, 0.9]
purple_light = [.5,0,.5, 0.3]
green = [0, 1, 0, 0.9]
green_light = [0,1,0,.3]
blue = [1,.3,0, .9]
blue_light = [1,.3,0,.3]
grey = [.5,.5,.5, .9]
grey_light = [.5,.5,.5, .3]
colors = [purple,green, blue]
colors_light = [purple_light, green_light, blue_light]
names = ["untreated", "treated"]

fig,ax=plt.subplots()
untreated= [[0.7430620717720843, 0.1927905580937579, 0.5315183448345845],[0.6303341193062824, 0.07734204718307244, 0.3546971097209438]]
treated = [[0.6932474006726265, 0.21935202076276508, 1.7986680791597873], [0.5203904269023357, 0.12982362527362076, 1.4763338740393293]]
untreated = np.array(untreated)
print(untreated)
treated = np.array(treated)
area_set = []
std_set = []
for index, i in enumerate([untreated, treated]):
	plt.cla()
	areas = np.mean(i, axis=0)
	stds = np.std(i, axis=0)
	ax.bar([0,1,2], areas, yerr=stds, tick_label=["IBM", "Crista Junction", "Crista Body"], color = colors_light, edgecolor=colors)
	ax.set_ylabel("Relative Surface Area")
	ax.set_xlabel("Component of IMM")
	area_set.append(areas)
	std_set.append(stds)
	ax.set_ylim([0,2])
	# ax.legend()
	filename = names[index]+".png"
	fig.savefig(filename)


colors_light.append(grey_light)
colors.append(grey)
plt.cla()
area_set = np.array(area_set)
std_set = np.array(std_set)
total_untreated = np.sum(untreated, axis=1)
total_treated = np.sum(treated, axis=1)

fold_change_total = np.mean(total_treated)/np.mean(total_untreated)
fold_change_std = np.sqrt((np.std(total_untreated)/np.mean(total_untreated))**2+(np.std(total_untreated)/np.mean(total_untreated))**2)

fold_changes = area_set[1]/area_set[0]
fold_changes = np.append(fold_changes, fold_change_total)
fold_std = np.sqrt((std_set[0]/area_set[0])**2 + (std_set[1]/area_set[1])**2)
fold_std = np.append(fold_std, fold_change_std)

ax.bar([0,1,2,3], fold_changes, yerr=fold_std, tick_label=["IBM", "Crista Junction", "Crista Body", "Total"], color = colors_light, edgecolor=colors)
ax.set_ylabel("Fold Change in Relative Surface Area")
ax.set_xlabel("Component of IMM")
ax.axhline(1, ls="-", lw="0.1")
# ax.set_yscale('log')
# ax.set_xscale("log")
# ax.set_xlim()
# ax.legend()
ax.set_ylim([0,4 ])
filename = "fold_change_area.png"
fig.savefig(filename)