import click
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

from cycler import cycler
purple = [.5, 0, .5, 0.9]
purple_light = [.5,0,.5, 0.3]
green = [0, 1, 0, 0.9]
green_light = [0,1,0,.3]
blue = [1,.3,0, .9]
blue_light = [1,.3,0,.3]

colors = [purple,green, blue]
colors_light = [purple_light, green_light, blue_light]
labels = ["untreated", "treated"]
@click.command()
@click.argument('imm', nargs=2)
@click.argument('omm', nargs=2)
def main(imm,omm):
	fig, ax = plt.subplots()
	
	for i in [0,1]:
		omm_df = pd.read_csv(omm[i])
		imm_df = pd.read_csv(imm[i])
		# imm_df = imm_df[imm_df.omm_dist > 19]
		imm_df = imm_df[imm_df.omm_dist > 40]
		print(imm_df)
		neighbor_x = omm_df.n_v_x.values[imm_df.omm_neighbor_index]
		neighbor_y = omm_df.n_v_y.values[imm_df.omm_neighbor_index]
		neighbor_z = omm_df.n_v_z.values[imm_df.omm_neighbor_index]

		vectors = np.array([imm_df.n_v_x.values,imm_df.n_v_y.values,imm_df.n_v_z.values]).transpose()
		neighbor_vectors = np.array([neighbor_x,neighbor_y,neighbor_z]).transpose()
		angles =[np.arccos(np.abs(np.dot(vectors[i], neighbor_vectors[i])))*180/np.pi for i in range(len(imm_df))]
		print(max(angles))
		weights = imm_df.area.values
		ax.hist(angles,weights = weights,label=labels[i], ec=colors[i], fc=colors_light[i],histtype="stepfilled", bins=30, density=True, range = [0,90])

	ax.set_ylabel("Relative Area")
	ax.set_xlabel("Angle Relative to the Nearest OMM (Degrees)")
	# ax.set_yscale('log')
	# ax.set_xscale("log")
	# ax.set_xlim()
	ax.set_xticks([0,15,30,45,60,75,90])
	# ax.set_ylim([0.001,100])
	ax.legend()
	filename="relative_angles.png"
	# plt.cla()
	# print(areas)
	# relative_areas = [areas[1][i]/areas[0][i] for i in [0,1,2]]
	# ax.bar([0,1,2], relative_areas, tick_label=["IBM", "Crista Junction", "Crista Body"], color = colors_light, edgecolor=colors)
	# ax.set_ylabel("Relative Area")
	# ax.set_xlabel("IMM-Nearest OMM Angle")
	# ax.axhline(1, ls="-", lw="0.1")
	# # ax.set_yscale('log')
	# # ax.set_xscale("log")
	# # ax.set_xlim()
	# # ax.legend()
	# ax.set_ylim([0,4 ])
	# filename = "fold_change_area.png"
	fig.savefig(filename)

	


if __name__ == "__main__":
	main()