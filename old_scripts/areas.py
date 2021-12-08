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

@click.command()
@click.argument('imm', nargs=2)
@click.argument('omm', nargs=2)
def main(imm,omm):
	fig, ax = plt.subplots()
	areas = []
	for i in [0,1]:
		plt.cla()
		omm_data = pd.read_csv(omm[i])
		omm_area = sum(omm_data.area)
		imm_data = pd.read_csv(imm[i])
		ibm = imm_data[imm_data.omm_dist < 19]
		junction = imm_data[imm_data.omm_dist > 19]
		junction = junction[junction.omm_dist < 40]
		crista = imm_data[imm_data.omm_dist > 40]

		ibm_area = sum(ibm.area)/omm_area
		junction_area = sum(junction.area)/omm_area
		crista_area = sum(crista.area)/omm_area
		print(imm[i], ibm_area, junction_area, crista_area)
		areas.append([ibm_area, junction_area, crista_area])
		ax.bar([0,1,2], [ibm_area, junction_area, crista_area], tick_label=["IBM", "Crista Junction", "Crista Body"], color = colors_light, edgecolor=colors)
		ax.set_ylabel("Relative Surface Area")
		ax.set_xlabel("Component of IMM")
		# ax.set_yscale('log')
		# ax.set_xscale("log")
		# ax.set_xlim()
		ax.set_ylim([0,2])
		# ax.legend()
		filename = imm[i]+".png"
		fig.savefig(filename)

	plt.cla()
	print(areas)
	relative_areas = [areas[1][i]/areas[0][i] for i in [0,1,2]]
	ax.bar([0,1,2], relative_areas, tick_label=["IBM", "Crista Junction", "Crista Body"], color = colors_light, edgecolor=colors)
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

	


if __name__ == "__main__":
	main()