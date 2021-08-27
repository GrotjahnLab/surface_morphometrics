import click
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

from cycler import cycler
purple = [.5, 0, .5, 0.9]
purple_light = [.5,0,.5, 0.3]
green = [0, 1, 0, 0.9]
green_light = [0,1,0,.3]

colors = [purple,green]
colors_light = [purple_light, green_light]

@click.command()
@click.argument('inputs', nargs=-1)
def main(inputs):
	if len(inputs) == 0:
		exit(0)
	fig, ax = plt.subplots()
	curves=[]
	weights=[]
	for input_csv in inputs:
		df = pd.read_csv(input_csv, sep=",", header=0)
		df = df[df['omm_dist'] > 40] 
		# df = df[df['omm_dist'] < 40] 
		verticality = [np.arccos(np.abs(i))*180/np.pi for i in df.n_v_z.values]
		print(df)

		# curve = df.curvedness_VV
		curves.append(list(verticality))
		weight = df.area.values
		weights.append(weight)
	if len(curves) == 1:
		curves = curves[0]
	if len(weights) == 1:
		weights = weights[0]
	names = [i.split("_")[0] for i in inputs]
	for i in range(len(inputs)):
		ax.hist(curves[i],weights = weights[i],label=names[i], ec=colors[i], fc=colors_light[i],histtype="stepfilled", bins=30, density=True, range = [0,90])
	ax.set_ylabel("Relative Area")
	# ax.set_xlabel("Angle Relative to the Growth Plane (Degrees)")
	ax.set_xlabel("IMM curvedness (AU)")
	# ax.set_yscale('log')
	# ax.set_xscale("log")
	ax.set_xlim()
	ax.set_xticks([0,30,60,90])
	# ax.set_ylim([0.001,100])
	ax.legend()
	filename = "_vs_".join(names)+"_crista_verticality.png"
	fig.savefig(filename)


if __name__ == "__main__":
	main()