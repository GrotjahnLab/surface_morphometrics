import click
import glob
import mrcfile
import numpy as np
import pandas


@click.command()
@click.argument('input', type=str)
@click.argument('output', type=str)
@click.option('-l', '--label', type=int, default=1, help="Label for feature of interest")
def convert(input, output, label):
	mrc = mrcfile.mmap(input, mode="r+", permissive=True)
	# voxel_size = mrc.voxel_size.x/10 # nm
	voxel_size = 1
	print(voxel_size)
	data = np.where(mrc.data == label)
	print(data)
	df = pandas.DataFrame(data={'x': data[2], 'y': data[1], 'z': data[0]})
	df = df * voxel_size 
	df.to_csv(output, sep=" ", index=False, header=False)


def convert_mitochondria(input_filename):
	labels = ["OMM", "IMM", "ER"]
	base = input_filename.split("_")[0]
	mrc = mrcfile.mmap(input_filename, mode="r", permissive=True)
	# voxel_size = mrc.voxel_size.x/10 # nm
	# voxel_size = 1
	# print(voxel_size)
	for index, label in enumerate(labels):
		data = np.where(mrc.data == index+1)
		# print(data)
		df = pandas.DataFrame(data={'x': data[2], 'y': data[1], 'z': data[0]})
		# df = df * voxel_size 
		output = base+"_{}.xyz".format(label)
		print(output)
		df.to_csv(output, sep=" ", index=False, header=False)

if __name__=="__main__":
	data = glob.glob("*_labels.mrc")
	# data = ["TT9_labels.mrc"]
	for datum in data:
		print(datum)
		convert_mitochondria(datum)
	# convert()