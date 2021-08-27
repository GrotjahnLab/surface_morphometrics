import click
import mrcfile
import numpy as np
import pandas


@click.command()
@click.argument('input', type=str)
@click.argument('output', type=str)
@click.option('-l', '--label', type=int, default=1)
def convert(input, output, label):
	mrc = mrcfile.mmap(input, mode="r+", permissive=True)
	voxel_size = mrc.voxel_size.x/10 # nm
	# voxel_size = 1
	print(voxel_size)
	data = np.where(mrc.data == label)
	print(data)
	df = pandas.DataFrame(data={'x': data[2], 'y': data[1], 'z': data[0]})
	df = df * voxel_size 
	df.to_csv(output, sep=" ", index=False, header=False)




if __name__=="__main__":
	convert()