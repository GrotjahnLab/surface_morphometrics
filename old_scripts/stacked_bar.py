from matplotlib import pyplot as plt                                                        

untreated = [21, 3, 30, 25, 31, 13]
treated = [2, 10, 1, 10, 17, 23]
labels = ["Sparse", "Swollen","Tubular", "Mixed", "Bent Lamellar", "Orthogonal Lamellar"]

untreated = [i/sum(untreated) for i in untreated]
treated = [i/sum(treated) for i in treated]

fig,ax = plt.subplots()

xaxis = ["Untreated", "Treated"]
floor = [0,0]
for index, label in enumerate(labels):
    values = [untreated[index], treated[index]]
    ax.bar(xaxis, values, label=label, bottom=floor)
    floor = [floor[i] + values[i] for i in [0,1]]


ax.legend()
fig.savefig("Stacked_bar.png")