import sys
import matplotlib.pyplot as plt
import numpy as np

def read_datafile(filename, x_col, y_col):
	data = np.loadtxt(filename)
	x = data[:, x_col]
	y = data[:, y_col]
	return x, y

def plot_from_instructions(instructions_file):
	with open(instructions_file, 'r') as file:
		lines = file.readlines()


	total_axes = 0
	layout = (1, 1)
	figsize = (5, 5)
	axes_instructions = []
	save_file = 'output_plot.pdf'

	current_axes = None
	for line in lines:
		line = line.strip()
		if line.startswith("total_axes:"):
			total_axes = int(line.split(":")[1].strip())
		elif line.startswith("layout:"):
			layout = tuple(map(int, line.split(":")[1].strip().split(',')))
		elif line.startswith("figsize:"):
			figsize = tuple(map(int, line.split(":")[1].strip().split(',')))
		elif line.startswith("axes"):
			if current_axes is not None:
				axes_instructions.append(current_axes)
			current_axes = {"plots": 0, "data": []}
		elif line.startswith("plots:"):
			current_axes["plots"] = int(line.split(":")[1].strip())
		elif line.startswith("file:"):
			data_info = {"file": line.split(":")[1].strip()}
			current_axes["data"].append(data_info)
		elif line.startswith("x_col:"):
			if current_axes["data"]:
				current_axes["data"][-1]["x_col"] = int(line.split(":")[1].strip())
			else:
				print("Error: 'x_col' specified without a preceding 'file'")
		elif line.startswith("y_col:"):
			if current_axes["data"]:
				current_axes["data"][-1]["y_col"] = int(line.split(":")[1].strip())
			else:
				print("Error: 'y_col' specified without a preceding 'file'")
		elif line.startswith("title:"):
			if current_axes["data"]:
				current_axes["data"][-1]["title"] = line.split(":")[1].strip()
			else:
				print("Error: 'title' specified without a preceding 'file'")
		elif line.startswith("xlabel:"):
			if current_axes["data"]:
				current_axes["data"][-1]["xlabel"] = line.split(":")[1].strip()
			else:
				print("Error: 'xlabel' specified without a preceding 'file'")
		elif line.startswith("ylabel:"):
			if current_axes["data"]:
				current_axes["data"][-1]["ylabel"] = line.split(":")[1].strip()
			else:
				print("Error: 'ylabel' specified without a preceding 'file'")
		elif line.startswith("color:"):
			if current_axes["data"]:
				current_axes["data"][-1]["color"] = line.split(":")[1].strip()
			else:
				print("Error: 'color' specified without a preceding 'file'")
		elif line.startswith("dot_size:"):
			if current_axes["data"]:
				current_axes["data"][-1]["dot_size"] = int(line.split(":")[1].strip())
			else:
				print("Error: 'dot_size' specified without a preceding 'file'")
		elif line.startswith("save_file:"):
			save_file = line.split(":")[1].strip()


	if current_axes is not None:
		axes_instructions.append(current_axes)


	fig, axes = plt.subplots(layout[0], layout[1], figsize=figsize)
	if layout[0] == 1 and layout[1] == 1:
		axes = [axes]
	elif layout[0] == 1 or layout[1] == 1:
		axes = axes.flatten()
	else:
		axes = axes.ravel()

	for ax_index, ax_info in enumerate(axes_instructions):
		ax = axes[ax_index]
		for plot_info in ax_info['data']:
			x, y = read_datafile(plot_info['file'], plot_info['x_col'], plot_info['y_col'])
			ax.scatter(x, y, color=plot_info['color'], s=plot_info['dot_size'], label=plot_info['title'])
			ax.set_title(plot_info['title'])
			ax.set_xlabel(plot_info['xlabel'])
			ax.set_ylabel(plot_info['ylabel'])
			ax.legend()

	plt.tight_layout()
	plt.savefig(save_file)	

inst_file = sys.argv[1]
print(inst_file)
plot_from_instructions(inst_file)
