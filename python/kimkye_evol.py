import sys
import numpy as np

def map_kimkye(xn, p):

	alpha, beta, n = p
	xn1 = np.array(xn, dtype=float)

	for _ in range(n):
		x_next = np.array([0, 0], dtype=float)
		x_next[0] = (4 * alpha * xn1[0] + beta * xn1[1]) * (1 - xn1[0])
		x_next[1] = (4 * alpha * xn1[1] + beta * xn1[0]) * (1 - xn1[1])
		xn1 = x_next

	return xn1

def round_row(row, decimal_places):
	return tuple(map(lambda x: round(x, decimal_places), row))

def reinjection(xn, xn1, xfixed, c):

	if (xn > xfixed + c or xn < xfixed - c):
		if (xfixed - c < xn1 < xfixed + c):
			return True

	return False

x0 = np.random.uniform(0, 1, 2)

writefiles = sys.argv[1]

with open(writefiles, 'r') as f:

	lines = f.readlines()

	for i, line in enumerate(lines):
		line = line.strip().replace("'$'\n'", "")
		if i == 0:
			filename_evol = line
		elif i == 1:
			filename_fixedpoints = line
		elif i == 2:
			filename_reinj = line
		elif i == 3:
			filename_rpd = line
		elif i == 4:
			filename_M = line


N = 50000
transient = 50000

xn = x0
xn1 = np.zeros(2)

alpha, beta, n = 0.689005, 0.5, 10

p = [alpha, beta, n]

xn1 = map_kimkye(xn, p)

for i in range(transient):
	xn1 = map_kimkye(xn, p)
	xn = xn1

xlast = np.zeros((N, 2))

with open(filename_evol, 'w') as f:

	for i in range(N):
		xn1 = map_kimkye(xn, p)
		f.write(f'{i} {xn[0]} {xn[1]} {xn1[0]} {xn1[1]}\n')
		xlast[i] = xn1
		xn = xn1

fixedpoints = set(round_row(row, 3) for row in xlast.tolist())
print(fixedpoints)

fp_arr = np.array(list(fixedpoints))

with open(filename_fixedpoints, 'w') as f:
	for fpi in fp_arr:
		f.write(f'{fpi[0]} {fpi[1]}')

fp = fp_arr[0]

N2 = 5000000
rtarget = 5000
rcount = 0
xreinjected = np.zeros(rtarget)
c = 5e-3

p = [0.689, beta, n]

with open(filename_reinj, 'w') as f:
	for i in range(N2):
		xn1 = map_kimkye(xn, p)
		
		if (reinjection(xn[0], xn1[0], fp[0], c)):
			xreinjected[rcount] = xn1[0]
			f.write(f'{xreinjected[rcount]}\n')
			rcount+=1
			if (rcount >= rtarget):
				break

		xn = xn1

xbins = np.linspace(fp[0] - c, fp[0] + c, 100)

hist, bins = np.histogram(xreinjected, bins=xbins)

with open(filename_rpd, 'w') as f:

	for i in range(len(hist)):
		f.write(f'{bins[i] - fp[0] + c} {hist[i]}\n')


xreinj_sort = np.sort(xreinjected)

with open(filename_M, 'w') as f:
	Mi = 0
	for i in range(1, len(xreinj_sort)):
		Mi += xreinj_sort[i]

		f.write(f'{xreinj_sort[i] - fp[0] + c} {Mi / i}\n')