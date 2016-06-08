import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
length = 94
path = "Friedman/table1.dat"
readmepath = "Friedman/ReadMe"
table = Table.read(path, readme = readmepath, format = "ascii.cds")
mybin = []
start = 0.0025
binSize = 0.005
edge = 0.08
while (start <= edge):
	mybin.append(start)
	start = start + binSize 
Type1AZ = []
def IsNormalSN1a(index):
	if (table['Type'][i] != "Ia" and table['Type'][i] != "Ia?"):
		return False
	elif (table['SN'][i] == "SN2006E" or table['SN'][i] == "SN2006mq"):
		return False
	else:
		return True
for i in range(length):
	if (IsNormalSN1a(i)):
		Type1AZ.append(table['z'][i])
median = np.median(Type1AZ)
maxZ = np.max(Type1AZ)
minZ = np.min(Type1AZ)
print median, maxZ,minZ
fig, ax = plt.subplots()
plt.hist(Type1AZ,bins = mybin,edgecolor = 'b')
plt.xlabel("Heliocentric Redshift Z")
plt.ylabel("Number of CfAIR2 SN 1a")
plt.axis([0.00,0.08,0,25])
majorLocator = MultipleLocator(0.02)
#majorFormatter = FormatStrFormatter("%d")
minorLocator = MultipleLocator(0.005)
ax.xaxis.set_major_locator(majorLocator)
#ax.xaxis.set_major_formatter(majorFormatter)
ax.xaxis.set_minor_locator(minorLocator)
ax.minorticks_on()
plt.show()
