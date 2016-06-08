from astropy.table import Table, Column, MaskedColumn
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
path1 = "Friedman/table5.dat"
readmepath1 = "Friedman/ReadMe"
table1 = Table.read(path1, readme = readmepath1, format = "ascii.cds")
#table1 comes from Table 5 of Friedman et al
print (table1.colnames)
table1.pprint()
path2 = "Contreras/table3.dat"
readmepath2 = "Contreras/Readme"
table2 = Table.read(path2, readme = readmepath2, format = "ascii.cds")
table2['Kmag'].fill_values = "invalid"
print (table2.filled())
table2 = table2.filled()
#(table2['Kmag'].filled()).pprint(max_lines = -1)
#table2 comes from Table 3 of Contreras et al 2010 
path3 = "Stritzinger/table3.dat"
readmepath3 = "Stritzinger/Readme"
table3 = Table.read(path3, readme = readmepath3, format = "ascii.cds")
print (table3.colnames)
#table1['SN'].pprint(max_lines = -1)
table3['Kmag'].fill_values = "invalid"
print (table3.filled())
table3 = table3.filled()

#table3 comes from Table 3 of Stritzinger et al 2011
def FindOverLappingSN(t1,ListofCSP,bandname):
    CSPlen = len(ListofCSP)
    diffD = {}
    #ListofCSP: a list of CSP paper(tables) that 
    #contains supernovaes in Friedman et al
    # dictionaries of supernovaes that exist in both CSP and CfAIR2
    for i in range(CSPlen):
        diffD.update(findoverlappingSN(t1,ListofCSP[i],bandname))
    return diffD
def findoverlappingSN(t1,t2,bandname):
	length = len(t2['SN'])
	diffD = {}
   	for i in range(length):
   		if (checkSNexist(t1, t2['SN'][i]) != -1):
   			t1Startingindex = checkSNexist(t1, t2['SN'][i])
   			if (find2MassStarindex(t1['Star'],t1Startingindex,t2['Seq'][i]) != -1):
   				t1index = find2MassStarindex(t1['Star'],t1Startingindex,t2['Seq'][i])
   				#the position of the SN and star in t1
   				if (validmagnitude(bandname,i,t1index,t1,t2)):
                                    #print 3
                                    key = (t1['SN'][t1index][2:],t1['Star'][t1index])
                                    diff = t2[bandname][i] - t1[bandname][t1index]
                                    ele = [t1[bandname][t1index],diff]
                                    new = {key:ele}
                                    diffD.update(new)
        return diffD
def checkSNexist(table,name):
    length = 3397
    for i in range(length):
        if (hasSamename(name,table['SN'][i])):
            return i
    return -1
def hasSamename(name1,name2):
    if ((ord(name1[-1]) != ord(name2[-1]))):
        return False
    elif (ord(name1[-1]) >= 97):
    	return (name1[-6:] == name2[-6:])
    else:
    	return (name1[-5:] == name2[-5:])
def find2MassStarindex(list1,index1,number):
	index = index1
	while (list1[index + 1] == list1[index] + 1):
		if (number == list1[index]):
			return index
		index = index + 1
	return -1
def validmagnitude(bandname,i,t1index,t1,t2):
   	if (t1[bandname][t1index] < 0):
            return False
   	elif (t2[bandname][i] > 25 ):
            return False
   	else:
            return True
#print(FindOverLappingSN(table1,[table2,table3],"Kmag"))
#print(FindOverLappingSN(table1,[table2,table3],"Hmag"))
#print(FindOverLappingSN(table1,[table2,table3],"Jmag"))
def Dict2List(d):
	new = []
	for key in d:
            new.append(d[key])
	return new
KDict = FindOverLappingSN(table1,[table2,table3],"Kmag")
KList = Dict2List(KDict)
#print len(KList)
print KList
XK = [x for [x,y] in KList]
print XK
YK = [y for [x,y] in KList]
plt.figure(1)
plt.scatter(XK,YK)
plt.xlabel("PAIRITEL[mag]")
plt.ylabel("CSP - PAIRITEL[mag]")
plt.title("K band")
plt.show()
JDict = FindOverLappingSN(table1,[table2,table3],"Jmag")
JList = Dict2List(JDict)
print JList
XJ = [x for [x,y] in JList]
print XK
YJ = [y for [x,y] in JList]
plt.figure(2)
plt.scatter(XJ,YJ)
plt.xlabel("PAIRITEL[mag]")
plt.ylabel("CSP - PAIRITEL[mag]")
plt.title("J band")
plt.show()
HDict = FindOverLappingSN(table1,[table2,table3],"Hmag")
HList = Dict2List(HDict)
print HList
XH = [x for [x,y] in HList]
print XK
YH = [y for [x,y] in HList]
plt.figure(3)
plt.scatter(XH,YH)
plt.xlabel("PAIRITEL[mag]")
plt.ylabel("CSP - PAIRITEL[mag]")
plt.title("H band")
plt.show()



