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
print (table2.colnames)
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
            t1Startingindex = checkSNexist(t1, t2['SN'][i])
            if (t1Startingindex > 0):
                t1index = find2MassStarindex(t1Startingindex,t1,i,t2)
   	        #the position of the SN and star in t1
   	        if (validmagnitude(bandname,t1index,t1,i,t2) and t1index > 0):
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
def find2MassStarindex(index1,t1,index2,t2):
	i1 = index1
	i2 = index2
	list1 = t1['Star']
	while (list1[i1+ 1] == list1[i1] + 1):
            a = lambda x,y: abs(t2['RAs'][i2] - t1['RAs'][i1]) <= 1
            b = lambda x,y: abs(t2['DEs'][i2] - t1['DEs'][i1]) <= 1
            if (a(i2,i1) and b(i2,i1)):
                g = lambda x,y: t2['RAm'][x] == t1['RAm'][y]
                h = lambda x,y: t2['DEm'][x] == t1['DEm'][y]
                if (g(i2,i1) and h(i2,i1)):
                    r = lambda x,y: t2['RAh'][x] == t1['RAh'][y]
                    k = lambda x,y: t2['DE-'][x] == t1['DE-'][y]
                    l = lambda x,y: t2['DEd'][x] == t1['DEd'][y]
                    if (r(i2,i1) and k(i2,i1) and l(i2,i1)):
                        return i1
            i1 = i1 +1
        return -1
def validmagnitude(bandname,t1index,t1,i,t2):
    if (t1[bandname][t1index] < 0):
        return False
    elif (t2[bandname][i] > 25 or t2[bandname][i] < 0 ):
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
#print KDict
KList = Dict2List(KDict)
#print KList
XK = [x for [x,y] in KList]
#print XK
YK = [y for [x,y] in KList]
print "K Band %d stars"%len(KList)
plt.figure(1)
plt.scatter(XK,YK,s=20)
plt.xlabel("PAIRITEL[mag]")
plt.ylabel("CSP - PAIRITEL[mag]")
plt.axis([11,17,-0.3,0.3])
plt.title("K band")
plt.show()
JDict = FindOverLappingSN(table1,[table2,table3],"Jmag")
JList = Dict2List(JDict)
print "J Band %d stars"%len(JList)
XJ = [x for [x,y] in JList]
#print XK
YJ = [y for [x,y] in JList]
plt.figure(2)
plt.scatter(XJ,YJ)
plt.xlabel("PAIRITEL[mag]")
plt.ylabel("CSP - PAIRITEL[mag]")
plt.axis([11,17,-0.3,0.3])
plt.title("J band")
plt.show()
HDict = FindOverLappingSN(table1,[table2,table3],"Hmag")
HList = Dict2List(HDict)
print "H Band %d stars"%len(HList)
XH = [x for [x,y] in HList]
#print XK
YH = [y for [x,y] in HList]
plt.figure(3)
plt.scatter(XH,YH)
plt.xlabel("PAIRITEL[mag]")
plt.ylabel("CSP - PAIRITEL[mag]")
plt.axis([11,17,-0.3,0.3])
plt.title("H band")
plt.show()



