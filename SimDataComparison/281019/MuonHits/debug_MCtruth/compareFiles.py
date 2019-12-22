myFile     = set(line.strip() for line in open('MCcomparison.txt')) 
bennetFile = set(line.strip() for line in open('events.txt')) 

f = open("outEvtsBenNoFran.txt","w")
for line in bennetFile.difference(myFile):
	f.write(line + "\n")
f.close()

f1 = open("outEvtsFranNoBen.txt","w")
for line in myFile.difference(bennetFile):
	f1.write(line + "\n")
f1.close()
#print(myFile.difference(bennetFile))
