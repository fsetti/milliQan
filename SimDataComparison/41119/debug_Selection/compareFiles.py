newFile     = set(line.strip() for line in open('newSelection.txt')) 
oldFile     = set(line.strip() for line in open('retrievePastScripts/oldSelection.txt')) 

print(len(newFile))
print(len(oldFile))
f = open("outEventsNewMinusOld.txt","w")
for line in newFile.difference(oldFile):
	f.write(line + "\n")
f.close()

f1 = open("outEventsOldMinusNew.txt","w")
for line in oldFile.difference(newFile):
	f1.write(line + "\n")
f1.close()
#print(myFile.difference(bennetFile))
