import os

snps = 0
insertions=0
deletions=0
directory="/lustre/nobackup/WUR/ABGC/moiti001/polish-assembly-shortreads/logs_rules/pilon"

for file in os.listdir(directory):
    with open(os.path.join(directory,file), "r") as infile:
        for line in infile:
            if line.startswith("Corrected"):
                extracted = [int(s) for s in line.split() if s.isdigit()]
                # print(extracted)
                # print(line)
                snps += extracted[0]
                insertions += extracted[1]
                deletions += extracted[3]
  
print("Corrected")              
print("snps: ", snps)
print("insertions: ", insertions)
print("deletions: ",deletions)

