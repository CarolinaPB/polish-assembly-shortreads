import os, sys

snps = 0
insertions=0
deletions=0
directory="logs_rules/pilon"
directory=os.path.join(os.getcwd(),"logs_rules/pilon")

# outfile = sys.argv[1]

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfile", help="File where summary will be saved")
args = parser.parse_args()

outfile = os.path.join(os.getcwd(),args.outfile)

for file in os.listdir(directory):
    with open(os.path.join(directory,file), "r") as infile:
        for line in infile:
            if line.startswith("Corrected"):
                extracted = [int(s) for s in line.split() if s.isdigit()]
                snps += extracted[0]
                insertions += extracted[1]
                deletions += extracted[3]
                
print("Corrected\n")
print(f"snps: {snps}\n")
print(f"insertions: {insertions}\n")
print(f"deletions: {deletions}\n")

with open(outfile, "w") as out:
    out.write("Corrected\n")
    out.write(f"snps: {snps}\n")
    out.write(f"insertions: {insertions}\n")
    out.write(f"deletions: {deletions}\n")