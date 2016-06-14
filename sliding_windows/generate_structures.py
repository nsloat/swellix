import rna
import os
import sys
import subprocess

for i in range(rna.SEQUENCE+1):

    os.system("echo \""+rna.SEQ[max(0, i-rna.WINDOW):i]+"\n"+rna.MODS[max(0, i-rna.WINDOW):i]+"\" | ./subopt_and_label "+sys.argv[1]+" "+sys.argv[2]+" "+(i-rna.WINDOW-1).__str__()+" "+rna.SEQ+" "+rna.MODS+" "+rna.WINDOW.__str__()+" "+rna.LENGTH.__str__()+" "+rna.TERMINAL_MISMATCHES.__str__()+" "+rna.ASYMMETRY.__str__()+" "+rna.MISMATCHES.__str__()+" "+sys.argv[3])
#    print "start: "+(i-rna.WINDOW-1).__str__()+"\nsubSeq: "+rna.SEQ[max(0,i-rna.WINDOW):i]+"\nsubMod: "+rna.MODS[max(0, i-rna.WINDOW):i]+"\n"
