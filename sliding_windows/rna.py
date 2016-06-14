from optparse import OptionParser
import os

def SCORE_MODIFIER(m):
    return 0

parser = OptionParser()
(options, args) = parser.parse_args()

filename = args[0]
if os.path.isfile(filename):
    execfile(filename)
else:
    print "ERROR! FAILED TO FIND FILE"
    exit();
SEQUENCE = len(SEQ)


def innermost_pair(struct):
    for i in range(len(struct)):
        if struct[i]==')':
            break
    end = i
    while(struct[i] != '('):
        i = i -1
    return i, end

def first_pair(struct):
    for i in range(len(struct)):
        if struct[i]=='(':
            return i


def write_to_file(prefix, end, width, struct, helices, helix_info,
                  wcpairs, gupairs, asymmetry, inner_mismatches, terminal_mismatches, chemmod, thermo):
    fil = open(prefix+end.__str__()+"x"+width.__str__()+".lab", 'a')
    out = [end, width, struct, helices, helix_info,
           wcpairs, gupairs, asymmetry, inner_mismatches, 
           terminal_mismatches, chemmod, thermo]
    out = [x.__str__() for x in out]
    out = ",".join(out) + "\n"
    fil.write(out)
    fil.close()

def score_helix(h):
    return -(1*h['wcpairs'] + 
             2*h['gupairs'] + 
             2*h['chemmod']+
             3*h['terminal_mismatches'] + 
             4*h['inner_mismatches'] +
             #5*h['asymmetry'] +
             SCORE_MODIFIER(h))

def best_in_file(i, j, pullfrom):
    fil = open(pullfrom+i.__str__()+"x"+j.__str__()+".lab", 'r')
    helices = fil.readlines()
    helices = [read_helix(x) for x in helices]
    # fix for ties.
    helicea = [(score_helix(helix), helix['thermo'], helix) for helix in helices if score_helix(helix) > -10000]
    if helicea:
        m = max(helicea)
        return m[2]

# TODO rm. COUNT
COUNT=0
def read_helix(string):
    struct = string.split(',');
    global COUNT
    COUNT +=1
    return dict({'end':int(struct[0]), 
                 'width':int(struct[1]),
                 'struct':struct[2], 
                 'wcpairs':int(struct[3]), 
                 'gupairs':int(struct[4]), 
                 'asymmetry':int(struct[5]), 
                 'inner_mismatches':int(struct[6]), 
                 'terminal_mismatches':int(struct[7]),
                 'chemmod':int(struct[8]),
                 'thermo':float(struct[9])}) 

def get_from_best_file(i, j):
    fil = open("best/"+i.__str__()+"x"+j.__str__()+".lab", 'r')
    helix = fil.readline()
    fil.close()
    ret = read_helix(helix)
    return ret;

def badchem(i, j, struct, offset):
    ri = i+offset;
    rj = j+offset;
    if MODS[ri] == 'M' or MODS[rj] == 'M':
        if (struct[i-1] == '(' and struct[i+1] == '('):
            if (not gupair(ri, rj) and 
                not gupair(ri-1, rj+1) and
                not gupair(ri+1, rj-1)):
                return 1000000
            else:
                return 2
        return 1
    return 0

def gupair(i, j):
    if ((SEQ[i] == 'G' and SEQ[j] == 'U') or
        (SEQ[i] == 'U' and SEQ[j] == 'G')):
        return 1
    return 0

def pairable(a, b):
    if ((a =='G' and b=='C') or 
        (a =='C' and b=='G') or 

        (a =='A' and b=='U') or 
        (a =='U' and b=='A') or 

        (a =='G' and b=='U') or 
        (a =='U' and b=='G')):
        return 1
    return 0

def ijpairs(m):
    ret = []
    s = m['struct']
    j = s.find(")")
    i = s.rfind("(", 0, j)
    while i >= 0:
        while (i >= 0 and s[i] != '('): i-=1
        while (j < len(s) and s[j] != ')'): j+=1
        ret.append([i+m['end']-m['width']+1, j+m['end']-m['width']+1]);
        i-=1
        j+=1
    return ret;
