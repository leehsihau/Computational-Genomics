import sys
import operator
import collections

def main():
    dic = {} # empty dictionary for storing codons and frequencies
    with open(sys.argv[1], 'r') as readfile:
        with open(sys.argv[2], 'w') as writefile:
            line = readfile.readline()
            while line:
                if line[0]=="A" or line[0]=="T" or line[0]=="G" or line[0]=="C": #if it doesn't start with DNA sequence ignore that line
                    lenth = len(line) 
                    lenth-=1 #Ignore the last new line charachter of every line

                    #if it's the whole genome file, ignore the last two or one genome
                    if lenth%3!=0: 
                        lenth=lenth-(lenth%3)
                    i=0
                    #start iterating through this line of genome
                    while i < lenth:
                        if lenth-i <3: break
                        curr = line[i:i+3] #take these three character as one codon
                        i+=3
                        if curr in dic: dic[curr]+=1
                        else: dic[curr]=1
                line=readfile.readline()
            sorted_x = sorted(dic.items(), key=operator.itemgetter(1), reverse=True) #sort this dictionary by key
            sorted_dict = collections.OrderedDict(sorted_x)

            #write to the output file
            for i in sorted_dict:
                writefile.write(""+i+","+str(sorted_dict[i])+"\n")
    

main()