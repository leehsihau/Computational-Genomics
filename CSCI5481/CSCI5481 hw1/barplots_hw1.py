import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
import operator
import collections

with open('output_separate.csv', mode='r') as sep_infile:
    with open('output_whole.csv', mode='r') as whole_infile:
        with open('codon.csv', mode='r') as codon_infile:
            sep_reader = csv.reader(sep_infile)
            sep_mydict = {rows[0]:rows[1] for rows in sep_reader}

            whole_reader = csv.reader(whole_infile)
            whole_mydict = {rows[0]:rows[1] for rows in whole_reader}

            labels = []
            men_means=[]
            women_means=[]

            amino_list=[]
            codon_to_amino={}
            amino_count_sep={}
            amino_count_whole={}

            codon_reader=csv.reader(codon_infile)
            codon_to_amino={rows[0]:rows[1] for rows in codon_reader}
            for i in codon_to_amino:
                amino_list.append(codon_to_amino[i])
                #print(i)
            
            for i in sep_mydict:
                if codon_to_amino[i] not in amino_count_sep: amino_count_sep[codon_to_amino[i]]=int(sep_mydict[i])
                else: amino_count_sep[codon_to_amino[i]]+=int(sep_mydict[i])
            for i in whole_mydict:
                if codon_to_amino[i] not in amino_count_whole: amino_count_whole[codon_to_amino[i]]=int(whole_mydict[i])
                else: amino_count_whole[codon_to_amino[i]]+=int(whole_mydict[i])
            sorted_x = sorted(amino_count_sep.items(), key=operator.itemgetter(1), reverse=True)
            sorted_dict = collections.OrderedDict(sorted_x)
            for i in sorted_dict:
                labels.append(i)
                men_means.append(sorted_dict[i])
                women_means.append(amino_count_whole[i])
            

            #print(amino_list)

            #for i in sep_mydict:
            #    labels.append(i)
            #    men_means.append(int(sep_mydict[i]))
            #    women_means.append(int(whole_mydict[i]))
            #print(labels)
            

            #labels = ['G1', 'G2', 'G3', 'G4', 'G5']
            #men_means = [20, 34, 30, 35, 27]
            #women_means = [25, 32, 34, 20, 25]

            x = np.arange(len(labels))  # the label locations
            width = 0.25  # the width of the bars

            fig, ax = plt.subplots()
            rects1 = ax.bar(x - width/2, men_means, width, label='Coding Sequences')
            rects2 = ax.bar(x + width/2, women_means, width, label='Whole genome')

            # Add some text for labels, title and custom x-axis tick labels, etc.
            ax.set_ylabel('Frequency')
            ax.set_title('Amino Acid')
            ax.set_xticks(x)
            ax.set_xticklabels(labels)
            ax.legend()


            def autolabel(rects):
                """Attach a text label above each bar in *rects*, displaying its height."""
                for rect in rects:
                    height = rect.get_height()
                    ax.annotate('{}'.format(height),
                                xy=(rect.get_x() + rect.get_width() / 2, height),
                                xytext=(0, 20),  # 3 points vertical offset
                                textcoords="offset points",
                                ha='center', va='bottom')


            #autolabel(rects1)
            #autolabel(rects2)

            fig.tight_layout()
            fig.autofmt_xdate()
            plt.tick_params(axis='x', which='major', labelsize=8)

            plt.show()