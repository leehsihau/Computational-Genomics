import sys
import random
import numpy as np
import matplotlib.pyplot as plot
MISMATCH=-3
MATCH=1
GAP=-2


def needleman_wunsch(ref, query):
    matrix=[]
    ref_len=len(ref)
    query_len=len(query)
    #initialize matrix as all 0s
    for i in range(0, query_len+1):
        placeholder=[]
        for j in range(0, ref_len+1):
            placeholder.append(0)
        matrix.append(placeholder)
    #print(str(matrix))
    #calculate gaps for the first row, which is reference
    for i in range(1, ref_len+1):
        matrix[0][i]=matrix[0][i-1]+GAP
    #Calculate gaps for the first column, which is query
    for i in range(1, query_len+1):
        matrix[i][0]=matrix[i-1][0]+GAP
    for i in range(1, query_len+ 1):
        for j in range(1, ref_len+ 1):
            #check if it is a match or not match
            value = MISMATCH
            if ref[j-1] == query[i-1]:
                value = MATCH
            #Choose the max value of ref with gap, query with gap and ref/query with match or mismatch
            matrix[i][j]=max([matrix[i-1][j-1] + value, matrix[i-1][j] + GAP, matrix[i][j-1] + GAP])
    #start from the down right corner of the matrix
    start_ref=ref_len
    start_que=query_len
    score=matrix[start_que][start_ref]
    #create empty placeholders for reconstruction
    query_=""
    ref_=""
   # print("start_ref: "+str(start_ref))
   # print("start_que: "+str(start_que))
    #execute until one of query/reference is done
    while start_ref > 0 and start_que > 0:
        #calculate the left, up and diagonal for comparison
        diagonal=matrix[start_que-1][start_ref-1]
        up=matrix[start_que-1][start_ref]
        left=matrix[start_que][start_ref-1]
        curr=matrix[start_que][start_ref]
        #compare the chosen max with diagonal
        if (curr == diagonal+MATCH and ref[start_ref-1]==query[start_que-1]) or (curr == diagonal+MISMATCH and ref[start_ref-1]!=query[start_que-1]) :
            #meaning diagonal is the max, doesn't matter if it's a match or mismatch
            ref_ = ref[start_ref-1] + ref_
            query_ = query[start_que-1] + query_
            #goes diagonal
            start_ref -= 1
            start_que -= 1
        #compare the chosen max with left
        elif curr == left+GAP:
            #coming from left, meaning query is taking a gap
            ref_ = ref[start_ref-1] + ref_ 
            query_ = "_" + query_
            #goes left
            start_ref -= 1
        #compare the chosen max with up 
        elif curr == up+GAP:
            #coming from up , meaning reference is taking a gap
            ref_ = "_" + ref_
            query_ = query[start_que-1] + query_
            # goes up
            start_que -= 1
    #print("start_ref: "+str(start_ref))
    #print("start_que: "+str(start_que))
    #when start_ref is bigger than 0, meaning that it didn't take diagonal and reference still has genes that are not aligned while query is already done

    while start_ref>0:
        ref_=ref[start_ref-1]+ref_
        query_="_"+query_
        start_ref-=1
    
    #when start_que is bigger than 0, meaning that it didn't take diagonal and query still has genes that are not aligned while reference is already done
    while start_que>0:
        ref_="_"+ref_
        query_=query[start_que-1]+query_
        start_que-=1
    return score, ref_, query_

def anchored_helper(ref, query):
    print("ref: "+ref)
    print("query: "+query)
    length=len(ref)
    score=0
    #for anchored NW, just simply counting the matches and unmatches
    for i in range(0,length):
        if ref[i]==query[i]:
            score+=MATCH
        else:
            score+=MISMATCH
    print("anchored score: "+str(score))
    print("\n")

    return score

def anchored_needleman_wunsch(ref, query, matches_ref, matches_que):
    total=0
    last_ref_end, last_que_end, ref_start,ref_end,que_start,que_end=-1,-1,0,0,0,0
    ref_=""
    que_=""
    #pop from the matches of reference and matches of query within every iteration
    while len(matches_que)!=0:
        # indices of both reference and query
        ref_indices=matches_ref.pop(0)
        que_indices=matches_que.pop(0)
        #starting index and ending index of both strings 
        ref_start=ref_indices[0]-1
        ref_end=ref_indices[1]-1
        que_start=que_indices[0]-1
        que_end=que_indices[1]-1
        #run NW with un-anchored sections
        score, ref__, que__=needleman_wunsch(ref[last_ref_end+1:ref_start], query[last_que_end+1:que_start])
        ref_+=ref__
        que_+=que__
        total+=score
        #count the matches and unmatches in anchored sections
        score=anchored_helper(ref[ref_start:ref_end+1], query[que_start:que_end+1])
        ref_+=ref[ref_start:ref_end+1]
        que_+=query[que_start:que_end+1]
        total+=score
        #set for next iteration
        last_ref_end=ref_end
        last_que_end=que_end
    #if there is more after the last anchored section, align that with NW
    score, ref__, que__=needleman_wunsch(ref[last_ref_end+1:], query[last_que_end+1:])
    ref_+=ref__
    que_+=que__
    total+=score
    #print("total: "+str(total))
    print("verify: "+str(verify_score(ref_, que_)))
    return total

def permute(sequence):
    #use random library to permute the sequence
    return ''.join(random.sample(sequence, len(sequence)))

def permute_100_times(ref, query):
    results = []
    rand_q=query
    rand_r=ref
    #permute sequences for 100 times
    for i in range(0, 100):
        print(str(i))
        score, ref_, que_ = needleman_wunsch(rand_r, rand_q)
        rand_r = permute(ref)
        results.append(score)
    #convert list to a numpy array
    arr = np.array(results)
    plot.hist(arr, bins='auto')
    plot.show()

def verify_score(ref, que):
    lenght=len(ref)
    lenght=len(que)
    score=0
    #verification of NW alignment score
    for i in range(0,lenght):
        if ref[i]==que[i]:
            #print("m")
            score+=MATCH
        elif ref[i]=="_" or que[i]=="_":
            score+=GAP
            #print("G")
        else:
            if ref[i]=="_" and que[i]=="_":
                print("CNM")
            score+=MISMATCH
            #print("MIS")
    #print(score)
    return score


#the regular NW
if len(sys.argv)==5:
    #get filenames from argv
    reference_fn=sys.argv[2]
    query_fn=sys.argv[4]
    #print(reference_fn)
    #print(query_fn)
    #open the file and read the second line ignoring the newline
    reference=open(reference_fn, "r").readlines()[1][:-1]
    query=open(query_fn, "r").readlines()[1][:-1]
    #reference="ATGTACAAAAAA"
    #query="ATGCTAGGTAC"
    #print(reference)
    #print("\n")
    #print(query)
    #run regular NW
    score, ref_, que_=needleman_wunsch(reference, query)
    print(score)
    print("alignment result: ")
    print(ref_)
    print("\n")
    print(que_)
    print("\n")
    #verification of the alignment results
    verify_score(ref_, que_)
    #permute_100_times(reference, query)
#run anchored NW
if len(sys.argv)==7:
    #get filenames from argv
    reference_fn=sys.argv[2]
    query_fn=sys.argv[4]
    matches_fn=sys.argv[6]
    print(reference_fn)
    print(query_fn)
    print(matches_fn)
    #open the file and read the second line ignoring the newline

    reference=open(reference_fn, "r").readlines()[1][:-1]
    query=open(query_fn, "r").readlines()[1][:-1]
    matches_1=[]
    matches_2=[]
    for line in open(matches_fn, "r").readlines()[1:]:
        placeholder=[]
        #split by space, tab
        arr=line.split()
        #matches for sars_cov_1
        placeholder.append(int(arr[0]))
        placeholder.append(int(arr[1]))
        matches_1.append(placeholder.copy())
        #clear the placeholder list
        placeholder.clear()
        #do the same for sars_cov_2
        placeholder.append(int(arr[2]))
        placeholder.append(int(arr[3]))
        matches_2.append(placeholder.copy())
    print(matches_1)
    print(matches_2)
    score=anchored_needleman_wunsch(reference, query, matches_1, matches_2)
    print(score)

