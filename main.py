from Bio import SeqIO

def extractSeq(name,type): #Extract reads from the file. Output will be list like["string","string"....]
    seqs = [str(fa.seq) for fa in SeqIO.parse(name,type)]
    return seqs

seqList1=extractSeq("10k_reads.fastq", "fastq")


def doubling(s):
    n = len(s) #n:The length of the string we are processing.
    sa = []    # sa[i]:Where is the ith suffix in the original string?
    rk = []    # rk[i]:What is the rank of original string's suffix(which starts from index "i") in array "sa"?
                #rk will not be used in this project.
    for i in range(n):#Initialize the sa and rk of the string.
        rk.append(ord(s[i])-ord('a')) #Every suffix is ranked according to its initial letter.
        sa.append(i) #While the ith suffix is the suffix strat from the position where index=i

    l = 0 # l is the length that already sorted. Then sort in length 2l.
    sig = len(s) #The number of unique ranks.
    while True:
        p = []
        #For the suffixs not longer than l, their rank of second key work should be the shortest for position after l is empty.
        for i in range(n-l,n):
            p.append(i)
        # 对于其它长度的后缀来说，起始位置在`sa[i]`的后缀排名第i，而它的前l个字符恰好是起始位置为`sa[i]-l`的后缀的第二关键字
        #For suffixs in length of other value, the suffix starts at sa[i] ranks ith place. And its first l letters are the second key work of the suffix starts from sa[i]-l
        for i in range(n):
            if sa[i]>=l:
                p.append(sa[i]-l)
        # 然后开始基数排序，先对第一关键字进行统计
        #Randix sort.
        # 先统计每个值都有多少
        #First count the number of every value.
        cnt = [0]*sig
        for i in range(n):
            cnt[rk[i]] += 1
        # 做个前缀和，方便基数排序
        for i in range(1,sig):
            cnt[i] += cnt[i-1]
        # 然后利用基数排序计算新sa #Use Randix sort to calculate new sa.
        for i in range(n-1,-1,-1):
            cnt[rk[p[i]]] -= 1
            sa[cnt[rk[p[i]]]] = p[i]
        # 然后利用新sa计算新rk Then use new sa to calculate new rk
        def equal(i,j,l):
            if rk[i]!=rk[j]:return False
            if i+l>=n and j+n>=n:
                return True
            if i+l<n and j+l<n:
                return rk[i+l]==rk[j+l]
            return False
        sig = -1
        tmp = [None]*n
        for i in range(n):
            # 直接通过判断第一关键字的排名和第二关键字的排名来确定它们的前2l个字符是否相同
            if i==0 or not equal(sa[i],sa[i-1],l):
                sig += 1
            tmp[sa[i]] = sig
        rk = tmp
        sig += 1
        if sig==n:
            break
        # 更新有效长度 #Update the step length
        l = l << 1 if l > 0 else 1
    # Finally calculate the height array.
    k = 0
    height = [0]*n
    for i in range(n):
        if rk[i]>0:
            j = sa[rk[i]-1]
            while i+k<n and j+k<n and s[i+k]==s[j+k]:
                k += 1
            height [rk[i]] = k
            k = max(0,k-1) # 下一个height的值至少从max(0,k-1)开始
    return sa,rk,height


#假设字符串A的开头是核苷酸C
#去另一个字符串（B）里面找所有C开头的后缀
#记录下当前最长匹配长度和相应的字串索引（是这个文件里的哪个字串）



def saContain(strA,sa,strB,k,height):#Find if the prefix of second string is the suffix of the first string
    left=0
    right=len(sa)-1

    while left<=right: #O(logN)
        mid=left+((right-left)>>1)
        suffix=strA[sa[mid]:]
        strBPre=strB[:k]
        comp=compare1(suffix,strBPre)
        if comp==0: # If it exactly matches.
            return len(suffix),True
        elif comp==-1:
            left = mid + 1
        elif comp==1:
            right=mid-1
        else:#i.e. comp==2,When the prefix is a part of the suffix
            return extendCompare(strA,strB,k,sa,mid,height)

    return 0,False

def compare1(s1,s2):#s1:suffix from strA, s2:prefix from strB
    if s1==s2:
        return 0
    elif s1<s2:
        return -1
    else:
        if startwith(s1,s2):
            return 2
        else:
            return 1


def startwith(s1,s2): # Judge if the suffix s1 is started with s2
    res=True
    if len(s2)>len(s1):
        return False
    for i in range (len(s2)):
        if s1[i]!=s2[i]:
            res=False
    return res

def extendCompare(strA,strB,k,sa,mid,height):# Search in all the suffixs of strA, start with strB, with height array
    #height array of suffix could show us how many letters of the current string is the same as the former string.
    #Therefore we just need to compare those suffixs with the same prefix.
    minlcp=k
    i=mid
    #print(f"height{height},minlcp={minlcp}")
    while height[i]>=minlcp and i<len(height)-1:# Only compare the suffixs with the same prefix
        curr=strA[sa[i]:]
        lenSuf=len(curr)
        if strB[:lenSuf]==curr:# Find where it is and then return the index.
            return len(curr),True
        i+=1
    return 0,False
def wash(totalMap): #O(N) All of the sequences will be sorted by the function provided by python library.
    #Then the same reads will be removed.
    i=0
    totalMap.sort()
    while i< len(totalMap)-1:
        if totalMap[i+1]==totalMap[i] :
            totalMap.pop(i+1)
            i-=1
        i+=1

    return totalMap,i+1
#Mapping function.
def mapping(totalMap,seed):
    reads=len(totalMap)
    visited=[False for i in range(reads)]
    for i in range(reads):#To align reads to current read as many as possible.
        if visited[i]==False:
            strA=totalMap[i]
            sa, rk, height = doubling(strA)
            for j in range(reads):#Find other reads if it could align to strA.
                if j!=i and visited[j]==False:
                    strB=totalMap[j]
                    ixd,flag= saContain(strA,sa,strB,seed,height) # saContain() will be introduced.
                    if flag==True:#If strB's prefix is really a suffix of strA
                        totalMap[i]+=strB[ixd:]#Then align them.
                        visited[j]=True#Avoid to visit strB again.
                        strA = totalMap[i]#Update the current read, after expansion, then continue to find other reads suitable.
                        sa, rk, height = doubling(strA)
    res=[]
    for i in range(len(totalMap)): #Reduce the size of list to be processed later on.
        if visited[i]==False:
            res.append(totalMap[i])

    return res
#Do mapping in different length of suffix again.
def repeat(totalMap,seed):
    #count=0
    #length=[len(totalMap)]
    while seed>2:
        totalMap=mapping(totalMap,seed)
        seed-=1
    i=0
    count=0#number of reads not mapped
    while i < len(totalMap):
        if len(totalMap[i]) == 36:
            count+=1
        i+=1
    return totalMap,count

#m=["abcdefghijklmnoadayo","adayokimiaaaaagadare","uhiquhibbbbbhbaabcde","abcdeusacccccchojonx","adarekansaaaahauhbha"]
#print(mapping(m,4))

def N_metric(alignment):#"alignment" should be the final result of mapping, list type.
    #Then let's calculate the N50 metric
    alignment.sort(key=lambda i: len(i), reverse=True)
    total=0
    N=0
    longest=len(alignment[0])
    for lines in alignment:
        total+=len(lines)
    for i in range(len(alignment)):
        if(N/total>0.5):
            return i-1,len(alignment[i-1]),longest,total
        N+=len(alignment[i])


kmer=8
seqList1,i=wash(seqList1)
print(f"{i} reads remain after filtering the similar sequences from 10000 reads.")
rawSeq,notMapped=repeat(seqList1,kmer)
#output=finalProc(rawSeq)
#print(output,f"\nlength of the mapped sequences {len(output)}")
print(f"Percentage of reads mapped = {i-notMapped}/{i} = {100*(i-notMapped)/i}%")
place,n50Len,longest,totalLen=N_metric(rawSeq)
print(f"N50={n50Len} \nThe longest one's length is{longest}\nTotal length of all of the contigs:{totalLen}")
def writeIn(output):#To write and store the results in a txt file.
    count=1
    res=[]
    import time
    localtime = time.asctime(time.localtime(time.time()))
    for reads in output:
        info=f"\nSEQ:{count}  ,  LENGTH={len(reads)}\n"
        info+=reads
        res.append(info)
        count+=1
    finale=f"Last Modified: {localtime}\n"
    str1="".join(res)
    finale+=str1
    f=open("FASTQ_output.txt", "w+")
    f.write(finale)
    f.close()

writeIn(rawSeq)
