f1=open('1.txt', encoding='utf-8', errors='ignore' )
f2=open('2.txt', encoding='utf-8', errors='ignore' )
line1=f1.readline()
line2=f2.readline()
re=[]
count=1
while line1 and line2:
   
    if line1!=line2:
        re.append(count)
    count+=1
    line1=f1.readline()
    line2=f2.readline()
print(re)
