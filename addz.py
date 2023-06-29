from readdata import readdata
a = readdata('C:\Users\hwyd\Desktop\MoSe2.vasp')
for i in range(len(a)):
    a[i][2] += 0.50
file1=open('C:\Users\hwyd\Desktop\3213.vasp','w')
for j in range(len(a)):
    print("",end="     ", file=file1)
    for k in range(3):        
        print("%.9f"%a[j][k],end="         ", file=file1)
    print("",file=file1)