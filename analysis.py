import re
import numpy as np
import pandas as pd
from sklearn import linear_model
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pandas import DataFrame
import math

# analyze the cleavage
class analyze_c:
    def __init__(self,filename,suf1,suf2,pre1,pre2):
        self.filename=filename
        self.suf1=suf1
        self.suf2=suf2
        self.pre1=pre1
        self.pre2=pre2
    #get sequence's data
    def get_data_seq(self):
        seq=[]
        with open(self.filename,'r') as pf:
            for line in pf.readlines():
                data=re.split(r"[ ]+", line)
                data[0]=int(data[0])
                data[1]=data[1].rstrip('\n')
                if data[0]>30:
                    seq.append(data)
        return seq

    #get ref's information ,return A,W,ref_r(array)
    def get_ref_read(self):

        #seq中seq[i][1]是读取次数，seq[i][2]是序列
        seq=analyze_c.get_data_seq(self)
        
        ref_x = []
        ref = []
        ref_r = []  #存储reference  序列，长度，W读数, Z读数

        #读取ref
        with open("ref.txt", "r") as f3:
            for line in f3.readlines():
                data=line.rstrip('\n')
                if len(data)>1:
                    ref_x.append(str(data))
        
        tmp_ll=[]
        #互补倒序
        for i in range(len(ref_x)):
            tmp_ll.append([])
            for x in ref_x[i]:
                if x=='A':
                    tmp_ll[i].append('T')
                elif x=='T':
                    tmp_ll[i].append('A')
                elif x=='C':
                    tmp_ll[i].append('G')
                elif x=='G':
                    tmp_ll[i].append('C')

        for i in range(len(tmp_ll)):
            tmp_ll[i].reverse()
            x=''.join(tmp_ll[i])
            ref.append(x)

        # reference  序列，长度，Z读数, W读数
        for i in range(4):
            ref_r.append([ref[i],len(ref[i]),0,0,0])

        # 读取ref的reads
        for x in ref_r:
            for y in seq:
                if y[1].find(x[0])!=-1:
                    if y[1].find(self.pre1)!=-1:
                        x[2]+=y[0]
                    if y[1].find(self.pre2)!=-1:
                        x[3]+=y[0]


        # 第一、二位分别存放Z,W
        ref_final=[]
        for x in ref_r:
            ref_final.append([x[0],x[2],x[3]])

        return ref_final

    #find s and its ref ,return read_final(the first is s,the second is s'reads, the last is s'prefixed)
    def get_s_read(self):
        seq=analyze_c.get_data_seq(self)

        read_sf=[]

        for inf_seq in seq:
            list_tmp=[]
            if inf_seq[1].find(self.pre1)!=-1:
                if inf_seq[1].find(self.suf1)!=-1:
                    x=inf_seq[1][:inf_seq[1].find(self.suf1)]
                elif inf_seq[1].find(self.suf2)!=-1:
                    x=inf_seq[1][:inf_seq[1].find(self.suf2)]
                else:
                    x=inf_seq[1]
                list_tmp.append(x[x.find(self.pre1)+len(self.pre1):])
                list_tmp.append(inf_seq[0])
                list_tmp.append('Z')
            if inf_seq[1].find(self.pre2)!=-1:
                if inf_seq[1].find(self.suf1)!=-1:
                    x=inf_seq[1][:inf_seq[1].find(self.suf1)]
                elif inf_seq[1].find(self.suf2)!=-1:
                    x=inf_seq[1][:inf_seq[1].find(self.suf2)]
                else:
                    x=inf_seq[1]
                list_tmp.append(x[x.find(self.pre2)+len(self.pre2):])
                list_tmp.append(inf_seq[0])
                list_tmp.append('W')
            if len(list_tmp)!=0:
                read_sf.append(list_tmp)

        read_sf.sort()

        read_final=[]
        key=0

        #merge the same sequence
        for i in range(len(read_sf)-1):
            if i==key:
                for j in range(key+1,len(read_sf)):
                    if read_sf[i][0]==read_sf[j][0] and read_sf[i][2]==read_sf[j][2]:
                        read_sf[i][1]+=read_sf[j][1]
                    else:
                        key=j
                        break
                read_final.append([read_sf[i][0],read_sf[i][1],read_sf[i][2]])
        
        return read_final

    #calculate cleavage
    def calculate_cleavage(self):
        read_final=analyze_c.get_s_read(self)
        ref_final=analyze_c.get_ref_read(self)
        read_all={}
        for s in read_final:
            read_num=[]  # 存储读取次数
            rW=0.0
            rZ=0.0
            if s[2]=='W':
                rW+=s[1]
            else:
                rZ+=s[1]
            
            key=read_final.index(s)

            for i in range(key+1,len(read_final)):
                if read_final[i][0]==s[0]:
                    if read_final[i][2]=='W':
                        rW+=read_final[i][1]
                    else:
                        rZ+=read_final[i][1]
            
            if abs(len(s[0])-len(ref_final[0][0]))>abs(len(s[0])-len(ref_final[1][0])):
                Z=ref_final[3][1]
                W=ref_final[1][2]
            else:
                Z=ref_final[2][1]
                W=ref_final[0][2]

            clv=(rZ/Z)/(rZ/Z+rW/W)
            read_num.append(rW)
            read_num.append(rZ)
            read_num.append(clv)

            if read_all.get(s[0])==None:
                if read_num[0]!=0 and read_num[1]!=0:
                    read_all[s[0]]=read_num
        return read_all

    #write cleavage into file
    def w_file(self,resfile):
        read_all=analyze_c.calculate_cleavage(self)
        # cleavage写入文件
        with open(resfile, "w") as pf:
            for x in read_all.keys():
                pf.write(str(read_all[x][2]))
                pf.write("\t")
                pf.write(str(read_all[x][0]))
                pf.write("\t")
                pf.write(str(read_all[x][1]))
                pf.write("\t")
                pf.write(str(x))
                pf.write("\n")


# calculate fold change
def calculate_f(filel,file2,res,tmp1,tmp2,xlsxfile):
    # 读取文件数据
    def rf(s,llist):
        with open(s, "r") as f:
            for line in f.readlines():
                data=re.split(r"[\t]+", line)
                data[3]=data[3].rstrip('\n')
                data[0]=float(data[0])
                llist.append(data)

    clv1=[]
    clv2=[]
    clv1_f=[]
    clv2_f=[]

    rf(filel,clv1)
    rf(file2,clv2)

    # calculate fold change
    fd_chg=[]
    fd_read=[]
    for i in range(len(clv2)):
        for x in clv1:
            if clv2[i][3]==x[3]:
                fd=(1-(clv2[i][0]))/(1-(x[0]))
                fd_chg.append(fd)
                fd_read.append([clv2[i][1],clv2[i][2],x[1],x[2]])
                clv2_f.append(clv2[i][0])
                clv1_f.append(x[0])
                break

    #rerange fold change
    me=np.median(np.array(fd_chg))
    k=1/me

    for x in fd_chg:
        x=x*k

    # 写入文件
    with open(res,"w") as pf:
        for x,y in zip(fd_chg,fd_read):
            pf.write(str(x))
            pf.write("\t")
            pf.write(str(y[0]))
            pf.write("\t")
            pf.write(str(y[1]))
            pf.write("\t")
            pf.write(str(y[2]))
            pf.write("\t")
            pf.write(str(y[3]))
            pf.write("\n")

    with open(tmp1, "w") as pf:
        for x in clv1_f:
            pf.write(str(x))
            pf.write("\n")

    with open(tmp2, "w") as pf:
        for x in clv2_f:
            pf.write(str(x))
            pf.write("\n")

    data={
    'clv1':clv1_f,
    'clv2':clv2_f,
    }
    df=DataFrame(data)
    df.to_excel(xlsxfile)

# LinearRegression
def regr(filename):
    df = pd.read_excel(filename) 
    df['clv1'].fillna(0,inplace = True) 
    df['clv2'].fillna(0,inplace = True) 
    cdf=df[['clv1','clv2']]

    msk =np.random.rand(len(df))<0.8
    train=cdf[msk]
    test=cdf[~msk]

    regr = linear_model.LinearRegression()
    train_x=np.asanyarray(train[['clv1']])
    train_y=np.asanyarray(train[['clv2']])
    regr.fit (train_x, train_y)

    coef=float(regr.coef_)
    intercept=float(regr.intercept_)

    # bootstrap B=1000
    nboot=1000
    Q=0
    for i in range(min(len(df['clv1']),len(df['clv2']))):
        Q+=(df['clv2'][i]-df['clv1'][i]*coef-intercept)**2

    e=np.sqrt(Q/(nboot-2))

    Sx=0
    sum_tmp=0
    for x in df['clv1']:
        sum_tmp+=x
    for x in df['clv1']:
        Sx+=(x-sum_tmp)**2
    Sx=Sx/nboot

    tt=3.291
    ci=float(tt*e/np.sqrt(Sx))
    intercept-=ci

    return regr.coef_,regr.intercept_,intercept

# plot
def plotpic(k,b0,b_,file1,file2,file3):
    # 加载数据
    clv1 = np.loadtxt(file1, encoding='utf-8')
    clv2 = np.loadtxt(file2, encoding='utf-8')
    fold = np.loadtxt(file3, encoding='utf-8')

    #matplotlib画图中中文显示会有问题，设置默认字体
    plt.rcParams['font.sans-serif']=['Times New Roman']
    plt.rcParams['axes.unicode_minus'] = False

    plt.xlabel('Fraction cleaved (-ligand)',size=14)
    plt.ylabel('Fraction cleaved (+ligand)',size=14)
    plt.xlim(xmax=1,xmin=0)
    plt.ylim(ymax=1,ymin=0)

    clv1=np.array(clv1)
    clv2=np.array(clv2)

    colors1 = '#0095ff' #blue
    colors2 = '#DC143C' #red
    area = np.pi   # points' area
    x=np.linspace(0,1,100)

    x1=[]
    x2=[]
    sp1=[]
    sp2=[]

    b1=b0+b_  #上界
    b2=b0-b_  #下界
    y1=x*k+b0
    y2=x*k+b1
    y3=x*k+b2

    plt.tick_params(labelsize=14)  #坐标轴字体大小
    plt.scatter(clv1, clv2, s=area, c=colors1,alpha=0.2)
    for i in range(len(clv2)):
        if clv2[i]-clv1[i]*k-b1>0 or clv2[i]-clv1[i]*k-b2<0:
            if fold[i][0]>3:
                print(fold[i][0])
                x1.append(clv1[i])
                x2.append(clv2[i])
            else:
                sp1.append(clv1[i])
                sp2.append(clv2[i])
    # for i in range(len(clv2)):
    #     if clv2[i]-clv1[i]*k-b1>0 or clv2[i]-clv1[i]*k-b2<0:
    #         if fold[i]>2:
    #             x1.append(clv1[i])
    #             x2.append(clv2[i])
    #         else:
    #             sp1.append(clv1[i])
    #             sp2.append(clv2[i])

    x1=np.array(x1)
    x2=np.array(x2)
    sp1=np.array(sp1)
    sp2=np.array(sp2)
    plt.scatter(sp1, sp2, s=area, c=colors1,alpha=0.6)
    plt.scatter(x1, x2, s=area, c=colors2)
    plt.plot(x,y1,c='black',linestyle='dashdot',alpha=0.4)
    plt.plot(x,y2,c='r',linestyle='--',alpha=0.5)
    plt.plot(x,y3,c='r',linestyle='--',alpha=0.5)

    plt.show()