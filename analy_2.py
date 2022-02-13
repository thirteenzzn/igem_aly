
from analysis import analyze_c
import analysis
import time

# prefixed
pre_W2="CCGGGAAACAAACAAA"
pre_W1="GGGAAACAAACAAA"
pre_Z2="AGATGCGTGGGACAAAACAAAAC"
pre_Z1="GCGTGGGACAAAACAAAAC"

# suffixed
suf1="AAAAAGAAA"
suf2="AAAAAGAAA"

def main():
    time1=time.time()
    #I1
    # c1_analyze=analyze_c("R21003417-16-16.txt",suf1,suf2,pre_Z1,pre_W1)
    # c1_analyze.w_file("R21003417-16-16_clv.txt")
    # ref=c1_analyze.get_ref_read()
    # print(ref)

    #I2
    # c2_analyze=analyze_c("R21003406-20-20.txt",suf1,suf2,pre_Z2,pre_W2)
    # c2_analyze.w_file("R21003406-20-20_clv.txt")
    # ref=c2_analyze.get_ref_read()
    # print(ref)
    # analysis.calculate_f("R21003417-16-16_clv.txt","R21003406-20-20_clv.txt","fold_1620.txt","c1_16.txt","c2_20.txt",'clv_1620.xlsx')
    coef,intercept,intercept_=analysis.regr("clv_1620.xlsx")
    analysis.plotpic(float(coef),float(intercept),float(intercept_),'c1_16.txt','c2_20.txt','fold_1620.txt')
    time2=time.time()
    print("running time is "+str(time2-time1))
   

if __name__ == '__main__':
    main()