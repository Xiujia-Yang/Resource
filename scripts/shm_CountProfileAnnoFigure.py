from __future__ import division
import csv, sys, os, glob, re
import numpy as np
import pandas as pd

def main():
    infile = sys.argv[1]#"Fig.5c_d.profile.mut.rate.txt"
    df=pd.read_csv(infile, sep="\t")
    Family=df["Family"].unique()
    mydata = []

    for i in Family:
        mydf=df[(df["Family"] == i)]
        Germline_id = mydf.pop("Germline_id")
        Family=mydf.pop("Family")
        Clone_num=mydf.pop("Clone_number")

        data = []
        data.append(i)
        for index, col in mydf.iteritems():
            value_col=col.tolist()
            new_list = [i for i in value_col if i != -100]
            if len(new_list) == 0:
                avg_num = 0
            else:
                avg_num = sum(new_list)/len(new_list)
            data.append(avg_num)

        mydata.append(data)

    Germline_id = df.pop("Germline_id")
    Family = df.pop("Family")
    Clone_num=df.pop("Clone_number")
    data=[]
    data.append("All")

    num = 0
    for index_all, col_all in df.iteritems():
        num +=1
        value_all = col_all.tolist()
        new_value = [i for i in value_all if i != -100]
        avg_all = sum(new_value)/len(new_value)
        #print (new_value)
        #print (avg_all)
        data.append(avg_all)

    #print (num)
    pos_list = [i for i in range(num)]
    #print (pos_list)
    mydata.append(data)
    new_df=pd.DataFrame(np.array(mydata), columns=['Family']+pos_list)
    new_df.to_csv("Fig.5c_d.profile.mut.avg.rate.txt",sep="\t",index=False)


if __name__ == "__main__":
    main()
