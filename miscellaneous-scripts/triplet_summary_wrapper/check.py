import os,sys
import numpy as np
import os.path
import matplotlib.pyplot as plt
import re
def run(command):
	p1 = os.popen(command)
	temp = p1.readlines()
	p1.close()
	return temp
#f1 = open("ourMethod.txt","w+")
#f2 = open("mpest.txt","w+")

if __name__=="__main__":
    lines1 = run("cat 1")
    lines2 = run("cat 2")
    score = 0
    for line1 in lines1:
        for line2 in lines2:
            cnt1 = len(line1)
            cnt2 = len(line2)
            t1 = line1[:cnt1-2].rstrip()
            t2 = line2[:cnt2-2].rstrip()     
            if t1 == t2:
                res1 = line1[ cnt1-2 : cnt1-1]
                res2 = line2[ cnt2-2 : cnt2-1]
    #            print(res1, res2)
                #print(res)
                score += int(res1) * int(res2)
                #print(t1,t2)
                #print(line.rstrip()[:13])
    print(score)