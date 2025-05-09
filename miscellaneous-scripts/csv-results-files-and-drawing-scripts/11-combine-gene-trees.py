import os,sys
import numpy as np
import os.path
import matplotlib.pyplot as plt
import re
def run(command):
      #print(command)
      p1 = os.popen(command)
      temp = p1.readlines()
      p1.close()
      #print(temp)
      return temp


def run1(command):
      p1 = os.popen(command)
      lines = p1.readlines()
      genetrees = ""
      p1.close()
      #print(temp)
      for line in lines:
            line = line.replace("e","").replace("-","")
            genetrees += line
            #run("echo "+line+" > gt.tre")
      return genetrees.rstrip()


chars = ['A','B','C','D','E','F','G','H','I','J','K']

def label(source, des): 
    f = open(source,"r")
    results=f.readlines()
    results = [x.strip() for x in results]
    f.close()
    
    f1 = open(des,"w")
    for res in results:
          cnt = 11
          for i in range(len(chars)-1,-1,-1):
                res = res.replace(str(cnt)+"."+str(cnt), chars[i])
                res = res.replace(str(cnt),chars[i])
                cnt-=1
          f1.write(res+"\n")
    f1.close()

folders = ["estimated_Xgenes_strongILS", "simulated_Xgenes_strongILS"]
genes = ["5", "15", "25", "50", "100"]
if __name__=="__main__": 
      #folders = next(os.walk('.'))[1] # get the folder list
      for  folder in folders:
            for gene in genes:
                        gt = folder + "/" + folder.replace("X", gene)
                        print(gt)
                        for i in range(1,21):
                              location = gt + "/Rep"+str(i)+"_"+gene+"genes"
                              print(location)
                              genetrees = run1("cat "+location)
                              genetrees.rstrip()
                              #print(genetrees)
                              run('echo "'+ genetrees + '" > gt.tre')
                              run ("./strip_edge_support2.pl -i gt.tre -o striped.gt.tre")
                              label("striped.gt.tre","label.striped.gt.tre")
                              run("cp label.striped.gt.tre "+location+".strip.labeled")
