import os,sys
import numpy as np
import os.path
import matplotlib.pyplot as plt
import re
def run(command):
	p1 = os.popen(command)
	temp = p1.readline()
	p1.close()
	return temp.rstrip()
chars = ['A','B','C','D','E','F','G','H','I','J','K']
def label(source, des):
    res=run("cat "+source)
    cnt=10
    for i in range(len(chars)-1,-1,-1):
        res=res.replace(str(cnt),chars[i])
        cnt-=1
        #print(cnt)
    command = 'echo "'+res +'" > ' + des  
    run(command)

if __name__=="__main__": 
      paths = ["higher-ILS", "lower-ILS"]
      subpaths = ["estimated-genetrees", "true-genetrees"]
      for  path in paths:
            for subpath in subpaths:
                  for i in range(1,21):
                        st = path+ "/true-speciestrees/"+"R"+str(i)+".true.tre"
                        run ("./strip_edge_support2.pl -i "+st+" -o striped_true_tree.tree")
                        label("striped_true_tree.tree","rooted_striped_true_tree.tree")
                        run ("./reroot_tree.pl -t rooted_striped_true_tree.tree -r A -o  "+ path+ "/true-speciestrees/"+"R"+str(i)+".label.true.tre")
                        
                        f3 = open("combined_gene_tree.tree","w+")
                        for j in range(1,201):
                              J = str(j)
                              while(len(J)<4):
                                    J = "0"+J
                              if subpath == "true-genetrees" : gt = path + "/"+subpath+"/" + "R" + str(i) + "/"+J+"/true.gt"
                              else: gt = path + "/"+subpath+"/" + "R" + str(i) + "/"+J+"/RAxML_bipartitions.final.f200"
                              res=run("cat "+gt)
                              res=res.replace("e-","").replace(" ","")
                              run(' echo "'+res+'" > temp.txt')      
                              run ("./strip_edge_support2.pl -i temp.txt -o striped_combined_gene_tree.tree")
                              label("striped_combined_gene_tree.tree", "label.striped_combined_gene_tree.tree")      
                              res=run("cat label.striped_combined_gene_tree.tree")
                              f3.write(res+"\n")
                              #break
                        f3.close()
                        run ("./reroot_tree.pl -t combined_gene_tree.tree -r A -o  rooted_striped_combined_gene_tree.tree")
                        run("cp rooted_striped_combined_gene_tree.tree "+path+ "/"+subpath+"/" + "R" + str(i)+"/all_gt.tre")
                        #run("java -jar dynadup-sp-from-triplets.jar -i "+gt+"/all_gt.tre -o "+gt+"/st_from_triplets.txt")
                        #res = run("python2 getFpFn.py -e  st_from_triplets.txt -t rooted_strpped_true-species.tre | sed 's/.//; s/,//' |awk '{print $1}'")
                        #f1.write(res+" ") 
                        #break
            #f1.write("\n")
                  #break

      #f1.close()
            


            
