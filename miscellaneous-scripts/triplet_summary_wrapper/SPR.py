# From Ishrat Jahan of UIU

from pylogeny.tree import tree     # import tree class definition


t = tree("(GAL,(ORN,((MAC,MON),(((CHO,DAS),(ECH,(LOX,PRO))),(((TUP,((ORY,OCH),(SPE,(CAV,(DIP,(MUS,RAT)))))),((OTO,MIC),(TAR,(CAL,(NEW,(PON,(GOR,(PAN,HOM)))))))),((SOR,ERI),((MYO,PTE),(((FEL,CAN),(VIC,(SUS,(TUR,BOS)))),EQU))))))));") # instantiate with this Newick string
#print t
topo = t.toTopology()       # converts a tree to a rich object of branches and nodes

rearrList = topo.allSPR()   # a list of rearrangement objects
for rearr in rearrList:
  move = rearr.doMove()     # perform the rearrangement and get a topology object
  print(move)                # prints the Newick string for this topology
  
  #foo(move)                # do something with this topology
