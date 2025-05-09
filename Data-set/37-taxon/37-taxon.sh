#!/bin/bash
set -ex

#timeFile="37-taxa-time.csv" # will keep the records of running time
#RFRateFile="37-RFRateFile.csv" # will keep the records of RF rates


folders=( noscale.200g.500b  noscale.200g.250b noscale.200g.1000b noscale.200g.1500b noscale.50g.500b noscale.100g.500b noscale.200g.true scale2d.200g.500b scale2u.200g.500b noscale.25g.500b noscale.400g.500b noscale.800g.500b)
methodNames=(mpest astral.5.6.1 stelar SuperTriplets_v1.1)
#methodNames=(SuperTriplets_v1.1)
# Headers of the csv files
#echo "ilslevel, genetrees, basepair, methodName, modelConditionName,Robinson-Foulds distance" > $RFRateFile
#echo "name,time" > $timeFile

for folder in ${folders[@]}
do
	echo $folder
	#break
	for method in ${methodNames[@]}
	do
		for j in {1..20}
		do
			 #if [ $folder == "noscale.200g.true" ]; then
			  #gt=$folder/R$j/BS.1
  		 	#else
			#  gt=$folder/R$j/R$j.genetrees
			  #break
			# fi
			#echo $folder

			#sed 's/e-//g' $gt > $gt.temp #to replace e- with nothing
            #./strip_edge_support2.pl -i $gt.temp -o $gt.stripped
            #./reroot_tree.pl -t $gt.stripped -r GAL -o $gt.stripped.rooted
			speciesTree=$folder/R$j/R$j.$method
			gt=$folder/R$j/R$j.genetrees
			if [[ -f "$speciesTree.stripped" ]];
			then 
				continue
			fi

			if [ $method == "mpest" ]; then 
					cp $gt gt
					START=$(date +%s.%N)
					./mpest ./gt ./specieslist37 1
					END=$(date +%s.%N)
					DIFF=$(echo "$END - $START" | bc)
					cp gt.best.of.10.tre  $speciesTree
					rm -r gt*
			
			else 		
				if [ $method == "SuperTriplets_v1.1" ]; then
					if [[ -f "$speciesTree" ]]; 
					then 
						rm -r $speciesTree
					fi
					START=$(date +%s.%N)
			 		java -jar $method.jar  $gt  $speciesTree
					END=$(date +%s.%N)
        	        DIFF=$(echo "$END - $START" | bc)

				else
					
					START=$(date +%s.%N)
					java -jar $method.jar -i $gt -o $speciesTree
					END=$(date +%s.%N)
					DIFF=$(echo "$END - $START" | bc)
				fi
			fi

			./strip_edge_support2.pl -i $speciesTree -o $speciesTree.stripped
			RFdistance=$(python2 getFpFn.py -e $speciesTree.stripped -t true.tree.strip | sed 's/.//; s/,//' |awk '{print $3}')

			#printf $folder | tr '.' ',' >> $RFRateFile; echo ,$method,R$j,"${RFdistance//)}" >> $RFRateFile
			#echo $folder,R$j,$method,$DIFF >> $timeFile
			#break
		done
		#break
	done
	#break
done
