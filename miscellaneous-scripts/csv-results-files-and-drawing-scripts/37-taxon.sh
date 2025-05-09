#!/bin/bash
set -ex

timeFile="37-taxa-time.csv" # will keep the records of running time
RFRateFile="37-RFRateFile.csv" # will keep the records of RF rates


folders=( noscale.200g.500b noscale.50g.500b noscale.100g.500b noscale.200g.true scale2d.200g.500b scale2u.200g.500b )
methodNames=(mpest SuperTriplets astral.5.5.6 stelar )

# Headers of the csv files
echo "ilslevel, genetrees, basepair, methodName, modelConditionName,Robinson-Foulds distance" > $RFRateFile
echo "name,time" > $timeFile

for folder in ${folders[@]}
do
	echo $folder
	#break
	for method in ${methodNames[@]}
	do
		for j in {1..20}
		do
			 if [ $folder == "noscale.200g.true" ]; then
			  gt=$folder/R$j/BS.1
  		 else
			  gt=$folder/R$j/Best.1
			  #break
			 fi
			#echo $folder

			sed 's/e-//g' $gt > $gt.temp #to replace e- with nothing
                        ./strip_edge_support2.pl -i $gt.temp -o $gt.stripped
                        ./reroot_tree.pl -t $gt.stripped -r GAL -o $gt.stripped.rooted

			if [ $method == "mpest" ]; then
					cp $gt.stripped.rooted gt
					START=$(date +%s.%N)
					./mpest ./gt ./specieslist37 1
					END=$(date +%s.%N)
					DIFF=$(echo "$END - $START" | bc)
					speciesTree=$folder/R$j/Best.1.tre
					cp gt.best.of.10.tre  $speciesTree
					rm -r gt*
			fi


			if [ $method == "SuperTriplets" ]; then
					rm -r $speciesTree
					speciesTree=$folder/R$j/$method
					START=$(date +%s.%N)
			 		java -jar $method.jar  $gt.stripped.rooted  $speciesTree
					END=$(date +%s.%N)
        	                        DIFF=$(echo "$END - $START" | bc)

			else
				speciesTree=$folder/R$j/$method
				START=$(date +%s.%N)
				java -jar $method.jar -i $gt.stripped.rooted -o $speciesTree
				END=$(date +%s.%N)
				DIFF=$(echo "$END - $START" | bc)
			fi

			./strip_edge_support2.pl -i $speciesTree -o $speciesTree.stripped
			RFdistance=$(python2 getFpFn.py -e $speciesTree.stripped -t true.tree.strip | sed 's/.//; s/,//' |awk '{print $3}')

			printf $folder | tr '.' ',' >> $RFRateFile; echo ,$method,R$j,"${RFdistance//)}" >> $RFRateFile
			echo $folder,R$j,$method,$DIFF >> $timeFile
			#break
		done
		#break
	done
	#break
done
