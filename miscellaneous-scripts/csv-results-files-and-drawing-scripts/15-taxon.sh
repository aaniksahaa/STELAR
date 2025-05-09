 #!/usr/bin/bash
set -xe
outputFile="15-taxa-fp-fn-RF-rates.csv"
#outputFile="SuperTriplets.time.csv"


echo "modelCondition,rep,method,FP-rate,FN-rate,Robinson-Foulds-distance" > $outputFile

#echo "name,time" > $outputFile

locations="$(ls -d */)"

for location in $locations
do
	if [ $location == "bin/" ]; then
		continue
	elif [ $location == "lib/" ]; then
		continue
	fi
	   for i in {1..10}
	    do
		 if [[ "$location" == "100gene-true/" ]]; then
		       cplocation=$location"R"$i
		 elif [ $location == "1000gene-true/" ]; then
		 	   cplocation=$location"R"$i
		 else
		       #continue
		       cplocation=$location"estimated-genetrees/R"$i
		fi

		#Stripping and rerooting after replacing e- and spaces
		#sed "s/ //g"  $cplocation/genetrees.tre >> gt.temp.1
		#sed "s/e-//g" gt.temp.1 >> gt.temp.2
		#./strip_edge_support2.pl -i gt.temp.2 -o  $cplocation/genetrees.tre.stripped
		#./reroot_tree.pl -t $cplocation/genetrees.tre.stripped -r A -o $cplocation/genetrees.tre.stripped.rooted


		#Running MPEST
		#cat $cplocation/genetrees.tre.stripped.rooted > gt.tre
		#./mpest ./gt.tre ./specieslist15 1
		#./strip_edge_support2.pl -i gt.tre.best.of.10.tre -o $cplocation".mpest"

		RFdistance=$(python2 getFpFn.py -e $cplocation.mpest -t true-species.tre.strip | sed 's/.//; s/,//' |awk '{print $1,$2,$3}')
		echo "${location///}",R$i,mpest,"${RFdistance//)}" >> $outputFile


		#Running SuperTriplets
	  #  START=$(date +%s.%N)
		#java -jar SuperTriplets.jar $cplocation/genetrees.tre.stripped.rooted sp.tre
		#END=$(date +%s.%N)
    #    DIFF=$(echo "$END - $START" | bc)
		#./strip_edge_support2.pl -i sp.tre -o $cplocation.SuperTriplets
		RFdistance=$(python2 getFpFn.py -e $cplocation.SuperTriplets -t true-species.tre.strip | sed 's/.//; s/,//' |awk '{print $1,$2,$3}')
		echo "${location///}",R$i,SuperTriplets,"${RFdistance//)}" >> $outputFile
		#echo "${location///},$DIFF" >> $outputFile

		#Running Astral
	  #java -jar astral.5.5.6.jar -i gt.tre -o sp.tre
		#./strip_edge_support2.pl -i sp.tre -o $cplocation.Astral
		RFdistance=$(python2 getFpFn.py -e $cplocation.Astral -t true-species.tre.strip | sed 's/.//; s/,//' |awk '{print $1,$2,$3}')
		echo "${location///}",R$i,Astral3,"${RFdistance//)}" >> $outputFile

		#Running Stelar
	    #java -jar stelar.jar -i gt.tre -o sp.tre
		#./strip_edge_support2.pl -i sp.tre -o $cplocation.stelar
		RFdistance=$(python2 getFpFn.py -e $cplocation.stelar -t true-species.tre.strip | sed 's/.//; s/,//' |awk '{print $1,$2,$3}')
		echo "${location///}",R$i,stelar,"${RFdistance//)}" >> $outputFile
		#rm  -r gt.*
		#break
	   done
	   #break
done
