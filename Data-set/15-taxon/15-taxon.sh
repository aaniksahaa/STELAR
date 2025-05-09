 #!/usr/bin/bash
set -xe
outputFileRFDistance="data.csv"
#outputFileTime="SuperTriplets.time.csv"


echo "modelCondition,rep,method,RobinsonFouldsdistance" > $outputFileRFDistance

#echo "modelCondition,rep,method,time" > $outputFileTime

locations=(100gene-100bp 100gene-1000bp 100gene-true 1000gene-1000bp 1000gene-100bp 1000gene-true )
#locations=(1000gene-1000bp )
for location in "${locations[@]}"
do
	   for i in {1..10}
	   do
		cplocation=15-taxa/$location/estimated-genetrees/R$i.genetrees.tree
		sp=15-taxa/$location/estimated-genetrees/R$i

		#Running MPEST

		cat $cplocation > gt.tre
		#START=$(date +%s.%N)
		./mpest ./gt.tre ./specieslist15 1
		#END=$(date +%s.%N)
		#DIFF=$(echo "$END - $START" | bc)
		./strip_edge_support2.pl -i gt.tre.best.of.10.tre -o $sp".mpest"
		# Outputing results
		#RFdistance=$(python2 getFpFn.py -e $cplocation.mpest -t true-species.tre.strip | sed 's/.//; s/,//' |awk '{print $3}')
		#echo "${location///}",R$i,mpest,"${RFdistance//)}" >> $outputFileRFDistance
		#echo "${location///}",R$i,mpest,"$DIFF" >> $outputFileTime
		rm  -r gt.*

		#Running SuperTriplets
	    
		#START=$(date +%s.%N)
		#if [[ -f "sp.tre" ]]; then 
		#	rm sp.tre
		#fi
		#java -jar SuperTriplets_v1.1.jar $cplocation sp.tre
		#END=$(date +%s.%N)
        #DIFF=$(echo "$END - $START" | bc)
		#./strip_edge_support2.pl -i sp.tre -o $sp.SuperTriplets
		#rm sp.tre

		#if [[ -f "$cplocation.SuperTriplets" ]]; then 
		#	rm $cplocation.SuperTriplets 
		#fi
		
		#RFdistance=$(python2 getFpFn.py -e $sp.SuperTriplets -t true-species.tre.strip | sed 's/.//; s/,//' |awk '{print $3}')
		#echo "${location///}",R$i,SuperTriplets,"${RFdistance//)}" >> $outputFileRFDistance
		#echo "${location///}",R$i,SuperTriplets,"$DIFF" >> $outputFileTime
		
		

		#Running Astral3
		#START=$(date +%s.%N)
	  	#java -jar astral.5.6.1.jar -i $cplocation -o sp.tre -x
		#END=$(date +%s.%N)
        #DIFF=$(echo "$END - $START" | bc)
		#./strip_edge_support2.pl -i sp.tre -o $sp.astral.5.6.1
		#rm sp.tre 


		#RFdistance=$(python2 getFpFn.py -e $cplocation.Astral -t true-species.tre.strip | sed 's/.//; s/,//' |awk '{print $3}')
		#echo "${location///}",R$i,Astral3,"${RFdistance//)}" >> $outputFileRFDistance
		#echo "${location///}",R$i,Astral3,"$DIFF" >> $outputFileTime
		

		#Running Stelar
		#START=$(date +%s.%N)
	    #java -jar STELAR.jar 	-i gt.tre -o sp.tre -xt
		#END=$(date +%s.%N)
		#./strip_edge_support2.pl -i sp.tre -o $cplocation.stelar
		#rm sp.tre 

        #RFdistance=$(python2 getFpFn.py -e $cplocation.stelar -t true-species.tre.strip | sed 's/.//; s/,//' |awk '{print $3}')
		#echo "${location///}",R$i,stelar,"${RFdistance//)}" >> $outputFileRFDistance
		#echo "${location///}",R$i,stelar,"$DIFF" >> $outputFileTime
		#break
		done
	#break
done
