for age in $(seq 8000 1000 80000);
	do
	echo "dupe_prob_output.py /Users/joshuaschraiber/Desktop/msmsdir/msms3.2rc-b163.jar -Y 110 -fi 0.001 -tSel ${age} -nrep 1000000 -sAA 0.0 -sAa 0.0"
done
