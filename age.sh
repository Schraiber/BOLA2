for age in $(seq 8000 500 50000);
	do
	python dupe_prob_sel.py ~/Desktop/msmsdir/msms3.2rc-b163.jar -Y 110 -fi 0.001 -tSel ${age} -nrep 200000 -sAA 0.0 -sAa 0.0 > test_freq_age_${age}.txt
done
