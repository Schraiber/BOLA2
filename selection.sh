for s in $(seq 0.0009 0.0001 0.0024);
	do
	s1=$(bc -l <<< "${s}/2.0")
	python dupe_prob_sel.py ~/Desktop/msmsdir/msms3.2rc-b163.jar -Y 110 -fi 0.0001 -tSel 11280 -nrep 1000000 -sAA ${s} -sAa ${s1} > test_freq_sel_${s}.txt
done
