java -cp lib/*:bin/ util.GenRunnable ../FlexSC differentialAbundance.DifferentialAbundance data/kenya_case.txt data/kenya_control.txt $1&

java -cp lib/*:bin/ util.EvaRunnable ../FlexSC differentialAbundance.DifferentialAbundance data/gambia_case.txt data/gambia_control.txt $1
