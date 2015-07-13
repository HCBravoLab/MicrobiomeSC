java -cp lib/*:bin/:../FlexSC/lib/*:../FlexSC/bin util.GenRunnable differentialAbundance.DifferentialAbundance data/kenya_case.txt data/kenya_control.txt $1&

java -cp lib/*:bin/:../FlexSC/lib/*:../FlexSC/bin util.EvaRunnable differentialAbundance.DifferentialAbundance data/gambia_case.txt data/gambia_control.txt $1
