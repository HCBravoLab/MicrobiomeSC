java -cp lib/*:bin/:../FlexSC/lib/*:../FlexSC/bin util.GenRunnable ../FlexSC differentialAbundance.DifferentialAbundance data/hmp_species_case_part1.txt data/hmp_species_control_part1.txt $1&

java -cp lib/*:bin/:../FlexSC/lib/*:../FlexSC/bin util.EvaRunnable ../FlexSC differentialAbundance.DifferentialAbundance data/hmp_species_case_part2.txt data/hmp_species_control_part2.txt $1
