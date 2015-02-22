java -cp bin/:../FlexSC/lib/*:../FlexSC/bin:lib/* util.GenRunnable ../FlexSC differentialAbundance.DifferentialAbundance data/msd_species_case_part1.txt data/msd_species_control_part1.txt  $1&

java -cp bin/:../FlexSC/lib/*:../FlexSC/bin:lib/* util.EvaRunnable ../FlexSC differentialAbundance.DifferentialAbundance data/msd_species_case_part2.txt data/msd_species_control_part2.txt $1
