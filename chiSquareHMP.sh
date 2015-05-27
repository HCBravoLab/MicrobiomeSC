java -cp lib/*:bin/:../FlexSC/lib/*:../FlexSC/bin util.GenRunnable chiSquare.ChiSquare -s data/hmp_species_case_part1_nozeros.txt -t data/hmp_species_control_part1_nozeros.txt  $1&

java -cp lib/*:bin/:../FlexSC/lib/*:../FlexSC/bin util.EvaRunnable chiSquare.ChiSquare -s data/hmp_species_case_part2_nozeros.txt -t data/hmp_species_control_part2_nozeros.txt $1
