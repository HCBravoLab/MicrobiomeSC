java -cp bin/:../FlexSC/lib/*:../FlexSC/bin:lib/* util.GenRunnable ../FlexSC chiSquare.ChiSquare -s data/pgp_species_case_part1.txt -t data/pgp_species_control_part1.txt  $1&

java -cp bin/:../FlexSC/lib/*:../FlexSC/bin:lib/* util.EvaRunnable ../FlexSC chiSquare.ChiSquare -s data/pgp_species_case_part2.txt -t data/pgp_species_control_part2.txt $1
