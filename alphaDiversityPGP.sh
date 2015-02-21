java -cp bin/:../FlexSC/lib/*:../FlexSC/bin util.GenRunnable ../FlexSC alphaDiversity.AlphaDiversity data/pgp_species_case_part1_transpose.txt data/pgp_species_control_part1_transpose.txt  $1&

java -cp bin/:../FlexSC/lib/*:../FlexSC/bin util.EvaRunnable ../FlexSC alphaDiversity.AlphaDiversity data/pgp_species_case_part2_transpose.txt data/pgp_species_control_part2_transpose.txt $1
