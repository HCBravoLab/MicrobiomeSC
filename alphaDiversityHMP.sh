java -cp bin/:../FlexSC/lib/*:../FlexSC/bin:lib/* util.GenRunnable ../FlexSC alphaDiversity.AlphaDiversity data/hmp_species_case_part1_transpose.txt data/hmp_species_control_part1_transpose.txt  $1&

java -cp bin/:../FlexSC/lib/*:../FlexSC/bin:lib/* util.EvaRunnable ../FlexSC alphaDiversity.AlphaDiversity data/hmp_species_case_part2_transpose.txt data/hmp_species_control_part2_transpose.txt $1
