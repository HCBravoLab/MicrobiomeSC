java -cp lib/*:bin/:../FlexSC/lib/*:../FlexSC/bin util.GenRunnable ../FlexSC alphaDiversity.AlphaDiversity data/hmp_species_case_part1_nozeros_transpose.txt data/hmp_species_control_part1_nozeros_transpose.txt  $1&

java -cp lib/*:bin/:../FlexSC/lib/*:../FlexSC/bin util.EvaRunnable ../FlexSC alphaDiversity.AlphaDiversity data/hmp_species_case_part2_nozeros_transpose.txt data/hmp_species_control_part2_nozeros_transpose.txt $1
