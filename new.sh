java -cp lib/*:bin/: util.GenRunnable ../FlexSC alphaDiversity.AlphaDiversity data/hmp_species_case_part1.txt data/hmp_species_control_part1.txt $1&

java -cp lib/*:bin/: util.EvaRunnable ../FlexSC alphaDiversity.AlphaDiversity data/hmp_species_case_part2.txt data/hmp_species_control_part2.txt $1
