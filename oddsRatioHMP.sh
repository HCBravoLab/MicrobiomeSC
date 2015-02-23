java -cp lib/*:bin/:../FlexSC/lib/*:../FlexSC/bin util.GenRunnable ../FlexSC oddsRatio.OddsRatio data/hmp_species_case_part1.txt data/hmp_species_control_part1.txt  $1&

java -cp lib/*:bin/:../FlexSC/lib/*:../FlexSC/bin util.EvaRunnable ../FlexSC oddsRatio.OddsRatio data/hmp_species_case_part2.txt data/hmp_species_control_part2.txt $1
