java -cp lib/*:bin/:../FlexSC/lib/*:../FlexSC/bin/ util.GenRunnable ./ oddsRatio.OddsRatio data/hmp_species_case_part1_nozeros.txt data/hmp_species_control_part1_nozeros.txt  $1&

java -cp lib/*:bin/:../FlexSC/lib/*:../FlexSC/bin/ util.EvaRunnable ./ oddsRatio.OddsRatio data/hmp_species_case_part2_nozeros.txt data/hmp_species_control_part2_nozeros.txt $1
