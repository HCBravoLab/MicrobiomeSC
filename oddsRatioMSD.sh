java -cp lib/*:bin/:../FlexSC/lib/*:../FlexSC/bin util.GenRunnable ../FlexSC oddsRatio.OddsRatio data/pgp_species_case_part1_nozeros.txt data/pgp_species_control_part1_nozeros.txt  $1&

java -cp lib/*:bin/:../FlexSC/lib/*:../FlexSC/bin util.EvaRunnable ../FlexSC oddsRatio.OddsRatio data/pgp_species_case_part2_nozeros.txt data/pgp_species_control_part2_nozeros.txt $1
