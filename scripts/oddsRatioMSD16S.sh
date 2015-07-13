java -cp bin:lib/* util.GenRunnable ../FlexSC oddsRatio.OddsRatio data/msd_species_case_part1.txt data/msd_species_control_part1.txt  &

java -cp bin:lib/* util.EvaRunnable ../FlexSC oddsRatio.OddsRatio data/msd_species_case_part2.txt data/msd_species_control_part2.txt  
