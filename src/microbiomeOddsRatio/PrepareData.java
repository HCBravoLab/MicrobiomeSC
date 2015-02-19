package microbiomeOddsRatio;


import java.io.File;
import java.io.FileNotFoundException;
import java.util.LinkedList;
import java.util.Scanner;

import microbiomeOddsRatio.Statistics;
import microbiomeOddsRatio.PrepareData.StatisticsData;

public class PrepareData {
	
	public static void CountFre(Statistics sta, String[] content) {
		sta.totalNum = content.length-1;
		for(int i = 1; i< content.length; i++){
			if(Double.parseDouble(content[i]) > 0.0){
				sta.numOfPresent++;
			}
		}
	}
	
	public static class StatisticsData {
		public Statistics[] data;
		public int numberOftuples;
	}
	
	public static StatisticsData readFile(String filename) {
		File file = new File(filename);
		Scanner scanner; 
		StatisticsData d = new StatisticsData(); 
		LinkedList<Statistics> lsta = new LinkedList<Statistics>();
		try {
			scanner = new Scanner(file);
			String line = scanner.nextLine();
			while(scanner.hasNextLine()) {
				Statistics sta = new Statistics();
				line = scanner.nextLine();
				String[] counts = line.split(" ");
				CountFre(sta, counts);
				lsta.add(sta);
			}
			scanner.close();			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		d.data = lsta.toArray(new Statistics[0]);
		d.numberOftuples = d.data.length;
		return d;
	}
	
	public static void main(String[] args) {
		readFile("data/case_chr2_29504091_30044866.txt");
	}
}
