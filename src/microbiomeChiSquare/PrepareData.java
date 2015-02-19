package microbiomeChiSquare;


import java.io.File;
import java.io.FileNotFoundException;
import java.util.LinkedList;
import java.util.Scanner;

public class PrepareData {
	
	public static void CountFre(Statistics sta, String[] content) {
		sta.totalNum = content.length;
		for(int i = 0; i< content.length; ++i){
			if(content[i].charAt(0) == '1'){
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
			System.out.println(line);
			d.numberOftuples = line.split(" ").length;
			System.out.println(d.numberOftuples);
			while(scanner.hasNextLine()) {
				line = scanner.nextLine();
				System.out.println("next line in for loop " + line);
				if (line.contains("OTU")){
				Statistics sta = new Statistics();
				line = scanner.nextLine();
				System.out.println(line);
				String[] snp = line.split(" ");
				CountFre(sta, snp);
				lsta.add(sta);
				}
				else{
					continue;
				}
			}
			scanner.close();			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		d.data = lsta.toArray(new Statistics[0]); 
		return d;
	}
	
	public static void main(String[] args) {
		readFile("data/case_chr2_29504091_30044866.txt");
	}
}
