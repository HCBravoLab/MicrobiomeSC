package alphaDiversity;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

import alphaDiversity.Statistics;

public class PrepareData {

	public static void CountFre(Statistics dataCounts, String[] content) {
		dataCounts.rowSum = 0.0;
		dataCounts.simpsonsIndex = 0.0;
		for(int i = 1; i < content.length; i++){
			dataCounts.rowSum += Double.parseDouble((content[i]));
		}
		for(int i = 1; i < content.length; i++){
			dataCounts.simpsonsIndex += (Math.pow((Double.parseDouble(content[i])/dataCounts.rowSum),2));
		}
		dataCounts.simpsonsSumOfSquares += Math.pow(dataCounts.simpsonsIndex, 2);
		dataCounts.simpsonsIndexTotalSum += dataCounts.simpsonsIndex;
		dataCounts.numOfSamples = content.length-1;
		dataCounts.numOfSimpsonsIndices++;
	}

	public static Statistics readFile(String filename) {
		File file = new File(filename);
		Scanner scanner; 
		Statistics sta = new Statistics();
		sta.simpsonsSumOfSquares = 0.0;
		sta.simpsonsIndexTotalSum = 0.0;
		//StatisticsData d = new StatisticsData(); 
		//LinkedList<Statistics> lsta = new LinkedList<Statistics>();
		try {
			scanner = new Scanner(file);
			String line = scanner.nextLine();
			//System.out.println(line);
			while(scanner.hasNextLine()) {
				line = scanner.nextLine();
				String[] counts = line.split(" ");
				CountFre(sta, counts);
			}
			scanner.close();			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return sta;
	}
}
