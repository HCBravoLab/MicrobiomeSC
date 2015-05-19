package naive;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class PrepareDataDANaive {
	
	public static double[][] readFile(String filename) {
		File file = new File(filename);
		Scanner scanner; 
		int numRows = 0;
		int numColumns = 0;
		try {
			scanner = new Scanner(file);
			String line = scanner.nextLine();
			while(scanner.hasNextLine()) {
				numRows++;
				line = scanner.nextLine();
				String[] counts = line.split(" ");
				numColumns = counts.length-1;
			}
			scanner.close();			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		double[][] mat = new double[numRows][numColumns];
		int row = 0;
		try {
			scanner = new Scanner(file);
			String line = scanner.nextLine();
			while(scanner.hasNextLine()) {
				line = scanner.nextLine();
				String[] counts = line.split(" ");
				for(int i = 0; i < counts.length-1; i++){
					mat[row][i] = Double.parseDouble(counts[i+1]);
				}
				row++;
			}
			scanner.close();			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

    	return mat;
	}
}
