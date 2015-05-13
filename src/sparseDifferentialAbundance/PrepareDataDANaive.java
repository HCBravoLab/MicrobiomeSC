
package sparseDifferentialAbundance;


import java.io.File;
import java.io.FileNotFoundException;
import java.util.LinkedList;
import java.util.Scanner;

public class PrepareDataDANaive {
	
	public static double[][] readFile(String filename) {
		File file = new File(filename);
		Scanner scanner; 
		LinkedList<LinkedList<Double>> indices = new LinkedList<LinkedList<Double>>();

        int numLines = 0;
        int numFeatures = 0;
        double[][] matrix;
		try {
			scanner = new Scanner(file);
			String line = scanner.nextLine();
			
			while(scanner.hasNextLine()) {
				line = scanner.nextLine();
				numLines++;
				String[] counts = line.split(" ");
				numFeatures = counts.length-1;
			}
			scanner.close();			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		matrix = new double[numLines][numFeatures];
		int i = 0;
		try {
			scanner = new Scanner(file);
			String line = scanner.nextLine();
			
			while(scanner.hasNextLine()) {
				line = scanner.nextLine();
				String[] counts = line.split(" ");
				for(int j =1; j<counts.length; j++){
				matrix[i][j] = Double.parseDouble(counts[j]);
				}
			}
			scanner.close();			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return matrix;
	}
}

