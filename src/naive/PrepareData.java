/*Original implementation of chi-square performed by xiao wang,
 *  slightly by justin wagner to handle microbiome count data
 */

package naive;


import java.io.File;
import java.io.FileNotFoundException;
import java.util.LinkedList;
import java.util.Scanner;

public class PrepareData {
	public static void CountFre(LinkedList<Integer> indices, String[] content) {
		for(int i = 1; i< content.length; i++){
			if(Double.parseDouble(content[i]) > 0.0){
				indices.add(new Integer(i));
			}
		}
	}
	
	public static int[] readFile(String filename) {
		File file = new File(filename);
		Scanner scanner; 
		LinkedList<Integer> indices = new LinkedList<Integer>();
        int numFeatures = 0;
        int numSamples = 0;
		try {
			scanner = new Scanner(file);
			String line = scanner.nextLine();
			while(scanner.hasNextLine()) {
				line = scanner.nextLine();
				String[] counts = line.split(" ");
				numFeatures = counts.length-1;
				numSamples++;
				CountFre(indices, counts);
				//lsta.add(indices);
			}
			scanner.close();			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		Integer[] ar = indices.toArray(new Integer[0]);
		int[]intAr = new int[ar.length+2];
		intAr[0] = numSamples;
		intAr[1] = numFeatures;
		for(int i = 0; i < ar.length; i++){
			intAr[i+2] = ar[i].intValue();
		}
		return intAr;
	}
}
