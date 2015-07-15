/*
The MIT License (MIT)

Copyright (c) <2015> <Justin Wagner>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/*Original implementation of chi-square performed by xiao wang,
 *  modified by justin wagner to handle microbiome count data
 */

package perFeatureNaive;


import java.io.File;
import java.io.FileNotFoundException;
import java.util.LinkedList;
import java.util.Scanner;

public class PrepareDataDA {
	public static void CountFre(LinkedList<Double> indices, LinkedList<Double> values, String[] content) {
		for(int i = 1; i< content.length; i++){
			if(Double.parseDouble(content[i]) > 0.0){
				indices.add(new Double(i));
				values.add(new Double(content[i]));
			}
		}
	}
	
	public static double[][] readFile(String filename) {
		File file = new File(filename);
		Scanner scanner; 
		LinkedList<Double> indices = new LinkedList<Double>();
		LinkedList<Double> values = new LinkedList<Double>();

        float numFeatures = 0;
        float numSamples = 0;
		try {
			scanner = new Scanner(file);
			String line = scanner.nextLine();
			while(scanner.hasNextLine()) {
				line = scanner.nextLine();
				String[] counts = line.split(" ");
				numFeatures = (float)counts.length-1;
				numSamples++;
				CountFre(indices, values, counts);
				//lsta.add(indices);
			}
			scanner.close();			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		Double[] indicesAr = indices.toArray(new Double[0]);		
		Double[] valuesAr = values.toArray(new Double[0]);
		double[][] out = new double[2][indicesAr.length+2];

		out[0][0] = numSamples;
		out[1][0] = numSamples;

		System.out.println("num samples PrepareData :" + numSamples);
		out[0][1] = numFeatures;
		out[1][1] = numFeatures;

		System.out.println("num features PrepareData :" + numFeatures);

		for(int i = 0; i < indicesAr.length; i++){
			out[0][i+2] = indicesAr[i].doubleValue();
			out[1][i+2] =valuesAr[i].doubleValue();
		}
	
		return out;
	}
}
