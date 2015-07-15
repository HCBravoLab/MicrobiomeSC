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

package precompute;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

import precompute.StatisticsAlphaDiversity;

public class PrepareDataAlphaDiversity {

	public static void CountFre(StatisticsAlphaDiversity dataCounts, String[] content) {
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

	public static StatisticsAlphaDiversity readFile(String filename) {
		File file = new File(filename);
		Scanner scanner; 
		StatisticsAlphaDiversity sta = new StatisticsAlphaDiversity();
		sta.simpsonsSumOfSquares = 0.0;
		sta.simpsonsIndexTotalSum = 0.0;
		try {
			scanner = new Scanner(file);
			String line = scanner.nextLine();
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
