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

package precomputeFilter;


import java.io.File;
import java.io.FileNotFoundException;
import java.util.LinkedList;
import java.util.Scanner;

import precomputeFilter.StatisticsDifferentialAbundance;


public class PrepareDataDifferentialAbundance {
	
	public static void CountFre(StatisticsDifferentialAbundance dataCounts, String[] content) {
		dataCounts.sumOfSquares = 0.0;
		for(int i = 1; i < content.length; i++){
			dataCounts.sumOfSquares += (Double.parseDouble((content[i]))*Double.parseDouble((content[i])));
			dataCounts.totalSum += Double.parseDouble((content[i]));
			if(Double.parseDouble(content[i]) > 0.0){
				dataCounts.numOfPresent++;
			}
		}
		dataCounts.numOfSamples = content.length-1;
	}
	
	public static class StatisticsData {
		public StatisticsDifferentialAbundance[] data;
		public int numOfTuples;
	}
	
	public static StatisticsData readFile(String filename) {
		File file = new File(filename);
		Scanner scanner; 
		StatisticsData d = new StatisticsData(); 
		LinkedList<StatisticsDifferentialAbundance> lsta = new LinkedList<StatisticsDifferentialAbundance>();
		try {
			scanner = new Scanner(file);
			String line = scanner.nextLine();
			while(scanner.hasNextLine()) {
				StatisticsDifferentialAbundance sta = new StatisticsDifferentialAbundance();
				line = scanner.nextLine();
				
				String[] counts = line.split(" ");
				CountFre(sta, counts);
				lsta.add(sta);
			}
			scanner.close();			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		d.data = lsta.toArray(new StatisticsDifferentialAbundance[0]);
		d.numOfTuples = d.data.length;
		return d;
	}
}
