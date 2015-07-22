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


import java.util.HashSet;

import util.EvaRunnable;
import util.GenRunnable;
import util.Utils;
import circuits.CircuitLib;
import circuits.arithmetic.FixedPointLib;
import circuits.arithmetic.FloatLib;
import circuits.arithmetic.IntegerLib;
import flexsc.CompEnv;

import org.apache.commons.math3.distribution.TDistribution;

import precomputeFilter.PrepareDataDifferentialAbundance;
import precomputeFilter.PrepareDataDifferentialAbundance.StatisticsData;

public class DifferentialAbundance {
	static public int Width = 32;
	static public int FWidth = 54;
	static public int FOffset = 11;
	static int filterThreshold = 5;

	static public<T> T[] filter(CompEnv<T> gen, T[][][] aliceCasePresent,
			T[][][] bobCasePresent,
			T[][][] aliceControlPresent,
			T[][][] bobControlPresent, int numOfTests){

		FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
		IntegerLib<T> ilib = new IntegerLib<T>(gen, Width);
		T[] threshold = flib.publicValue(filterThreshold);
		T[] zero = flib.publicValue(0.0);

		T[] filterResults = gen.newTArray(numOfTests);
		T[] caseNum;
		T[] controlNum;
		T[] totalPresent;
		T aboveThreshold;
		for(int i = 0; i < numOfTests; i++){			
			caseNum = flib.add(aliceCasePresent[i][3], bobCasePresent[i][3]);
			controlNum = flib.add(aliceControlPresent[i][3], bobControlPresent[i][3]);
			totalPresent = flib.add(caseNum, controlNum);
			aboveThreshold = ilib.not(flib.leq(totalPresent, threshold));
			filterResults[i] = aboveThreshold;
		}
		
		return filterResults;
	}
	
	public static<T> T[][][] compute(CompEnv<T> gen, T[][][] aliceCase,
			T[][][] bobCase,
			T[][][] aliceControl,
			T[][][] bobControl, int numOfTests) {	

		T[][][] res = gen.newTArray(numOfTests, 2, 0);
		FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
		IntegerLib<T> ilib = new IntegerLib<T>(gen, Width);
		T[] zero = flib.publicValue(0.0);
		T[] one = flib.publicValue(1.0);
		T[] negOne = flib.publicValue(-1.0);

		CircuitLib<T> cl = new CircuitLib<T>(gen);

		for(int i = 0; i < numOfTests; i++) {

			T[] caseSum = flib.add(aliceCase[i][0], bobCase[i][0]);
			T[] controlSum = flib.add(aliceControl[i][0], bobControl[i][0]);
			
			T[] caseSumOfSquares = flib.add(aliceCase[i][1], bobCase[i][1]);
			T[] controlSumOfSquares = flib.add(aliceControl[i][1], bobControl[i][1]);

			T[] caseNum = flib.add(aliceCase[i][2], bobCase[i][2]);
			T[] controlNum  = flib.add(aliceControl[i][2], bobControl[i][2]);
			
			T[] caseNumPresent = flib.add(aliceCase[i][3], bobCase[i][3]);
			T[] controlNumPresent= flib.add(aliceControl[i][3], bobControl[i][3]);

			T[] zeroInt = flib.publicValue(0.0);
			T caseIsZero = flib.eq(caseNumPresent, zeroInt);
			T controlIsZero = flib.eq(controlNumPresent, zeroInt);
			T oneIsZero = ilib.or(caseIsZero, controlIsZero);
			T sendOneOrNegOne = ilib.not(flib.leq(caseNumPresent, controlNumPresent));
			T[] tStatZero = cl.mux(one, negOne, sendOneOrNegOne);
			
			T[] caseMeanAbundance = flib.div(caseSum, caseNum);
			T[] caseVarianceSecondTerm = flib.div(flib.multiply(caseSum, caseSum), caseNum);
			T[] caseVariance = flib.div(flib.sub(caseSumOfSquares, caseVarianceSecondTerm), caseNum);
			
			T[] controlMeanAbundance = flib.div(controlSum, controlNum);		    
			T[] controlVarianceSecondTerm = flib.div(flib.multiply(controlSum, controlSum), controlNum);
			T[] controlVariance = flib.div(flib.sub(controlSumOfSquares, controlVarianceSecondTerm), controlNum);

			T[] tUpper = flib.sub(controlMeanAbundance, caseMeanAbundance);
			T tLowerFirstIsZero = flib.eq(caseNumPresent, zero);
			T[] tLowerFirst = flib.div(caseVariance, caseNum);
			tLowerFirst = cl.mux(tLowerFirst, zero, tLowerFirstIsZero);
			
			T tLowerSecondIsZero = flib.eq(controlNumPresent, zero);
			T[] tLowerSecond = flib.div(controlVariance, controlNum);
			tLowerSecond = cl.mux(tLowerSecond, zero, tLowerSecondIsZero);
			
			T[] tLowerSqrt = flib.sqrt(flib.add(tLowerFirst,tLowerSecond));
			T[] tStatNonZero = flib.div(tUpper, tLowerSqrt);
			
			T[] degreesOfFreedomTop = flib.add(tLowerFirst, tLowerSecond);
			degreesOfFreedomTop = flib.multiply(degreesOfFreedomTop, degreesOfFreedomTop);

			T[] degreesOfFreedomBottomFirst = flib.div(caseVariance, flib.sub(caseNum, one));
			degreesOfFreedomBottomFirst = flib.multiply(degreesOfFreedomBottomFirst,degreesOfFreedomBottomFirst);
			degreesOfFreedomBottomFirst = flib.div(degreesOfFreedomBottomFirst, flib.sub(caseNum, one));

			T[] degreesOfFreedomBottomSecond = flib.div(controlVariance, flib.sub(controlNum, one));
			degreesOfFreedomBottomSecond = flib.multiply(degreesOfFreedomBottomSecond,degreesOfFreedomBottomSecond);
			degreesOfFreedomBottomSecond = flib.div(degreesOfFreedomBottomSecond, flib.sub(controlNum, one));

			T[] degreesOfFreedom = flib.div(degreesOfFreedomTop, flib.add(degreesOfFreedomBottomFirst,degreesOfFreedomBottomSecond));
			
			T[] oneInt = flib.publicValue(1.0);
			
			T caseIsOne = flib.eq(caseNumPresent, oneInt);
			T controlIsOne = flib.eq(controlNumPresent, oneInt);
			T greaterThanOne = ilib.and(ilib.xor(caseIsOne, controlIsOne), oneIsZero);
			
			T returnCase = ilib.not(flib.leq(caseNumPresent, controlNumPresent));
				
			T[] tStat = cl.mux(tStatNonZero, tStatZero,  greaterThanOne);
			T[] dfZero = cl.mux(flib.sub(controlNum, one), flib.sub(caseNum, one), returnCase);
			T[] df = cl.mux(degreesOfFreedom, dfZero, oneIsZero);

			res[i][0] = tStat;
			res[i][1] = df;
		}
		return res;
	}

	public static class Generator<T> extends GenRunnable<T> {
		T[][][] aliceCase;
		T[][][] bobCase;
		T[][][] aliceControl;
		T[][][] bobControl;

		T[][][] aliceCaseFiltered;
		T[][][] bobCaseFiltered;
		T[][][] aliceControlFiltered;
		T[][][] bobControlFiltered;

		T[][][] aliceCaseHolder;
		T[][][] aliceControlHolder;
		T[][][] bobCaseHolder;
		T[][][] bobControlHolder;
		T[] l;
		int numOfTests;
		@Override
		public void prepareInput(CompEnv<T> gen) {
			StatisticsData caseSta = PrepareDataDifferentialAbundance.readFile(args[0]);
			StatisticsData controlSta = PrepareDataDifferentialAbundance.readFile(args[1]);
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
			l = flib.publicValue(0.0);
			System.out.println(caseSta.numOfTuples);
			System.out.println(controlSta.numOfTuples);

			boolean[][][] caseData = new boolean[caseSta.numOfTuples][4][l.length];
			boolean[][][] caseData2 = new boolean[caseSta.numOfTuples][4][l.length];
			boolean[][][] caseData3 = new boolean[caseSta.numOfTuples][4][l.length];

			for(int i = 0; i < caseSta.numOfTuples; i++) {
				caseData[i][0] = Utils.fromFloat(caseSta.data[i].totalSum, FWidth, FOffset);
				caseData[i][1] = Utils.fromFloat(caseSta.data[i].sumOfSquares, FWidth, FOffset);
				caseData[i][2] = Utils.fromFloat(caseSta.data[i].numOfSamples, FWidth, FOffset);
				caseData[i][3] = Utils.fromFloat(caseSta.data[i].numOfPresent, FWidth, FOffset);

				caseData2[i][0] = Utils.fromFloat(caseSta.data[i].totalSum, FWidth, FOffset);
				caseData2[i][1] = Utils.fromFloat(caseSta.data[i].sumOfSquares, FWidth, FOffset);
				caseData2[i][2] = Utils.fromFloat(caseSta.data[i].numOfSamples, FWidth, FOffset);
				caseData2[i][3] = Utils.fromFloat(caseSta.data[i].numOfPresent, FWidth, FOffset);

				caseData3[i][0] = Utils.fromFloat(caseSta.data[i].totalSum, FWidth, FOffset);
				caseData3[i][1] = Utils.fromFloat(caseSta.data[i].sumOfSquares, FWidth, FOffset);
				caseData3[i][2] = Utils.fromFloat(caseSta.data[i].numOfSamples, FWidth, FOffset);
				caseData3[i][3] = Utils.fromFloat(caseSta.data[i].numOfPresent, FWidth, FOffset);
			}

			boolean[][][] controlData = new boolean[controlSta.numOfTuples][4][l.length];
			boolean[][][] controlData2 = new boolean[controlSta.numOfTuples][4][l.length];
			boolean[][][] controlData3 = new boolean[controlSta.numOfTuples][4][l.length];

			for(int i = 0; i < controlSta.numOfTuples; i++) {
				controlData[i][0] = Utils.fromFloat(controlSta.data[i].totalSum, FWidth, FOffset);
				controlData[i][1] = Utils.fromFloat(controlSta.data[i].sumOfSquares, FWidth, FOffset);
				controlData[i][2] = Utils.fromFloat(controlSta.data[i].numOfSamples, FWidth, FOffset);
				controlData[i][3] = Utils.fromFloat(controlSta.data[i].numOfPresent, FWidth, FOffset);
				
				controlData2[i][0] = Utils.fromFloat(controlSta.data[i].totalSum, FWidth, FOffset);
				controlData2[i][1] = Utils.fromFloat(controlSta.data[i].sumOfSquares, FWidth, FOffset);
				controlData2[i][2] = Utils.fromFloat(controlSta.data[i].numOfSamples, FWidth, FOffset);
				controlData2[i][3] = Utils.fromFloat(controlSta.data[i].numOfPresent, FWidth, FOffset);

				controlData3[i][0] = Utils.fromFloat(controlSta.data[i].totalSum, FWidth, FOffset);
				controlData3[i][1] = Utils.fromFloat(controlSta.data[i].sumOfSquares, FWidth, FOffset);
				controlData3[i][2] = Utils.fromFloat(controlSta.data[i].numOfSamples, FWidth, FOffset);
				controlData3[i][3] = Utils.fromFloat(controlSta.data[i].numOfPresent, FWidth, FOffset);
			}
			aliceCase = gen.inputOfAlice(caseData);
			aliceControl = gen.inputOfAlice(controlData);
			bobCase = gen.inputOfBob(new boolean[caseSta.numOfTuples][4][l.length]);
			bobControl = gen.inputOfBob(new boolean[controlSta.numOfTuples][4][l.length]);

			aliceCaseFiltered = gen.inputOfAlice(caseData2);
			aliceControlFiltered = gen.inputOfAlice(controlData2);
			bobCaseFiltered = gen.inputOfBob(new boolean[caseSta.numOfTuples][4][l.length]);
			bobControlFiltered = gen.inputOfBob(new boolean[controlSta.numOfTuples][4][l.length]);

			aliceCaseHolder = gen.inputOfAlice(caseData3);
			aliceControlHolder = gen.inputOfAlice(controlData3);
			bobCaseHolder = gen.inputOfBob(new boolean[caseSta.numOfTuples][4][l.length]);
			bobControlHolder = gen.inputOfBob(new boolean[controlSta.numOfTuples][4][l.length]);
			numOfTests = caseSta.numOfTuples;
		}

		T[][][] res;
		int[] indices;
		int numFiltered;
		T[]filterRes;

		@Override
		public void secureCompute(CompEnv<T> gen) {
			numFiltered = 0;
			CircuitLib<T> cl = new CircuitLib<T>(gen);
			filterRes = filter(gen, aliceCase, bobCase, aliceControl, bobControl, numOfTests);
			boolean[] filResOut = gen.outputToAlice(filterRes);
			for(int i =0; i < numOfTests; i++){
				gen.channel.writeBoolean(filResOut[i]);
				gen.channel.flush();
			}

			indices = new int[numOfTests];
			
			for(int i = 0; i < numOfTests; i++){
				indices[i] = -1;
				if(!filResOut[i]){
					indices[numFiltered] = i;
					numFiltered++;
				}
				
			}
			gen.flush();

			System.out.println(numFiltered);
			System.out.println(" ");
			for(int i = 0; i < numFiltered; i++){				
				System.arraycopy(aliceCaseFiltered[i][0], 0, aliceCaseHolder[indices[i]][0], 0, l.length);
				System.arraycopy(aliceCaseFiltered[i][1], 0, aliceCaseHolder[indices[i]][1], 0, l.length);

				System.arraycopy(aliceControlFiltered[i][0], 0, aliceControlHolder[indices[i]][0], 0, l.length);
				System.arraycopy(aliceControlFiltered[i][1],0, aliceControlHolder[indices[i]][1], 0, l.length);
			}

			res = compute(gen, aliceCaseFiltered, bobCaseFiltered, aliceControlFiltered, bobControlFiltered, numFiltered);
			gen.flush();

		}

		@Override
		public void prepareOutput(CompEnv<T> gen) {
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
			System.out.println("tstat,df,pval");
			boolean[] outTstat;
			boolean[] outDF;
			int counter = 0;
			HashSet<Integer> indicesArL = new HashSet<Integer>();
			for(int i =0; i < numOfTests; i++){
				indicesArL.add(indices[i]);
			}
			CircuitLib<T> cl = new CircuitLib<T>(gen);
			for(int i = 0; i < numOfTests; i++){
				double tStat;
				double df;

				if(!indicesArL.contains(i)){
					tStat = 0.0;
					df = 0.0;
				}
				else{
					outTstat = gen.outputToAlice(res[counter][0]);
					outDF = gen.outputToAlice(res[counter][1]);
					tStat = Utils.toFloat(outTstat, FWidth, FOffset);
					df = Utils.toFloat(outDF, FWidth, FOffset);
					counter++;
				}
				System.out.println(tStat + "," + df);
			}
		}
	}

	public static class Evaluator<T> extends EvaRunnable<T> {
		T[][][] aliceCase;
		T[][][] bobCase;
		T[][][] aliceControl;
		T[][][] bobControl;
		
		T[][][] aliceCaseFiltered;
		T[][][] bobCaseFiltered;
		T[][][] aliceControlFiltered;
		T[][][] bobControlFiltered;

		T[][][] aliceCaseHolder;
		T[][][] aliceControlHolder;
		T[][][] bobCaseHolder;
		T[][][] bobControlHolder;
		int numOfTests;
		T[] l;
		@Override
		public void prepareInput(CompEnv<T> gen) {
			StatisticsData caseSta = PrepareDataDifferentialAbundance.readFile(args[0]);
			StatisticsData controlSta = PrepareDataDifferentialAbundance.readFile(args[1]);			
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
			l = flib.publicValue(0.0);

			System.out.println(caseSta.numOfTuples);
			System.out.println(controlSta.numOfTuples);
			boolean[][][] caseData = new boolean[caseSta.numOfTuples][4][l.length];
			boolean[][][] caseData2 = new boolean[caseSta.numOfTuples][4][l.length];
			boolean[][][] caseData3 = new boolean[caseSta.numOfTuples][4][l.length];

			for(int i = 0; i < caseSta.numOfTuples; i++) {
				caseData[i][0] = Utils.fromFloat(caseSta.data[i].totalSum, FWidth, FOffset);
				caseData[i][1] = Utils.fromFloat(caseSta.data[i].sumOfSquares, FWidth, FOffset);
				caseData[i][2] = Utils.fromFloat(caseSta.data[i].numOfSamples, FWidth, FOffset);
				caseData[i][3] = Utils.fromFloat(caseSta.data[i].numOfPresent, FWidth, FOffset);

				caseData2[i][0] = Utils.fromFloat(caseSta.data[i].totalSum, FWidth, FOffset);
				caseData2[i][1] = Utils.fromFloat(caseSta.data[i].sumOfSquares, FWidth, FOffset);
				caseData2[i][2] = Utils.fromFloat(caseSta.data[i].numOfSamples, FWidth, FOffset);
				caseData2[i][3] = Utils.fromFloat(caseSta.data[i].numOfPresent, FWidth, FOffset);

				caseData3[i][0] = Utils.fromFloat(caseSta.data[i].totalSum, FWidth, FOffset);
				caseData3[i][1] = Utils.fromFloat(caseSta.data[i].sumOfSquares, FWidth, FOffset);
				caseData3[i][2] = Utils.fromFloat(caseSta.data[i].numOfSamples, FWidth, FOffset);
				caseData3[i][3] = Utils.fromFloat(caseSta.data[i].numOfPresent, FWidth, FOffset);
			}

			boolean[][][] controlData = new boolean[controlSta.numOfTuples][4][l.length];
			boolean[][][] controlData2 = new boolean[controlSta.numOfTuples][4][l.length];
			boolean[][][] controlData3 = new boolean[controlSta.numOfTuples][4][l.length];

			for(int i = 0; i < controlSta.numOfTuples; i++) {
				controlData[i][0] = Utils.fromFloat(controlSta.data[i].totalSum, FWidth, FOffset);
				controlData[i][1] = Utils.fromFloat(controlSta.data[i].sumOfSquares, FWidth, FOffset);
				controlData[i][2] = Utils.fromFloat(controlSta.data[i].numOfSamples, FWidth, FOffset);
				controlData[i][3] = Utils.fromFloat(controlSta.data[i].numOfPresent, FWidth, FOffset);
				
				controlData2[i][0] = Utils.fromFloat(controlSta.data[i].totalSum, FWidth, FOffset);
				controlData2[i][1] = Utils.fromFloat(controlSta.data[i].sumOfSquares, FWidth, FOffset);
				controlData2[i][2] = Utils.fromFloat(controlSta.data[i].numOfSamples, FWidth, FOffset);
				controlData2[i][3] = Utils.fromFloat(controlSta.data[i].numOfPresent, FWidth, FOffset);

				controlData3[i][0] = Utils.fromFloat(controlSta.data[i].totalSum, FWidth, FOffset);
				controlData3[i][1] = Utils.fromFloat(controlSta.data[i].sumOfSquares, FWidth, FOffset);
				controlData3[i][2] = Utils.fromFloat(controlSta.data[i].numOfSamples, FWidth, FOffset);
				controlData3[i][3] = Utils.fromFloat(controlSta.data[i].numOfPresent, FWidth, FOffset);
			}			
			aliceCase = gen.inputOfAlice(caseData);
			aliceControl = gen.inputOfAlice(controlData);
			bobCase = gen.inputOfBob(new boolean[caseSta.numOfTuples][4][l.length]);
			bobControl = gen.inputOfBob(new boolean[controlSta.numOfTuples][4][l.length]);

			aliceCaseFiltered = gen.inputOfAlice(caseData2);
			aliceControlFiltered = gen.inputOfAlice(controlData2);
			bobCaseFiltered = gen.inputOfBob(new boolean[caseSta.numOfTuples][4][l.length]);
			bobControlFiltered = gen.inputOfBob(new boolean[controlSta.numOfTuples][4][l.length]);

			aliceCaseHolder = gen.inputOfAlice(caseData3);
			aliceControlHolder = gen.inputOfAlice(controlData3);
			bobCaseHolder = gen.inputOfBob(new boolean[caseSta.numOfTuples][4][l.length]);
			bobControlHolder = gen.inputOfBob(new boolean[controlSta.numOfTuples][4][l.length]);
			numOfTests = caseSta.numOfTuples;
		}
		T[][][] res;
		int[] indices;
		T[]filterRes;

		int numFiltered;
		@Override
		public void secureCompute(CompEnv<T> gen) {
			numFiltered = 0;
			CircuitLib<T> cl = new CircuitLib<T>(gen);
			filterRes = filter(gen, aliceCase, bobCase, aliceControl, bobControl, numOfTests);
			gen.outputToAlice(filterRes);
			boolean[] filResOut = new boolean[numOfTests];
			for(int i = 0; i < numOfTests; i++){
				filResOut[i] = gen.channel.readBoolean();
				gen.channel.flush();
			}
			indices = new int[numOfTests];
			for(int i = 0; i < numOfTests; i++){
				indices[i] = -1;
				if(filResOut[i]){
					indices[numFiltered] = i;
					numFiltered++;
				}
			}
			gen.flush();
			System.out.println(numFiltered);
			System.out.println(" ");

			System.out.println("d");

			for(int i = 0; i < numFiltered; i++){				
				System.arraycopy(bobCaseFiltered[i][0], 0, bobCaseHolder[indices[i]][0], 0, l.length);
				System.arraycopy(bobCaseFiltered[i][1], 0, bobCaseHolder[indices[i]][1], 0, l.length);

				System.arraycopy(bobControlFiltered[i][0], 0, bobControlHolder[indices[i]][0], 0, l.length);
				System.arraycopy(bobControlFiltered[i][1], 0, bobControlHolder[indices[i]][1], 0, l.length);
			}
			System.out.println("d2 ");

			res = compute(gen, aliceCaseFiltered, bobCaseFiltered, aliceControlFiltered, bobControlFiltered, numFiltered);
			gen.flush();

		}


		@Override
		public void prepareOutput(CompEnv<T> gen) {
			FixedPointLib<T> flib = new FixedPointLib<T>(gen, FWidth, FOffset);
			int counter = 0;
			HashSet<Integer> indicesArL = new HashSet<Integer>();
			for(int i =0; i < numFiltered; i++){
				indicesArL.add(indices[i]);
			}
			for(int i = 0; i < numOfTests; i++){
				if(!indicesArL.contains(i)){
					continue;
				}
				else{
					Utils.toFloat(gen.outputToAlice(res[counter][0]), FWidth, FOffset);
					Utils.toFloat(gen.outputToAlice(res[counter][1]), FWidth, FOffset);
					counter++;
				}
			}
		}
	}
}
