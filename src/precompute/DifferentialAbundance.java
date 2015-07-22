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


import util.EvaRunnable;
import util.GenRunnable;
import util.Utils;
import circuits.CircuitLib;
import circuits.arithmetic.FixedPointLib;
import circuits.arithmetic.FloatLib;
import circuits.arithmetic.IntegerLib;
import flexsc.CompEnv;

import org.apache.commons.math3.distribution.TDistribution;

import precompute.PrepareDataDifferentialAbundance;
import precompute.PrepareDataDifferentialAbundance.StatisticsData;


public class DifferentialAbundance {
	static public int Width = 32;
	static public int FWidth = 54;
	static public int FOffset = 11;

	public static<T> T[][][] compute(CompEnv<T> gen, T[][][] aliceCase,
			T[][][] bobCase,
			T[][][] aliceControl,
			T[][][] bobControl, T[][] numAliceCase, T[][] numBobCase, 
			T[][] numAliceControl, T[][] numBobControl, int numOfTests) {	

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
			
			T[] caseNumPresent = ilib.add(numAliceCase[i], numBobCase[i]);
			T[] controlNumPresent= ilib.add(numAliceControl[i], numBobControl[i]);

			T[] zeroInt = ilib.publicValue(0.0);
			T caseIsZero = ilib.eq(caseNumPresent, zeroInt);
			T controlIsZero = ilib.eq(controlNumPresent, zeroInt);
			T oneIsZero = ilib.or(caseIsZero, controlIsZero);
			T sendOneOrNegOne = ilib.geq(caseNumPresent, controlNumPresent);
			T[] tStatZero = cl.mux(one, negOne, sendOneOrNegOne);
			
			T[] caseMeanAbundance = flib.div(caseSum, caseNum);
			T[] caseVarianceSecondTerm = flib.div(flib.multiply(caseSum, caseSum), caseNum);
			T[] caseVariance = flib.div(flib.sub(caseSumOfSquares, caseVarianceSecondTerm), caseNum);
			
			T[] controlMeanAbundance = flib.div(controlSum, controlNum);		    
			T[] controlVarianceSecondTerm = flib.div(flib.multiply(controlSum, controlSum), controlNum);
			T[] controlVariance = flib.div(flib.sub(controlSumOfSquares, controlVarianceSecondTerm), controlNum);
			
			
			T[] tUpper = flib.sub(controlMeanAbundance, caseMeanAbundance);
			T tLowerFirstIsZero = ilib.eq(caseNumPresent, zero);
			T[] tLowerFirst = flib.div(caseVariance, caseNum);
			tLowerFirst = cl.mux(tLowerFirst, zero, tLowerFirstIsZero);
			
			T tLowerSecondIsZero = ilib.eq(controlNumPresent, zero);
			
			T[] tLowerSecond = flib.div(controlVariance, controlNum);
			tLowerSecond = cl.mux(tLowerSecond, zero, tLowerSecondIsZero);
			
			T[] tLowerSqrt = flib.sqrt(flib.add(flib.add(zero,tLowerFirst), flib.add(zero,tLowerSecond)));
			T[] tStatNonZero = flib.div(tUpper, tLowerSqrt);
			
			T[] degreesOfFreedomTop = flib.add(tLowerFirst, tLowerSecond);
			degreesOfFreedomTop = flib.multiply(flib.add(zero,degreesOfFreedomTop), flib.add(zero,degreesOfFreedomTop));

			T[] degreesOfFreedomBottomFirst = flib.div(caseVariance, flib.sub(caseNum, one));
			degreesOfFreedomBottomFirst = flib.multiply(flib.add(zero,degreesOfFreedomBottomFirst), flib.add(zero,degreesOfFreedomBottomFirst));
			degreesOfFreedomBottomFirst = flib.div(degreesOfFreedomBottomFirst, flib.sub(caseNum, one));

			T[] degreesOfFreedomBottomSecond = flib.div(controlVariance, flib.sub(controlNum, one));
			degreesOfFreedomBottomSecond = flib.multiply(flib.add(zero,degreesOfFreedomBottomSecond), flib.add(zero,degreesOfFreedomBottomSecond));
			degreesOfFreedomBottomSecond = flib.div(degreesOfFreedomBottomSecond, flib.sub(controlNum, one));

			T[] degreesOfFreedom = flib.div(degreesOfFreedomTop, flib.add(degreesOfFreedomBottomFirst, degreesOfFreedomBottomSecond));
			T[] oneInt = ilib.publicValue(1.0);
			
			T caseIsOne = ilib.eq(caseNumPresent, oneInt);
			T controlIsOne = ilib.eq(controlNumPresent, oneInt);
			T greaterThanOne = ilib.and(ilib.xor(caseIsOne, controlIsOne), oneIsZero);
			
			T returnCase = ilib.not(ilib.leq(caseNumPresent, controlNumPresent));
				
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
		T[][] numAliceCase;
		T[][] numBobCase;
		T[][] numAliceControl;
		T[][] numBobControl;

		int numOfTests;
		@Override
		public void prepareInput(CompEnv<T> gen) {
			StatisticsData caseSta = PrepareDataDifferentialAbundance.readFile(args[0]);
			StatisticsData controlSta = PrepareDataDifferentialAbundance.readFile(args[1]);
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
			T[] l = flib.publicValue(0.0);

			boolean[][][] caseData = new boolean[caseSta.numOfTuples][3][l.length];
			boolean[][] numCaseData = new boolean[caseSta.numOfTuples][Width];
			for(int i = 0; i < caseSta.numOfTuples; i++) {
				caseData[i][0] = Utils.fromFloat(caseSta.data[i].totalSum, FWidth, FOffset);
				caseData[i][1] = Utils.fromFloat(caseSta.data[i].sumOfSquares, FWidth, FOffset);
				caseData[i][2] = Utils.fromFloat(caseSta.data[i].numOfSamples, FWidth, FOffset);
				numCaseData[i] = Utils.fromInt(caseSta.data[i].numOfPresent, Width);
			}

			boolean[][][] controlData = new boolean[controlSta.numOfTuples][3][l.length];
			boolean[][] numControlData = new boolean[controlSta.numOfTuples][Width];

			for(int i = 0; i < controlSta.numOfTuples; i++) {
				controlData[i][0] = Utils.fromFloat(controlSta.data[i].totalSum, FWidth, FOffset);
				controlData[i][1] = Utils.fromFloat(controlSta.data[i].sumOfSquares, FWidth, FOffset);
				controlData[i][2] = Utils.fromFloat(controlSta.data[i].numOfSamples, FWidth, FOffset);
				numControlData[i] = Utils.fromInt(controlSta.data[i].numOfPresent, Width);
			}
			aliceCase = gen.inputOfAlice(caseData);
			aliceControl = gen.inputOfAlice(controlData);
			bobCase = gen.inputOfBob(caseData);
			bobControl = gen.inputOfBob(controlData);

			numAliceCase = gen.inputOfAlice(numCaseData);
			numAliceControl = gen.inputOfAlice(numControlData);
			numBobCase = gen.inputOfBob(numCaseData);
			numBobControl = gen.inputOfBob(numControlData);
			numOfTests = caseSta.numOfTuples;
		}

		T[][][] res;
		@Override
		public void secureCompute(CompEnv<T> gen) {
			res = compute(gen, aliceCase, bobCase, aliceControl, bobControl, numAliceCase, numBobCase,
					numAliceControl, numBobControl,numOfTests);
		}

		@Override
		public void prepareOutput(CompEnv<T> gen) {
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);

			System.out.println("t-stat,Degrees Of Freedom,p-value");
			for(int i = 0; i < numOfTests; i++){
				double tStat = flib.outputToAlice(res[i][0]);
				double df = flib.outputToAlice(res[i][1]);
				if (tStat == 0.0){
					System.out.println("NA,NA,NA");
					continue;
				}
				if (df <= 0.0){
					System.out.println(tStat +",NA,NA");
					continue;
				}
				TDistribution tDistribution = new TDistribution(df);
				if(tStat > 0.0)
					System.out.println(tStat + "," + df + "," + (1-tDistribution.cumulativeProbability(tStat))*2.0);
				else
					System.out.println(tStat + "," + df + "," +  tDistribution.cumulativeProbability(tStat)*2.0);
			}
		}
	}

	public static class Evaluator<T> extends EvaRunnable<T> {
		T[][][] aliceCase;
		T[][][] bobCase;
		T[][][] aliceControl;
		T[][][] bobControl;

		T[][] numAliceCase;
		T[][] numBobCase;
		T[][] numAliceControl;
		T[][] numBobControl;
		int numOfTests;
		@Override
		public void prepareInput(CompEnv<T> gen) {
			StatisticsData caseSta = PrepareDataDifferentialAbundance.readFile(args[0]);
			StatisticsData controlSta = PrepareDataDifferentialAbundance.readFile(args[1]);			
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
			T[] l = flib.publicValue(0.0);

			boolean[][][] caseData = new boolean[caseSta.numOfTuples][3][l.length];
			boolean[][] numCaseData = new boolean[caseSta.numOfTuples][l.length];
			for(int i = 0; i < caseSta.numOfTuples; i++) {								
				caseData[i][0] = Utils.fromFloat(caseSta.data[i].totalSum, FWidth, FOffset);
				caseData[i][1] = Utils.fromFloat(caseSta.data[i].sumOfSquares, FWidth, FOffset);
				caseData[i][2] = Utils.fromFloat(caseSta.data[i].numOfSamples, FWidth, FOffset);
				numCaseData[i] = Utils.fromInt(caseSta.data[i].numOfPresent, Width);
			}

			boolean[][][] controlData = new boolean[controlSta.numOfTuples][3][l.length];
			boolean[][] numControlData = new boolean[controlSta.numOfTuples][l.length];

			for(int i = 0; i < controlSta.numOfTuples; i++) {
				controlData[i][0] = Utils.fromFloat(controlSta.data[i].totalSum, FWidth, FOffset);
				controlData[i][1] = Utils.fromFloat(controlSta.data[i].sumOfSquares, FWidth, FOffset);
				controlData[i][2] = Utils.fromFloat(controlSta.data[i].numOfSamples, FWidth, FOffset);
				numControlData[i] = Utils.fromInt(controlSta.data[i].numOfPresent, Width);
			}
			aliceCase = gen.inputOfAlice(caseData);
			aliceControl = gen.inputOfAlice(controlData);
			bobCase = gen.inputOfBob(caseData);
			bobControl = gen.inputOfBob(controlData);
			
			numAliceCase = gen.inputOfAlice(numCaseData);
			numAliceControl = gen.inputOfAlice(numControlData);
			numBobCase = gen.inputOfBob(numCaseData);
			numBobControl = gen.inputOfBob(numControlData);
			numOfTests = caseSta.numOfTuples;
		}
		T[][][] res;

		@Override
		public void secureCompute(CompEnv<T> gen) {
			res = compute(gen, aliceCase, bobCase, aliceControl, bobControl, numAliceCase, numBobCase,
					numAliceControl, numBobControl,numOfTests);
		}

		@Override
		public void prepareOutput(CompEnv<T> gen) {
			FixedPointLib<T> flib = new FixedPointLib<T>(gen, FWidth, FOffset);
			for(int i = 0; i < numOfTests; i++){
				gen.outputToAlice(res[i][0]);
				gen.outputToAlice(res[i][1]);
			}
		}
	}
}
