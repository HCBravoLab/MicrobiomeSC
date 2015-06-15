package perFeatureDifferentialAbundance;


import util.EvaRunnable;
import util.GenRunnable;
import util.Utils;
import circuits.CircuitLib;
import circuits.arithmetic.FixedPointLib;
import circuits.arithmetic.FloatLib;
import circuits.arithmetic.IntegerLib;
import perFeatureDifferentialAbundance.PrepareData;
import perFeatureDifferentialAbundance.PrepareData.StatisticsData;
import flexsc.CompEnv;

import org.apache.commons.math3.distribution.TDistribution;

public class DifferentialAbundance {
	static public int Width = 32;
	static public int FWidth = 54;
	static public int FOffset = 11;
	static int filterThreshold = 5;


	static public<T> T[][] dummyVariable(CompEnv<T> gen){
		FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
		T[] zero = flib.publicValue(0.0);
		T[][] res = gen.newTArray(2, 0);
		res[0] = zero;
		res[1] = zero;
		return res;
	}
	
	static public<T> T filter(CompEnv<T> gen, T[] numAliceCase, T[] aliceCaseTotal,
			T[] numBobCase, T[] bobCaseTotal, T[] numAliceControl, T[] aliceControlTotal,
			T[] numBobControl, T[] bobControlTotal){

		IntegerLib<T> ilib = new IntegerLib<T>(gen, 32);
		T[] threshold = ilib.publicValue(filterThreshold);
		T[] caseNum = ilib.add(numAliceCase, numBobCase);
		T[] controlNum = ilib.add(numAliceControl, numBobControl);
		T[] totalPresent = ilib.add(caseNum, controlNum);
		T aboveThreshold = ilib.geq(totalPresent, threshold);
		return aboveThreshold;
	}
	
	public static<T> T[][] compute(CompEnv<T> gen, T[] aliceCaseTotalSum, T[] aliceCaseSumOfSquares, 
			T[] aliceCaseNumSamples,
			T[] bobCaseTotalSum, T[] bobCaseSumOfSquares, T[] bobCaseNumSamples, 
			T[] aliceControlTotalSum, T[] aliceControlSumOfSquares, T[] aliceControlNumSamples,
			T[] bobControlTotalSum, T[] bobControlSumOfSquares, T[] bobControlNumSamples, T[] numAliceCase,
			T[] numBobCase, T[] numAliceControl, 
			T[] numBobControl) {	

		FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
		IntegerLib<T> ilib = new IntegerLib<T>(gen, 32);
		CircuitLib<T> cl = new CircuitLib<T>(gen);

		T[] zero = flib.publicValue(0.0);
		T[] one = flib.publicValue(1.0);
		T[] negOne = flib.publicValue(-1.0);
		
		T[] caseSumOfSquares;
		T[] controlSumOfSquares;

		T[] caseSum;
		T[] controlSum;

		T[] caseNum;
		T[] controlNum;

		caseSumOfSquares = flib.add(aliceCaseSumOfSquares,bobCaseSumOfSquares);
		controlSumOfSquares = flib.add(aliceControlSumOfSquares, bobControlSumOfSquares);

		caseSum = flib.add(aliceCaseTotalSum, bobCaseTotalSum);
		controlSum = flib.add(aliceControlTotalSum, bobControlTotalSum);

		T[] caseNumPresent = ilib.add(numAliceCase, numBobCase);
		T[] controlNumPresent= ilib.add(numAliceControl,numBobControl);

		caseNum = flib.add(aliceCaseNumSamples, bobCaseNumSamples);
		controlNum = flib.add(aliceControlNumSamples, bobControlNumSamples);
		
		T[][] res = gen.newTArray(2, 0);

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
		T tLowerFirstIsZero = flib.eq(caseNumPresent, zero);
		T[] tLowerFirst = flib.div(caseVariance, caseNum);
		tLowerFirst = cl.mux(tLowerFirst, zero, tLowerFirstIsZero);
		
		T tLowerSecondIsZero = flib.eq(controlNumPresent, zero);
		
		T[] tLowerSecond = flib.div(controlVariance, controlNum);
		tLowerSecond = cl.mux(tLowerSecond, zero, tLowerSecondIsZero);
		
		T[] tLowerSqrt = flib.sqrt(flib.add(flib.add(zero,tLowerFirst), flib.add(zero,tLowerSecond)));
		T[] tStatNonZero = flib.div(tUpper, tLowerSqrt);
		
		T[] degreesOfFreedomTop = flib.add(tLowerFirst, tLowerSecond);
		degreesOfFreedomTop = flib.multiply(flib.add(zero,degreesOfFreedomTop), flib.add(zero,degreesOfFreedomTop));

		T[] degreesOfFreedomBottomFirst = flib.div(caseVariance, flib.sub(caseNum, one));
		degreesOfFreedomBottomFirst = flib.multiply(flib.add(zero,degreesOfFreedomBottomFirst), flib.add(zero,degreesOfFreedomBottomFirst));
		degreesOfFreedomBottomFirst = flib.div(degreesOfFreedomBottomFirst, flib.sub(caseNum, one));

		T[] degreesOfFreedomBottomSecond = flib.div(controlVariance, flib.sub(caseNum, one));
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

		res[0] = tStat;
		res[1] = df;
		return res;
	}

	public static class Generator<T> extends GenRunnable<T> {
		T[] aliceCaseTotalSum;
		T[] aliceCaseSumOfSquares;
		T[] aliceCaseNumOfSamples;
		T[] aliceCaseNumPresent;

		T[] bobCaseTotalSum;
		T[] bobCaseSumOfSquares;
		T[] bobCaseNumOfSamples;
		T[] bobCaseNumPresent;

		T[] aliceControlTotalSum;
		T[] aliceControlSumOfSquares;
		T[] aliceControlNumOfSamples;
		T[] aliceControlNumPresent;

		T[] bobControlTotalSum;
		T[] bobControlSumOfSquares;
		T[] bobControlNumOfSamples;
		T[] bobControlNumPresent;

		double[] caseTotalSum;
		double[] caseSumOfSquares;
		int[] caseNumOfSamples;
		int[] caseNumPresent;

		double[] controlTotalSum;
		double[] controlSumOfSquares;
		int[] controlNumOfSamples;
		int[] controlNumPresent;

		int numOfTests;
		T[][][] res;
		T[] l;

		@Override
		public void prepareInput(CompEnv<T> gen) {
			StatisticsData caseSta = PrepareData.readFile(args[0]);
			StatisticsData controlSta = PrepareData.readFile(args[1]);
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
			l = flib.publicValue(0.0);
			caseTotalSum = new double[caseSta.numOfTuples];
			caseSumOfSquares = new double[caseSta.numOfTuples];
			caseNumOfSamples = new int[caseSta.numOfTuples];
			caseNumPresent = new int[caseSta.numOfTuples];

			for(int i = 0; i < caseSta.numOfTuples; i++) {
				caseTotalSum[i] = caseSta.data[i].totalSum;
				caseSumOfSquares[i] = caseSta.data[i].sumOfSquares;
				caseNumOfSamples[i] = caseSta.data[i].numOfSamples;
				caseNumPresent[i] = caseSta.data[i].numOfPresent;
			}

			controlTotalSum = new double[controlSta.numOfTuples];
			controlSumOfSquares = new double[controlSta.numOfTuples];
			controlNumOfSamples = new int[controlSta.numOfTuples];
			controlNumPresent = new int[controlSta.numOfTuples];

			for(int i = 0; i < controlSta.numOfTuples; i++) {
				controlTotalSum[i] = controlSta.data[i].totalSum;
				controlSumOfSquares[i] = controlSta.data[i].sumOfSquares;
				controlNumOfSamples[i] = controlSta.data[i].numOfSamples;
				controlNumPresent[i] = controlSta.data[i].numOfPresent;
			}
			numOfTests = caseSta.numOfTuples;
			res = gen.newTArray(numOfTests, 2, 0);
		}

		@Override
		public void secureCompute(CompEnv<T> gen) {
			T result;
			
			boolean[] filteredFeatures = new boolean[numOfTests];
			for(int i =0; i < numOfTests; i++){				
				aliceCaseNumPresent =gen.newTArray(Width);
				aliceCaseNumOfSamples = gen.newTArray(Width);
				bobCaseNumPresent =gen.newTArray(Width);
				bobCaseNumOfSamples = gen.newTArray(Width);	
				aliceControlNumPresent =gen.newTArray(Width);
				aliceControlNumOfSamples = gen.newTArray(Width);
				bobControlNumPresent =gen.newTArray(Width);
				bobControlNumOfSamples = gen.newTArray(Width);
				
				aliceCaseNumPresent = gen.inputOfAlice(Utils.fromInt(caseNumPresent[i], Width));
				aliceCaseNumOfSamples = gen.inputOfAlice(Utils.fromInt(caseNumOfSamples[i], Width));
				bobCaseNumPresent = gen.inputOfBob(new boolean[Width]);
				bobCaseNumOfSamples = gen.inputOfBob(new boolean[Width]);
				aliceControlNumPresent = gen.inputOfAlice(Utils.fromInt(controlNumPresent[i], Width));
				aliceControlNumOfSamples = gen.inputOfAlice(Utils.fromInt(controlNumOfSamples[i], Width));
				bobControlNumPresent = gen.inputOfBob(new boolean[Width]);
				bobControlNumOfSamples = gen.inputOfBob(new boolean[Width]);
				
				result = filter(gen, aliceCaseNumPresent, aliceCaseNumOfSamples, 
					bobCaseNumPresent, bobCaseNumOfSamples,
					aliceControlNumPresent, aliceControlNumOfSamples,
					bobControlNumPresent, bobControlNumOfSamples);
				filteredFeatures[i] = gen.outputToAlice(result);
			}
				
			for(int i = 0; i < numOfTests; i++){
				if(!(filteredFeatures[i])){
					res[i] = dummyVariable(gen);
					continue;
				}
				aliceCaseTotalSum = gen.newTArray(Width);
				aliceCaseSumOfSquares = gen.newTArray(Width);
				aliceCaseNumOfSamples = gen.newTArray(Width);
				bobCaseTotalSum = gen.newTArray(Width);
				bobCaseSumOfSquares = gen.newTArray(Width);
				bobCaseNumOfSamples = gen.newTArray(Width);
				aliceControlTotalSum = gen.newTArray(Width);
				aliceControlSumOfSquares = gen.newTArray(Width);
				aliceControlNumOfSamples = gen.newTArray(Width);
				bobControlTotalSum = gen.newTArray(Width);
				bobControlSumOfSquares = gen.newTArray(Width);
				bobControlNumOfSamples = gen.newTArray(Width);
				
				aliceCaseTotalSum = gen.inputOfAlice(Utils.fromFloat(caseTotalSum[i], FWidth, FOffset));
				aliceCaseSumOfSquares = gen.inputOfAlice(Utils.fromFloat(caseSumOfSquares[i], FWidth, FOffset));
				aliceCaseNumOfSamples = gen.inputOfAlice(Utils.fromFloat(caseNumOfSamples[i], FWidth, FOffset));
				bobCaseTotalSum = gen.inputOfBob(new boolean[l.length]);
				bobCaseSumOfSquares = gen.inputOfBob(new boolean[l.length]);
				bobCaseNumOfSamples = gen.inputOfBob(new boolean[l.length]);

				aliceControlTotalSum = gen.inputOfAlice(Utils.fromFloat(controlTotalSum[i], FWidth, FOffset));
				aliceControlSumOfSquares = gen.inputOfAlice(Utils.fromFloat(controlSumOfSquares[i], FWidth, FOffset));
				aliceControlNumOfSamples = gen.inputOfAlice(Utils.fromFloat(controlNumOfSamples[i], FWidth, FOffset));
				bobControlTotalSum = gen.inputOfBob(new boolean[l.length]);
				bobControlSumOfSquares = gen.inputOfBob(new boolean[l.length]);
				bobControlNumOfSamples = gen.inputOfBob(new boolean[l.length]);
				
				aliceCaseNumPresent = gen.inputOfAlice(Utils.fromInt(caseNumPresent[i], Width));
				bobCaseNumPresent = gen.inputOfBob(new boolean[Width]);
				aliceControlNumPresent = gen.inputOfAlice(Utils.fromInt(controlNumPresent[i], Width));
				bobControlNumPresent = gen.inputOfBob(new boolean[Width]);
				
				res[i] = compute(gen, aliceCaseTotalSum, aliceCaseSumOfSquares, aliceCaseNumOfSamples,
						bobCaseTotalSum, bobCaseSumOfSquares, bobCaseNumOfSamples,
						aliceControlTotalSum, aliceControlSumOfSquares, aliceControlNumOfSamples,
						bobControlTotalSum, bobControlSumOfSquares, bobControlNumOfSamples,
						aliceCaseNumPresent, bobCaseNumPresent, aliceControlNumPresent, bobControlNumPresent);
				}		}

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
		T[] aliceCaseTotalSum;
		T[] aliceCaseSumOfSquares;
		T[] aliceCaseNumOfSamples;
		T[] aliceCaseNumPresent;

		T[] bobCaseTotalSum;
		T[] bobCaseSumOfSquares;
		T[] bobCaseNumOfSamples;
		T[] bobCaseNumPresent;

		T[] aliceControlTotalSum;
		T[] aliceControlSumOfSquares;
		T[] aliceControlNumOfSamples;
		T[] aliceControlNumPresent;

		T[] bobControlTotalSum;
		T[] bobControlSumOfSquares;
		T[] bobControlNumOfSamples;
		T[] bobControlNumPresent;

		double[] caseTotalSum;
		double[] caseSumOfSquares;
		int[] caseNumOfSamples;
		int[] caseNumPresent;

		double[] controlTotalSum;
		double[] controlSumOfSquares;
		int[] controlNumOfSamples;
		int[] controlNumPresent;

		int numOfTests;
		T[][][] res;
		T[] l;
		@Override
		public void prepareInput(CompEnv<T> gen) {
			StatisticsData caseSta = PrepareData.readFile(args[0]);
			StatisticsData controlSta = PrepareData.readFile(args[1]);			
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
			l = flib.publicValue(0.0);

			caseTotalSum = new double[caseSta.numOfTuples];
			caseSumOfSquares = new double[caseSta.numOfTuples];
			caseNumOfSamples = new int[caseSta.numOfTuples];
			caseNumPresent = new int[caseSta.numOfTuples];
			for(int i = 0; i < caseSta.numOfTuples; i++) {
				caseTotalSum[i] = caseSta.data[i].totalSum;
				caseSumOfSquares[i] = caseSta.data[i].sumOfSquares;
				caseNumOfSamples[i] = caseSta.data[i].numOfSamples;
				caseNumPresent[i] = caseSta.data[i].numOfPresent;
			}

			controlTotalSum = new double[controlSta.numOfTuples];
			controlSumOfSquares = new double[controlSta.numOfTuples];
			controlNumOfSamples = new int[controlSta.numOfTuples];
			controlNumPresent = new int[controlSta.numOfTuples];
			for(int i = 0; i < controlSta.numOfTuples; i++) {
				controlTotalSum[i] = controlSta.data[i].totalSum;
				controlSumOfSquares[i] = controlSta.data[i].sumOfSquares;
				controlNumOfSamples[i] = controlSta.data[i].numOfSamples;
				controlNumPresent[i] = controlSta.data[i].numOfPresent;
			}
			numOfTests = caseSta.numOfTuples;
			res = gen.newTArray(numOfTests, 2, 0);
		}

		@Override
		public void secureCompute(CompEnv<T> gen) {			
		T result;
		
		boolean[] filteredFeatures = new boolean[numOfTests];
		for(int i =0; i < numOfTests; i++){				
			aliceCaseNumPresent =gen.newTArray(Width);
			aliceCaseNumOfSamples = gen.newTArray(Width);
			bobCaseNumPresent =gen.newTArray(Width);
			bobCaseNumOfSamples = gen.newTArray(Width);	
			aliceControlNumPresent =gen.newTArray(Width);
			aliceControlNumOfSamples = gen.newTArray(Width);
			bobControlNumPresent =gen.newTArray(Width);
			bobControlNumOfSamples = gen.newTArray(Width);
			
			aliceCaseNumPresent = gen.inputOfAlice(Utils.fromInt(caseNumPresent[i], Width));
			aliceCaseNumOfSamples = gen.inputOfAlice(Utils.fromInt(caseNumOfSamples[i], Width));
			bobCaseNumPresent = gen.inputOfBob(new boolean[Width]);
			bobCaseNumOfSamples = gen.inputOfBob(new boolean[Width]);
			aliceControlNumPresent = gen.inputOfAlice(Utils.fromInt(controlNumPresent[i], Width));
			aliceControlNumOfSamples = gen.inputOfAlice(Utils.fromInt(controlNumOfSamples[i], Width));
			bobControlNumPresent = gen.inputOfBob(new boolean[Width]);
			bobControlNumOfSamples = gen.inputOfBob(new boolean[Width]);
			
			result = filter(gen, aliceCaseNumPresent, aliceCaseNumOfSamples, 
				bobCaseNumPresent, bobCaseNumOfSamples,
				aliceControlNumPresent, aliceControlNumOfSamples,
				bobControlNumPresent, bobControlNumOfSamples);
			filteredFeatures[i] = gen.outputToAlice(result);
		}
			
		for(int i = 0; i < numOfTests; i++){
			if(!(filteredFeatures[i])){
				res[i] = dummyVariable(gen);
				continue;
			}
			aliceCaseTotalSum = gen.newTArray(Width);
			aliceCaseSumOfSquares = gen.newTArray(Width);
			aliceCaseNumOfSamples = gen.newTArray(Width);
			bobCaseTotalSum = gen.newTArray(Width);
			bobCaseSumOfSquares = gen.newTArray(Width);
			bobCaseNumOfSamples = gen.newTArray(Width);
			aliceControlTotalSum = gen.newTArray(Width);
			aliceControlSumOfSquares = gen.newTArray(Width);
			aliceControlNumOfSamples = gen.newTArray(Width);
			bobControlTotalSum = gen.newTArray(Width);
			bobControlSumOfSquares = gen.newTArray(Width);
			bobControlNumOfSamples = gen.newTArray(Width);
			
			aliceCaseTotalSum = gen.inputOfAlice(Utils.fromFloat(caseTotalSum[i], FWidth, FOffset));
			aliceCaseSumOfSquares = gen.inputOfAlice(Utils.fromFloat(caseSumOfSquares[i], FWidth, FOffset));
			aliceCaseNumOfSamples = gen.inputOfAlice(Utils.fromFloat(caseNumOfSamples[i], FWidth, FOffset));
			aliceCaseNumPresent = gen.inputOfAlice(Utils.fromInt(caseNumPresent[i], Width));

			bobCaseTotalSum = gen.inputOfBob(new boolean[l.length]);
			bobCaseSumOfSquares = gen.inputOfBob(new boolean[l.length]);
			bobCaseNumOfSamples = gen.inputOfBob(new boolean[l.length]);
			bobCaseNumPresent = gen.inputOfBob(new boolean[Width]);

			aliceControlTotalSum = gen.inputOfAlice(Utils.fromFloat(controlTotalSum[i], FWidth, FOffset));
			aliceControlSumOfSquares = gen.inputOfAlice(Utils.fromFloat(controlSumOfSquares[i], FWidth, FOffset));
			aliceControlNumOfSamples = gen.inputOfAlice(Utils.fromFloat(controlNumOfSamples[i], FWidth, FOffset));
			aliceControlNumPresent = gen.inputOfAlice(Utils.fromInt(controlNumPresent[i], Width));

			bobControlTotalSum = gen.inputOfBob(new boolean[l.length]);
			bobControlSumOfSquares = gen.inputOfBob(new boolean[l.length]);
			bobControlNumOfSamples = gen.inputOfBob(new boolean[l.length]);
			bobControlNumPresent = gen.inputOfBob(new boolean[Width]);
			
			res[i] = compute(gen, aliceCaseTotalSum, aliceCaseSumOfSquares, aliceCaseNumOfSamples,
					bobCaseTotalSum, bobCaseSumOfSquares, bobCaseNumOfSamples,
					aliceControlTotalSum, aliceControlSumOfSquares, aliceControlNumOfSamples,
					bobControlTotalSum, bobControlSumOfSquares, bobControlNumOfSamples, 
					aliceCaseNumPresent, bobCaseNumPresent, aliceControlNumPresent, bobControlNumPresent);
			}
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
