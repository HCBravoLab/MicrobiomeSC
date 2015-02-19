package differentialAbundanceBatch;


import org.apache.commons.math3.distribution.TDistribution;

import util.EvaRunnable;
import util.GenRunnable;
import util.Utils;
import circuits.arithmetic.FixedPointLib;
import circuits.arithmetic.FloatLib;
import circuits.arithmetic.IntegerLib;
import differentialAbundanceBatch.PrepareData;
import differentialAbundanceBatch.PrepareData.StatisticsData;
import flexsc.CompEnv;

public class DifferentialAbundance {
	static public int Width = 9;
	static public int FWidth = 48;
	static public int FOffset = 14;

	public static<T> T[][][] compute(CompEnv<T> gen, T[][][] aliceCase,
			T[][][] bobCase,
			T[][][] aliceControl,
			T[][][] bobControl, int numOfTests) {	
		
		T[][][] res = gen.newTArray(numOfTests, 2, 0);
		FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
		T[] tStat;
		
		for(int i = 0; i < numOfTests; i++) {
			T[] caseSumOfSquares;
			T[] controlSumOfSquares;
	
			T[] caseTotalSum;
			T[] controlTotalSum;
	
			T[] caseNum;
			T[] controlNum;
			
			T[] caseVariance;
			T[] controlVariance;
			T[] caseVarianceSecondTerm;
			T[] controlVarianceSecondTerm;
			T[] caseMeanAbundance;
			T[] controlMeanAbundance;
			T[] tUpper;
			T[] tLowerFirst;
			T[] tLowerSecond;
			T[] tLowerSqrt;
		    
		    caseSumOfSquares = flib.add(aliceCase[i][1], bobCase[i][1]);
		    controlSumOfSquares = flib.add(aliceControl[i][1], bobControl[i][1]);
		    
		    caseTotalSum = flib.add(aliceCase[i][0], bobCase[i][0]);
		    controlTotalSum = flib.add(aliceControl[i][0], bobControl[i][0]);
		    
		    caseNum = flib.add(aliceCase[i][2], bobCase[i][2]);
		    controlNum = flib.add(aliceControl[i][2], bobControl[i][2]);
		    
		    caseMeanAbundance = flib.div(caseTotalSum, caseNum);
		    caseVarianceSecondTerm = flib.div(flib.multiply(caseTotalSum, caseTotalSum), caseNum);
		    caseVariance = flib.div(flib.sub(caseSumOfSquares, caseVarianceSecondTerm), caseNum);
		    
		    controlMeanAbundance = flib.div(controlTotalSum, controlNum);
		    controlVarianceSecondTerm = flib.div(flib.multiply(controlTotalSum, controlTotalSum), controlNum);
		    controlVariance = flib.div(flib.sub(controlSumOfSquares, controlVarianceSecondTerm), controlNum);
		 
		    tUpper = flib.sub(controlMeanAbundance, caseMeanAbundance);
		    tLowerFirst = flib.div(caseVariance, caseNum);
		    tLowerSecond = flib.div(controlVariance, controlNum);
		    tLowerSqrt = flib.sqrt(flib.add(tLowerFirst, tLowerSecond));
		    tStat = flib.div(tUpper, tLowerSqrt);

		    T[] degreesOfFreedomTop;
		    T[] degreesOfFreedomBottomFirst;
		    T[] degreesOfFreedomBottomSecond;
		    T[] caseVarianceSquared = flib.multiply(caseVariance, caseVariance);
		    T[] controlVarianceSquared = flib.multiply(controlVariance, controlVariance);
		    degreesOfFreedomTop = flib.add(flib.div(caseVariance, caseNum), flib.div(controlVariance, controlNum));
		    degreesOfFreedomTop = flib.multiply(degreesOfFreedomTop, degreesOfFreedomTop);
		    degreesOfFreedomBottomFirst = flib.div(caseVarianceSquared, flib.multiply(caseNum,flib.multiply(caseNum,(flib.sub(caseNum, flib.publicValue(1.0))))));
		    degreesOfFreedomBottomSecond = flib.div(controlVarianceSquared, flib.multiply(controlNum,flib.multiply(controlNum,(flib.sub(controlNum, flib.publicValue(1.0))))));
		    T[] degreesOfFreedom = flib.div(degreesOfFreedomTop, flib.add(degreesOfFreedomBottomFirst, degreesOfFreedomBottomSecond));
		    res[i][0] = tStat;
		    res[i][1] = degreesOfFreedom;
		}
		return res;
	}
	
	public static class Generator<T> extends GenRunnable<T> {
		T[][][] aliceCase;
		T[][][] bobCase;
		T[][][] aliceControl;
		T[][][] bobControl;

		int numOfTests;
		@Override
		public void prepareInput(CompEnv<T> gen) {
			StatisticsData caseSta = PrepareData.readFile(args[0]);
			StatisticsData controlSta = PrepareData.readFile(args[1]);
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
			T[] l = flib.publicValue(0.0);

			boolean[][][] caseData = new boolean[caseSta.numOfTuples][3][l.length];
			for(int i = 0; i < caseSta.numOfTuples; i++) {
				caseData[i][0] = Utils.fromFloat(caseSta.data[i].totalSum, FWidth, FOffset);
				caseData[i][1] = Utils.fromFloat(caseSta.data[i].sumOfSquares, FWidth, FOffset);
				caseData[i][2] = Utils.fromFloat(caseSta.data[i].numOfSamples, FWidth, FOffset);
			}
			
			boolean[][][] controlData = new boolean[controlSta.numOfTuples][3][l.length];
			for(int i = 0; i < controlSta.numOfTuples; i++) {
				controlData[i][0] = Utils.fromFloat(controlSta.data[i].totalSum, FWidth, FOffset);
				controlData[i][1] = Utils.fromFloat(controlSta.data[i].sumOfSquares, FWidth, FOffset);
				controlData[i][2] = Utils.fromFloat(controlSta.data[i].numOfSamples, FWidth, FOffset);
			}
			aliceCase = gen.inputOfAlice(caseData);
			aliceControl = gen.inputOfAlice(controlData);
			bobCase = gen.inputOfBob(caseData);
			bobControl = gen.inputOfBob(controlData);
			numOfTests = caseSta.numOfTuples;
		}

		T[][][] res;
		@Override
		public void secureCompute(CompEnv<T> gen) {
			res = compute(gen, aliceCase, bobCase, aliceControl, bobControl, numOfTests);
		}

		@Override
		public void prepareOutput(CompEnv<T> gen) {
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
			for(int i = 0; i < numOfTests; i++){
				double tStat = flib.outputToAlice(res[i][0]);
				double df = flib.outputToAlice(res[i][1]);
				System.out.print("t-stat = " + tStat + " ");
				System.out.print("degrees of freedom = " + df + " ");
				if (df < 1.0){
					System.out.print("\n");
					continue;
				}
				TDistribution tDistribution = new TDistribution(df);
				if(tStat > 0.0)
					System.out.println("p-value " + (1-tDistribution.cumulativeProbability(tStat))*2.0 + "\n");
				else
					System.out.println("p-value " + tDistribution.cumulativeProbability(tStat)*2.0 + "\n");
			}
		}
	}

	public static class Evaluator<T> extends EvaRunnable<T> {
		T[][][] aliceCase;
		T[][][] bobCase;
		T[][][] aliceControl;
		T[][][] bobControl;
		int numOfTests;
		@Override
		public void prepareInput(CompEnv<T> gen) {
			StatisticsData caseSta = PrepareData.readFile(args[0]);
			StatisticsData controlSta = PrepareData.readFile(args[1]);			
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
			T[] l = flib.publicValue(0.0);
			
			boolean[][][] caseData = new boolean[caseSta.numOfTuples][3][l.length];
			for(int i = 0; i < caseSta.numOfTuples; i++) {								
				caseData[i][0] = Utils.fromFloat(caseSta.data[i].totalSum, FWidth, FOffset);
				caseData[i][1] = Utils.fromFloat(caseSta.data[i].sumOfSquares, FWidth, FOffset);
				caseData[i][2] = Utils.fromFloat(caseSta.data[i].numOfSamples, FWidth, FOffset);
			}

			boolean[][][] controlData = new boolean[controlSta.numOfTuples][3][l.length];
			for(int i = 0; i < controlSta.numOfTuples; i++) {
				controlData[i][0] = Utils.fromFloat(controlSta.data[i].totalSum, FWidth, FOffset);
				controlData[i][1] = Utils.fromFloat(controlSta.data[i].sumOfSquares, FWidth, FOffset);
				controlData[i][2] = Utils.fromFloat(controlSta.data[i].numOfSamples, FWidth, FOffset);
			}
			
			aliceCase = gen.inputOfAlice(caseData);
			aliceControl = gen.inputOfAlice(controlData);
			bobCase = gen.inputOfBob(caseData);
			bobControl = gen.inputOfBob(controlData);
			numOfTests = caseSta.numOfTuples;
		}
		T[][][] res;

		@Override
		public void secureCompute(CompEnv<T> gen) {
			res = compute(gen, aliceCase, bobCase, aliceControl, bobControl, numOfTests);
		}

		@Override
		public void prepareOutput(CompEnv<T> gen) {
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
			for(int i = 0; i < numOfTests; i++){
				flib.outputToAlice(res[i][0]);
				flib.outputToAlice(res[i][1]);
			}
		}
	}
}
