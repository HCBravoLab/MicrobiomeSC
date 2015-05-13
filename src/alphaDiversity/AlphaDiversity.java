package alphaDiversity;

import alphaDiversity.PrepareData;
import util.EvaRunnable;
import util.GenRunnable;
import util.Utils;
import circuits.arithmetic.FloatLib;
import flexsc.CompEnv;
import org.apache.commons.math3.distribution.*;

public class AlphaDiversity {
	static public int Width = 9;
	static public int FWidth = 49;
	static public int FOffset = 8;

	public static<T> T[][] compute(CompEnv<T> gen, T[][][] aliceCase,
			T[][][] bobCase,
			T[][][] aliceControl,
			T[][][] bobControl, int numOfTests) {

		T[][] res = gen.newTArray(numOfTests, 0);
		FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
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
		T[] tStat;

		caseSumOfSquares = flib.add(aliceCase[0][1], bobCase[0][1]);
		controlSumOfSquares = flib.add(aliceControl[0][1], bobControl[0][1]);

		caseTotalSum = flib.add(aliceCase[0][0], bobCase[0][0]);
		controlTotalSum = flib.add(aliceControl[0][0], bobControl[0][0]);

		caseNum = flib.add(aliceCase[0][2], bobCase[0][2]);
		controlNum = flib.add(aliceControl[0][2], bobControl[0][2]);

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

		T[] caseVarianceSquared = flib.multiply(caseVariance, caseVariance);
		T[] controlVarianceSquared = flib.multiply(controlVariance, controlVariance);

		T[] degreesOfFreedomTop = flib.add(tLowerFirst, tLowerSecond);
		degreesOfFreedomTop = flib.multiply(degreesOfFreedomTop, degreesOfFreedomTop);
		T[] degreesOfFreedomBottomFirst = flib.div(caseVarianceSquared, flib.multiply(caseNum,flib.multiply(caseNum,(flib.sub(caseNum, flib.publicValue(1.0))))));
		T[] degreesOfFreedomBottomSecond = flib.div(controlVarianceSquared, flib.multiply(controlNum,flib.multiply(controlNum,(flib.sub(controlNum, flib.publicValue(1.0))))));
		T[] degreesOfFreedom = flib.div(degreesOfFreedomTop, flib.add(degreesOfFreedomBottomFirst, degreesOfFreedomBottomSecond));

		res[0] = tStat;
		res[1] = degreesOfFreedom;
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
			Statistics caseSta = PrepareData.readFile(args[0]);
			Statistics controlSta = PrepareData.readFile(args[1]);
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
			T[] l = flib.publicValue(0.0);
			boolean[][][] caseData = new boolean[1][4][l.length];
			caseData[0][0] = Utils.fromFixPoint(caseSta.simpsonsIndexTotalSum, FWidth, FOffset);
			caseData[0][1] = Utils.fromFixPoint(caseSta.simpsonsSumOfSquares, FWidth, FOffset);
			caseData[0][2] = Utils.fromFixPoint(caseSta.numOfSimpsonsIndices, FWidth, FOffset);
			caseData[0][3] = Utils.fromFixPoint(caseSta.numOfSamples, FWidth, FOffset);

			boolean[][][] controlData = new boolean[1][4][l.length];
			controlData[0][0] = Utils.fromFixPoint(controlSta.simpsonsIndexTotalSum, FWidth, FOffset);
			controlData[0][1] = Utils.fromFixPoint(controlSta.simpsonsSumOfSquares, FWidth, FOffset);
			controlData[0][2] = Utils.fromFixPoint(controlSta.numOfSimpsonsIndices, FWidth, FOffset);
			controlData[0][3] = Utils.fromFixPoint(controlSta.numOfSamples, FWidth, FOffset);

			aliceCase = gen.inputOfAlice(caseData);
			aliceControl = gen.inputOfAlice(controlData);
			bobCase = gen.inputOfBob(caseData);
			bobControl = gen.inputOfBob(controlData);
			numOfTests = caseSta.numOfSimpsonsIndices;
		}

		T[][] res;
		@Override
		public void secureCompute(CompEnv<T> gen) {
			res = compute(gen, aliceCase, bobCase, aliceControl, bobControl, numOfTests);
		}

		@Override
		public void prepareOutput(CompEnv<T> gen) {
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
			double tStat = flib.outputToAlice(res[0]);
			double df = flib.outputToAlice(res[1]);
			System.out.print("t-stat = " + tStat + " ");
			System.out.print("degrees of freedom = " + df + " ");
			TDistribution tDistribution = new TDistribution(df);
			if(tStat > 0.0)
				System.out.println("p-value " + (1-tDistribution.cumulativeProbability(tStat))*2.0 + "\n");
			else
				System.out.println("p-value " + tDistribution.cumulativeProbability(tStat)*2.0 + "\n");
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
			Statistics caseSta = PrepareData.readFile(args[0]);
			Statistics controlSta = PrepareData.readFile(args[1]);			
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
			T[] l = flib.publicValue(0.0);
			boolean[][][] caseData = new boolean[1][4][l.length];
			caseData[0][0] = Utils.fromFixPoint(caseSta.simpsonsIndexTotalSum, FWidth, FOffset);
			caseData[0][1] = Utils.fromFixPoint(caseSta.simpsonsSumOfSquares, FWidth, FOffset);
			caseData[0][2] = Utils.fromFixPoint(caseSta.numOfSimpsonsIndices, FWidth, FOffset);
			caseData[0][3] = Utils.fromFixPoint(caseSta.numOfSamples, FWidth, FOffset);

			boolean[][][] controlData = new boolean[1][4][l.length];
			controlData[0][0] = Utils.fromFixPoint(controlSta.simpsonsIndexTotalSum, FWidth, FOffset);
			controlData[0][1] = Utils.fromFixPoint(controlSta.simpsonsSumOfSquares, FWidth, FOffset);
			controlData[0][2] = Utils.fromFixPoint(controlSta.numOfSimpsonsIndices, FWidth, FOffset);
			controlData[0][3] = Utils.fromFixPoint(controlSta.numOfSamples, FWidth, FOffset);

			aliceCase = gen.inputOfAlice(caseData);
			aliceControl = gen.inputOfAlice(controlData);
			bobCase = gen.inputOfBob(caseData);
			bobControl = gen.inputOfBob(controlData);
			numOfTests = caseSta.numOfSimpsonsIndices;
		}
		T[][] res;

		@Override
		public void secureCompute(CompEnv<T> gen) {
			res = compute(gen, aliceCase, bobCase, aliceControl, bobControl, numOfTests);
		}

		@Override
		public void prepareOutput(CompEnv<T> gen) {
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
			flib.outputToAlice(res[0]);
			flib.outputToAlice(res[1]);
		}
	}
}
