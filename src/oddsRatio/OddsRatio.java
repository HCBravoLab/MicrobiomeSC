package oddsRatio;

import oddsRatio.PrepareData;
import oddsRatio.PrepareData.StatisticsData;
import util.EvaRunnable;
import util.GenRunnable;
import util.Utils;
import circuits.arithmetic.FloatLib;
import circuits.arithmetic.IntegerLib;
import flexsc.CompEnv;

public class OddsRatio {
	static public int Width = 9;
	static public int FWidth = 50;
	static public int FOffset = 11;

	public static<T> T[][] compute(CompEnv<T> env, T[][][] aliceCase,
			T[][][] bobCase,
			T[][][] aliceControl,
			T[][][] bobControl, int numOfTests, boolean precise) {
		T[][] res = env.newTArray(numOfTests, 0);
		IntegerLib<T> lib = new IntegerLib<T>(env);
		FloatLib<T> flib = new FloatLib<T>(env, FWidth, FOffset);
		for(int i = 0; i < numOfTests; ++i) {
			T[] a = lib.add(aliceCase[i][0], bobCase[i][0]);
			T[] b = lib.add(aliceCase[i][1], bobCase[i][1]);
			T[] c = lib.add(aliceControl[i][0], bobControl[i][0]);
			T[] d = lib.add(aliceControl[i][1], bobControl[i][1]);

			T[] fa = lib.toSecureFloat(a, flib);
			T[] fb = lib.toSecureFloat(b, flib);
			T[] fc = lib.toSecureFloat(c, flib);
			T[] fd = lib.toSecureFloat(d, flib);

			res[i] = flib.div(flib.multiply(fa, fd), flib.multiply(fb, fc));
		}
		return res;
	}
	
	public static class Generator<T> extends GenRunnable<T> {
		T[][][] aliceCase;
		T[][][] bobCase;
		T[][][] aliceControl;
		T[][][] bobControl;

		int numOfTests;
		double extraFactor;
		boolean precise;
		@Override
		public void prepareInput(CompEnv<T> env) throws Exception {
			StatisticsData caseSta = PrepareData.readFile(args[0]);
			StatisticsData controlSta = PrepareData.readFile(args[1]);
			FloatLib<T> flib = new FloatLib<T>(env, FWidth, FOffset);
			T[] l = flib.publicValue(0.0);

			boolean[][][] caseData = new boolean[caseSta.numberOftuples][2][l.length];
			for(int i = 0; i < caseSta.numberOftuples; ++i) {
				caseData[i][0] = Utils.fromInt(caseSta.data[i].numOfPresent, Width);
				caseData[i][1] = Utils.fromInt(caseSta.data[i].totalNum - caseSta.data[i].numOfPresent, Width);
			}

			boolean[][][] controlData = new boolean[controlSta.numberOftuples][2][Width];
			for(int i = 0; i < controlSta.numberOftuples; ++i) {
				controlData[i][0] = Utils.fromInt(controlSta.data[i].numOfPresent, Width);
				controlData[i][1] = Utils.fromInt(controlSta.data[i].totalNum - controlSta.data[i].numOfPresent, Width);
			}
			aliceCase = env.inputOfAlice(caseData);
			aliceControl = env.inputOfAlice(controlData);
			bobCase = env.inputOfBob(caseData);
			bobControl = env.inputOfBob(controlData);
			numOfTests = caseSta.numberOftuples;
		}

		T[][] res;
		@Override
		public void secureCompute(CompEnv<T> env) {
			res = compute(env, aliceCase, bobCase, aliceControl, bobControl, numOfTests, precise);
		}

		@Override
		public void prepareOutput(CompEnv<T> env) {
			FloatLib<T> flib = new FloatLib<T>(env, FWidth, FOffset);
			for(int i = 0; i < res.length; ++i)
				System.out.println("odds ratio:  " + flib.outputToAlice(res[i]));
		}

	}

	public static class Evaluator<T> extends EvaRunnable<T> {
		T[][][] aliceCase;
		T[][][] bobCase;
		T[][][] aliceControl;
		T[][][] bobControl;
		int numOfTests;
		boolean precise;
		@Override
		public void prepareInput(CompEnv<T> env) throws Exception {
			StatisticsData caseSta = PrepareData.readFile(args[0]);
			StatisticsData controlSta = PrepareData.readFile(args[1]);
			FloatLib<T> flib = new FloatLib<T>(env, FWidth, FOffset);
			T[] l = flib.publicValue(0.0);

			boolean[][][] caseData = new boolean[caseSta.numberOftuples][2][l.length];
			for(int i = 0; i < caseSta.numberOftuples; ++i) {
				caseData[i][0] = Utils.fromInt(caseSta.data[i].numOfPresent, Width);
				caseData[i][1] = Utils.fromInt(caseSta.data[i].totalNum - caseSta.data[i].numOfPresent, Width);
			}

			boolean[][][] controlData = new boolean[controlSta.numberOftuples][2][Width];
			for(int i = 0; i < controlSta.numberOftuples; ++i) {
				controlData[i][0] = Utils.fromInt(controlSta.data[i].numOfPresent, Width);
				controlData[i][1] = Utils.fromInt(controlSta.data[i].totalNum - controlSta.data[i].numOfPresent, Width);
			}
			aliceCase = env.inputOfAlice(caseData);
			aliceControl = env.inputOfAlice(controlData);
			bobCase = env.inputOfBob(caseData);
			bobControl = env.inputOfBob(controlData);
			numOfTests = caseSta.numberOftuples;
		}
		
		T[][] res;

		@Override
		public void secureCompute(CompEnv<T> env) {
			res = compute(env, aliceCase, bobCase, aliceControl, bobControl, numOfTests, precise);
		}

		@Override
		public void prepareOutput(CompEnv<T> env) {
			FloatLib<T> flib = new FloatLib<T>(env, FWidth, FOffset);
			for(int i = 0; i < numOfTests; ++i)
				flib.outputToAlice(res[i]);
		}
	}
}
