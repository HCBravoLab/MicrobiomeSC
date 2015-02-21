package chiSquare;

import chiSquare.PrepareData;
import chiSquare.PrepareData.StatisticsData;


import util.EvaRunnable;
import util.GenRunnable;
import util.Utils;
import circuits.arithmetic.FixedPointLib;
import circuits.arithmetic.FloatLib;
import circuits.arithmetic.IntegerLib;
import flexsc.CompEnv;

public class ChiSquare {
	static public int Width = 9;
	static public int FWidth = 50;
	static public int FOffset = 11;

	public static<T> T[][] compute(CompEnv<T> gen, T[][][] aliceCase,
			T[][][] bobCase,
			T[][][] aliceControl,
			T[][][] bobControl, int numOfTests) {
		
		T[][] res = gen.newTArray(numOfTests, 0);

		IntegerLib<T> lib = new IntegerLib<T>(gen);
		FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
		
		for(int i = 0; i < numOfTests; i++){			
			T[] a = lib.add(aliceCase[i][0], bobCase[i][0]);
			T[] b = lib.add(aliceCase[i][1], bobCase[i][1]);
			T[] c = lib.add(aliceControl[i][0], bobControl[i][0]);
			T[] d = lib.add(aliceControl[i][1], bobControl[i][1]);
			T[] g = lib.add(a, c);
			T[] k = lib.add(b, d);
	
			T[] fa = lib.toSecureFloat(a, flib);
			T[] fb = lib.toSecureFloat(b, flib);
			T[] fc = lib.toSecureFloat(c, flib);
			T[] fd = lib.toSecureFloat(d, flib);
			T[] fg = lib.toSecureFloat(g, flib);
			T[] fk = lib.toSecureFloat(k, flib);
	
			T[] tmp = flib.sub(flib.multiply(fa, fd), flib.multiply(fb, fc));
			tmp = flib.multiply(tmp, tmp);
			res[i] = flib.div(tmp, flib.multiply(fg, fk));
			return res;
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

			boolean[][][] caseData = new boolean[caseSta.numberOftuples][2][l.length];
			
			for(int i = 0; i < caseSta.numberOftuples; ++i) {
				System.out.println(caseSta.data[i].numOfPresent);
				caseData[i][0] = Utils.fromInt(caseSta.data[i].numOfPresent, Width);
				System.out.println(caseSta.data[i].totalNum - caseSta.data[i].numOfPresent);
				caseData[i][1] = Utils.fromInt(caseSta.data[i].totalNum - caseSta.data[i].numOfPresent, Width);
			}

			boolean[][][] controlData = new boolean[controlSta.numberOftuples][2][Width];
			for(int i = 0; i < controlSta.numberOftuples; ++i) {
				caseData[i][0] = Utils.fromInt(controlSta.data[i].numOfPresent, Width);
				controlData[i][0] = Utils.fromInt(controlSta.data[i].numOfPresent, Width);
				System.out.println(controlSta.data[i].totalNum - controlSta.data[i].numOfPresent);
				controlData[i][1] = Utils.fromInt(controlSta.data[i].totalNum - controlSta.data[i].numOfPresent, Width);
			}
			aliceCase = gen.inputOfAlice(caseData);
			aliceControl = gen.inputOfAlice(controlData);
			bobCase = gen.inputOfBob(caseData);
			bobControl = gen.inputOfBob(controlData);
			numOfTests = caseSta.numberOftuples;
		}

		T[][] res;
		@Override
		public void secureCompute(CompEnv<T> gen) {
			res = compute(gen, aliceCase, bobCase, aliceControl, bobControl, numOfTests);
		}

		@Override
		public void prepareOutput(CompEnv<T> gen) {
			FixedPointLib<T> flib = new FixedPointLib<T>(gen, FWidth, FOffset);
			for(int i = 0; i < numOfTests; ++i)
				System.out.println(flib.outputToAlice(res[i]));
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

			boolean[][][] caseData = new boolean[caseSta.numberOftuples][2][l.length];
			
			for(int i = 0; i < caseSta.numberOftuples; ++i) {
				System.out.println(caseSta.data[i].numOfPresent);
				caseData[i][0] = Utils.fromInt(caseSta.data[i].numOfPresent, Width);
				System.out.println(caseSta.data[i].totalNum - caseSta.data[i].numOfPresent);
				caseData[i][1] = Utils.fromInt(caseSta.data[i].totalNum - caseSta.data[i].numOfPresent, Width);
			}

			boolean[][][] controlData = new boolean[controlSta.numberOftuples][2][Width];
			for(int i = 0; i < controlSta.numberOftuples; ++i) {
				caseData[i][0] = Utils.fromInt(controlSta.data[i].numOfPresent, Width);
				controlData[i][0] = Utils.fromInt(controlSta.data[i].numOfPresent, Width);
				System.out.println(controlSta.data[i].totalNum - controlSta.data[i].numOfPresent);
				controlData[i][1] = Utils.fromInt(controlSta.data[i].totalNum - controlSta.data[i].numOfPresent, Width);
			}
			aliceCase = gen.inputOfAlice(caseData);
			aliceControl = gen.inputOfAlice(controlData);
			bobCase = gen.inputOfBob(caseData);
			bobControl = gen.inputOfBob(controlData);
			numOfTests = caseSta.numberOftuples;
		}
		T[][] res;

		@Override
		public void secureCompute(CompEnv<T> gen) {
			res = compute(gen, aliceCase, bobCase, aliceControl, bobControl, numOfTests);
		}

		@Override
		public void prepareOutput(CompEnv<T> gen) {
			FixedPointLib<T> flib = new FixedPointLib<T>(gen, FWidth, FOffset);
			for(int i = 0; i < numOfTests; ++i)
				flib.outputToAlice(res[i]);
		}
	}
}
