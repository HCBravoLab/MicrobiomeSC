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

import precompute.PrepareDataOddsRatio;
import precompute.PrepareDataOddsRatio.StatisticsData;
import util.EvaRunnable;
import util.GenRunnable;
import util.Utils;
import circuits.arithmetic.FixedPointLib;
import circuits.arithmetic.FloatLib;
import circuits.arithmetic.IntegerLib;
import flexsc.CompEnv;

public class OddsRatio {
	static public int Width = 32;
	static public int FWidth = 54;
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
			StatisticsData caseSta = PrepareDataOddsRatio.readFile(args[0]);
			StatisticsData controlSta = PrepareDataOddsRatio.readFile(args[1]);
			FloatLib<T> flib = new FloatLib<T>(env, FWidth, FOffset);
			T[] l = flib.publicValue(0.0);

			boolean[][][] caseData = new boolean[caseSta.numberOftuples][2][l.length];
			for(int i = 0; i < caseSta.numberOftuples; ++i) {
				caseData[i][0] = Utils.fromInt(caseSta.data[i].numOfPresent, Width);
				caseData[i][1] = Utils.fromInt(caseSta.data[i].totalNum - caseSta.data[i].numOfPresent, Width);
			}

			boolean[][][] controlData = new boolean[controlSta.numberOftuples][2][l.length];
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
			System.out.println("Odds Ratio");
			for(int i = 0; i < res.length; ++i)
				System.out.println(flib.outputToAlice(res[i]));
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
			StatisticsData caseSta = PrepareDataOddsRatio.readFile(args[0]);
			StatisticsData controlSta = PrepareDataOddsRatio.readFile(args[1]);
			FloatLib<T> flib = new FloatLib<T>(env, FWidth, FOffset);
			T[] l = flib.publicValue(0.0);

			boolean[][][] caseData = new boolean[caseSta.numberOftuples][2][l.length];
			for(int i = 0; i < caseSta.numberOftuples; ++i) {
				caseData[i][0] = Utils.fromInt(caseSta.data[i].numOfPresent, Width);
				caseData[i][1] = Utils.fromInt(caseSta.data[i].totalNum - caseSta.data[i].numOfPresent, Width);
			}

			boolean[][][] controlData = new boolean[controlSta.numberOftuples][2][l.length];
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
