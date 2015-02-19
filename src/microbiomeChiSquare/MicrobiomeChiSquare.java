package microbiomeChiSquare;

import java.nio.ByteBuffer;

import network.Server;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;

import microbiomeChiSquare.PrepareData;
import microbiomeChiSquare.PrepareData.StatisticsData;
import microbiomeChiSquare.Statistics;
import util.EvaRunnable;
import util.GenRunnable;
import util.Utils;
import circuits.arithmetic.FloatLib;
import circuits.arithmetic.IntegerLib;
import flexsc.CompEnv;

public class MicrobiomeChiSquare {
	static public int Width = 12;

	public static<T> FloatLib<T> getfloatLib(CompEnv<T> env, boolean precise){
		if(precise)
			return new FloatLib<T>(env, 30, 10);
		else return new FloatLib<T>(env, 20, 10);
	}
	
	public static<T> T[][] compute(CompEnv<T> env, T[][][] aliceCase,
			T[][][] bobCase,
			T[][][] aliceControl,
			T[][][] bobControl, int numOfTests, boolean precise) {
		T[][] res = env.newTArray(numOfTests, 0);
		IntegerLib<T> lib = new IntegerLib<T>(env);
		FloatLib<T> flib = getfloatLib(env, precise);
		for(int i = 0; i < numOfTests; ++i) {
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
			Options options = new Options();
			options.addOption("h", false, "high precision");
			options.addOption("s", "case", true, "case");
			options.addOption("t", "control", true, "control");

			CommandLineParser parser = new BasicParser();
			CommandLine cmd = parser.parse(options, args);

			precise = cmd.hasOption("t");
			if(!cmd.hasOption("s") || !cmd.hasOption("t")) {
				throw new Exception("wrong input");
			}
			StatisticsData caseInput = PrepareData.readFile(cmd.getOptionValue("s"));
			StatisticsData controlInput = PrepareData.readFile(cmd.getOptionValue("t"));
			Statistics[] caseSta = caseInput.data;
			Statistics[] controlSta = controlInput.data;
			boolean[][][] caseData = new boolean[caseSta.length][2][Width];

			int caseLength = ByteBuffer.wrap(Server.readBytes(env.is, 4)).getInt();
			int controlLength = ByteBuffer.wrap(Server.readBytes(env.is, 4)).getInt();

			//extraFactor = n/(r*s)
			
			extraFactor = 2*(caseLength + controlLength + caseInput.numberOftuples+ controlInput.numberOftuples);
			extraFactor /= 2*(caseLength + caseInput.numberOftuples);
			extraFactor /= 2*(controlLength + controlInput.numberOftuples);
			
			for(int i = 0; i < caseSta.length; ++i) {
				caseData[i][0] = Utils.fromInt(caseSta[i].numOfPresent, Width);
				caseData[i][1] = Utils.fromInt(caseSta[i].totalNum - caseSta[i].numOfPresent, Width);
			}

			boolean[][][] controlData = new boolean[controlSta.length][2][Width];
			for(int i = 0; i < caseSta.length; ++i) {
				controlData[i][0] = Utils.fromInt(controlSta[i].numOfPresent, Width);
				controlData[i][1] = Utils.fromInt(controlSta[i].totalNum - controlSta[i].numOfPresent, Width);
			}
			aliceCase = env.inputOfAlice(caseData);
			aliceControl = env.inputOfAlice(controlData);
			bobCase = env.inputOfBob(caseData);
			bobControl = env.inputOfBob(controlData);
			numOfTests = caseSta.length;
		}

		T[][] res;
		@Override
		public void secureCompute(CompEnv<T> env) {
			res = compute(env, aliceCase, bobCase, aliceControl, bobControl, numOfTests, precise);
		}

		@Override
		public void prepareOutput(CompEnv<T> env) {
			FloatLib<T> flib = getfloatLib(env, precise);
			for(int i = 0; i < res.length; ++i)
				System.out.println(flib.outputToAlice(res[i])*extraFactor);
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
			Options options = new Options();
			options.addOption("h", "high_precision", false, "high precision");
			options.addOption("s", "case", true, "case");
			options.addOption("t", "control", true, "control");

			CommandLineParser parser = new BasicParser();
			CommandLine cmd = parser.parse(options, args);

			precise = cmd.hasOption("h");
			if(!cmd.hasOption("s") || !cmd.hasOption("t")) {
				throw new Exception("wrong input");
			}

			StatisticsData caseInput = PrepareData.readFile(cmd.getOptionValue("s"));
			StatisticsData controlInput = PrepareData.readFile(cmd.getOptionValue("t"));
			Statistics[] caseSta = caseInput.data;
			Statistics[] controlSta = controlInput.data;
			boolean[][][] caseData = new boolean[caseSta.length][2][Width];

			env.os.write(ByteBuffer.allocate(4).putInt(caseInput.numberOftuples).array());
			env.os.write(ByteBuffer.allocate(4).putInt(controlInput.numberOftuples).array());
			env.os.flush();

			for(int i = 0; i < caseSta.length; ++i) {
				caseData[i][0] = Utils.fromInt(caseSta[i].numOfPresent, Width);
				caseData[i][1] = Utils.fromInt(caseSta[i].totalNum - caseSta[i].numOfPresent, Width);
			}

			boolean[][][] controlData = new boolean[controlSta.length][2][Width];
			for(int i = 0; i < caseSta.length; ++i) {
				controlData[i][0] = Utils.fromInt(controlSta[i].numOfPresent, Width);
				controlData[i][1] = Utils.fromInt(controlSta[i].totalNum - controlSta[i].numOfPresent, Width);
			}
			aliceCase = env.inputOfAlice(caseData);
			aliceControl = env.inputOfAlice(controlData);
			bobCase = env.inputOfBob(caseData);
			bobControl = env.inputOfBob(controlData);
			numOfTests = caseSta.length;
		}
		
		T[][] res;

		@Override
		public void secureCompute(CompEnv<T> env) {
			res = compute(env, aliceCase, bobCase, aliceControl, bobControl, numOfTests, precise);
		}

		@Override
		public void prepareOutput(CompEnv<T> env) {
			FloatLib<T> flib = getfloatLib(env, precise);
			for(int i = 0; i < numOfTests; ++i)
				flib.outputToAlice(res[i]);
		}
	}
}
