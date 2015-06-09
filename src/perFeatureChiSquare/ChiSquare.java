/*Original implementation of chi-square performed by xiao wang,
 *  slightly by justin wagner to handle microbiome count data
 */

package perFeatureChiSquare;

import java.nio.ByteBuffer;

import network.Server;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;

import perFeatureChiSquare.PrepareData;
import perFeatureChiSquare.PrepareData.StatisticsData;

import util.EvaRunnable;
import util.GenRunnable;
import util.Utils;
import circuits.arithmetic.FloatLib;
import circuits.arithmetic.IntegerLib;
import flexsc.CompEnv;

public class ChiSquare {
	static public int Width = 32;
	static public int FWidth = 54;
	static public int FOffset = 11;

	public static<T> T[] compute(CompEnv<T> gen, T[][] aliceCase,
			T[][] bobCase,
			T[][] aliceControl,
			T[][] bobControl) {

		IntegerLib<T> lib = new IntegerLib<T>(gen);
		FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
		T[] res = flib.publicValue(0.0);

		T[] a = lib.add(aliceCase[0], bobCase[0]);
		T[] b = lib.add(aliceCase[1], bobCase[1]);
		T[] c = lib.add(aliceControl[0], bobControl[0]);
		T[] d = lib.add(aliceControl[1], bobControl[1]);

		T[] fa = lib.toSecureFloat(a, flib);
		T[] fb = lib.toSecureFloat(b, flib);
		T[] fc = lib.toSecureFloat(c, flib);
		T[] fd = lib.toSecureFloat(d, flib);

		T[] upperFirst = flib.add(fa, flib.add(fb, flib.add(fc, fd)));
		T[] upperSecond = flib.sub(flib.multiply(fb, fc), flib.multiply(fa, fd));
		upperSecond = flib.multiply(upperSecond, upperSecond);
		T[] upper = flib.multiply(upperFirst, upperSecond);
		T[] lower = flib.multiply(flib.multiply(flib.add(fa, fb), flib.add(fa, fc)), flib.multiply(flib.add(fb, fd), flib.add(fc, fd)));
		res = flib.div(upper, lower);

		return res;

	}
	public static class Generator<T> extends GenRunnable<T> {
		T[][] aliceCase;
		T[][] bobCase;
		T[][] aliceControl;
		T[][] bobControl;

		Statistics[] caseSta;
		Statistics[] controlSta;

		int numOfTests;
		T[][] res;

		double extraFactor;
		boolean precise;
		@Override
		public void prepareInput(CompEnv<T> gen) throws Exception {			
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
			caseSta = caseInput.data;
			controlSta = controlInput.data;

			numOfTests = caseSta.length;
			res = gen.newTArray(numOfTests, 0);
		}

		@Override
		public void secureCompute(CompEnv<T> gen) {
			aliceCase = gen.newTArray(2, 0);
			bobCase = gen.newTArray(2, 0);
			aliceControl = gen.newTArray(2, 0);
			bobControl = gen.newTArray(2, 0);
			for(int i = 1; i < 2; i++){
				aliceCase[0] = gen.inputOfAlice(Utils.fromInt(caseSta[i].numOfPresent, Width));
				aliceCase[1] = gen.inputOfAlice(Utils.fromInt(caseSta[i].totalNum - caseSta[i].numOfPresent, Width));
				bobCase[0] = gen.inputOfBob(new boolean[Width]);
				bobCase[1] = gen.inputOfBob(new boolean[Width]);
				aliceControl[0] = gen.inputOfAlice(Utils.fromInt(controlSta[i].numOfPresent, Width));
				aliceControl[1] = gen.inputOfAlice(Utils.fromInt(controlSta[i].totalNum - controlSta[i].numOfPresent, Width));
				bobControl[0] = gen.inputOfBob(new boolean[Width]);
				bobControl[1] = gen.inputOfBob(new boolean[Width]);
				res[i] = compute(gen, aliceCase, bobCase, aliceControl, bobControl);
			}
		}

		@Override
		public void prepareOutput(CompEnv<T> gen) {
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
			ChiSquaredDistribution chiDistribution = new ChiSquaredDistribution(1.0);
			System.out.println("chi,p-value");
			for(int i = 1; i < numOfTests; ++i){
				double chi = flib.outputToAlice(res[i]);// * extraFactor;
				if(chi == 0.0){
					System.out.println("NA,NA");
					continue;
				}
				System.out.println(chi + "," + (1-chiDistribution.cumulativeProbability(chi)));
			}
		}
	}

	public static class Evaluator<T> extends EvaRunnable<T> {
		T[][] aliceCase;
		T[][] bobCase;
		T[][] aliceControl;
		T[][] bobControl;
		int numOfTests;
		boolean precise;
		T[][] res;

		Statistics[] caseSta;
		Statistics[] controlSta;
		@Override
		public void prepareInput(CompEnv<T> gen) throws Exception {
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
			caseSta = caseInput.data;
			controlSta = controlInput.data;

			gen.channel.writeInt(caseInput.numberOftuples);
			gen.channel.writeInt(controlInput.numberOftuples);
			gen.channel.flush();
			
			numOfTests = caseSta.length;
			res = gen.newTArray(numOfTests, 0);
		}

		@Override
		public void secureCompute(CompEnv<T> gen) {
			aliceCase = gen.newTArray(2, 0);
			aliceControl = gen.newTArray(2, 0);
			bobCase = gen.newTArray(2, 0);
			bobControl = gen.newTArray(2, 0);
			for(int i = 1; i < 2; i++){
				aliceCase[0] = gen.inputOfAlice(new boolean[Width]);
				aliceCase[1] = gen.inputOfAlice(new boolean[Width]);
				bobCase[0] = gen.inputOfBob(Utils.fromInt(caseSta[i].numOfPresent, Width));
				bobCase[1] = gen.inputOfBob(Utils.fromInt(caseSta[i].totalNum - caseSta[i].numOfPresent, Width));
				aliceControl[0] = gen.inputOfAlice(new boolean[Width]);
				aliceControl[1] = gen.inputOfAlice(new boolean[Width]);
				bobControl[0] = gen.inputOfBob(Utils.fromInt(controlSta[i].numOfPresent, Width));
				bobControl[1] = gen.inputOfBob(Utils.fromInt(controlSta[i].totalNum - 
						controlSta[i].numOfPresent, Width));
				res[i] = compute(gen, aliceCase, bobCase, aliceControl, bobControl);
			}
		}

		@Override
		public void prepareOutput(CompEnv<T> gen) {
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
			for(int i = 1; i < numOfTests; ++i)
				flib.outputToAlice(res[i]);
		}
	}
}
