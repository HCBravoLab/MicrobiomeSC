/*Original implementation of chi-square performed by xiao wang,
 *  slightly by justin wagner to handle microbiome count data
 */

package perFeatureChiSquare;

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
import circuits.CircuitLib;
import circuits.arithmetic.FloatLib;
import circuits.arithmetic.IntegerLib;
import flexsc.CompEnv;

public class ChiSquare {
	static public int Width = 32;
	static public int FWidth = 54;
	static public int FOffset = 11;
	static int filterThreshold = 5;

	static public<T> T[] dummyVariable(CompEnv<T> gen){
		FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
		T[] zero = flib.publicValue(0.0);
		return zero;
	}
	
	static public<T> T[] filter(CompEnv<T> gen, T[][] aliceCasePresent, T[][] aliceCaseNum,
			T[][] bobCasePresent, T[][] bobCaseNum,
			T[][] aliceControlPresent, T[][] aliceControlNum,
			T[][] bobControlPresent, T[][] bobControlNum, int numOfTests){

		IntegerLib<T> ilib = new IntegerLib<T>(gen, 32);
		CircuitLib<T> cl = new CircuitLib<T>(gen);
		T[] threshold = ilib.publicValue(filterThreshold);
		T[] filterResults = gen.newTArray(numOfTests);
		T[] caseNum;
		T[] controlNum;
		T[] totalPresent;
		T aboveThreshold;
		for(int i = 0; i < numOfTests; i++){			
			caseNum = ilib.add(aliceCasePresent[i], bobCasePresent[i]);
			controlNum = ilib.add(aliceControlPresent[i], bobControlPresent[i]);
			totalPresent = ilib.add(caseNum, controlNum);
			aboveThreshold = ilib.geq(totalPresent, threshold);
			filterResults[i] = aboveThreshold;
		}
		
		return filterResults;
	}
	
	public static<T> T[][] compute2(CompEnv<T> gen, T[][] aliceCasePresent, T[][] aliceCaseNum,
			T[][] bobCasePresent, T[][] bobCaseNum,
			T[][] aliceControlPresent, T[][] aliceControlNum,
			T[][] bobControlPresent, T[][] bobControlNum, int numOfTests) {

		T[][] res = gen.newTArray(numOfTests, 0);
		CircuitLib<T> cl = new CircuitLib<T>(gen);
		T[] aboveThreshold = gen.newTArray(1);
		IntegerLib<T> lib = new IntegerLib<T>(gen);
		FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
		T[] threshold = lib.publicValue(filterThreshold);
		T[] zero = flib.publicValue(0.0);
		for(int i = 0; i < numOfTests; i++){			
			T[] a = lib.add(aliceCasePresent[i], bobCasePresent[i]);
			T[] b = lib.add(aliceCaseNum[i], bobCaseNum[i]);
			T[] c = lib.add(aliceControlPresent[i], bobControlPresent[i]);
			T[] d = lib.add(aliceControlNum[i], bobControlNum[i]);

			T[] fa = lib.toSecureFloat(a, flib);
			T[] fb = lib.toSecureFloat(b, flib);
			T[] fc = lib.toSecureFloat(c, flib);
			T[] fd = lib.toSecureFloat(d, flib);

			T[] upperFirst = flib.add(fa, flib.add(fb, flib.add(fc, fd)));
			T[] upperSecond = flib.sub(flib.multiply(fb, fc), flib.multiply(fa, fd));
			upperSecond = flib.multiply(upperSecond, upperSecond);
			T[] upper = flib.multiply(upperFirst, upperSecond);
			T[] lower = flib.multiply(flib.multiply(flib.add(fa, fb), flib.add(fa, fc)), flib.multiply(flib.add(fb, fd), flib.add(fc, fd)));
			res[i] = flib.div(upper, lower);
		}
		return res;
	}
	
	public static<T> T[] compute(CompEnv<T> gen, T[] aliceCaseNumPresent, T[] aliceCaseTotalNum,
			T[] bobCaseNumPresent, T[] bobCaseTotalNum, T[] aliceControlNumPresent, T[] aliceControlTotalNum,
			T[] bobControlNumPresent, T[] bobControlTotalNum) {

		IntegerLib<T> lib = new IntegerLib<T>(gen);
		FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
		T[] res = flib.publicValue(0.0);

		T[] a = lib.add(aliceCaseNumPresent, bobCaseNumPresent);
		T[] b = lib.sub(lib.add(aliceCaseTotalNum, bobCaseTotalNum), a);
		T[] c = lib.add(aliceControlNumPresent, bobControlNumPresent);
		T[] d = lib.sub(lib.add(aliceControlTotalNum, bobControlTotalNum), c);

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
		T[][] aliceCaseNumPresentAr;
		T[][] aliceCaseTotalNumAr;

		T[][] bobCaseNumPresentAr;
		T[][] bobCaseTotalNumAr;

		T[][] aliceControlNumPresentAr;
		T[][] aliceControlTotalNumAr;

		T[][] bobControlNumPresentAr;
		T[][] bobControlTotalNumAr;

		T[][] aliceCaseNumPresentAr2;
		T[][] aliceCaseTotalNumAr2;

		T[][] bobCaseNumPresentAr2;
		T[][] bobCaseTotalNumAr2;

		T[][] aliceControlNumPresentAr2;
		T[][] aliceControlTotalNumAr2;

		T[][] bobControlNumPresentAr2;
		T[][] bobControlTotalNumAr2;
		
		T[] aliceCaseNumPresent;
		T[] aliceCaseTotalNum;

		T[] bobCaseNumPresent;
		T[] bobCaseTotalNum;

		T[] aliceControlNumPresent;
		T[] aliceControlTotalNum;

		T[] bobControlNumPresent;
		T[] bobControlTotalNum;
		
		Statistics[] caseSta;
		Statistics[] controlSta;
		int[] caseNumPresent;
		int[] caseTotalNum;
		int[] controlNumPresent;
		int[] controlTotalNum;
		T[] numAliceCase;
		T[] numBobCase;
		T[] numAliceControl;
		T[] numBobControl;
		
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
			caseNumPresent = new int[numOfTests];
			caseTotalNum = new int[numOfTests];
			controlNumPresent = new int[numOfTests];
			controlTotalNum = new int[numOfTests];

			for(int i =0; i < numOfTests; i++){
				caseNumPresent[i] = caseSta[i].numOfPresent;
				caseTotalNum[i] = caseSta[i].totalNum;
				controlNumPresent[i] = controlSta[i].numOfPresent;
				controlTotalNum[i] = controlSta[i].totalNum;
			}
			res = gen.newTArray(numOfTests, 0);
		}

		@Override
		public void secureCompute(CompEnv<T> gen) {
			T[] result = gen.newTArray(numOfTests);
			CircuitLib<T> cl = new CircuitLib<T>(gen);
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
			T[] zero = flib.publicValue(0.0);
			int numFiltered = 0;
			aliceCaseNumPresentAr = gen.newTArray(numOfTests, 0);
			aliceCaseTotalNumAr = gen.newTArray(numOfTests,0);
			bobCaseNumPresentAr = gen.newTArray(numOfTests,0);
			bobCaseTotalNumAr = gen.newTArray(numOfTests,0);
			aliceControlNumPresentAr = gen.newTArray(numOfTests,0);
			aliceControlTotalNumAr = gen.newTArray(numOfTests,0);
			bobControlNumPresentAr = gen.newTArray(numOfTests,0);
			bobControlTotalNumAr = gen.newTArray(numOfTests,0);

			for(int i = 0; i < numOfTests; i++){				
				aliceCaseNumPresentAr[i] = gen.inputOfAlice(Utils.fromInt(caseNumPresent[i], Width));
				aliceCaseTotalNumAr[i] = gen.inputOfAlice(Utils.fromInt(caseTotalNum[i], Width));
				bobCaseNumPresentAr[i] = gen.inputOfBob(new boolean[Width]);
				bobCaseTotalNumAr[i] = gen.inputOfBob(new boolean[Width]);
				aliceControlNumPresentAr[i] = gen.inputOfAlice(Utils.fromInt(controlNumPresent[i], Width));
				aliceControlTotalNumAr[i] = gen.inputOfAlice(Utils.fromInt(controlTotalNum[i], Width));
				bobControlNumPresentAr[i] = gen.inputOfBob(new boolean[Width]);
				bobControlTotalNumAr[i] = gen.inputOfBob(new boolean[Width]);
			}
			result = filter(gen, aliceCaseNumPresentAr, aliceCaseTotalNumAr, 
					bobCaseNumPresentAr, bobCaseTotalNumAr,
					aliceControlNumPresentAr, aliceControlTotalNumAr,
					bobControlNumPresentAr, bobControlTotalNumAr, numOfTests);
			
			boolean[] filResOut = cl.declassifyToBoth(result);
		
			//System.out.println(filRes[0]);
			int[] indices = new int[numOfTests];
			for(int i = 0; i < numOfTests; i++){
				if(filResOut[i]){
					indices[numFiltered] = i;
					numFiltered++;
				}
			}
			aliceCaseNumPresentAr2 = gen.newTArray(numFiltered, 0);
			aliceCaseTotalNumAr2 = gen.newTArray(numFiltered,0);
			bobCaseNumPresentAr2 = gen.newTArray(numFiltered,0);
			bobCaseTotalNumAr2 = gen.newTArray(numFiltered,0);
			aliceControlNumPresentAr2 = gen.newTArray(numFiltered,0);
			aliceControlTotalNumAr2 = gen.newTArray(numFiltered,0);
			bobControlNumPresentAr2 = gen.newTArray(numFiltered,0);
			bobControlTotalNumAr2 = gen.newTArray(numFiltered,0);
			for(int i = 0; i < numFiltered; i++){				
					aliceCaseNumPresentAr2[i] = aliceCaseNumPresentAr[indices[i]];
					aliceCaseTotalNumAr2[i] = aliceCaseTotalNumAr[indices[i]];
					bobCaseNumPresentAr2[i] = bobCaseNumPresentAr[indices[i]];
					bobCaseTotalNumAr2[i] = bobCaseTotalNumAr[indices[i]];
					aliceControlNumPresentAr2[i] = aliceControlNumPresentAr[indices[i]];
					aliceControlTotalNumAr2[i] = aliceControlTotalNumAr[indices[i]];
					bobControlNumPresentAr2[i] = bobControlNumPresentAr[indices[i]];
					bobControlTotalNumAr2[i] = bobControlTotalNumAr[indices[i]];
			}
				
			res = compute2(gen, aliceCaseNumPresentAr2, aliceCaseTotalNumAr2, 
					bobCaseNumPresentAr2, bobCaseTotalNumAr2,
					aliceControlNumPresentAr2, aliceControlTotalNumAr2,
					bobControlNumPresentAr2, bobControlTotalNumAr2, numFiltered);
		}

		@Override
		public void prepareOutput(CompEnv<T> gen) {
			ChiSquaredDistribution chiDistribution = new ChiSquaredDistribution(1.0);
			System.out.println("chi,p-value");
			for(int i = 0; i < res.length; i++){
				double chi = Utils.toFloat(gen.outputToAlice(res[i]), FWidth, FOffset);
				if(chi == 0.0){
					System.out.println("NA,NA");
					continue;
				}
				System.out.println(chi + "," + (1-chiDistribution.cumulativeProbability(chi)));
			}	
		}
	}

	public static class Evaluator<T> extends EvaRunnable<T> {
		T[][] aliceCaseNumPresentAr;
		T[][] aliceCaseTotalNumAr;

		T[][] bobCaseNumPresentAr;
		T[][] bobCaseTotalNumAr;

		T[][] aliceControlNumPresentAr;
		T[][] aliceControlTotalNumAr;

		T[][] bobControlNumPresentAr;
		T[][] bobControlTotalNumAr;

		T[][] aliceCaseNumPresentAr2;
		T[][] aliceCaseTotalNumAr2;

		T[][] bobCaseNumPresentAr2;
		T[][] bobCaseTotalNumAr2;

		T[][] aliceControlNumPresentAr2;
		T[][] aliceControlTotalNumAr2;

		T[][] bobControlNumPresentAr2;
		T[][] bobControlTotalNumAr2;

		T[] aliceCaseNumPresent;
		T[] aliceCaseTotalNum;

		T[] bobCaseNumPresent;
		T[] bobCaseTotalNum;

		T[] aliceControlNumPresent;
		T[] aliceControlTotalNum;

		T[] bobControlNumPresent;
		T[] bobControlTotalNum;

		int[] caseNumPresent;
		int[] caseTotalNum;
		int[] controlNumPresent;
		int[] controlTotalNum;
		int numOfTests;
		boolean precise;
		T[][] res;

		T[] numAliceCase;
		T[] numBobCase;
		T[] numAliceControl;
		T[] numBobControl;
		
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
			
			numOfTests = caseSta.length;

			caseNumPresent = new int[numOfTests];
			caseTotalNum = new int[numOfTests];
			controlNumPresent = new int[numOfTests];
			controlTotalNum = new int[numOfTests];
			for(int i =0; i < numOfTests; i++){
				caseNumPresent[i] = caseSta[i].numOfPresent;
				caseTotalNum[i] = caseSta[i].totalNum;
				controlNumPresent[i] = controlSta[i].numOfPresent;
				controlTotalNum[i] = controlSta[i].totalNum;
			}
			res = gen.newTArray(numOfTests, 0);
		}

		@Override
		public void secureCompute(CompEnv<T> gen) {
			T[] result = gen.newTArray(numOfTests);
			CircuitLib<T> cl = new CircuitLib<T>(gen);
			FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
			
			T[] zero = flib.publicValue(0.0);
			aliceCaseNumPresentAr = gen.newTArray(numOfTests, 0);
			aliceCaseTotalNumAr = gen.newTArray(numOfTests,0);
			bobCaseNumPresentAr = gen.newTArray(numOfTests,0);
			bobCaseTotalNumAr = gen.newTArray(numOfTests,0);
			aliceControlNumPresentAr = gen.newTArray(numOfTests,0);
			aliceControlTotalNumAr = gen.newTArray(numOfTests,0);
			bobControlNumPresentAr = gen.newTArray(numOfTests,0);
			bobControlTotalNumAr = gen.newTArray(numOfTests,0);
			int numFiltered = 0;
			for(int i = 0; i < numOfTests; i++){
				aliceCaseNumPresentAr[i] = gen.inputOfAlice(new boolean[Width]);
				aliceCaseTotalNumAr[i] = gen.inputOfAlice(new boolean[Width]);
				bobCaseNumPresentAr[i] = gen.inputOfBob(Utils.fromInt(caseNumPresent[i], Width));
				bobCaseTotalNumAr[i] = gen.inputOfBob(Utils.fromInt(caseTotalNum[i], Width));
				aliceControlNumPresentAr[i] = gen.inputOfAlice(new boolean[Width]);
				aliceControlTotalNumAr[i] = gen.inputOfAlice(new boolean[Width]);
				bobControlNumPresentAr[i] = gen.inputOfBob(Utils.fromInt(controlNumPresent[i], Width));
				bobControlTotalNumAr[i] = gen.inputOfBob(Utils.fromInt(controlTotalNum[i], Width));
			}
			result = filter(gen, aliceCaseNumPresentAr, aliceCaseTotalNumAr, 
					bobCaseNumPresentAr, bobCaseTotalNumAr,
					aliceControlNumPresentAr, aliceControlTotalNumAr,
					bobControlNumPresentAr, bobControlTotalNumAr, numOfTests);
			
			boolean[] filResOut = cl.declassifyToBoth(result);
			//System.out.println(filRes[0]);
			int[] indices = new int[numOfTests];
			for(int i = 0; i < numOfTests; i++){
				if(filResOut[i]){
					indices[numFiltered] = i;
					numFiltered++;
				}
			}

			aliceCaseNumPresentAr2 = gen.newTArray(numFiltered, 0);
			aliceCaseTotalNumAr2 = gen.newTArray(numFiltered,0);
			bobCaseNumPresentAr2 = gen.newTArray(numFiltered,0);
			bobCaseTotalNumAr2 = gen.newTArray(numFiltered,0);
			aliceControlNumPresentAr2 = gen.newTArray(numFiltered,0);
			aliceControlTotalNumAr2 = gen.newTArray(numFiltered,0);
			bobControlNumPresentAr2 = gen.newTArray(numFiltered,0);
			bobControlTotalNumAr2 = gen.newTArray(numFiltered,0);
			for(int i = 0; i < numFiltered; i++){				
					aliceCaseNumPresentAr2[i] = aliceCaseNumPresentAr[indices[i]];
					aliceCaseTotalNumAr2[i] = aliceCaseTotalNumAr[indices[i]];
					bobCaseNumPresentAr2[i] = bobCaseNumPresentAr[indices[i]];
					bobCaseTotalNumAr2[i] = bobCaseTotalNumAr[indices[i]];
					aliceControlNumPresentAr2[i] = aliceControlNumPresentAr[indices[i]];
					aliceControlTotalNumAr2[i] = aliceControlTotalNumAr[indices[i]];
					bobControlNumPresentAr2[i] = bobControlNumPresentAr[indices[i]];
					bobControlTotalNumAr2[i] = bobControlTotalNumAr[indices[i]];
			}
	
			res = compute2(gen, aliceCaseNumPresentAr2, aliceCaseTotalNumAr2, 
					bobCaseNumPresentAr2, bobCaseTotalNumAr2,
					aliceControlNumPresentAr2, aliceControlTotalNumAr2,
					bobControlNumPresentAr2, bobControlTotalNumAr2, numFiltered);
		}
		@Override
		public void prepareOutput(CompEnv<T> gen) {
			for(int j =0; j<res.length; j++){
				gen.outputToAlice(res[j]);
			}
		}
	}
}
