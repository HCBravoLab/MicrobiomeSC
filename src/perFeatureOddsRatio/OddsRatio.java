/*Original implementation of chi-square performed by xiao wang,
 *  slightly by justin wagner to handle microbiome count data
 */

package perFeatureOddsRatio;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;

import perFeatureOddsRatio.PrepareData;
import perFeatureOddsRatio.PrepareData.StatisticsData;

import util.EvaRunnable;
import util.GenRunnable;
import util.Utils;
import circuits.CircuitLib;
import circuits.arithmetic.FloatLib;
import circuits.arithmetic.IntegerLib;
import flexsc.CompEnv;

public class OddsRatio {
	static public int Width = 32;
	static public int FWidth = 54;
	static public int FOffset = 11;
	static int filterThreshold = 1;

	static public<T> T[] dummyVariable(CompEnv<T> gen){
		FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
		T[] zero = flib.publicValue(0.0);
		return zero;
	}
	

	static public<T> T filter(CompEnv<T> gen, T[] numAliceCase, T[] aliceCaseTotal,
			T[] numBobCase, T[] bobCaseTotal, T[] numAliceControl, T[] aliceControlTotal,
			T[] numBobControl, T[] bobControlTotal){

		//FloatLib<T> flib = new FloatLib<T>(gen, FWidth, FOffset);
		IntegerLib<T> ilib = new IntegerLib<T>(gen, 32);
		CircuitLib<T> cl = new CircuitLib<T>(gen);
		
		//T[] threshold = flib.publicValue(filterThreshold);
		T[] caseTotal = ilib.add(aliceCaseTotal, bobCaseTotal);
		T[] controlTotal = ilib.add(aliceControlTotal, bobControlTotal);
		T greaterCaseOrControl = ilib.geq(caseTotal, controlTotal);
		T[] threshold = cl.mux(caseTotal, controlTotal, greaterCaseOrControl);
		T[] caseNum = ilib.add(numAliceCase, numBobCase);
		T[] controlNum = ilib.add(numAliceControl, numBobControl);
		T[] totalPresent = ilib.add(caseNum, controlNum);
		T aboveThreshold = ilib.geq(totalPresent, threshold);
		return aboveThreshold;
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

		res = flib.div(flib.multiply(fa, fd), flib.multiply(fb, fc));
		return res;
	}
	
	public static class Generator<T> extends GenRunnable<T> {
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
			T result;
			
			boolean[] filteredFeatures = new boolean[numOfTests];
			for(int i =0; i < numOfTests; i++){				

				aliceCaseNumPresent = gen.newTArray(Width);
				aliceCaseTotalNum = gen.newTArray(Width);
				bobCaseNumPresent = gen.newTArray(Width);
				bobCaseTotalNum = gen.newTArray(Width);
				aliceControlNumPresent = gen.newTArray(Width);
				aliceControlTotalNum = gen.newTArray(Width);
				bobControlNumPresent = gen.newTArray(Width);
				bobControlTotalNum = gen.newTArray(Width);
				aliceCaseNumPresent = gen.inputOfAlice(Utils.fromInt(caseNumPresent[i], Width));
				aliceCaseTotalNum = gen.inputOfAlice(Utils.fromInt(caseTotalNum[i], Width));
				bobCaseNumPresent = gen.inputOfBob(new boolean[Width]);
				bobCaseTotalNum = gen.inputOfBob(new boolean[Width]);
				aliceControlNumPresent = gen.inputOfAlice(Utils.fromInt(controlNumPresent[i], Width));
				aliceControlTotalNum = gen.inputOfAlice(Utils.fromInt(controlTotalNum[i], Width));
				bobControlNumPresent = gen.inputOfBob(new boolean[Width]);
				bobControlTotalNum = gen.inputOfBob(new boolean[Width]);
				result = filter(gen, aliceCaseNumPresent, aliceCaseTotalNum, 
					bobCaseNumPresent, bobCaseTotalNum,
					aliceControlNumPresent, aliceControlTotalNum,
					bobControlNumPresent, bobControlTotalNum);
				filteredFeatures[i] = gen.outputToAlice(result);
			}
				
			for(int i = 0; i < numOfTests; i++){
				if(!(filteredFeatures[i])){
					res[i] = dummyVariable(gen);
					continue;
				}
				aliceCaseNumPresent = gen.newTArray(Width);
				aliceCaseTotalNum = gen.newTArray(Width);
				bobCaseNumPresent = gen.newTArray(Width);
				bobCaseTotalNum = gen.newTArray(Width);
				aliceControlNumPresent = gen.newTArray(Width);
				aliceControlTotalNum = gen.newTArray(Width);
				bobControlNumPresent = gen.newTArray(Width);
				bobControlTotalNum = gen.newTArray(Width);
				
				aliceCaseNumPresent = gen.inputOfAlice(Utils.fromInt(caseNumPresent[i], Width));
				aliceCaseTotalNum = gen.inputOfAlice(Utils.fromInt(caseTotalNum[i], Width));
				bobCaseNumPresent = gen.inputOfBob(new boolean[Width]);
				bobCaseTotalNum = gen.inputOfBob(new boolean[Width]);
				aliceControlNumPresent = gen.inputOfAlice(Utils.fromInt(controlNumPresent[i], Width));
				aliceControlTotalNum = gen.inputOfAlice(Utils.fromInt(controlTotalNum[i], Width));
				bobControlNumPresent = gen.inputOfBob(new boolean[Width]);
				bobControlTotalNum = gen.inputOfBob(new boolean[Width]);
				res[i] = compute(gen, aliceCaseNumPresent, aliceCaseTotalNum, 
						bobCaseNumPresent, bobCaseTotalNum,
						aliceControlNumPresent, aliceControlTotalNum,
						bobControlNumPresent, bobControlTotalNum);
			}
		}

		@Override
		public void prepareOutput(CompEnv<T> gen) {
			System.out.println("Odds Ratio");
			for(int i = 0; i < res.length; i++){
				double OR = Utils.toFloat(gen.outputToAlice(res[i]), FWidth, FOffset);
				if(OR == 0.0){
					System.out.println("NA");
					continue;
				}
				System.out.println(OR);
			}	
		}
	}

	public static class Evaluator<T> extends EvaRunnable<T> {
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
			T result;
		
			boolean[] filteredFeatures = new boolean[numOfTests];
			for(int i =0; i < numOfTests; i++){

				aliceCaseNumPresent = gen.newTArray(Width);
				aliceCaseTotalNum = gen.newTArray(Width);
				bobCaseNumPresent = gen.newTArray(Width);
				bobCaseTotalNum = gen.newTArray(Width);
				aliceControlNumPresent = gen.newTArray(Width);
				aliceControlTotalNum = gen.newTArray(Width);
				bobControlNumPresent = gen.newTArray(Width);
				bobControlTotalNum = gen.newTArray(Width);
				
				aliceCaseNumPresent = gen.inputOfAlice(new boolean[Width]);
				aliceCaseTotalNum = gen.inputOfAlice(new boolean[Width]);
				bobCaseNumPresent = gen.inputOfBob(Utils.fromInt(caseNumPresent[i], Width));
				bobCaseTotalNum = gen.inputOfBob(Utils.fromInt(caseTotalNum[i], Width));
				aliceControlNumPresent = gen.inputOfAlice(new boolean[Width]);
				aliceControlTotalNum = gen.inputOfAlice(new boolean[Width]);
				bobControlNumPresent = gen.inputOfBob(Utils.fromInt(controlNumPresent[i], Width));
				bobControlTotalNum = gen.inputOfBob(Utils.fromInt(controlTotalNum[i], Width));
				result = filter(gen, aliceCaseNumPresent, aliceCaseTotalNum, 
						bobCaseNumPresent, bobCaseTotalNum,
						aliceControlNumPresent, aliceControlTotalNum,
						bobControlNumPresent, bobControlTotalNum);
				filteredFeatures[i] = gen.outputToAlice(result);
			}
			
			for(int i = 0; i < numOfTests; i++){
				if(!(filteredFeatures[i])){
					res[i] = dummyVariable(gen);
					continue;
				}
				aliceCaseNumPresent = gen.newTArray(Width);
				aliceCaseTotalNum = gen.newTArray(Width);
				bobCaseNumPresent = gen.newTArray(Width);
				bobCaseTotalNum = gen.newTArray(Width);
				aliceControlNumPresent = gen.newTArray(Width);
				aliceControlTotalNum = gen.newTArray(Width);
				bobControlNumPresent = gen.newTArray(Width);
				bobControlTotalNum = gen.newTArray(Width);
				
				aliceCaseNumPresent = gen.inputOfAlice(new boolean[Width]);
				aliceCaseTotalNum = gen.inputOfAlice(new boolean[Width]);
				bobCaseNumPresent = gen.inputOfBob(Utils.fromInt(caseNumPresent[i], Width));
				bobCaseTotalNum = gen.inputOfBob(Utils.fromInt(caseTotalNum[i], Width));
				aliceControlNumPresent = gen.inputOfAlice(new boolean[Width]);
				aliceControlTotalNum = gen.inputOfAlice(new boolean[Width]);
				bobControlNumPresent = gen.inputOfBob(Utils.fromInt(controlNumPresent[i], Width));
				bobControlTotalNum = gen.inputOfBob(Utils.fromInt(controlTotalNum[i], Width));
				
				res[i] = compute(gen, aliceCaseNumPresent, aliceCaseTotalNum, 
						bobCaseNumPresent, bobCaseTotalNum,
						aliceControlNumPresent, aliceControlTotalNum,
						bobControlNumPresent, bobControlTotalNum);
			}
		}

		@Override
		public void prepareOutput(CompEnv<T> gen) {
			for(int j =0; j<res.length; j++){
				gen.outputToAlice(res[j]);
			}
		}
	}
}
