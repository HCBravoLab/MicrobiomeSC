package sparseDifferentialAbundance;


import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;

import util.EvaRunnable;
import util.GenRunnable;
import util.Utils;
import circuits.BitonicSortLib;
import circuits.arithmetic.FixedPointLib;
import circuits.arithmetic.FloatLib;
import circuits.arithmetic.IntegerLib;
import flexsc.CompEnv;

public class DifferentialAbundanceNaive {
	static int PLength = 30;
	static int VLength = 10;
	static public<T> T[][][] compute(CompEnv<T> gen, T[][][] inputAliceCase, 
			T[][][] inputBobCase, T[][][] inputAliceControl, T[][][] inputBobControl){//, T[][][] inputAliceControl, T[][][] inputBobControl){
		T[][][] inCase = gen.newTArray((inputAliceCase.length + inputBobCase.length)/2 , inputAliceCase[0].length + inputBobCase[0].length, 0);
		FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
		for(int i = 0; i < (inputAliceCase.length + inputBobCase.length)/2; i++){
			System.arraycopy(inputAliceCase[i], 0, inCase[i], 0, inputAliceCase[i].length);
			System.arraycopy(inputBobCase[i], 0, inCase[i], inputAliceCase[i].length, inputBobCase[i].length);
		}

		T[][][] inControl = gen.newTArray((inputAliceControl.length + inputBobControl.length)/2, inputAliceControl[0].length + inputBobControl[0].length, 0);
		for(int i = 0; i < (inputAliceControl.length + inputBobControl.length)/2; i++){
			System.arraycopy(inputAliceControl[i], 0, inControl[i], 0, inputAliceControl[i].length);
			System.arraycopy(inputBobControl[i], 0, inControl[i], inputAliceControl[i].length, inputBobControl[i].length);
		}

		T[][][] caseSum = gen.newTArray((inputAliceCase.length + inputBobCase.length)/2, 2,0);
		T[][][] caseSumOfSquares = gen.newTArray((inputAliceCase.length + inputBobCase.length)/2, 2,0);

		for(int i = 0; i < (inputAliceCase.length + inputBobCase.length)/2; i++){
			for (int j = 0; j < inCase[0].length; j++){
				caseSum[i][0] = flib.add(caseSum[i][0], inCase[i][j]);
				caseSumOfSquares[i][1] = flib.add(caseSumOfSquares[i][0], flib.multiply(inCase[i][j], inCase[i][j]));
			}
		}
		T[] tStat;
/*
		for(int i = 0; i < inputAliceCase[0].length; i++){	
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

			caseSumOfSquares = outCase[3][i];
			controlSumOfSquares = outControl[3][i];

			caseTotalSum = outCase[2][i];
			controlTotalSum = outControl[2][i];

			caseNum = flib.add(ilib.toSecureFixPoint(aliceCaseNum, flib), ilib.toSecureFixPoint(bobCaseNum, flib));
			controlNum = flib.add(ilib.toSecureFixPoint(aliceControlNum, flib), ilib.toSecureFixPoint(bobControlNum, flib));

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

			T[] degreesOfFreedomTop = flib.add(tLowerFirst, tLowerSecond);
			degreesOfFreedomTop = flib.multiply(degreesOfFreedomTop, degreesOfFreedomTop);

			T[] degreesOfFreedomBottomFirst = flib.div(caseVariance, flib.sub(caseNum, flib.publicValue(1.0)));
			degreesOfFreedomBottomFirst = flib.multiply(degreesOfFreedomBottomFirst, degreesOfFreedomBottomFirst);
			degreesOfFreedomBottomFirst = flib.div(degreesOfFreedomBottomFirst, flib.sub(caseNum, flib.publicValue(1.0)));

			T[] degreesOfFreedomBottomSecond = flib.div(controlVariance, flib.sub(caseNum, flib.publicValue(1.0)));
			degreesOfFreedomBottomSecond = flib.multiply(degreesOfFreedomBottomSecond, degreesOfFreedomBottomSecond);
			degreesOfFreedomBottomSecond = flib.div(degreesOfFreedomBottomSecond, flib.sub(controlNum, flib.publicValue(1.0)));

			T[] degreesOfFreedom = flib.div(degreesOfFreedomTop, flib.add(degreesOfFreedomBottomFirst, degreesOfFreedomBottomSecond));
			res[0][i] = tStat;
			res[1][i] = degreesOfFreedom;
		}
		*/
		//return res;
		
		//System.arraycopy(outControl[2], 0, outCase[3], 0, inputCounters[0].length);

		//System.arraycopy(res[0], 0, outCase[4], 0, inputCounters[0].length);

		//return outCase;
		return inCase;
	}
	
	public static class Generator<T> extends GenRunnable<T> {
		T[][][] inputBobCase;
		T[][][] inputAliceCase;
		T[][][] inputAliceControl;
		T[][][] inputBobControl;
		T[][][] inputCounters;
		T[][][] in;
		
		T[] aliceCaseNum;
		T[] bobCaseNum;
		T[] aliceControlNum;
		T[] bobControlNum;

		@Override
		public void prepareInput(CompEnv<T> gen) throws Exception {
			Options options = new Options();
			options.addOption("s", "case", true, "case");
			options.addOption("t", "control", true, "control");

			CommandLineParser parser = new BasicParser();
			CommandLine cmd = parser.parse(options, args);

			if(!cmd.hasOption("s") || !cmd.hasOption("t")) {
			  throw new Exception("wrong input");
			}
			FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
			T[] l = flib.publicValue(0.0);
			double[][] caseInput = PrepareDataDA.readFile(cmd.getOptionValue("s"));
			double[][] controlInput = PrepareDataDA.readFile(cmd.getOptionValue("t"));

			int EvaCaseDim1 = gen.channel.readInt();
			gen.channel.flush();
			int EvaCaseDim2 = gen.channel.readInt();
			gen.channel.flush();
			int GenCaseDim1 = caseInput.length;
			gen.channel.writeInt(GenCaseDim1);
			gen.channel.flush();
			int GenCaseDim2 = caseInput[0].length;
			gen.channel.writeInt(GenCaseDim2);
			gen.channel.flush();
			int EvaControlDim1 = gen.channel.readInt();
			gen.channel.flush();
			int EvaControlDim2 = gen.channel.readInt();
			gen.channel.flush();
			int GenControlDim1 = controlInput.length;
			gen.channel.writeInt(GenControlDim1);
			int GenControlDim2 = controlInput[0].length;
			gen.channel.writeInt(GenControlDim2);

			inputAliceCase = gen.newTArray(EvaCaseDim1, EvaCaseDim2, 0);
			for(int i = 0; i < EvaCaseDim1; i++){
				for(int j = 0; j < EvaCaseDim2; j++){
					inputCounters[i][j] = gen.inputOfAlice(new boolean[l.length]);
				}
			}
			gen.flush();

			inputBobCase = gen.newTArray(GenCaseDim1, GenCaseDim2, 0);
			for (int i = 0; i < GenCaseDim1; i++) {
				for (int j = 0; j < GenCaseDim2; j++) {
					inputBobCase[i][j] = gen.inputOfBob(Utils.fromFloat(caseInput[i][j], PLength, VLength));
				}
			}
			gen.flush();
			System.out.println("Done with inputBobCase gen");
			
			inputAliceControl = gen.newTArray(EvaControlDim1, EvaControlDim2, 0);
			for(int i = 0; i < EvaControlDim1; i++){
				for(int j = 0; j < EvaControlDim2; j++){
					inputAliceControl[i][j] = gen.inputOfAlice(new boolean[l.length]);
				}
			}
			gen.flush();

			inputBobControl = gen.newTArray(GenControlDim1, GenControlDim2, 0);
			for (int i = 0; i < GenControlDim1; i++) {
				for (int j = 0; j < GenControlDim2; j++) {
					inputBobControl[i][j] = gen.inputOfBob(Utils.fromFloat(controlInput[i][j], PLength, VLength));
				}
			}
			System.out.println("Done with inputBobControl gen");
			gen.flush();
		}
		
		@Override
		public void secureCompute(CompEnv<T> gen) {
			in = compute(gen, inputAliceCase, inputBobCase, inputAliceControl, inputBobControl);
		}
		@Override
		public void prepareOutput(CompEnv<T> gen) {
			for(int i = 0; i < in.length; i++){
				for(int j = 0; j < in[0].length; j++){
					System.out.print(Utils.toFixPoint(gen.outputToAlice(in[i][j]), 20) + " ");
				}
				System.out.println();
			}			
		}
	}
	
	public static class Evaluator<T> extends EvaRunnable<T> {
		T[][][] inputBobCase;
		T[][][] inputAliceCase;
		T[][][] inputCounters;
		T[][][] inputAliceControl;
		T[][][] inputBobControl;
		T[] scResult;
		T[][][] in;
		@Override
		public void prepareInput(CompEnv<T> gen) throws Exception {
			Options options = new Options();
			options.addOption("s", "case", true, "case");
			options.addOption("t", "control", true, "control");

			CommandLineParser parser = new BasicParser();
			CommandLine cmd = parser.parse(options, args);

			if(!cmd.hasOption("s") || !cmd.hasOption("t")) {
			  throw new Exception("wrong input");
			}

			FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
			T[] l = flib.publicValue(0.0);
			double[][] caseInput = PrepareDataDA.readFile(cmd.getOptionValue("s"));
			double[][] controlInput = PrepareDataDA.readFile(cmd.getOptionValue("t"));

			int EvaCaseDim1 = caseInput.length;
			gen.channel.writeInt(EvaCaseDim1);
			gen.channel.flush();
			int EvaCaseDim2 = caseInput[0].length;
			gen.channel.writeInt(EvaCaseDim2);
			gen.channel.flush();			
			int GenCaseDim1 = gen.channel.readInt();
			gen.channel.flush();
			int GenCaseDim2 = gen.channel.readInt();
			gen.channel.flush();
			int EvaControlDim1 = controlInput.length;
			gen.channel.writeInt(EvaControlDim1);
			gen.channel.flush();
			int EvaControlDim2 = controlInput[0].length;
			gen.channel.writeInt(EvaControlDim2);
			gen.channel.flush();
			int GenControlDim1 = gen.channel.readInt();
			gen.channel.flush();
			int GenControlDim2 = gen.channel.readInt();
			gen.channel.flush();

			gen.flush();
			inputAliceCase = gen.newTArray(EvaCaseDim1, EvaCaseDim2, 0);
			for (int i = 0; i < EvaCaseDim1; ++i) {
				for (int j = 0; j < EvaCaseDim2; j++) {
					inputAliceCase[i][j] = gen.inputOfAlice(Utils.fromFloat(caseInput[i][j], PLength, VLength));
				}
			}

			gen.flush();
			System.out.println("Done with inputCounters eva");

			inputBobCase = gen.newTArray(GenCaseDim1, GenCaseDim2, 0);
			for (int i = 0; i < GenCaseDim1; ++i) {
				for (int j = 0; j < GenCaseDim2; j++) {
					inputBobCase[i][j] = gen.inputOfBob(new boolean[l.length]);
				}
			}

			gen.flush();
			System.out.println("Done with inputAliceCase eva");

			inputAliceCase = gen.newTArray(EvaControlDim1, EvaControlDim1, 0);		
			for (int i = 0; i < EvaControlDim1; ++i) {
				for (int j = 0; j < EvaControlDim1; j++) {
					inputAliceCase[i][j] = gen.inputOfAlice(Utils.fromFloat(controlInput[i][j], PLength, VLength));
				}
			}
			gen.flush();

			inputBobControl = gen.newTArray(GenControlDim1, GenControlDim2, 0);
			for (int i = 0; i < 4; ++i) {
				for (int j = 0; j < GenControlDim2; j++) {
					inputBobControl[i][j] = gen.inputOfAlice(new boolean[l.length]);
				}
			}
			gen.flush();

			System.out.println("Done with inputBobControl eva");

		}
		
		@Override
		public void secureCompute(CompEnv<T> gen) {
			in = compute(gen, inputAliceCase, inputBobCase, inputAliceControl, inputBobControl);
		}
		
		@Override
		public void prepareOutput(CompEnv<T> gen) {
			for(int i = 0; i < in.length; i++){
				for(int j =0; j<in[0].length; j++){
					gen.outputToAlice(in[i][j]);
				}
			}
		}
				
	}
	
}