package naive;

import naiveF.PrepareDataDANaive;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;

import util.EvaRunnable;
import util.GenRunnable;
import util.Utils;
import circuits.arithmetic.FloatLib;
import flexsc.CompEnv;

public class DifferentialAbundance {
	static int PLength = 54;
	static int VLength = 11;
	static public<T> T[][][] compute(CompEnv<T> gen, T[][][] inputAliceCase, 
			T[][][] inputBobCase, T[][][] inputAliceControl, T[][][] inputBobControl){//, T[][][] inputAliceControl, T[][][] inputBobControl){
		T[][][] inCase = gen.newTArray((inputAliceCase.length + inputBobCase.length)/2 , inputAliceCase[0].length + inputBobCase[0].length, 0);
		FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
		T[] zero = flib.publicValue(0.0);

		for(int i = 0; i < (inputAliceCase.length + inputBobCase.length)/2; i++){
			System.arraycopy(inputAliceCase[i], 0, inCase[i], 0, inputAliceCase[i].length);
			System.arraycopy(inputBobCase[i], 0, inCase[i], inputAliceCase[i].length, inputBobCase[i].length);
		}

		T[][][] inControl = gen.newTArray((inputAliceControl.length + inputBobControl.length)/2, inputAliceControl[0].length + inputBobControl[0].length, 0);
		for(int i = 0; i < (inputAliceControl.length + inputBobControl.length)/2; i++){
			System.arraycopy(inputAliceControl[i], 0, inControl[i], 0, inputAliceControl[i].length);
			System.arraycopy(inputBobControl[i], 0, inControl[i], inputAliceControl[i].length, inputBobControl[i].length);
		}

		T[][] caseSum = gen.newTArray(inCase.length,0);
		T[][] caseSumOfSquares = gen.newTArray(inCase.length,0);

		for(int i = 0; i < inCase.length; i++){
			caseSum[i] = inCase[i][0];
			caseSumOfSquares[i] = inCase[i][0];
		}
		
		for(int i = 0; i < inCase.length; i++){
			for (int j = 0; j < inCase[0].length; j++){
				caseSum[i] = flib.add(caseSum[i], inCase[i][j]);
				caseSumOfSquares[i] = flib.add(caseSumOfSquares[i], 
						flib.multiply(flib.add(zero, inCase[i][j]), flib.add(zero,inCase[i][j])));
			}
		}
		
		T[][] controlSum = gen.newTArray(inControl.length,0);
		T[][] controlSumOfSquares = gen.newTArray(inControl.length,0);

		for(int i = 0; i < inControl.length; i++){
			controlSum[i] = inControl[i][0];
			controlSumOfSquares[i] = inControl[i][0];
		}
		
		for(int i = 0; i < inControl.length; i++){
			for (int j = 0; j < inControl[0].length; j++){
				controlSum[i] = flib.add(controlSum[i], inControl[i][j]);
				controlSumOfSquares[i] = flib.add(controlSumOfSquares[i], 
						flib.multiply(flib.add(zero, inControl[i][j]), flib.add(zero,inControl[i][j])));
			}
		}
		
		T[] caseNum = flib.publicValue(inCase[0].length);
		T[] controlNum = flib.publicValue(inControl[0].length);
		T[] tStat;
		T[][][] res = gen.newTArray(4, inputAliceControl.length,0);
		
		for(int i = 0; i < inCase.length; i++){	
			T[] caseTotalSum;
			T[] controlTotalSum;

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

			caseTotalSum = caseSum[i];
			controlTotalSum = controlSum[i];

			caseMeanAbundance = flib.div(caseTotalSum, caseNum);
			caseVarianceSecondTerm = flib.div(flib.multiply(flib.add(zero, caseTotalSum), flib.add(zero,caseTotalSum)), caseNum);
			caseVariance = flib.div(flib.sub(caseSumOfSquares[i], caseVarianceSecondTerm), caseNum);
			controlMeanAbundance = flib.div(controlTotalSum, controlNum);		    
			controlVarianceSecondTerm = flib.div(flib.multiply(flib.add(zero, controlTotalSum), flib.add(zero, controlTotalSum)), controlNum);
			controlVariance = flib.div(flib.sub(controlSumOfSquares[i], controlVarianceSecondTerm), controlNum);

			tUpper = flib.sub(controlMeanAbundance, caseMeanAbundance);
			tLowerFirst = flib.div(caseVariance, caseNum);
			tLowerSecond = flib.div(controlVariance, controlNum);
			tLowerSqrt = flib.sqrt(flib.add(tLowerFirst, tLowerSecond));
			tStat = flib.div(tUpper, tLowerSqrt);

			T[] degreesOfFreedomTop = flib.add(tLowerFirst, tLowerSecond);
			degreesOfFreedomTop = flib.multiply(flib.add(zero, degreesOfFreedomTop), flib.add(zero, degreesOfFreedomTop));

			T[] degreesOfFreedomBottomFirst = flib.div(caseVariance, flib.sub(caseNum, flib.publicValue(1.0)));
			degreesOfFreedomBottomFirst = flib.multiply(flib.add(zero, degreesOfFreedomBottomFirst), flib.add(zero, degreesOfFreedomBottomFirst));
			degreesOfFreedomBottomFirst = flib.div(degreesOfFreedomBottomFirst, flib.sub(caseNum, flib.publicValue(1.0)));

			T[] degreesOfFreedomBottomSecond = flib.div(controlVariance, flib.sub(caseNum, flib.publicValue(1.0)));
			degreesOfFreedomBottomSecond = flib.multiply(flib.add(zero, degreesOfFreedomBottomSecond), flib.add(zero, degreesOfFreedomBottomSecond));
			degreesOfFreedomBottomSecond = flib.div(degreesOfFreedomBottomSecond, flib.sub(controlNum, flib.publicValue(1.0)));

			T[] degreesOfFreedom = flib.div(degreesOfFreedomTop, flib.add(degreesOfFreedomBottomFirst, degreesOfFreedomBottomSecond));
			res[0][i] = tStat;
			res[1][i] = degreesOfFreedom;
		}
		return res;
		
		//System.arraycopy(outControl[2], 0, outCase[3], 0, inputCounters[0].length);

		//System.arraycopy(res[0], 0, outCase[4], 0, inputCounters[0].length);

		//return outCase;
		//return controlSumOfSquares;
		//return inControl;
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
			double[][] caseInput = PrepareDataDANaive.readFile(cmd.getOptionValue("s"));
			double[][] controlInput = PrepareDataDANaive.readFile(cmd.getOptionValue("t"));
			for(int i = 0; i < caseInput.length; i++){
				for(int j =0; j < caseInput[0].length; j++){
					System.out.print(caseInput[i][j] + " ");
				}
				System.out.println();
			}

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
			gen.channel.flush();
			int GenControlDim2 = controlInput[0].length;
			gen.channel.writeInt(GenControlDim2);
			gen.channel.flush();

			System.out.println("Eva Case Dim 1 " + EvaCaseDim1);
			System.out.println("Eva Case Dim 2 " + EvaCaseDim2);
			System.out.println("Eva Control Dim 1 " + EvaControlDim1);
			System.out.println("Eva Control Dim 2 " + EvaControlDim2);
			System.out.println("Gen Case Dim 1 " + GenCaseDim1);
			System.out.println("Gen Case Dim 2 " + GenCaseDim2);
			System.out.println("Gen Control Dim 1 " + GenControlDim1);
			System.out.println("Gen Control Dim 2 " + GenControlDim2);


			inputAliceCase = gen.newTArray(GenCaseDim1, GenCaseDim2, 0);
			for (int i = 0; i < GenCaseDim1; i++) {
				for (int j = 0; j < GenCaseDim2; j++) {
					inputAliceCase[i][j] = gen.inputOfAlice(Utils.fromFloat(caseInput[i][j], PLength, VLength));
				}
			}
			gen.flush();
			System.out.println("Done with inputBobCase gen");
			
			inputBobCase = gen.newTArray(EvaCaseDim1, EvaCaseDim2, 0);
			for(int i = 0; i < EvaCaseDim1; i++){
				for(int j = 0; j < EvaCaseDim2; j++){
					inputBobCase[i][j] = gen.inputOfBob(new boolean[l.length]);
				}
			}
			gen.flush();
			System.out.println("Done with inputAliceCase gen");
		
			inputAliceControl = gen.newTArray(GenControlDim1, GenControlDim2, 0);
			for (int i = 0; i < GenControlDim1; i++) {
				for (int j = 0; j < GenControlDim2; j++) {
					inputAliceControl[i][j] = gen.inputOfAlice(Utils.fromFloat(controlInput[i][j], PLength, VLength));
				}
			}

			gen.flush();
			System.out.println("Done with inputAliceControl gen");
			inputBobControl = gen.newTArray(EvaControlDim1, EvaControlDim2, 0);
			for(int i = 0; i < EvaControlDim1; i++){
				for(int j = 0; j < EvaControlDim2; j++){
					inputBobControl[i][j] = gen.inputOfBob(new boolean[l.length]);
				}
			}
			gen.flush();
			System.out.println("Done with inputBobControl gen");

		}
		
		@Override
		public void secureCompute(CompEnv<T> gen) {
			in = compute(gen, inputAliceCase, inputBobCase, inputAliceControl, inputBobControl);
		}
		@Override
		public void prepareOutput(CompEnv<T> gen) {
			for(int i = 0; i < in.length; i++){
				for(int j = 0; j < in[0].length; j++){
					System.out.print(Utils.toFloat(gen.outputToAlice(in[0][j]), PLength, VLength) + " ");
					System.out.print(Utils.toFloat(gen.outputToAlice(in[1][j]), PLength, VLength) + " ");
					//System.out.print(Utils.toFloat(gen.outputToAlice(in[j]), PLength, VLength) + " ");
					System.out.println();
				}
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
			double[][] caseInput = PrepareDataDANaive.readFile(cmd.getOptionValue("s"));
			double[][] controlInput = PrepareDataDANaive.readFile(cmd.getOptionValue("t"));
			for(int i = 0; i < caseInput.length; i++){
				for(int j =0; j < caseInput[0].length; j++){
					System.out.print(caseInput[i][j] + " ");
				}
				System.out.println();
			}
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


			System.out.println("Eva Case Dim 1 " + EvaCaseDim1);
			System.out.println("Eva Case Dim 2 " + EvaCaseDim2);
			System.out.println("Eva Control Dim 1 " + EvaControlDim1);
			System.out.println("Eva Control Dim 2 " + EvaControlDim2);
			System.out.println("Gen Case Dim 1 " + GenCaseDim1);
			System.out.println("Gen Case Dim 2 " + GenCaseDim2);
			System.out.println("Gen Control Dim 1 " + GenControlDim1);
			System.out.println("Gen Control Dim 2 " + GenControlDim2);
			


			inputAliceCase = gen.newTArray(GenCaseDim1, GenCaseDim2, 0);
			for (int i = 0; i < GenCaseDim1; i++) {
				for (int j = 0; j < GenCaseDim2; j++) {
					inputAliceCase[i][j] = gen.inputOfAlice(new boolean[l.length]);
				}
			}


			gen.flush();
			System.out.println("Done with inputAliceCase eva");
			inputBobCase = gen.newTArray(EvaCaseDim1, EvaCaseDim2, 0);
			for (int i = 0; i < EvaCaseDim1; i++) {
				for (int j = 0; j < EvaCaseDim2; j++) {
					inputBobCase[i][j] = gen.inputOfBob(Utils.fromFloat(caseInput[i][j], PLength, VLength));
				}
			}
			
			gen.flush();
			System.out.println("Done with inputBobCase eva");

			inputAliceControl = gen.newTArray(GenControlDim1, GenControlDim2, 0);
			for (int i = 0; i < GenControlDim1; i++) {
				for (int j = 0; j < GenControlDim2; j++) {
					inputAliceControl[i][j] = gen.inputOfAlice(new boolean[l.length]);
				}
			}

			gen.flush();
			System.out.println("Done with inputAliceControl eva");
			
			inputBobControl = gen.newTArray(EvaControlDim1, EvaControlDim2, 0);		
			for (int i = 0; i < EvaControlDim1; i++) {
				for (int j = 0; j < EvaControlDim2; j++) {
					inputBobControl[i][j] = gen.inputOfBob(Utils.fromFloat(controlInput[i][j], PLength, VLength));
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