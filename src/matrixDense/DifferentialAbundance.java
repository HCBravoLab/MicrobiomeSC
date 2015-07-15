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

package matrixDense;

import naiveF.PrepareDataDANaive;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.distribution.TDistribution;

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
		T [] pointOne = flib.publicValue(0.000001);
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
						flib.multiply(flib.add(pointOne, inCase[i][j]), flib.add(pointOne,inCase[i][j])));
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
						flib.multiply(flib.add(pointOne, inControl[i][j]), flib.add(pointOne,inControl[i][j])));
			}
		}
		
		T[] caseNum = flib.publicValue(inCase[0].length);
		T[] controlNum = flib.publicValue(inControl[0].length);
		T[] tStat;
		T[][][] res = gen.newTArray(inCase.length,2,0);
		
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
			caseVarianceSecondTerm = flib.div(flib.multiply(flib.add(pointOne, caseTotalSum), flib.add(pointOne,caseTotalSum)), caseNum);
			caseVariance = flib.div(flib.sub(caseSumOfSquares[i], caseVarianceSecondTerm), caseNum);
			controlMeanAbundance = flib.div(controlTotalSum, controlNum);		    
			controlVarianceSecondTerm = flib.div(flib.multiply(flib.add(pointOne, controlTotalSum), flib.add(pointOne, controlTotalSum)), controlNum);
			controlVariance = flib.div(flib.sub(controlSumOfSquares[i], controlVarianceSecondTerm), controlNum);

			tUpper = flib.sub(controlMeanAbundance, caseMeanAbundance);
			tLowerFirst = flib.div(caseVariance, caseNum);
			tLowerSecond = flib.div(controlVariance, controlNum);
			tLowerSqrt = flib.sqrt(flib.add(tLowerFirst, tLowerSecond));
			tStat = flib.div(tUpper, tLowerSqrt);

			T[] degreesOfFreedomTop = flib.add(tLowerFirst, tLowerSecond);
			degreesOfFreedomTop = flib.multiply(flib.add(pointOne, degreesOfFreedomTop), flib.add(pointOne, degreesOfFreedomTop));

			T[] degreesOfFreedomBottomFirst = flib.div(caseVariance, flib.sub(caseNum, flib.publicValue(1.0)));
			degreesOfFreedomBottomFirst = flib.multiply(flib.add(pointOne, degreesOfFreedomBottomFirst), flib.add(pointOne, degreesOfFreedomBottomFirst));
			degreesOfFreedomBottomFirst = flib.div(degreesOfFreedomBottomFirst, flib.sub(caseNum, flib.publicValue(1.0)));

			T[] degreesOfFreedomBottomSecond = flib.div(controlVariance, flib.sub(caseNum, flib.publicValue(1.0)));
			degreesOfFreedomBottomSecond = flib.multiply(flib.add(pointOne, degreesOfFreedomBottomSecond), flib.add(pointOne, degreesOfFreedomBottomSecond));
			degreesOfFreedomBottomSecond = flib.div(degreesOfFreedomBottomSecond, flib.sub(controlNum, flib.publicValue(1.0)));

			T[] degreesOfFreedom = flib.div(degreesOfFreedomTop, flib.add(degreesOfFreedomBottomFirst, degreesOfFreedomBottomSecond));
			res[i][0] = tStat;
			res[i][1] = degreesOfFreedom;
		}
		return res;
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
			FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
			for(int j = 0; j < in.length; j++){
				double tStat = flib.outputToAlice(in[j][0]);
				double df = flib.outputToAlice(in[j][1]);
				if (tStat == 0.0){
					System.out.println("NA,NA,NA");
					continue;
				}
				if (df <= 0.0){
					System.out.println(tStat +",NA,NA");
					continue;
				}
				TDistribution tDistribution = new TDistribution(df);
				if(tStat > 0.0)
					System.out.println(tStat + "," + df + "," + (1-tDistribution.cumulativeProbability(tStat))*2.0);
				else
					System.out.println(tStat + "," + df + "," +  tDistribution.cumulativeProbability(tStat)*2.0);

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
			FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
			for(int i = 0; i < in.length; i++){
				flib.outputToAlice(in[i][0]);
				flib.outputToAlice(in[i][1]);
			}
		}
				
	}
	
}
