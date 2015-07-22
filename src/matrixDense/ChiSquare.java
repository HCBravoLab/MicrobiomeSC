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

import naiveF.PrepareDataNaive;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;

import util.EvaRunnable;
import util.GenRunnable;
import util.Utils;
import circuits.arithmetic.FloatLib;
import circuits.arithmetic.IntegerLib;
import flexsc.CompEnv;

public class ChiSquare {
	static int PLength = 54;
	static int VLength = 11;
	static public<T> T[][] compute(CompEnv<T> gen, T[][][] inputAliceCase, 
			T[][][] inputBobCase, T[][][] inputAliceControl, T[][][] inputBobControl){
		T[][][] inCase = gen.newTArray((inputAliceCase.length + inputBobCase.length)/2 , inputAliceCase[0].length + inputBobCase[0].length, 0);
		FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
		IntegerLib<T> ilib = new IntegerLib<T>(gen, 32);
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

		for(int i = 0; i < inCase.length; i++){
			caseSum[i] = ilib.publicValue(0.0);
		}
		
		for(int i = 0; i < inCase.length; i++){
			for (int j = 0; j < inCase[0].length; j++){
				caseSum[i] = ilib.add(caseSum[i], inCase[i][j]);
			}
		}
		
		T[][] controlSum = gen.newTArray(inControl.length,0);
		for(int i = 0; i < inControl.length; i++){
			controlSum[i] = ilib.publicValue(0.0);
		}
		
		for(int i = 0; i < inControl.length; i++){
			for (int j = 0; j < inControl[0].length; j++){
				controlSum[i] = ilib.add(controlSum[i], inControl[i][j]);
			}
		}
		
		T[] caseNum = ilib.publicValue(inCase[0].length);
		T[] controlNum = ilib.publicValue(inControl[0].length);
		T[][] res = gen.newTArray(inCase.length, 0);
		
		for(int i = 0; i < inCase.length; i++){
			T[] a = caseSum[i];
			T[] b = ilib.sub(caseNum, a);
			T[] c = controlSum[i];
			T[] d = ilib.sub(controlNum, c);
			
			T[] fa = ilib.toSecureFloat(a, flib);
			T[] fb = ilib.toSecureFloat(b, flib);
			T[] fc = ilib.toSecureFloat(c, flib);
			T[] fd = ilib.toSecureFloat(d, flib);
			
			T[] upperFirst = flib.add(fa, flib.add(fb, flib.add(fc, fd)));
			T[] upperSecond = flib.sub(flib.multiply(fb, fc), flib.multiply(fa, fd));
			upperSecond = flib.multiply(upperSecond, upperSecond);
			T[] upper = flib.multiply(upperFirst, upperSecond);
			T[] lower = flib.multiply(flib.multiply(flib.add(fa, fb), flib.add(fa, fc)), flib.multiply(flib.add(fb, fd), flib.add(fc, fd)));
			res[i] = flib.div(upper, lower);
		}
		
		return res;
	}
	
	public static class Generator<T> extends GenRunnable<T> {
		T[][][] inputBobCase;
		T[][][] inputAliceCase;
		T[][][] inputAliceControl;
		T[][][] inputBobControl;
		T[][][] inputCounters;
		T[][] in;
		
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
			IntegerLib<T> flib = new IntegerLib<T>(gen, 32);
			T[] l = flib.publicValue(0.0);
			int[][] caseInput = PrepareDataNaive.readFile(cmd.getOptionValue("s"));
			int[][] controlInput = PrepareDataNaive.readFile(cmd.getOptionValue("t"));
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
					inputAliceCase[i][j] = gen.inputOfAlice(Utils.fromInt(caseInput[i][j], 32));
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
					inputAliceControl[i][j] = gen.inputOfAlice(Utils.fromInt(controlInput[i][j], 32));
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
			ChiSquaredDistribution chiDistribution = new ChiSquaredDistribution(1.0);
			System.out.println("chi,p-value");
			for(int i = 0; i < in.length; i++){
				double chi = Utils.toFloat(gen.outputToAlice(in[i]), PLength, VLength);
				if(chi == 0.0){
					System.out.println("NA,NA");
					continue;
				}
				System.out.println(chi + "," + (1-chiDistribution.cumulativeProbability(chi)));
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
		T[][] in;
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

			IntegerLib<T> flib = new IntegerLib<T>(gen, 32);
			T[] l = flib.publicValue(0.0);
			int[][] caseInput = PrepareDataNaive.readFile(cmd.getOptionValue("s"));
			int[][] controlInput = PrepareDataNaive.readFile(cmd.getOptionValue("t"));
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
					inputBobCase[i][j] = gen.inputOfBob(Utils.fromInt(caseInput[i][j], 32));
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
					inputBobControl[i][j] = gen.inputOfBob(Utils.fromInt(controlInput[i][j], 32));
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
			for(int j =0; j<in.length; j++){
				gen.outputToAlice(in[j]);
			}
		}
				
	}
	
}
