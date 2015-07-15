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

package perFeatureNaive;

import perFeatureNaive.PrepareDataNaive;

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
	static int filterThreshold = 0;
	
	static public<T> T[] compute(CompEnv<T> gen, T[][] inputAliceCase, 
			T[][] inputBobCase, T[][] inputAliceControl, T[][] inputBobControl){
		FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
		IntegerLib<T> ilib = new IntegerLib<T>(gen, 32);
		
		T[] aliceCaseSum = ilib.publicValue(0.0);
		for (int i = 0; i < inputAliceCase.length; i++){
			aliceCaseSum = ilib.add(aliceCaseSum, inputAliceCase[i]);
		}
		
		T[] bobCaseSum = ilib.publicValue(0.0);
		for(int i = 0; i < inputBobCase.length; i++){
			bobCaseSum = ilib.add(bobCaseSum, inputBobCase[i]);
		}
		
		T[] aliceControlSum = ilib.publicValue(0.0);
		for (int i = 0; i < inputAliceControl.length; i++){
			aliceControlSum = ilib.add(aliceControlSum, inputAliceControl[i]);
		}
		
		T[] bobControlSum = ilib.publicValue(0.0);
		for(int i = 0; i < inputBobControl.length; i++){
			bobControlSum = ilib.add(bobControlSum, inputBobControl[i]);
		}
		
		T[] caseNum = ilib.publicValue(inputAliceCase.length + inputBobCase.length);
		T[] controlNum = ilib.publicValue(inputAliceControl.length + inputBobControl.length);
		T[] res = flib.publicValue(0.0);
		
		T[] a = ilib.add(aliceCaseSum, bobCaseSum);
		T[] b = ilib.sub(caseNum, a);
		T[] c = ilib.add(aliceControlSum, bobControlSum);
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
		res = flib.div(upper, lower);
		
		return res;
	}
	
	public static class Generator<T> extends GenRunnable<T> {
		T[][] inputBobCase;
		T[][] inputAliceCase;
		T[][] inputAliceControl;
		T[][] inputBobControl;
		T[][] inputCounters;
		T[][] in;
		int[][] caseInput;
		int[][] controlInput;
		T[] aliceCaseNum;
		T[] bobCaseNum;
		T[] aliceControlNum;
		T[] bobControlNum;
		
		int EvaCaseFeatures;
		int EvaCaseSamples;
		int GenCaseFeatures;
		int GenCaseSamples;
		int EvaControlFeatures;
		int EvaControlSamples;
		int GenControlFeatures;
		int GenControlSamples;

		T[] numAliceCase;
		T[] numBobCase;
		T[] numAliceControl;
		T[] numBobControl;
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
			caseInput = PrepareDataNaive.readFile(cmd.getOptionValue("s"));
			controlInput = PrepareDataNaive.readFile(cmd.getOptionValue("t"));

			EvaCaseFeatures = gen.channel.readInt();
			gen.channel.flush();
			EvaCaseSamples = gen.channel.readInt();
			gen.channel.flush();
			GenCaseFeatures = caseInput.length;
			gen.channel.writeInt(GenCaseFeatures);
			gen.channel.flush();
			GenCaseSamples = caseInput[0].length;
			gen.channel.writeInt(GenCaseSamples);
			gen.channel.flush();
			EvaControlFeatures = gen.channel.readInt();
			gen.channel.flush();
			EvaControlSamples = gen.channel.readInt();
			gen.channel.flush();
			GenControlFeatures = controlInput.length;
			gen.channel.writeInt(GenControlFeatures);
			gen.channel.flush();
			GenControlSamples = controlInput[0].length;
			gen.channel.writeInt(GenControlSamples);
			gen.channel.flush();

			in = gen.newTArray(EvaCaseFeatures, 0);

			System.out.println("Eva Case Dim 1 " + EvaCaseFeatures);
			System.out.println("Eva Case Dim 2 " + EvaCaseSamples);
			System.out.println("Eva Control Dim 1 " + EvaControlFeatures);
			System.out.println("Eva Control Dim 2 " + EvaControlSamples);
			System.out.println("Gen Case Dim 1 " + GenCaseFeatures);
			System.out.println("Gen Case Dim 2 " + GenCaseSamples);
			System.out.println("Gen Control Dim 1 " + GenControlFeatures);
			System.out.println("Gen Control Dim 2 " + GenControlSamples);
		}
		
		@Override
		public void secureCompute(CompEnv<T> gen) {
			
			for(int i =0; i < EvaCaseFeatures; i++){
				int numCase = 0;
				int numControl = 0;
				for(int j = 0; j < caseInput[0].length; j++){
					if(caseInput[i][j] > 0.0){
						numCase++;
					}
				}
				for(int j = 0; j < controlInput[0].length; j++){
					if(controlInput[i][j] > 0.0){
						numControl++;
					}
				}
				/*
				numAliceCase = gen.inputOfAlice(Utils.fromInt(numCase, 32));
				numBobCase = gen.inputOfBob(new boolean[32]);
				numAliceControl = gen.inputOfAlice(Utils.fromInt(numControl, 32));
				numBobControl = gen.inputOfBob(new boolean[32]);
*/
				inputAliceCase = gen.newTArray(GenCaseSamples, 0);
				inputBobCase = gen.newTArray(EvaCaseSamples, 0);
				inputAliceControl = gen.newTArray(GenControlSamples, 0);
				inputBobControl = gen.newTArray(EvaControlSamples, 0);

				for(int j =0; j < GenCaseSamples; j++){
					inputAliceCase[j] = gen.inputOfAlice(Utils.fromInt(caseInput[i][j], 32));
				}
				gen.flush();
				for(int j =0; j < EvaCaseSamples; j++){
					inputBobCase[j] = gen.inputOfBob(new boolean[32]);
				}
				gen.flush();
				for(int j =0; j < GenControlSamples; j++){
					inputAliceControl[j] = gen.inputOfAlice(Utils.fromInt(controlInput[i][j], 32));
				}
				gen.flush();
				for(int j =0; j < EvaControlSamples; j++){
					inputBobControl[j] = gen.inputOfBob(new boolean[32]);
				}
				gen.flush();
				in[i] = compute(gen, inputAliceCase, inputBobCase, inputAliceControl, inputBobControl);
			}
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
		T[][] inputBobCase;
		T[][] inputAliceCase;
		T[][] inputCounters;
		T[][] inputAliceControl;
		T[][] inputBobControl;
		T[] scResult;
		T[][] in;
		
		int[][] caseInput;
		int[][] controlInput;
		
		int EvaCaseFeatures;
		int EvaCaseSamples;
		int GenCaseFeatures;
		int GenCaseSamples;
		int EvaControlFeatures;
		int EvaControlSamples;
		int GenControlFeatures;
		int GenControlSamples;

		T[] numAliceCase;
		T[] numBobCase;
		T[] numAliceControl;
		T[] numBobControl;
		
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

			caseInput = PrepareDataNaive.readFile(cmd.getOptionValue("s"));
			controlInput = PrepareDataNaive.readFile(cmd.getOptionValue("t"));

			EvaCaseFeatures = caseInput.length;
			gen.channel.writeInt(EvaCaseFeatures);
			gen.channel.flush();
			EvaCaseSamples = caseInput[0].length;
			gen.channel.writeInt(EvaCaseSamples);
			gen.channel.flush();			
			GenCaseFeatures = gen.channel.readInt();
			gen.channel.flush();
			GenCaseSamples = gen.channel.readInt();
			gen.channel.flush();
			EvaControlFeatures = controlInput.length;
			gen.channel.writeInt(EvaControlFeatures);
			gen.channel.flush();
			EvaControlSamples = controlInput[0].length;
			gen.channel.writeInt(EvaControlSamples);
			gen.channel.flush();
			GenControlFeatures = gen.channel.readInt();
			gen.channel.flush();
			GenControlSamples = gen.channel.readInt();
			gen.channel.flush();

			in = gen.newTArray(EvaCaseFeatures, 0);

			System.out.println("Eva Case Dim 1 " + EvaCaseFeatures);
			System.out.println("Eva Case Dim 2 " + EvaCaseSamples);
			System.out.println("Eva Control Dim 1 " + EvaControlFeatures);
			System.out.println("Eva Control Dim 2 " + EvaControlSamples);
			System.out.println("Gen Case Dim 1 " + GenCaseFeatures);
			System.out.println("Gen Case Dim 2 " + GenCaseSamples);
			System.out.println("Gen Control Dim 1 " + GenControlFeatures);
			System.out.println("Gen Control Dim 2 " + GenControlSamples);
			
			System.out.println("Done with inputBobControl eva");
		}
		
		@Override
		public void secureCompute(CompEnv<T> gen) {
			
			for(int i =0; i < EvaCaseFeatures; i++){
				int numCase = 0;
				int numControl = 0;
				for(int j = 0; j < caseInput[0].length; j++){
					if(caseInput[i][j] > 0.0){
						numCase++;
					}
				}
				for(int j = 0; j < controlInput[0].length; j++){
					if(controlInput[i][j] > 0.0){
						numControl++;
					}
				}
				/*
				numAliceCase = gen.inputOfAlice(new boolean[32]);
				numBobCase = gen.inputOfBob(Utils.fromInt(numCase, 32));
				numAliceControl = gen.inputOfAlice(new boolean[32]);
				numBobControl = gen.inputOfBob(Utils.fromInt(numControl, 32));
				*/
				inputAliceCase = gen.newTArray(GenCaseSamples, 0);
				inputBobCase = gen.newTArray(EvaCaseSamples, 0);
				inputAliceControl = gen.newTArray(GenControlSamples, 0);
				inputBobControl = gen.newTArray(EvaControlSamples, 0);

				for(int j =0; j < GenCaseSamples; j++){
					inputAliceCase[j] = gen.inputOfAlice(new boolean[32]);
				}
				gen.flush();
				for(int j =0; j < EvaCaseSamples; j++){
					inputBobCase[j] = gen.inputOfBob(Utils.fromInt(caseInput[i][j], 32));
				}
				gen.flush();
				for(int j =0; j < GenControlSamples; j++){
					inputAliceControl[j] = gen.inputOfAlice(new boolean[32]);
				}
				gen.flush();
				for(int j =0; j < EvaControlSamples; j++){
					inputBobControl[j] = gen.inputOfBob(Utils.fromInt(controlInput[i][j], 32));
				}
				gen.flush();

				in[i] = compute(gen, inputAliceCase, inputBobCase, inputAliceControl, inputBobControl);
			}
		}
		
		@Override
		public void prepareOutput(CompEnv<T> gen) {
			for(int j =0; j<in.length; j++){
				gen.outputToAlice(in[j]);
			}
		}
				
	}
	
}
