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

package perFeatureSparse;

import perFeatureNaive.PrepareDataNaive;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;

import util.EvaRunnable;
import util.GenRunnable;
import util.Utils;
import circuits.CircuitLib;
import circuits.arithmetic.FloatLib;
import circuits.arithmetic.IntegerLib;
import flexsc.CompEnv;

public class OddsRatio {
	static int PLength = 54;
	static int VLength = 11;
	static public<T> T[] compute(CompEnv<T> gen, T[][] inputAliceCase, 
			T[][] inputBobCase, T[][] inputAliceControl, T[][] inputBobControl, T[] numAliceCase, T[] numBobCase,
			T[] numAliceControl, T[] numBobControl){
		FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
		IntegerLib<T> ilib = new IntegerLib<T>(gen, 32);
		CircuitLib<T> cl = new CircuitLib<T>(gen);

		T[] zero = flib.publicValue(0.0);
		T[] zeroInt = ilib.publicValue(0.0);
		T[] pointFive = flib.publicValue(0.5);

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
		
		T[] caseNum = ilib.add(numAliceCase, numBobCase);
		T[] controlNum = ilib.add(numAliceControl, numBobControl);
		T[] res = flib.publicValue(0.0);
		
		T[] a = ilib.add(aliceCaseSum, bobCaseSum);
		T[] b = ilib.sub(caseNum, a);
		T[] c = ilib.add(aliceControlSum, bobControlSum);
		T[] d = ilib.sub(controlNum, c);
			
		T[] fa = ilib.toSecureFloat(a, flib);
		T faIsZero = flib.eq(fa, zero);
		T[] fb = ilib.toSecureFloat(b, flib);
		T fbIsZero = flib.eq(fb, zero);
		T[] fc = ilib.toSecureFloat(c, flib);
		T fcIsZero = flib.eq(fc, zero);
		T[] fd = ilib.toSecureFloat(d, flib);
		T fdIsZero = flib.eq(fd, zero);
		T aOrDZero = ilib.or(faIsZero, fdIsZero);
		T bOrCZero = ilib.or(fbIsZero, fcIsZero);

		res = flib.div(flib.multiply(fa, fd), flib.multiply(fb, fc));
		res = cl.mux(res, flib.publicValue(Double.MAX_VALUE), bOrCZero);
		//res = fd;
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
		T[] numAliceCase;
		T[] numBobCase;
		T[] numAliceControl;
		T[] numBobControl;
		
		int NumFeatures;
		int[] EvaCaseSamples;
		int[] GenCaseSamples;
		int[] EvaControlSamples;
		int[] GenControlSamples;

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
			caseInput = PrepareDataNaive.readFile(cmd.getOptionValue("s"));
			controlInput = PrepareDataNaive.readFile(cmd.getOptionValue("t"));
			
			NumFeatures = caseInput.length;
			
			in = gen.newTArray(NumFeatures, 0);

			EvaCaseSamples = new int[NumFeatures];
			GenCaseSamples = new int[NumFeatures];
			EvaControlSamples = new int[NumFeatures];
			GenControlSamples = new int[NumFeatures];
			
			for(int i =0; i < NumFeatures; i++){
				int numCase = 0;
				int numControl = 0;
				for(int j = 0; j < caseInput[0].length; j++){
					if(caseInput[i][j] > 0){
						numCase++;
					}
				}
				for(int j = 0; j < controlInput[0].length; j++){
					if(controlInput[i][j] > 0){
						numControl++;
					}
				}
				GenCaseSamples[i] = numCase;
				GenControlSamples[i] = numControl;
			}
			for(int i = 0; i < NumFeatures; i++){
				EvaCaseSamples[i] = gen.channel.readInt();
				gen.channel.flush();
				gen.channel.writeInt(GenCaseSamples[i]);
				gen.channel.flush();
				EvaControlSamples[i] = gen.channel.readInt();
				gen.channel.flush();
				gen.channel.writeInt(GenControlSamples[i]);
				gen.channel.flush();
			}
		}
		
		@Override
		public void secureCompute(CompEnv<T> gen) {

			numAliceCase = gen.inputOfAlice(Utils.fromInt(caseInput[0].length, 32));
			numBobCase = gen.inputOfBob(new boolean[32]);
			numAliceControl = gen.inputOfAlice(Utils.fromInt(controlInput[0].length, 32));
			numBobControl = gen.inputOfBob(new boolean[32]);
			for(int i =0; i < NumFeatures; i++){
				inputAliceCase = gen.newTArray(GenCaseSamples[i], 0);
				inputBobCase = gen.newTArray(EvaCaseSamples[i], 0);
				inputAliceControl = gen.newTArray(GenControlSamples[i], 0);
				inputBobControl = gen.newTArray(EvaControlSamples[i], 0);

				for(int j =0; j < GenCaseSamples[i]; j++){
					inputAliceCase[j] = gen.inputOfAlice(Utils.fromInt(1, 32));
				}
				gen.flush();
				for(int j =0; j < EvaCaseSamples[i]; j++){
					inputBobCase[j] = gen.inputOfBob(new boolean[32]);
				}
				gen.flush();
				for(int j =0; j < GenControlSamples[i]; j++){
					inputAliceControl[j] = gen.inputOfAlice(Utils.fromInt(1, 32));
				}
				gen.flush();
				for(int j =0; j < EvaControlSamples[i]; j++){
					inputBobControl[j] = gen.inputOfBob(new boolean[32]);
				}
				gen.flush();

				in[i] = compute(gen, inputAliceCase, inputBobCase, inputAliceControl, inputBobControl, 
						numAliceCase, numBobCase, numAliceControl, numBobControl);
			}
		}
		@Override
		public void prepareOutput(CompEnv<T> gen) {
			System.out.println("odds ratio");
			for(int i = 0; i < in.length; i++){
				double OR = Utils.toFloat(gen.outputToAlice(in[i]), PLength, VLength);
				System.out.println(OR);
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

		T[] numAliceCase;
		T[] numBobCase;
		T[] numAliceControl;
		T[] numBobControl;
		int[][] caseInput;
		int[][] controlInput;

		int NumFeatures;
		int[] EvaCaseSamples;
		int[] GenCaseSamples;
		int[] EvaControlSamples;
		int[] GenControlSamples;
		
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
			caseInput = PrepareDataNaive.readFile(cmd.getOptionValue("s"));
			controlInput = PrepareDataNaive.readFile(cmd.getOptionValue("t"));

			NumFeatures = caseInput.length;
			in = gen.newTArray(NumFeatures, 0);
			
			EvaCaseSamples = new int[NumFeatures];
			GenCaseSamples = new int[NumFeatures];
			EvaControlSamples = new int[NumFeatures];
			GenControlSamples = new int[NumFeatures];
			
			for(int i =0; i < NumFeatures; i++){
				int numCase = 0;
				int numControl = 0;
				for(int j = 0; j < caseInput[0].length; j++){
					if(caseInput[i][j] > 0){
						numCase++;
					}
				}
				for(int j = 0; j < controlInput[0].length; j++){
					if(controlInput[i][j] > 0){
						numControl++;
					}
				}
				EvaCaseSamples[i] = numCase;
				EvaControlSamples[i] = numControl;
			}
			
			for(int i = 0; i < NumFeatures; i++){		
				gen.channel.writeInt(EvaCaseSamples[i]);
				gen.channel.flush();	
				GenCaseSamples[i] = gen.channel.readInt();
				gen.channel.flush();
				gen.channel.writeInt(EvaControlSamples[i]);
				gen.channel.flush();
				GenControlSamples[i] = gen.channel.readInt();
				gen.channel.flush();
			}
		}
		
		@Override
		public void secureCompute(CompEnv<T> gen) {


			numAliceCase = gen.inputOfAlice(new boolean[32]);
			numBobCase = gen.inputOfBob(Utils.fromInt(caseInput[0].length, 32));
			numAliceControl = gen.inputOfAlice(new boolean[32]);
			numBobControl = gen.inputOfBob(Utils.fromInt(controlInput[0].length, 32));
			for(int i =0; i < NumFeatures; i++){
				inputAliceCase = gen.newTArray(GenCaseSamples[i], 0);
				inputBobCase = gen.newTArray(EvaCaseSamples[i], 0);
				inputAliceControl = gen.newTArray(GenControlSamples[i], 0);
				inputBobControl = gen.newTArray(EvaControlSamples[i], 0);

				for(int j =0; j < GenCaseSamples[i]; j++){
					inputAliceCase[j] = gen.inputOfAlice(new boolean[32]);
				}
				gen.flush();
				for(int j =0; j < EvaCaseSamples[i]; j++){
					inputBobCase[j] = gen.inputOfBob(Utils.fromInt(1, 32));
				}
				gen.flush();
				for(int j =0; j < GenControlSamples[i]; j++){
					inputAliceControl[j] = gen.inputOfAlice(new boolean[32]);
				}
				gen.flush();
				for(int j =0; j < EvaControlSamples[i]; j++){
					inputBobControl[j] = gen.inputOfBob(Utils.fromInt(1, 32));
				}
				gen.flush();
				
				in[i] = compute(gen, inputAliceCase, inputBobCase, inputAliceControl, inputBobControl, 
						numAliceCase, numBobCase, numAliceControl, numBobControl);
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
