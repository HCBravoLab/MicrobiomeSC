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
import circuits.arithmetic.IntegerLib;
import flexsc.CompEnv;

public class AlphaDiversity {
	static int PLength = 54;
	static int VLength = 11;
	
	static public<T> T[] computeSimpsons(CompEnv<T> gen, T[][] input){
		FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
		T[] pointOne = flib.publicValue(0.000001);
		T[] one = flib.publicValue(1.0);
		T[] zero = flib.publicValue(0.0);
		T[] simpsons = flib.publicValue(0.0);
		T[] inputAdjusted = flib.publicValue(0.0);
		T[] simpsonsUpper = flib.publicValue(0.0);
		T[] simpsonsLower = flib.publicValue(0.0);
		
		for(int j = 0; j < input.length; j++){
			inputAdjusted = flib.sub(input[j], pointOne);//flib.sub(flib.add(input[j], pointOne), pointOne);
			simpsonsUpper = flib.add(simpsonsUpper, flib.multiply(inputAdjusted, flib.sub(inputAdjusted, one)));
			simpsonsLower = flib.add(simpsonsLower, inputAdjusted);
		}
		simpsonsLower = flib.add(pointOne, simpsonsLower);
		simpsons = flib.div(flib.add(simpsonsUpper,pointOne), flib.multiply(simpsonsLower,flib.sub(simpsonsLower,one)));

		return simpsons;
	}
	
	static public<T> T[][] compute(CompEnv<T> gen, T[][] aliceCaseSimpsons, 
			T[][] bobCaseSimpsons, T[][] aliceControlSimpsons, T[][] bobControlSimpsons){
		FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
		T[] pointOne = flib.publicValue(0.1);
		T[] zero = flib.publicValue(0.0);
		T[] caseSum = flib.publicValue(0.0);
		T[] caseSumOfSquares = flib.publicValue(0.0);
		T[] controlSum = flib.publicValue(0.0);
		T[] controlSumOfSquares = flib.publicValue(0.0);
		
		T[] caseShift = pointOne;
		for(int i = 0; i < aliceCaseSimpsons.length; i++){
			caseSum = flib.add(caseSum, flib.sub(aliceCaseSimpsons[i], caseShift));
			caseSumOfSquares = flib.add(caseSumOfSquares, 
					flib.multiply(flib.sub(aliceCaseSimpsons[i], caseShift), flib.sub(aliceCaseSimpsons[i], caseShift)));
		}
		
		for(int i = 0; i < bobCaseSimpsons.length; i++){
			caseSum = flib.add(caseSum, flib.sub(bobCaseSimpsons[i], caseShift));
			caseSumOfSquares = flib.add(caseSumOfSquares, 
					flib.multiply(flib.sub(bobCaseSimpsons[i], caseShift), flib.sub(bobCaseSimpsons[i], caseShift)));
		}
		T[] controlShift = pointOne;

		for(int i = 0; i < aliceControlSimpsons.length; i++){
			controlSum = flib.add(controlSum, flib.sub(aliceControlSimpsons[i], controlShift));
			T[] addSquare = flib.multiply(flib.sub(aliceControlSimpsons[i], controlShift), flib.sub(aliceControlSimpsons[i], controlShift));
			controlSumOfSquares = flib.add(controlSumOfSquares, addSquare);
		}
		
		for(int i = 0; i < bobControlSimpsons.length; i++){
			controlSum = flib.add(controlSum, flib.sub(bobControlSimpsons[i], controlShift));
			controlSumOfSquares = flib.add(controlSumOfSquares, 
					flib.multiply(flib.sub(bobControlSimpsons[i], controlShift), flib.sub(bobControlSimpsons[i], controlShift)));
		}

		T[] caseNum = flib.publicValue(aliceCaseSimpsons.length + bobCaseSimpsons.length);
		T[] controlNum = flib.publicValue(aliceControlSimpsons.length + bobControlSimpsons.length);
		T[] tStat;
		T[][] res = gen.newTArray(2, 0);

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

		caseTotalSum = caseSum;
		controlTotalSum = controlSum;

		caseMeanAbundance = flib.div(flib.add(caseTotalSum, caseShift), caseNum);
		caseVarianceSecondTerm = flib.div(flib.multiply(caseTotalSum, caseTotalSum), caseNum);
		caseVariance = flib.div(flib.sub(caseSumOfSquares, caseVarianceSecondTerm), caseNum);
		controlMeanAbundance = flib.div(flib.add(controlSum, controlShift), controlNum);		    
		controlVarianceSecondTerm = flib.div(flib.multiply(controlSum, controlTotalSum), controlNum);
		controlVariance = flib.div(flib.sub(controlSumOfSquares, controlVarianceSecondTerm), controlNum);

		tUpper = flib.sub(controlMeanAbundance, caseMeanAbundance);
		tLowerFirst = flib.div(caseVariance, caseNum);
		tLowerSecond = flib.div(controlVariance, controlNum);
		tLowerSqrt = flib.sqrt(flib.add(tLowerFirst, tLowerSecond));
		tStat = flib.div(tUpper, tLowerSqrt);

		T[] degreesOfFreedomTop = flib.add(tLowerFirst, tLowerSecond);
		degreesOfFreedomTop = flib.multiply(flib.add(zero,degreesOfFreedomTop), flib.add(zero,degreesOfFreedomTop));

		T[] degreesOfFreedomBottomFirst = flib.div(caseVariance, flib.sub(caseNum, flib.publicValue(1.0)));
		degreesOfFreedomBottomFirst = flib.multiply(flib.add(zero,degreesOfFreedomBottomFirst), flib.add(zero,degreesOfFreedomBottomFirst));
		degreesOfFreedomBottomFirst = flib.div(degreesOfFreedomBottomFirst, flib.sub(caseNum, flib.publicValue(1.0)));

		T[] degreesOfFreedomBottomSecond = flib.div(controlVariance, flib.sub(caseNum, flib.publicValue(1.0)));
		degreesOfFreedomBottomSecond = flib.multiply(flib.add(zero,degreesOfFreedomBottomSecond), flib.add(zero,degreesOfFreedomBottomSecond));
		degreesOfFreedomBottomSecond = flib.div(degreesOfFreedomBottomSecond, flib.sub(controlNum, flib.publicValue(1.0)));

		T[] degreesOfFreedom = flib.div(degreesOfFreedomTop, flib.add(degreesOfFreedomBottomFirst, degreesOfFreedomBottomSecond));
		res[0] = tStat;
		res[1] = degreesOfFreedom;

		return res;
	}
	
	public static class Generator<T> extends GenRunnable<T> {
		T[][] inputBobCase;
		T[][] inputAliceCase;
		T[][] inputAliceControl;
		T[][] inputBobControl;
		T[][] inputCounters;
		T[][] simpsonsResAliceCase;
		T[][] simpsonsResBobCase;
		T[][] simpsonsResAliceControl;
		T[][] simpsonsResBobControl;

		T[][] alphaDiversityRes;
		
		double[][] caseInput;
		double[][] controlInput;
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
			caseInput = PrepareDataDANaive.readFile(cmd.getOptionValue("s"));
			controlInput = PrepareDataDANaive.readFile(cmd.getOptionValue("t"));

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

			simpsonsResAliceCase  = gen.newTArray(GenCaseFeatures, 0);
			simpsonsResBobCase  = gen.newTArray(EvaCaseFeatures, 0);
			simpsonsResAliceControl = gen.newTArray(GenControlFeatures, 0);
			simpsonsResBobControl = gen.newTArray(EvaControlFeatures, 0);
			
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
			FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
			T[] l = flib.publicValue(0.0);
				inputAliceCase = gen.newTArray(GenCaseSamples, 0);
				inputBobCase = gen.newTArray(EvaCaseSamples, 0);
				inputAliceControl = gen.newTArray(GenControlSamples, 0);
				inputBobControl = gen.newTArray(EvaControlSamples, 0);
				for(int i = 0; i < GenCaseFeatures; i++){
					boolean[][] aliceCase = new boolean[GenCaseSamples][l.length];
					for(int j =0; j < caseInput[i].length; j++){
						aliceCase[j] = Utils.fromFloat(caseInput[i][j], PLength, VLength);
					}
					gen.flush();
					simpsonsResAliceCase[i] = computeSimpsons(gen, gen.inputOfAlice(aliceCase));
				}
				for(int i = 0; i < EvaCaseFeatures; i++){
					boolean[][] bobCase = new boolean[EvaCaseSamples][l.length];
					gen.flush();
					simpsonsResBobCase[i] = computeSimpsons(gen, gen.inputOfBob(bobCase));
				}
				for(int i = 0; i < GenControlFeatures; i++){
					boolean[][] aliceControl = new boolean[GenControlSamples][l.length];
					for(int j =0; j < controlInput[i].length; j++){
						aliceControl[j] = Utils.fromFloat(controlInput[i][j], PLength, VLength);
					}
					gen.flush();
					simpsonsResAliceControl[i] = computeSimpsons(gen, gen.inputOfAlice(aliceControl));
				}
				for(int i = 0; i < EvaControlFeatures; i++){
					boolean[][] bobControl = new boolean[EvaControlSamples][l.length];
					gen.flush();
					simpsonsResBobControl[i] = computeSimpsons(gen, gen.inputOfBob(bobControl));
				}

				alphaDiversityRes = compute(gen, simpsonsResAliceCase, simpsonsResBobCase, 
						simpsonsResAliceControl, simpsonsResBobControl);

		}
		@Override
		public void prepareOutput(CompEnv<T> gen) {
			FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
			double tStat = flib.outputToAlice(alphaDiversityRes[0]);
			double df = flib.outputToAlice(alphaDiversityRes[1]);
			if (tStat == 0.0){
				System.out.println("NA,NA,NA");
				return;
			}
			if (df <= 0.0){
				System.out.println(tStat +",NA,NA");
				return;
			}
			TDistribution tDistribution = new TDistribution(df);
			if(tStat > 0.0)
				System.out.println(tStat + "," + df + "," + (1-tDistribution.cumulativeProbability(tStat))*2.0);
			else
				System.out.println(tStat + "," + df + "," +  tDistribution.cumulativeProbability(tStat)*2.0);
			}			
		
	}
	
	public static class Evaluator<T> extends EvaRunnable<T> {
		T[][] inputBobCase;
		T[][] inputAliceCase;
		T[][] inputCounters;
		T[][] inputAliceControl;
		T[][] inputBobControl;
		T[] scResult;
		T[][] simpsonsResAliceCase;
		T[][] simpsonsResAliceControl;
		T[][] simpsonsResBobCase;
		T[][] simpsonsResBobControl;

		T[][] alphaDiversityRes;
		
		double[][] caseInput;
		double[][] controlInput;
		
		int EvaCaseFeatures;
		int EvaCaseSamples;
		int GenCaseFeatures;
		int GenCaseSamples;
		int EvaControlFeatures;
		int EvaControlSamples;
		int GenControlFeatures;
		int GenControlSamples;
		
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
			caseInput = PrepareDataDANaive.readFile(cmd.getOptionValue("s"));
			controlInput = PrepareDataDANaive.readFile(cmd.getOptionValue("t"));

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

			simpsonsResAliceCase  = gen.newTArray(GenCaseFeatures, 0);
			simpsonsResBobCase  = gen.newTArray(EvaCaseFeatures, 0);
			simpsonsResAliceControl = gen.newTArray(GenControlFeatures, 0);
			simpsonsResBobControl = gen.newTArray(EvaControlFeatures, 0);
			
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
			FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
			T[] l = flib.publicValue(0.0);
				inputAliceCase = gen.newTArray(GenCaseSamples, 0);
				inputBobCase = gen.newTArray(EvaCaseSamples, 0);
				inputAliceControl = gen.newTArray(GenControlSamples, 0);
				inputBobControl = gen.newTArray(EvaControlSamples, 0);
				for(int i =0; i < GenCaseFeatures; i++){
					boolean[][] aliceCase = new boolean[GenCaseSamples][l.length];
					gen.flush();
					simpsonsResAliceCase[i] = computeSimpsons(gen, gen.inputOfAlice(aliceCase));
				}
				for(int i =0; i < EvaCaseFeatures; i++){
					boolean[][]  bobCase = new boolean[EvaCaseSamples][l.length];
					for(int j =0; j < caseInput[i].length; j++){
						bobCase[j] = Utils.fromFloat(caseInput[i][j], PLength, VLength);
					}
					gen.flush();
					simpsonsResBobCase[i] = computeSimpsons(gen, gen.inputOfBob(bobCase));
				}
				for(int i =0; i < GenControlFeatures; i++){
					boolean[][] aliceControl = new boolean[GenControlSamples][l.length];
					gen.flush();
					simpsonsResAliceControl[i] = computeSimpsons(gen, gen.inputOfAlice(aliceControl));
				}
				for(int i =0; i < EvaControlFeatures; i++){
					boolean[][] bobControl = new boolean[EvaControlSamples][l.length];
					for(int j =0; j < controlInput[i].length; j++){
							bobControl[j] = Utils.fromFloat(controlInput[i][j], PLength, VLength);
					}
					gen.flush();
					simpsonsResBobControl[i] = computeSimpsons(gen, gen.inputOfBob(bobControl));
				}

				alphaDiversityRes = compute(gen, simpsonsResAliceCase, simpsonsResBobCase, 
						simpsonsResAliceControl, simpsonsResBobControl);
		}
		
		@Override
		public void prepareOutput(CompEnv<T> gen) {
			FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
			flib.outputToAlice(alphaDiversityRes[0]);
			flib.outputToAlice(alphaDiversityRes[1]);
		}
				
	}
}
