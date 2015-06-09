package perFeature;

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
		T[] simpsonsUpper = flib.publicValue(0.0);
		T[] simpsonsLower = flib.publicValue(0.0);

		for(int i = 0; i < input.length ; i++){
				simpsonsUpper = flib.add(simpsonsUpper, 
						flib.multiply(flib.add(zero, input[i]), flib.add(zero,flib.sub(input[i], one))));
				simpsonsLower = flib.add(simpsonsLower, flib.add(zero,input[i]));
		}
		
		simpsons = flib.div(simpsonsUpper, flib.multiply(flib.add(zero,
					simpsonsLower), flib.add(zero, flib.sub(simpsonsLower,one))));

		return simpsons;
	}
	
	static public<T> T[][] compute(CompEnv<T> gen, T[][] aliceCaseSimpsons, 
			T[][] bobCaseSimpsons, T[][] aliceControlSimpsons, T[][] bobControlSimpsons){
		FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
		T[] pointOne = flib.publicValue(0.000001);
		T[] zero = flib.publicValue(0.0);
		T[] caseSum = flib.publicValue(0.0);
		T[] caseSumOfSquares = flib.publicValue(0.0);
		
		for(int i = 0; i < aliceCaseSimpsons.length; i++){
			caseSum = flib.add(caseSum, flib.add(zero,aliceCaseSimpsons[i]));
			caseSumOfSquares = flib.add(caseSumOfSquares, 
					flib.multiply(flib.add(zero, aliceCaseSimpsons[i]), flib.add(zero,aliceCaseSimpsons[i])));
		}
		
		for(int i = 0; i < bobCaseSimpsons.length; i++){
			caseSum = flib.add(caseSum, flib.add(zero,bobCaseSimpsons[i]));
			caseSumOfSquares = flib.add(caseSumOfSquares, 
					flib.multiply(flib.add(zero, bobCaseSimpsons[i]), flib.add(zero,bobCaseSimpsons[i])));
		}
		
		T[] controlSum = flib.publicValue(0.0);
		T[] controlSumOfSquares = flib.publicValue(0.0);
		
		for(int i = 0; i < aliceControlSimpsons.length; i++){
			controlSum = flib.add(controlSum, flib.add(zero,aliceControlSimpsons[i]));
			T[] addSquare = flib.multiply(flib.add(zero, aliceControlSimpsons[i]), flib.add(zero,aliceControlSimpsons[i]));
			controlSumOfSquares = flib.add(controlSumOfSquares, addSquare);
		}
		
		for(int i = 0; i < bobControlSimpsons.length; i++){
			controlSum = flib.add(controlSum, flib.add(zero,bobControlSimpsons[i]));
			controlSumOfSquares = flib.add(controlSumOfSquares, 
					flib.multiply(flib.add(zero, bobControlSimpsons[i]), flib.add(zero,bobControlSimpsons[i])));
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

		caseMeanAbundance = flib.div(caseTotalSum, caseNum);
		caseVarianceSecondTerm = flib.div(flib.multiply(flib.add(pointOne, caseTotalSum), flib.add(pointOne,caseTotalSum)), caseNum);
		caseVariance = flib.div(flib.sub(caseSumOfSquares, caseVarianceSecondTerm), caseNum);
		controlMeanAbundance = flib.div(controlTotalSum, controlNum);		    
		controlVarianceSecondTerm = flib.div(flib.multiply(flib.add(pointOne, controlTotalSum), flib.add(pointOne, controlTotalSum)), controlNum);
		controlVariance = flib.div(flib.sub(controlSumOfSquares, controlVarianceSecondTerm), controlNum);

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
					for(int j =0; j < GenCaseSamples; j++){
						inputAliceCase[j] = gen.inputOfAlice(Utils.fromFloat(caseInput[i][j], PLength, VLength));
					}
					gen.flush();
					simpsonsResAliceCase[i] = computeSimpsons(gen, inputAliceCase);
				}
				for(int i = 0; i < EvaCaseFeatures; i++){
					for(int j =0; j < EvaCaseSamples; j++){
						inputBobCase[j] = gen.inputOfBob(new boolean[l.length]);
					}
					gen.flush();
					simpsonsResBobCase[i] = computeSimpsons(gen, inputBobCase);
				}
				for(int i = 0; i < GenControlFeatures; i++){
					for(int j =0; j < GenControlSamples; j++){
						inputAliceControl[j] = gen.inputOfAlice(Utils.fromFloat(controlInput[i][j], PLength, VLength));
					}
					gen.flush();
					simpsonsResAliceControl[i] = computeSimpsons(gen, inputAliceControl);
				}
				for(int i = 0; i < EvaControlFeatures; i++){
					for(int j =0; j < EvaControlSamples; j++){
						inputBobControl[j] = gen.inputOfBob(new boolean[l.length]);
					}
					gen.flush();
					simpsonsResBobControl[i] = computeSimpsons(gen, inputBobControl);
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
					for(int j =0; j < GenCaseSamples; j++){
						inputAliceCase[j] = gen.inputOfAlice(new boolean[l.length]);
					}
					gen.flush();
					simpsonsResAliceCase[i] = computeSimpsons(gen, inputAliceCase);
				}
				for(int i =0; i < EvaCaseFeatures; i++){
					for(int j =0; j < EvaCaseSamples; j++){
						inputBobCase[j] = gen.inputOfBob(Utils.fromFloat(caseInput[i][j], PLength, VLength));
					}
					gen.flush();
					simpsonsResBobCase[i] = computeSimpsons(gen, inputBobCase);
				}
				for(int i =0; i < GenControlFeatures; i++){
					for(int j =0; j < GenControlSamples; j++){
						inputAliceControl[j] = gen.inputOfAlice(new boolean[l.length]);
					}
					gen.flush();
					simpsonsResAliceControl[i] = computeSimpsons(gen, inputAliceControl);
				}
				for(int i =0; i < EvaControlFeatures; i++){
					for(int j =0; j < EvaControlSamples; j++){
						inputBobControl[j] = gen.inputOfBob(Utils.fromFloat(controlInput[i][j], PLength, VLength));
					}
					gen.flush();
					simpsonsResBobControl[i] = computeSimpsons(gen, inputBobControl);
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