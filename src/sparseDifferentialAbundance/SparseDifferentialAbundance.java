package sparseDifferentialAbundance;

import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;

import util.EvaRunnable;
import util.GenRunnable;
import util.Utils;
import circuits.BitonicSortLib;
import circuits.arithmetic.FloatLib;
import circuits.arithmetic.IntegerLib;
import flexsc.CompEnv;

public class SparseDifferentialAbundance {
	static public int VLength = 49;
	static public int PLength = 8;
	static public<T> T[][] compute(CompEnv<T> gen, T[][][] inputCounters, T[][][] inputAliceCase, T[][][] inputAliceControl, 
			T[][][] inputBobCase, T[][][] inputBobControl, T[] numAliceCase, T[] numAliceControl, T[] numBobCase, T[] numBobControl){
		
		FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
		T[] l = flib.publicValue(0.0);
		int d1 = inputAliceCase.length;
		int d2Case = inputCounters[0].length + inputAliceCase[0].length + inputBobCase[0].length;
		int d3 = l.length;
		int d2Control = inputCounters[0].length + inputAliceControl[0].length + inputBobControl[0].length;

		T[][][] inCase = gen.newTArray(d1, d2Case, d3);
		T[][][] inControl = gen.newTArray(d1, d2Control, d3);

		for(int i = 0; i < inputCounters.length; i++){
			System.arraycopy(inputCounters[i], 0, inCase[i], 0, inputCounters[i].length);
			System.arraycopy(inputAliceCase[i], 0, inCase[i], inputCounters[i].length, inputAliceCase[i].length);
			System.arraycopy(inputBobCase[i], 0, inCase[i], inputCounters[i].length + inputAliceCase[i].length, inputBobCase[i].length);
		}

		for(int i = 0; i < inputCounters.length; i++){
			System.arraycopy(inputCounters[i], 0, inControl[i], 0, inputCounters[i].length);
			System.arraycopy(inputAliceControl[i], 0, inControl[i], inputCounters[i].length, inputAliceControl[i].length);
			System.arraycopy(inputBobControl[i], 0, inControl[i], inputCounters[i].length + inputAliceControl[i].length, inputBobControl[i].length);
		}		

		BitonicSortLib<T> blib = new  BitonicSortLib<T>(gen);

		IntegerLib<T> ilib = new IntegerLib<T>(gen);

		int[] attachedRows = {1,2,3};
		//Sort with payload first row and third row as payload
		blib.sortWithPayloadM(inCase[0], inCase, attachedRows, 1, blib.SIGNAL_ONE);		
		
		//make accumulation pass
		// s[i][2] = s[i][1] + s[i+1][1]*s[i+1][2]
		
		for(int i = (inCase[0].length-2); i >= 0; i--){
			inCase[2][i] = flib.add(inCase[2][i], flib.multiply(inCase[1][i+1], inCase[2][i+1]));
			inCase[3][i] = flib.add(flib.multiply(inCase[2][i], inCase[2][i]), flib.multiply(inCase[1][i+1], inCase[3][i+1]));
		}
		int[] attachedRows2 = {0,1,3};
		blib.sortWithPayloadM(inCase[2], inCase, attachedRows2, 0, blib.SIGNAL_ONE);		

		//T[][] in2UnsortedCase = gen.newTArray(in2Case.length, 0);
		//System.arraycopy(in2Case, 0, in2UnsortedCase, 0, in1Case.length);
		
		//Sort with payload first row and second row as payload
		//lib.sortWithPayload(in2Case, in1Case, lib.SIGNAL_ONE);		
		
		//Sort with payload first row and third row as payload
		//lib.sortWithPayload(in2UnsortedCase, in3Case, lib.SIGNAL_ONE);
		
		T[][][] outCase = gen.newTArray(d1, inputCounters[0].length, d3);

		System.arraycopy(inCase[0], 0, outCase[0], 0, inputCounters[0].length);
		System.arraycopy(inCase[1], 0, outCase[1], 0, inputCounters[0].length);
		System.arraycopy(inCase[2], 0, outCase[2], 0, inputCounters[0].length);
		System.arraycopy(inCase[3], 0, outCase[3], 0, inputCounters[0].length);
		
		int[] attachedRows3 = {1,2,3};
		blib.sortWithPayloadM(outCase[0], outCase, attachedRows3, 0, blib.SIGNAL_ONE);		

		int[] attachedRows4 = {1,2,3};
		blib.sortWithPayloadM(inControl[0], inControl, attachedRows4, 1, blib.SIGNAL_ONE);		

		//make accumulation pass
		// s[i][2] = s[i][1] + s[i+1][1]*s[i+1][2]
		
		for(int i = (inControl[0].length-2); i >= 0; i--){
			inControl[2][i] = flib.add(inControl[2][i], flib.multiply(inControl[1][i+1], inControl[2][i+1]));
			inControl[3][i] = flib.add(flib.multiply(inControl[2][i], inControl[2][i]), flib.multiply(inControl[1][i+1], inControl[3][i+1]));
		}


		int[] attachedRows5 = {0,1,3};
		blib.sortWithPayloadM(inControl[2], inControl, attachedRows5, 0, blib.SIGNAL_ONE);		
				
		T[][][] outControl = gen.newTArray(d1, inputCounters[0].length, d3);

		System.arraycopy(inControl[0], 0, outControl[0], 0, inputCounters[0].length);
		System.arraycopy(inControl[1], 0, outControl[1], 0, inputCounters[0].length);
		System.arraycopy(inControl[2], 0, outControl[2], 0, inputCounters[0].length);
		System.arraycopy(inControl[3], 0, outControl[3], 0, inputCounters[0].length);


		int[] attachedRows6 = {1,2,3};
		blib.sortWithPayloadM(outControl[2], outControl, attachedRows6, 0, blib.SIGNAL_ONE);
		
		T[][][] res = gen.newTArray(2, inputCounters[0].length, l.length);

		T[] tStat;

		for(int i = 0; i < inputCounters[0].length; i++){	
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

			caseNum = flib.add(numAliceCase, numBobCase);
			controlNum = flib.add(numAliceControl, numBobControl);

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
		
		return res[0];
		//return in;
		
		//return inCase[0];
		//return outControl;

		//return outCase[2];
		//return res[0];
		//return outControl;
	}
	
	public static class Generator<T> extends GenRunnable<T> {
		T[][][] inputBobCase;
		T[][][] inputBobControl;
		T[][][] inputAliceCase;
		T[][][] inputAliceControl;
		T[][][] inputCounters;
		T[] numBobCase;
		T[] numBobControl;
		T[] numAliceCase;
		T[] numAliceControl;
		T[][] ret;
		int counterRows;

		@Override
		public void prepareInput(CompEnv<T> gen) throws Exception {
			//Take command line input from Generator and parse file
			counterRows = 4;
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

			// Store the number of samples, which will be used to compute b and d of contingency counts table
			numAliceCase = gen.inputOfAlice(Utils.fromFloat(caseInput[0][0], VLength, PLength));
			numAliceControl = gen.inputOfAlice(Utils.fromFloat(controlInput[0][0], VLength, PLength));
			gen.flush();

			numBobCase = gen.inputOfBob(new boolean[l.length]);
			numBobControl = gen.inputOfBob(new boolean[l.length]);
			gen.flush();

			System.out.println("Gen number of samples case: " + caseInput[0][0]);
			System.out.println("Gen number of samples control: " + controlInput[0][0]);

			System.out.println("Gen number of features case: " + caseInput[0][1]);
			System.out.println("Gen number of features control: " + controlInput[0][1]);

			int numCounters = (int)Math.round(caseInput[0][1]);
			System.out.println("Gen number of counts case: " + caseInput[0].length);
			System.out.println("Gen number of counts control: " + controlInput[0].length);
			
			int EvaCaseCounts = gen.channel.readInt();
			gen.channel.flush();
			int GenCaseCounts = caseInput[0].length-2;
			gen.channel.writeInt(GenCaseCounts);
			gen.channel.flush();

			int EvaControlCounts = gen.channel.readInt();
			gen.channel.flush();
			int GenControlCounts = controlInput[0].length-2;
			gen.channel.writeInt(GenControlCounts);
			gen.channel.flush();
			System.out.println("Passed counts nums gen");		
/*
			ArrayList<Float> alCase = new ArrayList<Float>();
			for(int i= 2; i < caseInput[0].length; i++){
				alCase.add(caseInput[0][i]);
			}
			Collections.sort(alCase);
			
			ArrayList<Float> alControl = new ArrayList<Float>();
			for(int i= 2; i < controlInput[0].length; i++){
				alControl.add(controlInput[0][i]);
			}
			Collections.sort(alControl);

			for(int i =2; i <caseInput[0].length;i++){
				caseInput[0][i] = alCase.get(i-2);
			}
			for(int i =2; i <controlInput[0].length;i++){
				controlInput[0][i] = alControl.get(i-2);
			}
			*/
			System.out.println("Passed sorting gen");		

			inputCounters = gen.newTArray(counterRows, numCounters, l.length);
			System.out.println("Num counters gen: " + numCounters);

			for(int i = 0; i < numCounters;i++){
				inputCounters[0][i] = gen.inputOfAlice(Utils.fromFloat(i+1, VLength, PLength));
			}
			for(int i = 0; i < numCounters;i++){
				inputCounters[1][i] = gen.inputOfAlice(Utils.fromFloat(0, VLength, PLength));
			}
			for(int i = 0; i < numCounters;i++){
				inputCounters[2][i] = gen.inputOfAlice(Utils.fromFloat(0, VLength, PLength));
			}
			for(int i = 0; i < numCounters;i++){
				inputCounters[3][i] = gen.inputOfAlice(Utils.fromFloat(0, VLength, PLength));
			}
			gen.flush();
			
			inputAliceCase = gen.newTArray(counterRows, GenCaseCounts, l.length);
			for(int i = 0; i < GenCaseCounts; i++){
				inputAliceCase[0][i] = gen.inputOfAlice(Utils.fromFloat(caseInput[0][i+2], VLength, PLength));
			}
			for(int i = 0; i < GenCaseCounts; i++){
				inputAliceCase[1][i] = gen.inputOfAlice(Utils.fromFloat(1, VLength, PLength));
			}
			for(int i = 0; i < GenCaseCounts; i++){
				inputAliceCase[2][i] = gen.inputOfAlice(Utils.fromFloat(caseInput[1][i+2], VLength, PLength));
			}
			for(int i = 0; i < GenCaseCounts; i++){
				inputAliceCase[3][i] = gen.inputOfAlice(Utils.fromFloat(caseInput[1][i+2], VLength, PLength));
			}
			gen.flush();
			
			inputAliceControl = gen.newTArray(counterRows, GenControlCounts, l.length);
			for(int i = 0; i < GenControlCounts; i++){
				inputAliceControl[0][i] = gen.inputOfAlice(Utils.fromFloat(controlInput[0][i+2], VLength, PLength));
			}
			for(int i = 0; i < GenControlCounts; i++){
				inputAliceControl[1][i] = gen.inputOfAlice(Utils.fromFloat(1, VLength, PLength));
			}
			for(int i = 0; i < GenControlCounts; i++){
				inputAliceControl[2][i] = gen.inputOfAlice(Utils.fromFloat(controlInput[1][i+2], VLength, PLength));
			}
			for(int i = 0; i < GenControlCounts; i++){
				inputAliceControl[3][i] = gen.inputOfAlice(Utils.fromFloat(controlInput[1][i+2], VLength, PLength));
			}
			gen.flush();
			
			inputBobCase = gen.newTArray(counterRows, EvaCaseCounts, l.length);
			for(int i = 0; i < counterRows; i++){
				for(int j = 0; j < EvaCaseCounts; j++){
					inputBobCase[i][j] = gen.inputOfBob(new boolean[l.length]);		
				}
			}
			gen.flush();

			inputBobControl = gen.newTArray(counterRows, EvaControlCounts, l.length);
			for(int i = 0; i < counterRows; i++){
				for(int j = 0; j < EvaControlCounts; j++){
					inputBobControl[i][j] = gen.inputOfBob(new boolean[l.length]);
				}
			}
			gen.flush();
			System.out.println("onto secure compute gen");		

		}
		
		@Override
		public void secureCompute(CompEnv<T> gen) {
			ret = compute(gen, inputCounters, inputAliceCase, inputAliceControl, inputBobCase, inputBobControl,
					 numAliceCase, numAliceControl, numBobCase, numBobControl);
		}
		@Override
		public void prepareOutput(CompEnv<T> gen) {
			//for(int i = 0; i < ret.length; i++){
				for(int j = 0; j < ret.length; j++){
					//System.out.print(Utils.toInt(gen.outputToAlice(ret[j])) + " ");
					System.out.print(Utils.toFloat(gen.outputToAlice(ret[j]), VLength, PLength) + " ");
				}
				System.out.println();
			//}	
		}
	}
	
	public static class Evaluator<T> extends EvaRunnable<T> {
		T[][][] inputBobCase;
		T[][][] inputBobControl;
		T[] scResult;
		T[][][] inputAliceCase;
		T[][][] inputAliceControl;

		T[][][] inputCounters;

		T[] numBobCase;
		T[] numBobControl;
		T[] numAliceCase;
		T[] numAliceControl;
		T[][] ret;
		int counterRows;
		
		@Override
		public void prepareInput(CompEnv<T> gen) throws Exception {
			counterRows = 4;
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

			numAliceCase = gen.inputOfAlice(new boolean[l.length]);
			numAliceControl = gen.inputOfAlice(new boolean[l.length]);
			gen.flush();

			numBobCase = gen.inputOfBob(Utils.fromFloat(caseInput[0][0], VLength, PLength));
			numBobControl = gen.inputOfBob(Utils.fromFloat(controlInput[0][0], VLength, PLength));
			gen.flush();

			System.out.println("Eva number of samples case: " + caseInput[0][0]);
			System.out.println("Eva number of samples control: " + controlInput[0][0]);

			System.out.println("Eva number of features case: " + caseInput[0][1]);
			System.out.println("Eva number of features control: " + controlInput[0][1]);
			int numCounters = (int)Math.round(caseInput[0][1]);

			System.out.println("Eva number of counts case: " + caseInput[0].length);
			System.out.println("Eva number of counts control: " + controlInput[0].length);		

			int EvaCaseCounts = caseInput[0].length-2;
			gen.channel.writeInt(EvaCaseCounts);
			gen.channel.flush();
			int GenCaseCounts = gen.channel.readInt();
			gen.channel.flush();

			int EvaControlCounts = controlInput[0].length-2;
			gen.channel.writeInt(EvaControlCounts);
			gen.channel.flush();
			int GenControlCounts = gen.channel.readInt();
			gen.channel.flush();
			System.out.println("Passed Counts Nums Eva");		
/*
			ArrayList<Float> alCase = new ArrayList<Float>();
			for(int i= 1; i < caseInput[0].length; i++){
				alCase.add(caseInput[0][i]);
			}
			Collections.sort(alCase);
			ArrayList<Float> alControl = new ArrayList<Float>();
			for(int i= 1; i < controlInput[0].length; i++){
				alControl.add(controlInput[0][i]);
			}
			Collections.sort(alControl);

			for(int i =1; i <caseInput.length;i++){
				caseInput[0][i] = alCase.get(i-1);
			}
			for(int i =1; i <controlInput.length;i++){
				controlInput[0][i] = alControl.get(i-1);
			}
			*/
			inputCounters = gen.newTArray(counterRows, numCounters, l.length);
			System.out.println("Num counters eva: " + numCounters);

			for(int i = 0; i < numCounters;i++){
				inputCounters[0][i] = gen.inputOfAlice(new boolean[l.length]);
			}
			for(int i = 0; i < numCounters;i++){
				inputCounters[1][i] = gen.inputOfAlice(new boolean[l.length]);
			}
			for(int i = 0; i < numCounters;i++){
				inputCounters[2][i] = gen.inputOfAlice(new boolean[l.length]);
			}
			for(int i = 0; i < numCounters;i++){
				inputCounters[3][i] = gen.inputOfAlice(new boolean[l.length]);
			}
			gen.flush();
			System.out.println("Passed Counts Inputs eva");		

			inputAliceCase = gen.newTArray(counterRows, GenCaseCounts, l.length);
			for(int i = 0; i < counterRows; i++){
				for(int j = 0; j < GenCaseCounts; j++){
					inputAliceCase[i][j] = gen.inputOfAlice(new boolean[l.length]);		
				}
			}
			gen.flush();

			inputAliceControl = gen.newTArray(counterRows, GenControlCounts, l.length);
			for(int i = 0; i < counterRows; i++){
				for(int j = 0; j < GenControlCounts; j++){
					inputAliceControl[i][j] = gen.inputOfAlice(new boolean[l.length]);
				}
			}
			gen.flush();

			inputBobCase = gen.newTArray(counterRows, EvaCaseCounts, l.length);
			for(int i = 0; i < EvaCaseCounts; i++){
				inputBobCase[0][i] = gen.inputOfBob(Utils.fromFloat(caseInput[0][i+2], VLength, PLength));
			}
			for(int i = 0; i < EvaCaseCounts; i++){
				inputBobCase[1][i] = gen.inputOfBob(Utils.fromFloat(1, VLength, PLength));
			}
			for(int i = 0; i < EvaCaseCounts; i++){
				inputBobCase[2][i] = gen.inputOfBob(Utils.fromFloat(caseInput[1][i+2], VLength, PLength));
			}
			for(int i = 0; i < EvaCaseCounts; i++){
				inputBobCase[3][i] = gen.inputOfBob(Utils.fromFloat(1, VLength, PLength));
			}
			gen.flush();
			
			inputBobControl = gen.newTArray(counterRows, EvaControlCounts, l.length);
			for(int i = 0; i < EvaControlCounts; i++){
				inputBobControl[0][i] = gen.inputOfBob(Utils.fromFloat(controlInput[0][i+2], VLength, PLength));
			}
			for(int i = 0; i < EvaControlCounts; i++){
				inputBobControl[1][i] = gen.inputOfBob(Utils.fromFloat(1, VLength, PLength));
			}
			for(int i = 0; i < EvaControlCounts; i++){
				inputBobControl[2][i] = gen.inputOfBob(Utils.fromFloat(controlInput[1][i+2], VLength, PLength));
			}
			for(int i = 0; i < EvaControlCounts; i++){
				inputBobControl[3][i] = gen.inputOfBob(Utils.fromFloat(1, VLength, PLength));
			}
			gen.flush();
			System.out.println("outto secure compute eva");		

		}
		
		@Override
		public void secureCompute(CompEnv<T> gen) {
			ret = compute(gen, inputCounters, inputAliceCase, inputAliceControl, inputBobCase, inputBobControl,
					 numAliceCase, numAliceControl, numBobCase, numBobControl);
		}
		
		@Override
		public void prepareOutput(CompEnv<T> gen) {
			for(int i = 0; i < ret.length; i++){
			//	for(int j = 0; j < ret[0].length; j++){
					gen.outputToAlice(ret[i]);
				//}
			}
		}
				
	}
	
}

