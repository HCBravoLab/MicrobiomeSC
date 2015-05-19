package sparse;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.distribution.TDistribution;

import util.EvaRunnable;
import util.GenRunnable;
import util.Utils;
import circuits.BitonicSortLib;
import circuits.arithmetic.FloatLib;
import circuits.arithmetic.IntegerLib;
import flexsc.CompEnv;

public class DifferentialAbundance {
	static int width = 54;
	static int offset = 11;
	static public<T> T[][][] compute(CompEnv<T> gen, T[][][] inputCounters, T[][][] inputAliceCase, 
			T[][][] inputBobCase, T[][][] inputAliceControl, T[][][] inputBobControl,
			T[] aliceCaseNum, T[] bobCaseNum, T[] aliceControlNum, T[] bobControlNum){//, T[][][] inputAliceControl, T[][][] inputBobControl){
		T[][][] inCase = gen.newTArray(inputAliceCase.length, inputAliceCase[0].length + inputBobCase[0].length + inputCounters[0].length, 0);
		
		for(int i = 0; i < 4; i++){
			System.arraycopy(inputCounters[i], 0, inCase[i], 0, inputCounters[i].length);
			System.arraycopy(inputAliceCase[i], 0, inCase[i], inputCounters[i].length, inputAliceCase[i].length);
			System.arraycopy(inputBobCase[i], 0, inCase[i], inputCounters[i].length + inputAliceCase[i].length, inputBobCase[i].length);
		}

		T[][][] inControl = gen.newTArray(inputAliceControl.length, inputAliceControl[0].length + inputBobControl[0].length + inputCounters[0].length, 0);
		for(int i = 0; i < 4; i++){
			System.arraycopy(inputCounters[i], 0, inControl[i], 0, inputCounters[i].length);
			System.arraycopy(inputAliceControl[i], 0, inControl[i], inputCounters[i].length, inputAliceControl[i].length);
			System.arraycopy(inputBobControl[i], 0, inControl[i], inputCounters[i].length + inputAliceControl[i].length, inputBobControl[i].length);
		}

		BitonicSortLib<T> lib = new BitonicSortLib<T>(gen);
		IntegerLib<T> ilib = new IntegerLib<T>(gen);
		FloatLib<T> flib = new FloatLib<T>(gen, width, offset);
		T [] zero = flib.publicValue(0.0);
		int [] rows = {1,2,3};
		lib.sortWithPayloadM(inCase[0], inCase, rows, lib.SIGNAL_ONE);
		for(int i = (inCase[0].length-2); i >= 0; i--){
			T keyEq = ilib.eq(inCase[0][i], inCase[0][i+1]);
			T dataGreater = ilib.not(ilib.leq(inCase[1][i], inCase[1][i+1]));
			T swap = ilib.and(keyEq, dataGreater);
			for(int k = 0; k < rows.length; k++){
				T[] s = ilib.mux(inCase[rows[k]][i+1], inCase[rows[k]][i], swap);
				s = ilib.xor(s, inCase[rows[k]][i]);
				T[] ki = ilib.xor(inCase[rows[k]][i+1], s);
				T[] kj = ilib.xor(inCase[rows[k]][i], s);
				inCase[rows[k]][i] = ki;
				inCase[rows[k]][i+1] = kj;
			}
		 }
		for(int i = (inCase[0].length-2); i >= 0; i--){
			T[] first = inCase[1][i+1];
			T[] second = inCase[2][i+1];
			T[] result = flib.multiply(flib.add(zero, first), flib.add(zero, second));
			inCase[2][i] = flib.add(result, inCase[2][i]);
			inCase[3][i] = flib.add(flib.multiply(flib.add(zero,inCase[1][i+1]), flib.add(zero,inCase[3][i+1])), 
					flib.multiply(flib.add(zero,inCase[3][i]), flib.add(zero,inCase[3][i])));
		}

		int [] rows2 = {0,2,3};
		lib.sortWithPayloadM(inCase[1], inCase, rows2, lib.SIGNAL_ONE);
		
		T[][][] outCase = gen.newTArray(4, inputCounters[0].length, 0);
		System.arraycopy(inCase[0], 0, outCase[0], 0, inputCounters[0].length);
		System.arraycopy(inCase[1], 0, outCase[1], 0, inputCounters[0].length);
		System.arraycopy(inCase[2], 0, outCase[2], 0, inputCounters[0].length);
		System.arraycopy(inCase[3], 0, outCase[3], 0, inputCounters[0].length);

		lib.sortWithPayloadM(outCase[0], outCase, rows, lib.SIGNAL_ONE);

	    lib.sortWithPayloadM(inControl[0], inControl, rows, lib.SIGNAL_ONE);
		for(int i = (inControl[0].length-2); i >= 0; i--){
			T keyEq = ilib.eq(inControl[0][i], inControl[0][i+1]);
			T dataGreater = ilib.not(ilib.leq(inControl[1][i], inControl[1][i+1]));
			T swap = ilib.and(keyEq, dataGreater);
			for(int k = 0; k < rows.length; k++){
				T[] s = ilib.mux(inControl[rows[k]][i+1], inControl[rows[k]][i], swap);
				s = ilib.xor(s, inControl[rows[k]][i]);
				T[] ki = ilib.xor(inControl[rows[k]][i+1], s);
				T[] kj = ilib.xor(inControl[rows[k]][i], s);
				inControl[rows[k]][i] = ki;
				inControl[rows[k]][i+1] = kj;
			}
		 }
		for(int i = (inControl[0].length-2); i >= 0; i--){
			T[] first = inControl[1][i+1];
			T[] second = inControl[2][i+1];
			T[] result = flib.multiply(flib.add(zero, first), flib.add(zero, second));
			inControl[2][i] = flib.add(result, inControl[2][i]);
			inControl[3][i] = flib.add(flib.multiply(flib.add(zero,inControl[1][i+1]), flib.add(zero,inControl[3][i+1])), 
					flib.multiply(flib.add(zero,inControl[3][i]), flib.add(zero,inControl[3][i])));}
		lib.sortWithPayloadM(inControl[1], inControl, rows2, lib.SIGNAL_ONE);
		
		T[][][] outControl = gen.newTArray(4, inputCounters[0].length, 0);
		System.arraycopy(inControl[0], 0, outControl[0], 0, inputCounters[0].length);
		System.arraycopy(inControl[1], 0, outControl[1], 0, inputCounters[0].length);
		System.arraycopy(inControl[2], 0, outControl[2], 0, inputCounters[0].length);
		System.arraycopy(inControl[3], 0, outControl[3], 0, inputCounters[0].length);

		lib.sortWithPayloadM(outControl[0], outControl, rows, lib.SIGNAL_ONE);
		
		T[][][] res = gen.newTArray(inputCounters[0].length, 2, 0);
		
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

			caseNum = flib.add(ilib.toSecureFloat(aliceCaseNum, flib), ilib.toSecureFloat(bobCaseNum, flib));
			controlNum = flib.add(ilib.toSecureFloat(aliceControlNum, flib), ilib.toSecureFloat(bobControlNum, flib));

			caseMeanAbundance = flib.div(caseTotalSum, caseNum);
			caseVarianceSecondTerm = flib.div(flib.multiply(flib.add(zero, caseTotalSum), flib.add(zero,caseTotalSum)), caseNum);
			caseVariance = flib.div(flib.sub(caseSumOfSquares, caseVarianceSecondTerm), caseNum);
			controlMeanAbundance = flib.div(controlTotalSum, controlNum);		    
			controlVarianceSecondTerm = flib.div(flib.multiply(flib.add(zero, controlTotalSum), flib.add(zero, controlTotalSum)), controlNum);
			controlVariance = flib.div(flib.sub(controlSumOfSquares, controlVarianceSecondTerm), controlNum);

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
			FloatLib<T> flib = new FloatLib<T>(gen, width, offset);
			T[] l = flib.publicValue(0.0);
			double[][] caseInput = PrepareDataDA.readFile(cmd.getOptionValue("s"));
			double[][] controlInput = PrepareDataDA.readFile(cmd.getOptionValue("t"));
			int numCounters = (int)caseInput[0][1];
			System.out.println("num counters: " + numCounters);
			int aliceCaseNumInt = (int)caseInput[0][0];			
			System.out.println("alice case num: " + aliceCaseNumInt);
			int aliceControlNumInt = (int)controlInput[0][0];
			System.out.println("alice control num: " + aliceControlNumInt);

			aliceCaseNum = gen.inputOfAlice(Utils.fromInt(aliceCaseNumInt, 32));
			bobCaseNum = gen.inputOfBob(new boolean[32]);
			aliceControlNum = gen.inputOfAlice(Utils.fromInt(aliceControlNumInt, 32));
			bobControlNum = gen.inputOfBob(new boolean[32]);

			int EvaCaseNum = gen.channel.readInt();
			gen.channel.flush();
			int GenCaseNum = caseInput[0].length-2;
			gen.channel.writeInt(GenCaseNum);
			gen.channel.flush();
			int EvaControlNum = gen.channel.readInt();
			gen.channel.flush();
			int GenControlNum = controlInput[0].length-2;
			gen.channel.writeInt(GenControlNum);

			inputCounters = gen.newTArray(4, numCounters, 0);
			for(int i = 0; i < numCounters; i++){
				inputCounters[0][i] = gen.inputOfAlice(Utils.fromFloat(i+1.0, width, offset));
			}
			for(int i = 1; i < 4; i++){
				for(int j = 0; j < numCounters; j++){
					inputCounters[i][j] = gen.inputOfAlice(Utils.fromFloat(0.0, width, offset));
				}
			}
			
			gen.flush();
			System.out.println("Done with inputCounters gen");
			inputAliceCase = gen.newTArray(4, GenCaseNum, 0);			
			for(int i = 0; i < GenCaseNum; i++){
				inputAliceCase[0][i] = gen.inputOfAlice(Utils.fromFloat(caseInput[0][i+2], width, offset));
			}
			for(int i = 0; i < GenCaseNum; i++){
				inputAliceCase[1][i] = gen.inputOfAlice(Utils.fromFloat(1.0, width, offset));
			}
			for(int i = 0; i < GenCaseNum; i++){
				inputAliceCase[2][i] = gen.inputOfAlice(Utils.fromFloat(caseInput[1][i+2], width, offset));
			}
			for(int i = 0; i < GenCaseNum; i++){
				inputAliceCase[3][i] = gen.inputOfAlice(Utils.fromFloat(caseInput[1][i+2], width, offset));
			}
			
			gen.flush();
			System.out.println("Done with inputAliceCase gen");

			inputBobCase = gen.newTArray(4, EvaCaseNum, 0);
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < EvaCaseNum; j++) {
					inputBobCase[i][j] = gen.inputOfBob(new boolean[l.length]);
				}
			}
			gen.flush();
			System.out.println("Done with inputBobCase gen");

			inputAliceControl = gen.newTArray(4,GenControlNum, 0);
			for(int i = 0; i < GenControlNum; i++){
				inputAliceControl[0][i] = gen.inputOfAlice(Utils.fromFloat(controlInput[0][i+2], width, offset));
			}
			for(int i = 0; i < GenControlNum; i++){
				inputAliceControl[1][i] = gen.inputOfAlice(Utils.fromFloat(1.0, width, offset));
			}
			for(int i = 0; i < GenControlNum; i++){
				inputAliceControl[2][i] = gen.inputOfAlice(Utils.fromFloat(controlInput[1][i+2], width, offset));
			}
			for(int i = 0; i < GenControlNum; i++){
				inputAliceControl[3][i] = gen.inputOfAlice(Utils.fromFloat(controlInput[1][i+2], width, offset));
			}
			
			System.out.println("Done with inputAliceControl gen");
			gen.flush();

			inputBobControl = gen.newTArray(4, EvaControlNum, 0);
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < EvaControlNum; j++) {
					inputBobControl[i][j] = gen.inputOfBob(new boolean[l.length]);
				}
			}
			System.out.println("Done with inputBobControl gen");
			gen.flush();
		}
		
		@Override
		public void secureCompute(CompEnv<T> gen) {
			in = compute(gen, inputCounters, inputAliceCase, inputBobCase, inputAliceControl, inputBobControl,
					aliceCaseNum, bobCaseNum, aliceControlNum, bobControlNum);
		}
		@Override
		public void prepareOutput(CompEnv<T> gen) {
			for(int j = 0; j < in.length; j++){
					double tStat = Utils.toFloat(gen.outputToAlice(in[j][0]), width, offset);
					double df = Utils.toFloat(gen.outputToAlice(in[j][1]), width, offset);
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

			FloatLib<T> flib = new FloatLib<T>(gen, width, offset);
			T[] l = flib.publicValue(0.0);
			double[][] caseInput = PrepareDataDA.readFile(cmd.getOptionValue("s"));
			double[][] controlInput = PrepareDataDA.readFile(cmd.getOptionValue("t"));
			int numCounters = (int)caseInput[0][1];
			int bobCaseNumInt = (int)caseInput[0][0];
			int bobControlNumInt = (int)controlInput[0][0];

			aliceCaseNum = gen.inputOfAlice(new boolean[32]);
			bobCaseNum = gen.inputOfBob(Utils.fromInt(bobCaseNumInt, 32));
			aliceControlNum = gen.inputOfAlice(new boolean[32]);
			bobControlNum = gen.inputOfBob(Utils.fromInt(bobControlNumInt, 32));

			int EvaCaseNum = caseInput[0].length-2;
			gen.channel.writeInt(EvaCaseNum);
			gen.channel.flush();
			int GenCaseNum = gen.channel.readInt();
			gen.channel.flush();
			int EvaControlNum = controlInput[0].length-2;
			gen.channel.writeInt(EvaControlNum);
			gen.channel.flush();
			int GenControlNum = gen.channel.readInt();
			gen.channel.flush();

			System.out.println(EvaCaseNum);
			System.out.println(EvaControlNum);
			gen.flush();
			inputCounters = gen.newTArray(4, numCounters, 0);
			for (int i = 0; i < 4; ++i) {
				for (int j = 0; j < numCounters; j++) {
					inputCounters[i][j] = gen.inputOfAlice(new boolean[l.length]);
				}
			}

			gen.flush();
			System.out.println("Done with inputCounters eva");

			inputAliceCase = gen.newTArray(4, GenCaseNum, 0);
			for (int i = 0; i < 4; ++i) {
				for (int j = 0; j < GenCaseNum; j++) {
					inputAliceCase[i][j] = gen.inputOfAlice(new boolean[l.length]);
				}
			}

			gen.flush();
			System.out.println("Done with inputAliceCase eva");

			inputBobCase = gen.newTArray(4, EvaCaseNum, 0);			
			for(int i = 0; i < EvaCaseNum; i++){
				inputBobCase[0][i] = gen.inputOfBob(Utils.fromFloat(caseInput[0][i+2], width, offset));
			}

			for(int i = 0; i < EvaCaseNum; i++){
				inputBobCase[1][i] = gen.inputOfBob(Utils.fromFloat(1.0, width, offset));
			}
			for(int i = 0; i < EvaCaseNum; i++){
				inputBobCase[2][i] = gen.inputOfBob(Utils.fromFloat(caseInput[1][i+2], width, offset));
			}

			for(int i = 0; i < EvaCaseNum; i++){
				inputBobCase[3][i] = gen.inputOfBob(Utils.fromFloat(caseInput[1][i+2],width, offset));
			}
			gen.flush();

			inputAliceControl = gen.newTArray(4, GenControlNum, 0);
			for (int i = 0; i < 4; ++i) {
				for (int j = 0; j < GenControlNum; j++) {
					inputAliceControl[i][j] = gen.inputOfAlice(new boolean[l.length]);
				}
			}
			gen.flush();
			System.out.println("Done with inputAliceControl eva");

			inputBobControl = gen.newTArray(4, EvaControlNum, 0);
			for(int i = 0; i < EvaControlNum; i++){
				inputBobControl[0][i] = gen.inputOfBob(Utils.fromFloat(controlInput[0][i+2], width, offset));
			}
			for(int i = 0; i < EvaControlNum; i++){
				inputBobControl[1][i] = gen.inputOfBob(Utils.fromFloat(1.0, width, offset));
			}
			for(int i = 0; i < EvaControlNum; i++){
				inputBobControl[2][i] = gen.inputOfBob(Utils.fromFloat(controlInput[1][i+2], width, offset));
			}
			for(int i = 0; i < EvaControlNum; i++){
				inputBobControl[3][i] = gen.inputOfBob(Utils.fromFloat(controlInput[1][i+2], width, offset));
			}
			
			gen.flush();

			System.out.println("Done with inputBobControl eva");

		}
		
		@Override
		public void secureCompute(CompEnv<T> gen) {
			in = compute(gen, inputCounters, inputAliceCase, inputBobCase, inputAliceControl, inputBobControl,
					aliceCaseNum, bobCaseNum, aliceControlNum, bobControlNum);
		}
		
		@Override
		public void prepareOutput(CompEnv<T> gen) {
			for(int j =0; j<in.length; j++){
				for(int i = 0; i < in[0].length; i++)
					gen.outputToAlice(in[j][i]);
			}
		}
				
	}
	
}