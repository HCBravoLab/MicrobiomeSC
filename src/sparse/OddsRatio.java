package sparse;

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

public class OddsRatio {
	static int width = 54;
	static int offset = 11;
	static public<T> T[][] compute(CompEnv<T> gen, T[][][] inputCounters, T[][][] inputAliceCase, 
			T[][][] inputBobCase, T[][][] inputAliceControl, T[][][] inputBobControl,
			T[] aliceCaseNum, T[] bobCaseNum, T[] aliceControlNum, T[] bobControlNum){//, T[][][] inputAliceControl, T[][][] inputBobControl){
		T[][][] inCase = gen.newTArray(inputAliceCase.length, inputAliceCase[0].length + inputBobCase[0].length + inputCounters[0].length, 0);
		
		for(int i = 0; i < 3; i++){
			System.arraycopy(inputCounters[i], 0, inCase[i], 0, inputCounters[i].length);
			System.arraycopy(inputAliceCase[i], 0, inCase[i], inputCounters[i].length, inputAliceCase[i].length);
			System.arraycopy(inputBobCase[i], 0, inCase[i], inputCounters[i].length + inputAliceCase[i].length, inputBobCase[i].length);
		}

		T[][][] inControl = gen.newTArray(inputAliceControl.length, inputAliceControl[0].length + inputBobControl[0].length + inputCounters[0].length, 0);
		for(int i = 0; i < 3; i++){
			System.arraycopy(inputCounters[i], 0, inControl[i], 0, inputCounters[i].length);
			System.arraycopy(inputAliceControl[i], 0, inControl[i], inputCounters[i].length, inputAliceControl[i].length);
			System.arraycopy(inputBobControl[i], 0, inControl[i], inputCounters[i].length + inputAliceControl[i].length, inputBobControl[i].length);
		}
		
		BitonicSortLib<T> lib = new BitonicSortLib<T>(gen);
		IntegerLib<T> ilib = new IntegerLib<T>(gen);
		int [] rows = {1,2};
		lib.sortWithPayloadM(inCase[0], inCase, rows, lib.SIGNAL_ONE);
		for(int i = (inCase[0].length-2); i >= 0; i--){
			T keyEq = ilib.eq(inCase[0][i], inCase[0][i+1]);
			T dataGreater = ilib.not(ilib.leq(inCase[1][i], inCase[1][i+1]));
			T swap = ilib.and(keyEq, dataGreater);
			//T swap = eq(greater, dir);
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
			inCase[2][i] = ilib.add(inCase[1][i], ilib.multiply(inCase[1][i+1], inCase[2][i+1]));
		}

		int [] rows2 = {0,2};
		lib.sortWithPayloadM(inCase[1], inCase, rows2, lib.SIGNAL_ONE);

		
		T[][][] outCase = gen.newTArray(3, inputCounters[0].length, 0);
		System.arraycopy(inCase[0], 0, outCase[0], 0, inputCounters[0].length);
		System.arraycopy(inCase[1], 0, outCase[1], 0, inputCounters[0].length);
		System.arraycopy(inCase[2], 0, outCase[2], 0, inputCounters[0].length);
		lib.sortWithPayloadM(outCase[0], outCase, rows, lib.SIGNAL_ONE);

	    lib.sortWithPayloadM(inControl[0], inControl, rows, lib.SIGNAL_ONE);
		for(int i = (inControl[0].length-2); i >= 0; i--){
			T keyEq = ilib.eq(inControl[0][i], inControl[0][i+1]);
			T dataGreater = ilib.not(ilib.leq(inControl[1][i], inControl[1][i+1]));
			T swap = ilib.and(keyEq, dataGreater);
			//T swap = eq(greater, dir);
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
			inControl[2][i] = ilib.add(inControl[1][i], ilib.multiply(inControl[1][i+1], inControl[2][i+1]));
		}
		lib.sortWithPayloadM(inControl[1], inControl, rows2, lib.SIGNAL_ONE);
		
		T[][][] outControl = gen.newTArray(3, inputCounters[0].length, 0);
		System.arraycopy(inControl[0], 0, outControl[0], 0, inputCounters[0].length);
		System.arraycopy(inControl[1], 0, outControl[1], 0, inputCounters[0].length);
		System.arraycopy(inControl[2], 0, outControl[2], 0, inputCounters[0].length);
		lib.sortWithPayloadM(outControl[0], outControl, rows, lib.SIGNAL_ONE);
		
		T[][] res = gen.newTArray(inputCounters[0].length, 0);
		FloatLib<T> flib = new FloatLib<T>(gen, width, offset);
		
		for(int i = 0; i < inputCounters[0].length; i++){
			T[] a = outCase[2][i];
			T[] b = ilib.sub(ilib.add(aliceCaseNum, bobCaseNum), a);
			T[] c = outControl[2][i];
			T[] d = ilib.sub(ilib.add(aliceControlNum, bobControlNum), c);

			T[] fa = ilib.toSecureFloat(a, flib);
			T[] fb = ilib.toSecureFloat(b, flib);
			T[] fc = ilib.toSecureFloat(c, flib);
			T[] fd = ilib.toSecureFloat(d, flib);

			res[i] = flib.div(flib.multiply(fa, fd), flib.multiply(fb, fc));
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
			
			int[] caseInput = PrepareData.readFile(cmd.getOptionValue("s"));
			int[] controlInput = PrepareData.readFile(cmd.getOptionValue("t"));
			int numCounters = caseInput[1];
			int aliceCaseNumInt = caseInput[0];
			int aliceControlNumInt = controlInput[0];
			aliceCaseNum = gen.inputOfAlice(Utils.fromInt(aliceCaseNumInt, 32));
			bobCaseNum = gen.inputOfBob(new boolean[32]);
			aliceControlNum = gen.inputOfAlice(Utils.fromInt(aliceControlNumInt, 32));
			bobControlNum = gen.inputOfBob(new boolean[32]);

			int EvaCaseNum = gen.channel.readInt();
			gen.channel.flush();
			int GenCaseNum = caseInput.length-2;
			gen.channel.writeInt(GenCaseNum);
			gen.channel.flush();
			int EvaControlNum = gen.channel.readInt();
			gen.channel.flush();
			int GenControlNum = controlInput.length-2;
			gen.channel.writeInt(GenControlNum);

			inputCounters = gen.newTArray(3, numCounters, 0);
			for(int i = 0; i < numCounters; i++){
				inputCounters[0][i] = gen.inputOfAlice(Utils.fromInt(i+1, 32));
			}
			for(int i = 1; i < 3; i++){
				for(int j = 0; j < numCounters; j++){
					inputCounters[i][j] = gen.inputOfAlice(Utils.fromInt(0, 32));
				}
			}
			gen.flush();
			System.out.println("Done with inputCounters gen");
			
			inputAliceCase = gen.newTArray(3, GenCaseNum, 0);			
			for(int i = 0; i < GenCaseNum; i++){
				inputAliceCase[0][i] = gen.inputOfAlice(Utils.fromInt(caseInput[i+2], 32));
			}
			for(int i = 1; i < 3; i++){
				for(int j = 0; j < GenCaseNum; j++){
					inputAliceCase[i][j] = gen.inputOfAlice(Utils.fromInt(1, 32));
				}
			}
			gen.flush();
			System.out.println("Done with inputAliceCase gen");

			inputBobCase = gen.newTArray(3, EvaCaseNum, 0);
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < EvaCaseNum; j++) {
					inputBobCase[i][j] = gen.inputOfBob(new boolean[32]);
				}
			}
			gen.flush();
			System.out.println("Done with inputBobCase gen");

			inputAliceControl = gen.newTArray(3,GenControlNum, 0);
			for(int i = 0; i < GenControlNum; i++){
				inputAliceControl[0][i] = gen.inputOfAlice(Utils.fromInt(controlInput[i+2], 32));
			}
			for(int i = 1; i < 3; i++){
				for(int j = 0; j < GenControlNum; j++){
					inputAliceControl[i][j] = gen.inputOfAlice(Utils.fromInt(1, 32));
				}
			}
			System.out.println("Done with inputAliceControl gen");
			gen.flush();

			inputBobControl = gen.newTArray(3, EvaControlNum, 0);
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < EvaControlNum; j++) {
					inputBobControl[i][j] = gen.inputOfBob(new boolean[32]);
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
					System.out.println(Utils.toFloat(gen.outputToAlice(in[j]), width, offset) + " ");
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
		T[] aliceCaseNum;
		T[] bobCaseNum;
		T[] aliceControlNum;
		T[] bobControlNum;
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
			
			int[] caseInput = PrepareData.readFile(cmd.getOptionValue("s"));
			int[] controlInput = PrepareData.readFile(cmd.getOptionValue("t"));
			int numCounters = caseInput[1];
			int bobCaseNumInt = caseInput[0];
			int bobControlNumInt = controlInput[0];

			aliceCaseNum = gen.inputOfAlice(new boolean[32]);
			bobCaseNum = gen.inputOfBob(Utils.fromInt(bobCaseNumInt, 32));
			aliceControlNum = gen.inputOfAlice(new boolean[32]);
			bobControlNum = gen.inputOfBob(Utils.fromInt(bobControlNumInt, 32));

			int EvaCaseNum = caseInput.length-2;
			gen.channel.writeInt(EvaCaseNum);
			gen.channel.flush();
			int GenCaseNum = gen.channel.readInt();
			gen.channel.flush();
			int EvaControlNum = controlInput.length-2;
			gen.channel.writeInt(EvaControlNum);
			gen.channel.flush();
			int GenControlNum = gen.channel.readInt();
			gen.channel.flush();

			inputCounters = gen.newTArray(3, numCounters, 0);
			for (int i = 0; i < 3; ++i) {
				for (int j = 0; j < numCounters; j++) {
					inputCounters[i][j] = gen.inputOfAlice(new boolean[32]);
				}
			}
			gen.flush();
			System.out.println("Done with inputCounters eva");

			inputAliceCase = gen.newTArray(3, GenCaseNum, 0);
			for (int i = 0; i < 3; ++i) {
				for (int j = 0; j < GenCaseNum; j++) {
					inputAliceCase[i][j] = gen.inputOfAlice(new boolean[32]);
				}
			}
			gen.flush();
			System.out.println("Done with inputAliceCase eva");

			inputBobCase = gen.newTArray(3, EvaCaseNum, 0);			
			for(int i = 0; i < EvaCaseNum; i++){
				inputBobCase[0][i] = gen.inputOfBob(Utils.fromInt(caseInput[i+2], 32));
			}
			for(int i = 1; i < 3; i++){
				for(int j = 0; j < EvaCaseNum; j++){
					inputBobCase[i][j] = gen.inputOfBob(Utils.fromInt(1, 32));
				}
			}
			gen.flush();
			System.out.println("Done with inputBobCase eva");

			inputAliceControl = gen.newTArray(3, GenControlNum, 0);
			for (int i = 0; i < 3; ++i) {
				for (int j = 0; j < GenControlNum; j++) {
					inputAliceControl[i][j] = gen.inputOfAlice(new boolean[32]);
				}
			}
			gen.flush();
			System.out.println("Done with inputAliceControl eva");

			inputBobControl = gen.newTArray(3, EvaControlNum, 0);
			for(int i = 0; i < EvaControlNum; i++){
				inputBobControl[0][i] = gen.inputOfBob(Utils.fromInt(controlInput[i+2], 32));
			}
			for(int i = 1; i < 3; i++){
				for(int j = 0; j < EvaControlNum; j++){
					inputBobControl[i][j] = gen.inputOfBob(Utils.fromInt(1, 32));
				}
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
					gen.outputToAlice(in[j]);
				}
		}	
	}
}
