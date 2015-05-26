package sparseOptimized;

import java.util.Arrays;
import java.util.Comparator;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;

import util.EvaRunnable;
import util.GenRunnable;
import util.Utils;
import circuits.BitonicSortLib;
import circuits.arithmetic.FloatLib;
import circuits.arithmetic.IntegerLib;
import flexsc.CompEnv;

public class OddsRatio {
	static int width= 54;
	static int offset = 11;
	static public<T> T[][] compute(CompEnv<T> gen, T[][][] aliceCase, 
			T[][][] bobCase, T[][][] aliceControl, T[][][] bobControl, T[][] inputCounters,
			T[] aliceCaseNum, T[] bobCaseNum, T[] aliceControlNum, T[] bobControlNum){//, T[][][] inputAliceControl, T[][][] inputBobControl){
		
		BitonicSortLib<T> lib = new BitonicSortLib<T>(gen);
		IntegerLib<T> ilib = new IntegerLib<T>(gen);
		FloatLib<T> flib = new FloatLib<T>(gen, width, offset);

		T[] zero = ilib.publicValue(0.0);
		T[] holder =  ilib.publicValue(0.0);
		T[] holder0 =  ilib.publicValue(0.0);
		
		for(int i = (aliceCase[0].length-2); i >= 0; i--){
			holder0 = aliceCase[0][i+1];
			T rowIdsEqual = ilib.not(ilib.eq(aliceCase[0][i], holder0));
			T[] addAndIncrement = ilib.add(aliceCase[1][i], aliceCase[1][i+1]);
			aliceCase[1][i] = ilib.mux(addAndIncrement, aliceCase[1][i], rowIdsEqual);
			aliceCase[1][i+1] = ilib.mux(zero, holder, rowIdsEqual);
			aliceCase[0][i+1] = ilib.mux(zero, holder0, rowIdsEqual);
			holder = aliceCase[1][i];
		}
		
		lib.sortWithPayload(aliceCase[0], aliceCase[1], lib.SIGNAL_ZERO);
		
		holder =  ilib.publicValue(0.0);
		holder0 =  ilib.publicValue(0.0);
		
		
		for(int i = (bobCase[0].length-2); i >= 0; i--){
			holder0 = bobCase[0][i+1];
			T rowIdsEqual = ilib.not(ilib.eq(bobCase[0][i], holder0));
			T[] addAndIncrement = ilib.add(bobCase[1][i], bobCase[1][i+1]);
			bobCase[1][i] = ilib.mux(addAndIncrement, bobCase[1][i], rowIdsEqual);
			bobCase[1][i+1] = ilib.mux(zero, holder, rowIdsEqual);
			bobCase[0][i+1] = ilib.mux(zero, holder0, rowIdsEqual);
			holder = bobCase[1][i];
		}
		
		lib.sortWithPayload(bobCase[0], bobCase[1], lib.SIGNAL_ZERO);

		holder =  ilib.publicValue(0.0);
		holder0 =  ilib.publicValue(0.0);

		for(int i = (bobControl[0].length-2); i >= 0; i--){
			holder0 = bobControl[0][i+1];
			T rowIdsEqual = ilib.not(ilib.eq(bobControl[0][i], holder0));
			T[] addAndIncrement = ilib.add(bobControl[1][i], bobControl[1][i+1]);
			bobControl[1][i] = ilib.mux(addAndIncrement, bobControl[1][i], rowIdsEqual);
			bobControl[1][i+1] = ilib.mux(zero, holder, rowIdsEqual);
			bobControl[0][i+1] = ilib.mux(zero, holder0, rowIdsEqual);
			holder = bobControl[1][i];
		}
		lib.sortWithPayload(bobControl[0], bobControl[1], lib.SIGNAL_ZERO);

		holder =  ilib.publicValue(0.0);
		holder0 =  ilib.publicValue(0.0);

		for(int i = (aliceControl[0].length-2); i >= 0; i--){
			holder0 = aliceControl[0][i+1];
			T rowIdsEqual = ilib.not(ilib.eq(aliceControl[0][i], holder0));
			T[] addAndIncrement = ilib.add(aliceControl[1][i], aliceControl[1][i+1]);
			aliceControl[1][i] = ilib.mux(addAndIncrement, aliceControl[1][i], rowIdsEqual);
			aliceControl[1][i+1] = ilib.mux(zero, holder, rowIdsEqual);
			aliceControl[0][i+1] = ilib.mux(zero, holder0, rowIdsEqual);
			holder = aliceControl[1][i];
		}
		lib.sortWithPayload(aliceControl[0], aliceControl[1], lib.SIGNAL_ZERO);
		
		T[][] res = gen.newTArray(inputCounters.length, 0);

		for(int i = 0; i < inputCounters.length; i++){
			T[] a = ilib.add(aliceCase[1][i],bobCase[1][i]);
			T[] b = ilib.sub(ilib.add(aliceCaseNum,bobCaseNum), a);
			T[] c = ilib.add(aliceControl[1][i], bobControl[1][i]);
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
		T[][] inputCounters;
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
			System.out.println("Alice case num " + aliceCaseNumInt);
			System.out.println("Alice control num " + aliceControlNumInt);

			int EvaCaseNum = gen.channel.readInt();
			gen.channel.flush();
			int GenCaseNum = caseInput.length-2;
			gen.channel.writeInt(GenCaseNum);
			gen.channel.flush();
			int EvaControlNum = gen.channel.readInt();
			gen.channel.flush();
			int GenControlNum = controlInput.length-2;
			gen.channel.writeInt(GenControlNum);

			inputCounters = gen.newTArray(numCounters, 0);
			for(int i = 0; i < numCounters; i++){
				inputCounters[i] = gen.inputOfAlice(Utils.fromInt(i+1, 32));
			}
			Comparator<Integer[]> comparator = new Comparator<Integer[]>(){
				@Override
				public int compare(Integer[] a, Integer[] b){
					return Integer.compare(a[0], b[0]);
				}
			};
			gen.flush();
			System.out.println("Done with inputCounters gen");
			Integer[][] caseIn = new Integer[GenCaseNum+numCounters][2];
			for(int i = 0; i < numCounters; i++){
				caseIn[i][0] = i+1;
				caseIn[i][1] = 0;
			}			
			for(int i = 0; i < GenCaseNum; i++){
				caseIn[i+numCounters][0] = caseInput[i+2];
				caseIn[i+numCounters][1] = 1;
			}

			Arrays.sort(caseIn, comparator);
			inputAliceCase = gen.newTArray(2, GenCaseNum+numCounters, 0);			
			for(int i = 0; i < GenCaseNum+numCounters; i++){
				inputAliceCase[0][i] = gen.inputOfAlice(Utils.fromInt(caseIn[i][0], 32));
				inputAliceCase[1][i] = gen.inputOfAlice(Utils.fromInt(caseIn[i][1], 32));
			}
			System.out.println();
			gen.flush();
			System.out.println("Done with inputAliceCase gen");

			inputBobCase = gen.newTArray(2, EvaCaseNum+numCounters, 0);
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < EvaCaseNum+numCounters; j++) {
					inputBobCase[i][j] = gen.inputOfBob(new boolean[32]);
				}
			}
			gen.flush();
			System.out.println("Done with inputBobCase gen");
			Integer[][] controlIn = new Integer[GenControlNum+numCounters][2];
			for(int i = 0; i < numCounters; i++){
				controlIn[i][0] = i+1;
				controlIn[i][1] = 0;
			}			
			for(int i = 0; i < GenControlNum; i++){
				controlIn[i+numCounters][0] = controlInput[i+2];
				controlIn[i+numCounters][1] = 1;
			}
			Arrays.sort(controlIn, comparator);

			inputAliceControl = gen.newTArray(2,GenControlNum+numCounters, 0);
			for(int i = 0; i < GenControlNum+numCounters; i++){
				inputAliceControl[0][i] = gen.inputOfAlice(Utils.fromInt(controlIn[i][0], 32));
				inputAliceControl[1][i] = gen.inputOfAlice(Utils.fromInt(controlIn[i][1], 32));
			}
			System.out.println("Done with inputAliceControl gen");
			gen.flush();

			inputBobControl = gen.newTArray(2, EvaControlNum+numCounters, 0);
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < EvaControlNum+numCounters; j++) {
					inputBobControl[i][j] = gen.inputOfBob(new boolean[32]);
				}
			}
			System.out.println("Done with inputBobControl gen");
			gen.flush();
		}
		
		@Override
		public void secureCompute(CompEnv<T> gen) {
			in = compute(gen, inputAliceCase, inputBobCase, inputAliceControl, inputBobControl, 
					inputCounters, aliceCaseNum, bobCaseNum, aliceControlNum, bobControlNum);
		}
		@Override
		public void prepareOutput(CompEnv<T> gen) {
			System.out.println("odds ratio");
			for(int i = in.length-1; i >= 0; i--){
				System.out.println(Utils.toFloat(gen.outputToAlice(in[i]), width, offset) + " ");
			}
		}
	}
	
	public static class Evaluator<T> extends EvaRunnable<T> {
		T[][][] inputBobCase;
		T[][][] inputAliceCase;
		T[][] inputCounters;
		T[][][] inputAliceControl;
		T[][][] inputBobControl;
		T[] scResult;
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
			int bobCaseNumInt = caseInput[0];
			int bobControlNumInt = controlInput[0];

			System.out.println("Bob case num " + bobCaseNumInt);
			System.out.println("Bob control num " + bobControlNumInt);

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

			inputCounters = gen.newTArray(numCounters, 0);
			for (int j = 0; j < numCounters; j++) {
					inputCounters[j] = gen.inputOfAlice(new boolean[32]);
			}
			gen.flush();
			System.out.println("Done with inputCounters eva");
			Comparator<Integer[]> comparator = new Comparator<Integer[]>(){
				@Override
				public int compare(Integer[] a, Integer[] b){
					return Integer.compare(a[0], b[0]);
				}
			};
			inputAliceCase = gen.newTArray(2, GenCaseNum+numCounters, 0);
			for (int i = 0; i < 2; ++i) {
				for (int j = 0; j < GenCaseNum+numCounters; j++) {
					inputAliceCase[i][j] = gen.inputOfAlice(new boolean[32]);
				}
			}

			gen.flush();
			System.out.println("Done with inputAliceCase eva");

			Integer[][] caseIn = new Integer[EvaCaseNum+numCounters][2];
			for(int i = 0; i < numCounters; i++){
				caseIn[i][0] = i+1;
				caseIn[i][1] = 0;
			}			
			for(int i = 0; i < EvaCaseNum; i++){
				caseIn[i+numCounters][0] = caseInput[i+2];
				caseIn[i+numCounters][1] = 1;
			}
			Arrays.sort(caseIn, comparator);

			inputBobCase = gen.newTArray(2, EvaCaseNum+numCounters, 0);			
			for(int i = 0; i < EvaCaseNum+numCounters; i++){
				inputBobCase[0][i] = gen.inputOfBob(Utils.fromInt(caseIn[i][0], 32));
			}
			for(int i = 0; i < EvaCaseNum+numCounters; i++){
				inputBobCase[1][i] = gen.inputOfBob(Utils.fromInt(caseIn[i][1], 32));
			}
			System.out.println();
			gen.flush();

			inputAliceControl = gen.newTArray(2, GenControlNum+numCounters, 0);
			for (int i = 0; i < 2; ++i) {
				for (int j = 0; j < GenControlNum+numCounters; j++) {
					inputAliceControl[i][j] = gen.inputOfAlice(new boolean[32]);
				}
			}
			gen.flush();
			System.out.println("Done with inputAliceControl eva");
			Integer[][] controlIn = new Integer[EvaControlNum+numCounters][2];
			for(int i = 0; i < numCounters; i++){
				controlIn[i][0] = i+1;
				controlIn[i][1] = 0;
			}			
			for(int i = 0; i < EvaControlNum; i++){
				controlIn[i+numCounters][0] = controlInput[i+2];
				controlIn[i+numCounters][1] = 1;
			}
			Arrays.sort(controlIn, comparator);

			inputBobControl = gen.newTArray(2,EvaControlNum+numCounters, 0);
			for(int i = 0; i < EvaControlNum+numCounters; i++){
				inputBobControl[0][i] = gen.inputOfBob(Utils.fromInt(controlIn[i][0], 32));
			}
			for(int i = 0; i < EvaControlNum+numCounters; i++){
				inputBobControl[1][i] = gen.inputOfBob(Utils.fromInt(controlIn[i][1], 32));
			}

			gen.flush();
			System.out.println("Done with inputBobControl eva");
		}
		
		@Override
		public void secureCompute(CompEnv<T> gen) {
			in = compute(gen, inputAliceCase, inputBobCase, inputAliceControl, inputBobControl, inputCounters,
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
