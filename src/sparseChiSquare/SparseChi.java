package sparseChiSquare;

import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import util.EvaRunnable;
import util.GenRunnable;
import util.Utils;
import circuits.BitonicSortLib;
import circuits.arithmetic.FloatLib;
import circuits.arithmetic.IntegerLib;
//import example.PrepareData.StatisticsData;
import flexsc.CompEnv;

public class SparseChi {
	static public<T> T[][] compute(CompEnv<T> gen, T[][][] inputCounters, T[][][] inputAliceCase, T[][][] inputAliceControl, 
			T[][][] inputBobCase, T[][][] inputBobControl, T[] numAliceCase, T[] numAliceControl, T[] numBobCase, T[] numBobControl){
		
		int d1 = inputAliceCase.length;
		int d2Case = inputCounters[0].length + inputAliceCase[0].length + inputBobCase[0].length;
		int d3 = 0;
		int d2Control = inputCounters[0].length + inputAliceControl[0].length + inputBobControl[0].length;

		T[][][] inCase = gen.newTArray(d1, d2Case, d3);
		T[][] in1Case = gen.newTArray(d2Case, 0);
		T[][] in2Case = gen.newTArray(d2Case, 0);
		T[][] in3Case = gen.newTArray(d2Case, 0);
		T[][][] inControl = gen.newTArray(d1, d2Control, d3);
		T[][] in1Control = gen.newTArray(d2Control, 0);
		T[][] in2Control = gen.newTArray(d2Control, 0);
		T[][] in3Control = gen.newTArray(d2Control, 0);

		for(int i = 0; i < 3; i++){
			System.arraycopy(inputCounters[i], 0, inCase[i], 0, inputCounters[i].length);
			System.arraycopy(inputAliceCase[i], 0, inCase[i], inputCounters[i].length, inputAliceCase[i].length);
			System.arraycopy(inputBobCase[i], 0, inCase[i], inputCounters[i].length + inputAliceCase[i].length, inputBobCase[i].length);
		}

		for(int i = 0; i < 3; i++){
			System.arraycopy(inputCounters[i], 0, inControl[i], 0, inputCounters[i].length);
			System.arraycopy(inputAliceControl[i], 0, inControl[i], inputCounters[i].length, inputAliceControl[i].length);
			System.arraycopy(inputBobControl[i], 0, inControl[i], inputCounters[i].length + inputAliceControl[i].length, inputBobControl[i].length);
		}

		System.arraycopy(inCase[0], 0, in1Case, 0, in1Case.length);
		System.arraycopy(inCase[1], 0, in2Case, 0, in2Case.length);
		System.arraycopy(inCase[2], 0, in3Case, 0, in3Case.length);

		BitonicSortLib<T> lib = new  BitonicSortLib<T>(gen);
		//Make two copies of first row
		T[][] in1UnsortedCase = gen.newTArray(in1Case.length, 0);
		System.arraycopy(in1Case, 0, in1UnsortedCase, 0, in1Case.length);
		
		//Sort with payload first row and second row as payload
		lib.sortWithPayloadS(in1Case, in2Case, lib.SIGNAL_ONE);		
		
		//Sort with payload first row and third row as payload
		lib.sortWithPayloadS(in1UnsortedCase, in3Case, lib.SIGNAL_ONE);		

		IntegerLib<T> lib1 = new IntegerLib<T>(gen);
		
		T[] one = lib1.publicValue(1.0);
		
		//make accumulation pass
		// s[i][2] = s[i][1] + s[i+1][1]*s[i+1][2]
		
		for(int i = (in1Case.length-2); i >= 0; i--){
			in3Case[i] = lib.add(in2Case[i], lib.multiply(in2Case[i+1], in3Case[i+1]));
		}

		T[][] in2UnsortedCase = gen.newTArray(in2Case.length, 0);
		System.arraycopy(in2Case, 0, in2UnsortedCase, 0, in1Case.length);
		
		//Sort with payload first row and second row as payload
		lib.sortWithPayload(in2Case, in1Case, lib.SIGNAL_ONE);		
		
		//Sort with payload first row and third row as payload
		lib.sortWithPayload(in2UnsortedCase, in3Case, lib.SIGNAL_ONE);
		
		T[][][] outCase = gen.newTArray(d1, inputCounters[0].length, d3);

		System.arraycopy(in1Case, 0, outCase[0], 0, inputCounters[0].length);
		System.arraycopy(in2Case, 0, outCase[1], 0, inputCounters[0].length);
		System.arraycopy(in3Case, 0, outCase[2], 0, inputCounters[0].length);

		T[][] out1UnsortedCase = gen.newTArray(outCase[0].length, 0);
		System.arraycopy(outCase[0], 0, out1UnsortedCase, 0, out1UnsortedCase.length);
		T[][] out1Case = gen.newTArray(outCase[0].length, 0);
		T[][] out2Case = gen.newTArray(outCase[0].length, 0);
		T[][] out3Case = gen.newTArray(outCase[0].length, 0);

		System.arraycopy(outCase[0], 0, out1Case, 0, out1UnsortedCase.length);
		System.arraycopy(outCase[1], 0, out2Case, 0, out1UnsortedCase.length);
		System.arraycopy(outCase[2], 0, out3Case, 0, out1UnsortedCase.length);

		//Sort with payload first row and second row as payload
		lib.sortWithPayload(out1Case, out2Case, lib.SIGNAL_ONE);		
		
		//Sort with payload first row and third row as payload
		lib.sortWithPayload(out1UnsortedCase, out3Case, lib.SIGNAL_ONE);

		System.arraycopy(out1Case, 0, outCase[0], 0, inputCounters[0].length);
		System.arraycopy(out2Case, 0, outCase[1], 0, inputCounters[0].length);
		System.arraycopy(out3Case, 0, outCase[2], 0, inputCounters[0].length);

		System.arraycopy(inControl[0], 0, in1Control, 0, in1Control.length);
		System.arraycopy(inControl[1], 0, in2Control, 0, in2Control.length);
		System.arraycopy(inControl[2], 0, in3Control, 0, in3Control.length);
		
		//Make two copies of first row
		T[][] in1UnsortedControl = gen.newTArray(in1Control.length, 0);
		System.arraycopy(in1Control, 0, in1UnsortedControl, 0, in1Control.length);
		
		//Sort with payload first row and second row as payload
		lib.sortWithPayloadS(in1Control, in2Control, lib.SIGNAL_ONE);		
		
		//Sort with payload first row and third row as payload
		lib.sortWithPayloadS(in1UnsortedControl, in3Control, lib.SIGNAL_ONE);		
		
		//make accumulation pass
		// s[i][2] = s[i][1] + s[i+1][1]*s[i+1][2]
		
		for(int i = (in1Control.length-2); i >= 0; i--){
			in3Control[i] = lib.add(in2Control[i], lib.multiply(in2Control[i+1], in3Control[i+1]));
		}

		T[][] in2UnsortedControl = gen.newTArray(in2Control.length, 0);
		System.arraycopy(in2Control, 0, in2UnsortedControl, 0, in1Control.length);
		
		//Sort with payload first row and second row as payload
		lib.sortWithPayload(in2Control, in1Control, lib.SIGNAL_ONE);		
		
		//Sort with payload first row and third row as payload
		lib.sortWithPayload(in2UnsortedControl, in3Control, lib.SIGNAL_ONE);
		
		T[][][] outControl = gen.newTArray(d1, inputCounters[0].length, d3);

		System.arraycopy(in1Control, 0, outControl[0], 0, inputCounters[0].length);
		System.arraycopy(in2Control, 0, outControl[1], 0, inputCounters[0].length);
		System.arraycopy(in3Control, 0, outControl[2], 0, inputCounters[0].length);

		T[][] out1UnsortedControl = gen.newTArray(outControl[0].length, 0);
		System.arraycopy(outControl[0], 0, out1UnsortedControl, 0, out1UnsortedControl.length);
		T[][] out1Control = gen.newTArray(outControl[0].length, 0);
		T[][] out2Control = gen.newTArray(outControl[0].length, 0);
		T[][] out3Control = gen.newTArray(outControl[0].length, 0);

		System.arraycopy(outControl[0], 0, out1Control, 0, out1UnsortedControl.length);
		System.arraycopy(outControl[1], 0, out2Control, 0, out1UnsortedControl.length);
		System.arraycopy(outControl[2], 0, out3Control, 0, out1UnsortedControl.length);

		//Sort with payload first row and second row as payload
		lib.sortWithPayload(out1Control, out2Control, lib.SIGNAL_ONE);		
		
		//Sort with payload first row and third row as payload
		lib.sortWithPayload(out1UnsortedControl, out3Control, lib.SIGNAL_ONE);

		System.arraycopy(out1Control, 0, outControl[0], 0, inputCounters[0].length);
		System.arraycopy(out2Control, 0, outControl[1], 0, inputCounters[0].length);
		System.arraycopy(out3Control, 0, outControl[2], 0, inputCounters[0].length);
		
		T[][] res = gen.newTArray(inputCounters[0].length, 0);

		IntegerLib<T> ilib = new IntegerLib<T>(gen);
		FloatLib<T> flib = new FloatLib<T>(gen, 49, 8);

		for(int i = 0; i < inputCounters[0].length; i++){			
			//T[] a = lib.add(aliceCase[i][0], bobCase[i][0]);
			T[] a = outCase[2][i];

			T[] b = ilib.sub(lib.add(numAliceCase, numBobCase), outCase[2][i]);
			
			//T[] c = lib.add(aliceControl[i][0], bobControl[i][0]);
			T[] c = outControl[2][i];

			T[] d = ilib.sub(ilib.add(numAliceControl, numBobControl), outControl[2][i]);
			
			T[] g = ilib.add(a, c);
			T[] k = ilib.add(b, d);

			T[] fa = ilib.toSecureFloat(a, flib);
			T[] fb = ilib.toSecureFloat(b, flib);
			T[] fc = ilib.toSecureFloat(c, flib);
			T[] fd = ilib.toSecureFloat(d, flib);
			T[] fg = ilib.toSecureFloat(g, flib);
			T[] fk = ilib.toSecureFloat(k, flib);

			T[] tmp = flib.sub(flib.multiply(fa, fd), flib.multiply(fb, fc));
			tmp = flib.multiply(tmp, tmp);
			res[i] = flib.div(tmp, flib.multiply(fg, fk));
		}
		//return res;
		//return in;
		
		//return outCase;

		return outControl[2];
		//return res;
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

		@Override
		public void prepareInput(CompEnv<T> gen) throws Exception {
			//Take command line input from Generator and parse file
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
			
			// Store the number of samples, which will be used to compute b and d of contingency counts table
			numAliceCase = gen.inputOfAlice(Utils.fromInt(caseInput[0], 32));
			numAliceControl = gen.inputOfAlice(Utils.fromInt(controlInput[0], 32));
			gen.flush();

			numBobCase = gen.inputOfBob(new boolean[32]);
			numBobControl = gen.inputOfBob(new boolean[32]);
			gen.flush();

			System.out.println("Gen number of samples case: " + caseInput[0]);
			System.out.println("Gen number of samples control: " + controlInput[0]);

			System.out.println("Gen number of features case: " + caseInput[1]);
			System.out.println("Gen number of features control: " + controlInput[1]);

			System.out.println("Gen number of counts case: " + caseInput.length);
			System.out.println("Gen number of counts control: " + controlInput.length);		
			int EvaCaseCounts = gen.channel.readInt();
			gen.channel.flush();
			int GenCaseCounts = caseInput.length-2;
			gen.channel.writeInt(GenCaseCounts);
			gen.channel.flush();

			int EvaControlCounts = gen.channel.readInt();
			gen.channel.flush();
			int GenControlCounts = controlInput.length-2;
			gen.channel.writeInt(GenControlCounts);
			gen.channel.flush();
			
			ArrayList<Integer> alCase = new ArrayList<Integer>();
			for(int i= 2; i < caseInput.length; i++){
				alCase.add(caseInput[i]);
			}
			Collections.sort(alCase);
			ArrayList<Integer> alControl = new ArrayList<Integer>();
			for(int i= 2; i < controlInput.length; i++){
				alControl.add(controlInput[i]);
			}
			for(int i =2; i <caseInput.length;i++){
				caseInput[i] = alCase.get(i-2);
			}
			for(int i =2; i <controlInput.length;i++){
				controlInput[i] = alControl.get(i-2);
			}
			
			inputCounters = gen.newTArray(3, caseInput[1], 0);
			for(int i = 0; i < caseInput[1];i++){
				inputCounters[0][i] = gen.inputOfAlice(Utils.fromInt(i+1, 32));
			}
			for(int i = 0; i < caseInput[1];i++){
				inputCounters[1][i] = gen.inputOfAlice(Utils.fromInt(0, 32));
			}
			for(int i = 0; i < caseInput[1];i++){
				inputCounters[2][i] = gen.inputOfAlice(Utils.fromInt(0, 32));
			}
			gen.flush();
			
			inputAliceCase = gen.newTArray(3, GenCaseCounts, 0);
			for(int i = 0; i < GenCaseCounts; i++){
				inputAliceCase[0][i] = gen.inputOfAlice(Utils.fromInt(caseInput[i+2], 32));
			}
			for(int i = 0; i < GenCaseCounts; i++){
				inputAliceCase[1][i] = gen.inputOfAlice(Utils.fromInt(1, 32));
			}
			for(int i = 0; i < GenCaseCounts; i++){
				inputAliceCase[2][i] = gen.inputOfAlice(Utils.fromInt(1, 32));
			}
			gen.flush();
			inputAliceControl = gen.newTArray(3, GenControlCounts, 0);
			for(int i = 0; i < GenControlCounts; i++){
				inputAliceControl[0][i] = gen.inputOfAlice(Utils.fromInt(controlInput[i+2], 32));
			}
			for(int i = 0; i < GenControlCounts; i++){
				inputAliceControl[1][i] = gen.inputOfAlice(Utils.fromInt(1, 32));
			}
			for(int i = 0; i < GenControlCounts; i++){
				inputAliceControl[2][i] = gen.inputOfAlice(Utils.fromInt(1, 32));
			}
			gen.flush();
			
			inputBobCase = gen.newTArray(3, EvaCaseCounts, 0);
			for(int i = 0; i < 3; i++){
				for(int j = 0; j < EvaCaseCounts; j++){
					inputBobCase[i][j] = gen.inputOfBob(new boolean[32]);		

				}
			}
			gen.flush();

			inputBobControl = gen.newTArray(3, EvaControlCounts, 0);
			for(int i = 0; i < 3; i++){
				for(int j = 0; j < EvaControlCounts; j++){
					inputBobControl[i][j] = gen.inputOfBob(new boolean[32]);
				}
			}
			gen.flush();
		}
		
		@Override
		public void secureCompute(CompEnv<T> gen) {
			ret = compute(gen, inputCounters, inputAliceCase, inputAliceControl, inputBobCase, inputBobControl,
					 numAliceCase, numAliceControl, numBobCase, numBobControl);
		}
		@Override
		public void prepareOutput(CompEnv<T> gen) {
			//for(int i = 0; i < ret.length; i++){
				for(int j = 0; j < inputCounters[0].length-1; j++){
					//System.out.print(Utils.toFloat(gen.outputToAlice(ret[j]), 49, 8) + " ");
					System.out.print(Utils.toInt(gen.outputToAlice(ret[j])) + " ");
				}
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

			numAliceCase = gen.inputOfAlice(new boolean[32]);
			numAliceControl = gen.inputOfAlice(new boolean[32]);
			gen.flush();

			numBobCase = gen.inputOfBob(Utils.fromInt(caseInput[0], 32));
			numBobControl = gen.inputOfBob(Utils.fromInt(controlInput[0], 32));
			gen.flush();

			System.out.println("Eva number of samples case: " + caseInput[0]);
			System.out.println("Eva number of samples control: " + controlInput[0]);

			System.out.println("Eva number of features case: " + caseInput[1]);
			System.out.println("Eva number of features control: " + controlInput[1]);

			System.out.println("Eva number of counts case: " + caseInput.length);
			System.out.println("Eva number of counts control: " + controlInput.length);		

			int EvaCaseCounts = caseInput.length-2;
			gen.channel.writeInt(EvaCaseCounts);
			gen.channel.flush();
			int GenCaseCounts = gen.channel.readInt();
			gen.channel.flush();

			int EvaControlCounts = controlInput.length-2;
			gen.channel.writeInt(EvaControlCounts);
			gen.channel.flush();
			int GenControlCounts = gen.channel.readInt();
			gen.channel.flush();
			
			ArrayList<Integer> alCase = new ArrayList<Integer>();
			for(int i= 1; i < caseInput.length; i++){
				alCase.add(caseInput[i]);
			}
			Collections.sort(alCase);
			ArrayList<Integer> alControl = new ArrayList<Integer>();
			for(int i= 1; i < controlInput.length; i++){
				alControl.add(controlInput[i]);
			}
			for(int i =1; i <caseInput.length;i++){
				caseInput[i] = alCase.get(i-1);
			}
			for(int i =1; i <controlInput.length;i++){
				controlInput[i] = alControl.get(i-1);
			}
			inputCounters = gen.newTArray(3, caseInput[1], 0);
			for(int i = 0; i < caseInput[1];i++){
				inputCounters[0][i] = gen.inputOfAlice(new boolean[32]);
			}
			for(int i = 0; i < caseInput[1];i++){
				inputCounters[1][i] = gen.inputOfAlice(new boolean[32]);
			}
			for(int i = 0; i < caseInput[1];i++){
				inputCounters[2][i] = gen.inputOfAlice(new boolean[32]);
			}
			gen.flush();
			
			inputAliceCase = gen.newTArray(3, GenControlCounts, 0);
			for(int i = 0; i < 3; i++){
				for(int j = 0; j < GenControlCounts; j++){
					inputAliceCase[i][j] = gen.inputOfAlice(new boolean[32]);		
				}
			}
			gen.flush();

			inputAliceControl = gen.newTArray(3, GenControlCounts, 0);
			for(int i = 0; i < 3; i++){
				for(int j = 0; j < GenControlCounts; j++){
					inputAliceControl[i][j] = gen.inputOfAlice(new boolean[32]);
				}
			}
			gen.flush();

			inputBobCase = gen.newTArray(3, EvaCaseCounts, 0);
			for(int i = 0; i < EvaCaseCounts; i++){
				inputBobCase[0][i] = gen.inputOfBob(Utils.fromInt(caseInput[i+2], 32));
			}
			for(int i = 0; i < EvaCaseCounts; i++){
				inputBobCase[1][i] = gen.inputOfBob(Utils.fromInt(1, 32));
			}
			for(int i = 0; i < EvaCaseCounts; i++){
				inputBobCase[2][i] = gen.inputOfBob(Utils.fromInt(1, 32));
			}
			gen.flush();
			
			inputBobControl = gen.newTArray(3, EvaControlCounts, 0);
			for(int i = 0; i < EvaControlCounts; i++){
				inputBobControl[0][i] = gen.inputOfBob(Utils.fromInt(controlInput[i+2], 32));
			}
			for(int i = 0; i < EvaControlCounts; i++){
				inputBobControl[1][i] = gen.inputOfBob(Utils.fromInt(1, 32));
			}
			for(int i = 0; i < EvaControlCounts; i++){
				inputBobControl[2][i] = gen.inputOfBob(Utils.fromInt(1, 32));
			}
			gen.flush();
		}
		
		@Override
		public void secureCompute(CompEnv<T> gen) {
			ret = compute(gen, inputCounters, inputAliceCase, inputAliceControl, inputBobCase, inputBobControl,
					 numAliceCase, numAliceControl, numBobCase, numBobControl);
		}
		
		@Override
		public void prepareOutput(CompEnv<T> gen) {
			//for(int i = 0; i < in.length; i++){
				for(int j = 0; j < inputCounters[0].length-1; j++){
					gen.outputToAlice(ret[j]);
				}
			//}
		}
				
	}
	
}
