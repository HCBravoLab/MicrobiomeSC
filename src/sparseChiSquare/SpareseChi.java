package sparseChiSquare;

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
//import example.PrepareData.StatisticsData;
import flexsc.CompEnv;

public class SpareseChi {
	static public<T> T[][] compute(CompEnv<T> gen, T[][][] inputCounters, T[][][] inputAliceCase, T[][][] inputAliceControl, 
			T[][][] inputBobCase, T[][][] inputBobControl, T[] numAliceCase, T[] numAliceControl, T[] numBobCase, T[] numBobControl){
		
		int d1 = inputAliceCase.length;
		int countersLength = inputCounters[0].length;
		int d2Case = inputCounters[0].length + inputAliceCase[0].length + inputBobCase[0].length;
		int d3 = 0;
		int d2Control = inputCounters[0].length + inputAliceControl[0].length + inputBobControl[0].length;

		T[][][] inCase = gen.newTArray(d1, d2Case, d3);
		T[][][] inControl = gen.newTArray(d1, d2Control, d3);

		T[][][] outCase = gen.newTArray(d1, countersLength, d3);
		T[][][] outControl = gen.newTArray(d1, countersLength, d3);
		
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

		BitonicSortLib<T> lib = new  BitonicSortLib<T>(gen);
		IntegerLib<T> ilib = new IntegerLib<T>(gen);

		//Sort with payload first row and second row as payload

		int[] attachedRows = {1,2};
		//Sort with payload first row and third row as payload
		lib.sortWithPayloadM(inCase[0], inCase, attachedRows, 1, lib.SIGNAL_ONE);	
		lib.sortWithPayloadM(inCase[0], inCase, attachedRows, 1, lib.SIGNAL_ONE);	

		//make accumulation pass
		// s[i][2] = s[i][1] + s[i+1][1]*s[i+1][2]
		
		for(int i = (inCase[0].length-2); i >= 0; i--){
			inCase[2][i] = ilib.add(inCase[1][i], ilib.multiply(inCase[1][i+1], inCase[2][i+1]));
		}


		int[] attachedRows2 = {0,2};
		//Sort with payload first row and third row as payload
		lib.sortWithPayloadM(inCase[1], inCase, attachedRows2, 0, lib.SIGNAL_ONE);	
		lib.sortWithPayloadM(inCase[1], inCase, attachedRows2, 0, lib.SIGNAL_ONE);	

		System.arraycopy(inCase[0], 0, outCase[0], 0, countersLength);
		System.arraycopy(inCase[1], 0, outCase[1], 0, countersLength);
		System.arraycopy(inCase[2], 0, outCase[2], 0, countersLength);

		int[] attachedRows3 = {1,2};
		lib.sortWithPayloadM(outCase[0], outCase, attachedRows3, 0, lib.SIGNAL_ONE);		
		lib.sortWithPayloadM(outCase[0], outCase, attachedRows3, 0, lib.SIGNAL_ONE);		

		//Sort with payload first row and second row as payload

		int[] attachedRows4 = {1,2};
		lib.sortWithPayloadM(inControl[0], inControl, attachedRows4, 1, lib.SIGNAL_ONE);	
		lib.sortWithPayloadM(inControl[0], inControl, attachedRows4, 1, lib.SIGNAL_ONE);	

		//make accumulation pass
		// s[i][2] = s[i][1] + s[i+1][1]*s[i+1][2]
		
		for(int i = (inControl[0].length-2); i >= 0; i--){
			inControl[2][i] = ilib.add(inControl[1][i], ilib.multiply(inControl[1][i+1], inControl[2][i+1]));
		}

		int[] attachedRows5 = {0,2};
		lib.sortWithPayloadM(inControl[1], inControl, attachedRows5, 0, lib.SIGNAL_ONE);	
		lib.sortWithPayloadM(inControl[1], inControl, attachedRows5, 0, lib.SIGNAL_ONE);	

		System.arraycopy(inControl[0], 0, outControl[0], 0, countersLength);
		System.arraycopy(inControl[1], 0, outControl[1], 0, countersLength);
		System.arraycopy(inControl[2], 0, outControl[2], 0, countersLength);
		
		int[] attachedRows6 = {1,2};
		lib.sortWithPayloadM(outControl[0], outControl, attachedRows6, 0, lib.SIGNAL_ONE);	
		lib.sortWithPayloadM(outControl[0], outControl, attachedRows6, 0, lib.SIGNAL_ONE);	


		int[] attachedRows7 = {1,2};
		lib.sortWithPayloadM(inControl[0], inControl, attachedRows7, 0, lib.SIGNAL_ONE);
		lib.sortWithPayloadM(inControl[0], inControl, attachedRows7, 0, lib.SIGNAL_ONE);

		T[][] res = gen.newTArray(inputCounters[0].length, 0);

		FloatLib<T> flib = new FloatLib<T>(gen, 30, 10);

		for(int i = 0; i < inputCounters[0].length; i++){			
			//T[] a = lib.add(aliceCase[i][0], bobCase[i][0]);
			T[] a = outCase[2][i];

			T[] b = ilib.sub(ilib.add(numAliceCase, numBobCase), outCase[2][i]);
			
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
		
		//return outCase[2];
		//lib.sortWithPayloadM(inControl[0], inControl, attachedRows, 0, lib.SIGNAL_ONE);
		//return inControl[0];
		//return outControl[2];
		return res;
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
				for(int j = 0; j < inputCounters[0].length; j++){
					System.out.print(Utils.toFloat(gen.outputToAlice(ret[j]), 30, 10) + " ");
					//System.out.print(Utils.toInt(gen.outputToAlice(ret[j])) + " ");

				}
				//System.out.println();
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
				for(int j = 0; j < inputCounters[0].length; j++){
					gen.outputToAlice(ret[j]);
				}
			//}
		}
				
	}
	
}
