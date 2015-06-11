package perFeatureSparse;

import perFeatureSparse.PrepareDataDANaive;

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
import circuits.CircuitLib;

import flexsc.CompEnv;

public class DifferentialAbundanceSafe {
	static int PLength = 54;
	static int VLength = 11;
	static int filterThreshold = 1;

	static public<T> T[][] dummyVariable(CompEnv<T> gen){
		FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
		T[] zero = flib.publicValue(0.0);
		T[][] res = gen.newTArray(2, 0);
		res[0] = zero;
		res[1] = zero;
		return res;
	}
	
	static public<T> T[][] onlyOne(CompEnv<T> gen, T[][] inputAliceCase, 
			T[][] inputBobCase, T[][] inputAliceControl, T[][] inputBobControl, T[] numAliceCase, T[] numBobCase,
			T[] numAliceControl, T[] numBobControl){
		FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
		IntegerLib<T> ilib = new IntegerLib<T>(gen, 32);

		CircuitLib<T> cl = new CircuitLib<T>(gen);

		T[] zero = flib.publicValue(0.0);
		T[] one = flib.publicValue(1.0);
		T[] negOne = flib.publicValue(-1.0);
		T[] pointOne = flib.publicValue(.01);

		T[][] res = gen.newTArray(2, 0);
		T[] caseSum = flib.publicValue(0.0);
		T[] caseSumOfSquares = flib.publicValue(0.0);
		T[] caseShift = pointOne;
		for(int i = 0; i < inputAliceCase.length; i++){
			caseSum = flib.add(caseSum, flib.sub(inputAliceCase[i], caseShift));
			T[] addSquare = flib.multiply(flib.sub(inputAliceCase[i], caseShift), flib.sub(inputAliceCase[i], caseShift));
			caseSumOfSquares = flib.add(caseSumOfSquares, addSquare);
		}
		
		for(int i = 0; i < inputBobCase.length; i++){
			caseSum = flib.add(caseSum, flib.sub(inputBobCase[i], caseShift));
			T[] addSquare = flib.multiply(flib.sub(inputBobCase[i], caseShift), flib.sub(inputBobCase[i], caseShift));
			caseSumOfSquares = flib.add(caseSumOfSquares, addSquare);
		}
		
		T[] controlSum = flib.publicValue(0.0);
		T[] controlSumOfSquares = flib.publicValue(0.0);
		T[] controlShift = pointOne;
		for(int i = 0; i < inputAliceControl.length; i++){
			controlSum = flib.add(controlSum, flib.sub(inputAliceControl[i], controlShift));
			T[] addSquare = flib.multiply(flib.sub(inputAliceControl[i], controlShift), flib.sub(inputAliceControl[i], controlShift));
			controlSumOfSquares = flib.add(controlSumOfSquares, addSquare);
		}
		
		for(int i = 0; i < inputBobControl.length; i++){
			controlSum = flib.add(controlSum, flib.sub(inputBobControl[i], controlShift));
			T[] addSquare = flib.multiply(flib.sub(inputBobControl[i], controlShift), flib.sub(inputBobControl[i], controlShift));
			controlSumOfSquares = flib.add(controlSumOfSquares, addSquare);
		}

		T[] caseNum = flib.add(ilib.toSecureFloat(numAliceCase, flib), ilib.toSecureFloat(numBobCase, flib));
		T[] controlNum = flib.add(ilib.toSecureFloat(numAliceControl, flib), ilib.toSecureFloat(numBobControl, flib));
		T[] caseNumPresent = flib.add(flib.publicValue(inputAliceCase.length), flib.publicValue(inputBobCase.length));
		T[] controlNumPresent = flib.add(flib.publicValue(inputAliceControl.length), flib.publicValue(inputBobControl.length));
		
		T sendOneOrZero = ilib.not(flib.leq(caseNumPresent, controlNumPresent));
		T[] tStatZero = cl.mux(negOne, one, sendOneOrZero);
		
		//T[] tStatCaseNonZero = flib.div(caseSum, caseNum);
		//T[] tStatControlNonZero = flib.div(controlSum, controlNum);
		T[] caseMeanAbundance = flib.div(flib.add(caseSum, caseShift), caseNum);
		T[] caseVarianceSecondTerm = flib.div(flib.multiply(caseSum, caseSum), caseNum);
		T[] caseVariance = flib.div(flib.sub(caseSumOfSquares, caseVarianceSecondTerm), caseNum);
		T[] controlMeanAbundance = flib.div(flib.add(controlSum, controlShift), controlNum);		    
		T[] controlVarianceSecondTerm = flib.div(flib.multiply(controlSum, controlSum), controlNum);
		T[] controlVariance = flib.div(flib.sub(controlSumOfSquares, controlVarianceSecondTerm), controlNum);

		T[] tUpper = flib.sub(controlMeanAbundance, caseMeanAbundance);
		T tLowerFirstIsZero = flib.eq(caseNum, zero);
		T[] tLowerFirst = flib.div(caseVariance, caseNum);
		tLowerFirst = cl.mux(tLowerFirst, zero, tLowerFirstIsZero);
		
		T tLowerFirstIsLessThanOne = ilib.and(flib.leq(tLowerFirst, pointOne), ilib.not(flib.eq(tLowerFirst, zero)));
		T tLowerSecondIsZero = flib.eq(controlNum, zero);
		
		T[] tLowerSecond = flib.div(controlVariance, controlNum);
		tLowerSecond = cl.mux(tLowerSecond, zero, tLowerSecondIsZero);
		T tLowerSecondIsLessThanOne = ilib.and(flib.leq(tLowerSecond, pointOne), ilib.not(flib.eq(tLowerSecond, zero)));
		T oneLessThanOne = ilib.or(tLowerFirstIsLessThanOne,tLowerSecondIsLessThanOne);
		
		//T[] tLowerSqrtLessThanOne = flib.sqrt(flib.add(flib.add(zero,tLowerFirst), flib.add(zero,tLowerSecond)));
		//T[] tLowerSqrt = flib.sqrt(flib.add(tLowerFirst,tLowerSecond));
		T[] tLowerSqrt = flib.sqrt(flib.add(flib.add(zero,tLowerFirst), flib.add(zero,tLowerSecond)));

		//tLowerSqrt = cl.mux(tLowerSqrt, tLowerSqrtLessThanOne, oneLessThanOne);
		T[] tStatNonZero = flib.div(tUpper, tLowerSqrt);
		//T[] tStatNonZero = flib.add(flib.add(zero,tLowerFirst), flib.add(zero,tLowerSecond));
		T caseIsGreaterThanOne = ilib.not(flib.leq(caseNumPresent, one));
		T controlIsGreaterThanOne = ilib.not(flib.leq(controlNumPresent, one));
		T caseIsOne = flib.eq(caseNumPresent, one);
		T controlIsOne = flib.eq(controlNumPresent, one);
		T isOne = ilib.or(caseIsOne, controlIsOne);
		T greaterThanOne = ilib.or(caseIsGreaterThanOne, controlIsGreaterThanOne);
		T nonOne = ilib.and(greaterThanOne, ilib.not(isOne));
		T returnCase = flib.leq(caseNumPresent, controlNumPresent);
		
		//T[] tStatNonZero = cl.mux(flib.sub(zero, tStatCaseNonZero), tStatControlNonZero, returnCase);
	
		T[] tStat = cl.mux(tStatZero, tStatNonZero, nonOne);
		T[] df = cl.mux(flib.sub(caseNum, one), flib.sub(controlNum, one), returnCase);
		
		res[0] = tStat;
		res[1] = df;
		return res;
	}
	
	static public<T> T filter(CompEnv<T> gen, T[] numAliceCase, 
			T[] numBobCase, T[] numAliceControl, T[] numBobControl){

		IntegerLib<T> ilib = new IntegerLib<T>(gen, 32);
		T[] zero = ilib.publicValue(0.0);
		T[] one = ilib.publicValue(0.0);

		T[] threshold = ilib.publicValue(filterThreshold);
		T[] caseNum = ilib.publicValue(0.0);
		T[] controlNum = ilib.publicValue(0.0);
		T caseAboveThreshold;
		T controlAboveThreshold;
		T bothAboveThreshold;
		caseNum = ilib.add(numAliceCase, numBobCase);
		controlNum = ilib.add(numAliceControl,numBobControl);
		caseAboveThreshold = ilib.geq(caseNum, threshold);
		controlAboveThreshold = ilib.geq(controlNum, threshold);
		T oneZero = ilib.not(ilib.xor(caseAboveThreshold, controlAboveThreshold));
		return oneZero;
	}
	
	static public<T> T[][] compute(CompEnv<T> gen, T[][] inputAliceCase, 
			T[][] inputBobCase, T[][] inputAliceControl, T[][] inputBobControl, T[] numAliceCase, T[] numBobCase,
			T[] numAliceControl, T[] numBobControl){
		FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
		IntegerLib<T> ilib = new IntegerLib<T>(gen, 32);
		CircuitLib<T> cl = new CircuitLib<T>(gen);

		T[] zero = flib.publicValue(0.0);
		T[] one = flib.publicValue(1.0);
		T[] negOne = flib.publicValue(-1.0);

		T[] pointOne = flib.publicValue(.01);
		//T[] zero = flib.publicValue(0.0);

		T[] caseSum = flib.publicValue(0.0);
		T[] caseSumOfSquares = flib.publicValue(0.0);
		T[] caseShift = pointOne;
		for(int i = 0; i < inputAliceCase.length; i++){
			caseSum = flib.add(caseSum, flib.sub(inputAliceCase[i], caseShift));
			T[] addSquare = flib.multiply(flib.sub(inputAliceCase[i], caseShift), flib.sub(inputAliceCase[i], caseShift));
			caseSumOfSquares = flib.add(caseSumOfSquares, addSquare);
		}
		
		for(int i = 0; i < inputBobCase.length; i++){
			caseSum = flib.add(caseSum, flib.sub(inputBobCase[i], caseShift));
			T[] addSquare = flib.multiply(flib.sub(inputBobCase[i], caseShift), flib.sub(inputBobCase[i], caseShift));
			caseSumOfSquares = flib.add(caseSumOfSquares, addSquare);
		}
		
		T[] controlSum = flib.publicValue(0.0);
		T[] controlSumOfSquares = flib.publicValue(0.0);
		T[] controlShift = pointOne;
		for(int i = 0; i < inputAliceControl.length; i++){
			controlSum = flib.add(controlSum, flib.sub(inputAliceControl[i], controlShift));
			T[] addSquare = flib.multiply(flib.sub(inputAliceControl[i], controlShift), flib.sub(inputAliceControl[i], controlShift));
			controlSumOfSquares = flib.add(controlSumOfSquares, addSquare);
		}
		
		for(int i = 0; i < inputBobControl.length; i++){
			controlSum = flib.add(controlSum, flib.sub(inputBobControl[i], controlShift));
			T[] addSquare = flib.multiply(flib.sub(inputBobControl[i], controlShift), flib.sub(inputBobControl[i], controlShift));
			controlSumOfSquares = flib.add(controlSumOfSquares, addSquare);
		}

		T[] caseNum = flib.add(ilib.toSecureFloat(numAliceCase, flib), ilib.toSecureFloat(numBobCase, flib));
		T[] controlNum = flib.add(ilib.toSecureFloat(numAliceControl, flib), ilib.toSecureFloat(numBobControl, flib));
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
		controlMeanAbundance = flib.div(flib.add(controlTotalSum, controlShift), controlNum);		    
		controlVarianceSecondTerm = flib.div(flib.multiply(flib.add(zero,controlTotalSum), flib.add(zero,controlTotalSum)), controlNum);
		controlVariance = flib.div(flib.sub(controlSumOfSquares, controlVarianceSecondTerm), controlNum);

		tUpper = flib.sub(controlMeanAbundance, caseMeanAbundance);
		//T divisorIsZero = flib.eq(caseNum, zero);
		tLowerFirst = flib.div(caseVariance, caseNum);
		//tLowerFirst = cl.mux(zero, tLowerFirst, divisorIsZero);
		//divisorIsZero = flib.eq(controlNum, zero);
		tLowerSecond = flib.div(controlVariance, controlNum);
		//tLowerSecond = cl.mux(zero, tLowerSecond, divisorIsZero);

		tLowerSqrt = flib.sqrt(flib.add(tLowerFirst, tLowerSecond));
		tStat = flib.div(tUpper, tLowerSqrt);

		T[] degreesOfFreedomTop = flib.add(tLowerFirst, tLowerSecond);
		degreesOfFreedomTop = flib.multiply(flib.add(zero,degreesOfFreedomTop), flib.add(zero,degreesOfFreedomTop));

		T[] degreesOfFreedomBottomFirst = flib.div(caseVariance, flib.sub(caseNum, one));
		degreesOfFreedomBottomFirst = flib.multiply(flib.add(zero,degreesOfFreedomBottomFirst), flib.add(zero,degreesOfFreedomBottomFirst));
		degreesOfFreedomBottomFirst = flib.div(degreesOfFreedomBottomFirst, flib.sub(caseNum, one));

		T[] degreesOfFreedomBottomSecond = flib.div(controlVariance, flib.sub(caseNum, one));
		degreesOfFreedomBottomSecond = flib.multiply(flib.add(zero,degreesOfFreedomBottomSecond), flib.add(zero,degreesOfFreedomBottomSecond));
		degreesOfFreedomBottomSecond = flib.div(degreesOfFreedomBottomSecond, flib.sub(controlNum, one));

		T[] degreesOfFreedom = flib.div(degreesOfFreedomTop, flib.add(degreesOfFreedomBottomFirst, degreesOfFreedomBottomSecond));
		/*
		T[] caseNumOnlyOne = ilib.add(numAliceCase, numBobCase);
		T[] controlNumOnlyOne = ilib.add(numAliceControl, numBobControl);
		T[] tStatOnlyOne = ilib.sub(controlSum, caseSum);
		T sendOneOrZero = ilib.eq(tStatOnlyOne, ilib.publicValue(1.0));
		tStatOnlyOne = cl.mux(ilib.publicValue(-1.0), ilib.publicValue(1.0), sendOneOrZero);
		T sendCaseDFOnlyOne = ilib.eq(tStatOnlyOne, ilib.publicValue(1.0));
		T[] dfOnlyOne = cl.mux(flib.sub(caseNumOnlyOne, ilib.publicValue(1.0)), flib.sub(controlNumOnlyOne, ilib.publicValue(1.0)), sendCaseDFOnlyOne);
		*/
		res[0] = tStat;
		res[1] = degreesOfFreedom;
		//res[0] = caseVariance;
		//res[1] = controlVariance;		

		//res[0] = caseNum;
		//res[1] = controlNum;		
		/*
		 	n    = 0
		    sum1 = 0
		    sum2 = 0
		 
		    for x in data:
		        n    = n + 1
		        sum1 = sum1 + x
		 
		    mean = sum1 / n
		 
		    for x in data:
		        sum2 = sum2 + (x - mean)*(x - mean)
		 
		    variance = sum2 / (n - 1)
		    return variance
		    
		        n = 0
    sum1 = 0
    for x in data:
        n = n + 1
        sum1 = sum1 + x
    mean = sum1/n
 
    sum2 = 0
    sum3 = 0
    for x in data:
        sum2 = sum2 + (x - mean)**2
        sum3 = sum3 + (x - mean)
    variance = (sum2 - sum3**2/n)/(n - 1)
    return variance
		 */
		/*
		CircuitLib<T> cl = new CircuitLib<T>(gen);
		T[] numCase = flib.publicValue(0.0);

		T[] caseMean = flib.publicValue(0.0);
		T[] caseSum1 = flib.publicValue(0.0);
		T[] caseSum2 = flib.publicValue(0.0);
		T[] caseSum3 = flib.publicValue(0.0);
		T[] xsubMean = flib.publicValue(0.0);

		T[] correctedX = flib.publicValue(0.0);
		T[] caseVariance = flib.publicValue(0.0);
		
		for(int i =0; i < inputAliceCase.length; i++){
			numCase = flib.add(numCase, one);
			correctedX = flib.add(zero, inputAliceCase[i]);
			caseSum1 = flib.add(caseSum1, correctedX);
		}
		for(int i =0; i < inputBobCase.length; i++){
			numCase = flib.add(numCase, one);
			correctedX = flib.add(zero, inputBobCase[i]);
			caseSum1 = flib.add(caseSum1, correctedX);
		}
		caseMean = flib.div(caseSum1, numCase);
		
		for(int i =0; i < inputAliceCase.length; i++){
			correctedX = flib.add(zero, inputAliceCase[i]);
			caseSum2 = flib.add(caseSum2, flib.multiply(flib.sub(correctedX, caseMean), flib.sub(correctedX, caseMean)));
		}
		for(int i =0; i < inputBobCase.length; i++){
			correctedX = flib.add(zero, inputBobCase[i]);
			caseSum2 = flib.add(caseSum2, flib.multiply(flib.sub(correctedX, caseMean), flib.sub(correctedX, caseMean)));
		}
	    caseVariance = flib.div(caseSum2, flib.sub(numCase, one));

		T[] numControl = flib.publicValue(0.0);

		T[] controlMean = flib.publicValue(0.0);
		T[] controlSum1 = flib.publicValue(0.0);
		T[] controlSum2 = flib.publicValue(0.0);

		T[] controlVariance = flib.publicValue(0.0);
		
		for(int i =0; i < inputAliceControl.length; i++){
			numControl = flib.add(numControl, one);
			correctedX = flib.add(zero, inputAliceControl[i]);
			controlSum1 = flib.add(controlSum1, correctedX);
		}
		for(int i =0; i < inputBobControl.length; i++){
			numControl = flib.add(numControl, one);
			correctedX = flib.add(zero, inputBobControl[i]);
			controlSum1 = flib.add(controlSum1, correctedX);
		}
		controlMean = flib.div(controlSum1, numControl);
		
		for(int i =0; i < inputAliceControl.length; i++){
			correctedX = flib.add(zero, inputAliceControl[i]);
			controlSum2 = flib.add(caseSum2, flib.multiply(flib.sub(correctedX, controlMean), flib.sub(correctedX, controlMean)));
		}
		for(int i =0; i < inputBobControl.length; i++){
			correctedX = flib.add(zero, inputBobControl[i]);
			controlSum2 = flib.add(controlSum2, flib.multiply(flib.sub(correctedX, controlMean), flib.sub(correctedX, controlMean)));
		}
	    controlVariance = flib.div(controlSum2, flib.sub(numControl, one));
		
		T[] differnceOfMeans = flib.sub(caseMean, controlMean);
		T[] caseVarianceDivNum = flib.div(caseVariance, numCase);
		T[] controlVarianceDivNum = flib.div(controlVariance, numControl);
		T[] varianceSum = flib.add(caseVarianceDivNum, controlVarianceDivNum);
		T[] varianceSumSqrt = flib.sqrt(varianceSum);
		T[] tStat = flib.div(differnceOfMeans, varianceSumSqrt);
		T[] varianceSumSquare = flib.multiply(varianceSum, varianceSum);
		T[] caseVarianceDivNumSq = flib.multiply(caseVariance, caseVariance);
		T[] dfBottomFirst = flib.div(caseVarianceDivNumSq, flib.sub(numCase, one));
		T[] controlVarianceDivNumSq = flib.multiply(controlVariance, controlVariance);
		T[] dfBottomSecond = flib.div(controlVarianceDivNumSq, flib.sub(numControl, one));
		T[] degreesOfFreedom = flib.div(varianceSumSquare, flib.add(dfBottomFirst, dfBottomSecond));
		
		T[][] res = gen.newTArray(2, 0);
		*/
		/*
		CircuitLib<T> cl = new CircuitLib<T>(gen);
		T[] numCase = flib.publicValue(0.0);
		T[] deltaCase;

		T[] caseMean = flib.publicValue(0.0);
		T[] caseM2 = flib.publicValue(0.0);
		T[] correctedX = flib.publicValue(0.0);
		T[] caseVariance = flib.publicValue(0.0);
		
		for(int i =0; i < inputAliceCase.length; i++){
			numCase = flib.add(numCase, one);
			correctedX = flib.add(zero, inputAliceCase[i]);
			deltaCase = flib.sub(correctedX, caseMean);
			caseMean = flib.add(caseMean, flib.div(deltaCase, numCase));
			caseM2 = flib.add(caseM2, flib.multiply(deltaCase,flib.sub(correctedX, caseMean)));
		}
		
		for(int i =0; i < inputBobCase.length; i++){
			numCase = flib.add(numCase, one);
			correctedX = flib.add(zero, inputBobCase[i]);
			deltaCase = flib.sub(correctedX, caseMean);
			caseMean = flib.add(caseMean, flib.div(deltaCase, numCase));
			caseM2 = flib.add(caseM2, flib.multiply(deltaCase,flib.sub(correctedX, caseMean)));
		}
		
		caseVariance =  flib.div(caseM2, flib.sub(numCase, one));
		T lessThan2 = flib.leq(numCase, flib.publicValue(1.0));
		caseVariance = cl.mux( caseVariance,zero,lessThan2);
		
		T[] numControl = flib.publicValue(0.0);
		T[] deltaControl;

		T[] controlMean = flib.publicValue(0.0);
		T[] controlM2 = flib.publicValue(0.0);
		T[] controlVariance = flib.publicValue(0.0);
		
		for(int i =0; i < inputAliceControl.length; i++){
			numControl = flib.add(numControl, one);
			correctedX = flib.add(zero, inputAliceControl[i]);
			deltaControl = flib.sub(correctedX, controlMean);
			controlMean = flib.add(controlMean, flib.div(deltaControl, numControl));
			controlM2 = flib.add(controlM2, flib.multiply(deltaControl,flib.sub(correctedX, controlMean)));
		}
		
		for(int i =0; i < inputBobControl.length; i++){
			numControl = flib.add(numControl, one);
			correctedX = flib.add(zero, inputBobControl[i]);
			deltaControl = flib.sub(correctedX, controlMean);
			controlMean = flib.add(controlMean, flib.div(deltaControl, numControl));
			controlM2 = flib.add(controlM2, flib.multiply(deltaControl,flib.sub(correctedX, controlMean)));
		}
		
		controlVariance =  flib.div(controlM2, flib.sub(numControl, one));
		lessThan2 = flib.leq(numControl, flib.publicValue(1.0));
		controlVariance = cl.mux( controlVariance,zero, lessThan2);
		
		T[] differnceOfMeans = flib.sub(caseMean, controlMean);
		T[] caseVarianceDivNum = flib.div(caseVariance, numCase);
		T[] controlVarianceDivNum = flib.div(controlVariance, numControl);
		T[] varianceSum = flib.add(caseVarianceDivNum, controlVarianceDivNum);
		T[] varianceSumSqrt = flib.sqrt(varianceSum);
		T[] tStat = flib.div(differnceOfMeans, varianceSumSqrt);
		T[] varianceSumSquare = flib.multiply(varianceSum, varianceSum);
		T[] caseVarianceDivNumSq = flib.multiply(caseVariance, caseVariance);
		T[] dfBottomFirst = flib.div(caseVarianceDivNumSq, flib.sub(numCase, one));
		T[] controlVarianceDivNumSq = flib.multiply(controlVariance, controlVariance);
		T[] dfBottomSecond = flib.div(controlVarianceDivNumSq, flib.sub(numControl, one));
		T[] degreesOfFreedom = flib.div(varianceSumSquare, flib.add(dfBottomFirst, dfBottomSecond));
		T[][] res = gen.newTArray(2, 0);
		*/
		//res[0] = caseVariance;
		//res[1] = controlVariance;
		
		return res;
	}
	
	public static class Generator<T> extends GenRunnable<T> {
		T[][] inputBobCase;
		T[][] inputAliceCase;
		T[][] inputAliceControl;
		T[][] inputBobControl;
		T[][] inputCounters;
		T[][][] in;
		double[][] caseInput;
		double[][] controlInput;
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
			FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
			T[] l = flib.publicValue(0.0);
			caseInput = PrepareDataDANaive.readFile(cmd.getOptionValue("s"));
			controlInput = PrepareDataDANaive.readFile(cmd.getOptionValue("t"));
			
			NumFeatures = caseInput.length;
			
			in = gen.newTArray(NumFeatures, 2, 0);

			EvaCaseSamples = new int[NumFeatures];
			GenCaseSamples = new int[NumFeatures];
			EvaControlSamples = new int[NumFeatures];
			GenControlSamples = new int[NumFeatures];
			
			for(int i =0; i < NumFeatures; i++){
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
			FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
			T[] l = flib.publicValue(0.0);
			T result;
			int caseCounter = 0;
			int controlCounter = 0;
			
			boolean[] filteredFeatures = new boolean[NumFeatures];
			for(int i =0; i < NumFeatures; i++){
				numAliceCase = gen.inputOfAlice(Utils.fromInt(GenCaseSamples[i], 32));
				numBobCase = gen.inputOfBob(new boolean[32]);
				numAliceControl = gen.inputOfAlice(Utils.fromInt(GenControlSamples[i], 32));
				numBobControl = gen.inputOfBob(new boolean[32]);
				result = filter(gen, numAliceCase, numBobCase, numAliceControl, numBobControl);
				filteredFeatures[i] = gen.outputToAlice(result);
			}

			
			numAliceCase = gen.inputOfAlice(Utils.fromInt(caseInput[0].length, 32));
			numBobCase = gen.inputOfBob(new boolean[32]);
			numAliceControl = gen.inputOfAlice(Utils.fromInt(controlInput[0].length, 32));
			numBobControl = gen.inputOfBob(new boolean[32]);
			for(int i =0; i < NumFeatures; i++){
				inputAliceCase = gen.newTArray(GenCaseSamples[i], 0);
				inputBobCase = gen.newTArray(EvaCaseSamples[i], 0);
				inputAliceControl = gen.newTArray(GenControlSamples[i], 0);
				inputBobControl = gen.newTArray(EvaControlSamples[i], 0);

				caseCounter = 0;
				controlCounter = 0;
				for(int j =0; j < caseInput[0].length; j++){
					if(caseInput[i][j] > 0.0){
						inputAliceCase[caseCounter++] = gen.inputOfAlice(Utils.fromFloat(caseInput[i][j], PLength, VLength));
					}
				}
				gen.flush();
				for(int j =0; j < EvaCaseSamples[i]; j++){
					inputBobCase[j] = gen.inputOfBob(new boolean[l.length]);
				}
				gen.flush();
				for(int j =0; j < controlInput[0].length; j++){
					if(controlInput[i][j] > 0.0){
						inputAliceControl[controlCounter++] = gen.inputOfAlice(Utils.fromFloat(controlInput[i][j], PLength, VLength));
					}
				}
				gen.flush();
				for(int j =0; j < EvaControlSamples[i]; j++){
					inputBobControl[j] = gen.inputOfBob(new boolean[l.length]);
				}
				gen.flush();
				if(!(filteredFeatures[i])){
					in[i] = onlyOne(gen, inputAliceCase, inputBobCase, inputAliceControl, inputBobControl, 
							numAliceCase, numBobCase, numAliceControl, numBobControl);
					continue;
				}
				in[i] = compute(gen, inputAliceCase, inputBobCase, inputAliceControl, inputBobControl, 
						numAliceCase, numBobCase, numAliceControl, numBobControl);
			}
		}
		@Override
		public void prepareOutput(CompEnv<T> gen) {
			FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
			for(int j = 0; j < in.length; j++){
				double tStat = flib.outputToAlice(in[j][0]);
				double df = flib.outputToAlice(in[j][1]);
				
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
		T[][] inputBobCase;
		T[][] inputAliceCase;
		T[][] inputCounters;
		T[][] inputAliceControl;
		T[][] inputBobControl;
		T[] scResult;
		T[][][] in;

		T[] numAliceCase;
		T[] numBobCase;
		T[] numAliceControl;
		T[] numBobControl;
		double[][] caseInput;
		double[][] controlInput;

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

			FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
			T[] l = flib.publicValue(0.0);
			caseInput = PrepareDataDANaive.readFile(cmd.getOptionValue("s"));
			controlInput = PrepareDataDANaive.readFile(cmd.getOptionValue("t"));

			NumFeatures = caseInput.length;
			in = gen.newTArray(NumFeatures, 2, 0);
			
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

			FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
			T[] l = flib.publicValue(0.0);
			int caseCounter = 0;
			int controlCounter = 0;
			T result;
			boolean[] filteredFeatures = new boolean[NumFeatures];
			for(int i =0; i < NumFeatures; i++){
				numAliceCase = gen.inputOfAlice(new boolean[32]);
				numBobCase = gen.inputOfBob(Utils.fromInt(EvaCaseSamples[i], 32));
				numAliceControl = gen.inputOfAlice(new boolean[32]);
				numBobControl = gen.inputOfBob(Utils.fromInt(EvaControlSamples[i], 32));
				result = filter(gen, numAliceCase, numBobCase, numAliceControl, numBobControl);
				filteredFeatures[i] = gen.outputToAlice(result);
			}
			numAliceCase = gen.inputOfAlice(new boolean[32]);
			numBobCase = gen.inputOfBob(Utils.fromInt(caseInput[0].length, 32));
			numAliceControl = gen.inputOfAlice(new boolean[32]);
			numBobControl = gen.inputOfBob(Utils.fromInt(controlInput[0].length, 32));
			for(int i =0; i < NumFeatures; i++){

				inputAliceCase = gen.newTArray(GenCaseSamples[i], 0);
				inputBobCase = gen.newTArray(EvaCaseSamples[i], 0);
				inputAliceControl = gen.newTArray(GenControlSamples[i], 0);
				inputBobControl = gen.newTArray(EvaControlSamples[i], 0);

				caseCounter = 0;
				controlCounter = 0;
				for(int j =0; j < GenCaseSamples[i]; j++){
					inputAliceCase[j] = gen.inputOfAlice(new boolean[l.length]);
				}
				gen.flush();
				for(int j =0; j < caseInput[0].length; j++){
					if(caseInput[i][j] > 0.0){
						inputBobCase[caseCounter++] = gen.inputOfBob(Utils.fromFloat(caseInput[i][j], PLength, VLength));
					}
				}
				gen.flush();
				for(int j =0; j < GenControlSamples[i]; j++){
					inputAliceControl[j] = gen.inputOfAlice(new boolean[l.length]);
				}
				gen.flush();
				for(int j =0; j < controlInput[0].length; j++){
					if(controlInput[i][j] > 0.0){
						inputBobControl[controlCounter++] = gen.inputOfBob(Utils.fromFloat(controlInput[i][j], PLength, VLength));
					}
				}
				gen.flush();
				if(!(filteredFeatures[i])){
					in[i] = onlyOne(gen, inputAliceCase, inputBobCase, inputAliceControl, inputBobControl, 
							numAliceCase, numBobCase, numAliceControl, numBobControl);
					continue;
				}
				in[i] = compute(gen, inputAliceCase, inputBobCase, inputAliceControl, inputBobControl, 
						numAliceCase, numBobCase, numAliceControl, numBobControl);
			}
		}
		
		@Override
		public void prepareOutput(CompEnv<T> gen) {
			FloatLib<T> flib = new FloatLib<T>(gen, PLength, VLength);
			for(int i = 0; i < in.length; i++){
				flib.outputToAlice(in[i][0]);
				flib.outputToAlice(in[i][1]);
			}
		}
				
	}
}