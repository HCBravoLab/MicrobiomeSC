package microbiomeOddsRatio;



public class Statistics {
	public String presenceAbsenceVector;
	public char g1;
	public char g2;
	public int numOfPresent;
	public int totalNum;
	
	public static char[] LookupTable = {'A', 'G', 'C', 'T'};
	
	public static int getindex(char c){
		for(int i = 0; i < 4; ++i)
			if(c == LookupTable[i])
				return i;
		try {
			throw new Exception("DNA type not supported!");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(1);
		}
		return -1;
	}
}
