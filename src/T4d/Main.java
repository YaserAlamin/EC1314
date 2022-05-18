package T4d;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.util.*;

public class Main {

	public static void main(String[] args) {
		
		Scanner sc = new Scanner(System.in);
		
		// Read input constants
		int N = 320;						// N - number of elements (must be a multiple of s)
		int l = 100;					// l - lenght of the string representation
		int s = 2;						// s - tournament size
		double pm = (double)0.2;		// Pm - probability of mutation
		double pc = 0.5;				// Pc - probability of crossover
		double r = 0.5;					// r - fraction of the worst elements from the 
										// old population that is going to be replaced.
										// The actual number of individuals should be round up.
		int g = 100;						// g - number of generations
		int experiments = 100;	
		
		ArrayList<Integer> Ns = new ArrayList<Integer>();
		Ns.add(10);
		Ns.add(40);
		Ns.add(320);
		Iterator<Integer> iterator_Ns = Ns.iterator();
		
		ArrayList<Double> rs = new ArrayList<Double>();
		rs.add(0.5);
		rs.add(1.0);
		Iterator<Double> iterator_rs = rs.iterator();
		
		ArrayList<Double> pms = new ArrayList<Double>();
		pms.add(0.0);
		pms.add(1.0/l);
		pms.add(10.0/l);
		Iterator<Double> iterator_pms = pms.iterator();
		
		ArrayList<Double> pcs = new ArrayList<Double>();
		pcs.add(0.0);
		pcs.add(0.5);
		pcs.add(1.0);
		Iterator<Double> iterator_pcs = pcs.iterator();
		
		/* SEL_method
		 * 0 -> Roulette Wheel
		 * 1 -> SUS
		 * 2 -> Tournament with s = 2
		 * 3 -> Tournament with s = 4
		 * 4 -> Truncation with 10% cutoff
		 * 5 -> Truncation with 50% cutoff
		 */
		
		/* CO_method
		 * 0 -> Uniform Crossover
		 * 1 -> One-point Crossover
		 * 2 -> Two-point Crossover
		 */
		DecimalFormat df = new DecimalFormat("#0.000");
		
		int total_iterations = Ns.size()*rs.size()*pms.size()*pcs.size()*6*3*experiments;
		int total_iterations_so_far = 0;
		
		try {

			// BufferedWriter to write on the constraint file
			File file =  new File("D:\\Universidade Mestrado\\Evolutionary Computing\\Lab4_results.txt");
			FileOutputStream fileOS = new FileOutputStream(file);
			OutputStreamWriter oSW = new OutputStreamWriter(fileOS);
			BufferedWriter out = new BufferedWriter(oSW);

			out.append("Elements with length of " + l + " on the ONEMAX problem." + System.getProperty("line.separator") +
					g + " generations per experience." + System.getProperty("line.separator") + 
					experiments + " experiments with 3 pre-runs." + System.getProperty("line.separator") + 
					"Total number of experiments: " + total_iterations + System.getProperty("line.separator")
					);
			
			out.append("Selection_Method " +
					"CrossOver_method " +
					"N " +
					//multirun_Data.getConfiguration().getL() + " " +
					//multirun_Data.getConfiguration().getS() + " " +
					"Pm " +
					"Pc " +
					"R " +

					"Avg_max " +
					"Avg_avg " +
					"Avg_min " +
					"Avg_stdDeviation " +
					"Max_max " +
					"Min_min " + 
			
					"Total_time_in_seconds " +
					"Avg_time_in_seconds " +
					"Min_time_in_seconds " +
					"Max_time_in_seconds " +
					"StdDeviation_time_in_seconds " +
			
					"Total_number_of_convergence " +
					"Average_number_of_convergence " +
					"Average_convergence_percentage " +
					"Min_iterations_to_converge " +
					"Avg_iterations_to_converge " +
					"Max_iterations_to_converge" + 	
					System.getProperty("line.separator")
					);

			
			while(iterator_Ns.hasNext()){
				N = iterator_Ns.next();
				while(iterator_rs.hasNext()){
					r = iterator_rs.next();
					while(iterator_pms.hasNext()){
						pm = iterator_pms.next();
						while(iterator_pcs.hasNext()){
							pc = iterator_pcs.next();
							for(int sE=0;sE<=5;sE++)
								for(int cO=0;cO<=2;cO++){
									//show how many iterations are 
									int iterations = 0;

									ArrayList<GA_Data> multirun_list = new ArrayList<GA_Data>();
									for(int i=0; i<experiments+2;i++){
										if(i<2)
											System.out.print("pre-caching..");
										else
											System.out.print(++total_iterations_so_far + " out of " + total_iterations + ". Locally " + ++iterations + "...");
										GA_Data data = new GA_Data();
										data = CustomGeneticAlghorithm.geneticONEMAX(sE, cO, N, l, s, pm, pc, r, g);
										if(i<2)
											continue;
										//data.print(false);
										multirun_list.add(data);

										System.out.println("done");
									}

									GA_multirun_Data multirun_Data = new GA_multirun_Data();
									multirun_Data.calculate_data(multirun_list);
									//multirun_Data.print(true);



									out.append(
											"\"" + multirun_Data.getConfiguration().getSEL_method_name() + "\" " +
													"\"" + multirun_Data.getConfiguration().getCO_method_name() + "\" " +
													multirun_Data.getConfiguration().getN() + " " +
													//multirun_Data.getConfiguration().getL() + " " +
													//multirun_Data.getConfiguration().getS() + " " +
													multirun_Data.getConfiguration().getPm() + " " +
													multirun_Data.getConfiguration().getPc() + " " +
													multirun_Data.getConfiguration().getR() + " " +

													df.format(multirun_Data.getAvg_max()) + " " +
													df.format(multirun_Data.getAvg_avg()) + " " +
													df.format(multirun_Data.getAvg_min()) + " " +
													df.format(multirun_Data.getAvg_stdDeviation()) + " " +
													df.format(multirun_Data.getMax_max()) + " " +
													df.format(multirun_Data.getMin_min()) + " " + 

													df.format(multirun_Data.getTotal_time_in_seconds()) + " " +
													df.format(multirun_Data.getAvg_time_in_seconds()) + " " +
													df.format(multirun_Data.getMin_time_in_seconds()) + " " +
													df.format(multirun_Data.getMax_time_in_seconds()) + " " +
													df.format(multirun_Data.getStdDeviation_time_in_seconds()) + " " +

													df.format(multirun_Data.getTotal_number_of_convergence()) + " " +
													df.format(multirun_Data.getAverage_number_of_convergence()) + " " +
													df.format(multirun_Data.getAverage_convergence_percentage()) + " " +
													df.format(multirun_Data.getMin_iterations_to_converge()) + " " +
													df.format(multirun_Data.getAvg_iterations_to_converge()) + " " +
													df.format(multirun_Data.getMax_iterations_to_converge())		
											);
									out.append(System.getProperty("line.separator"));
									out.flush();
								}	
							}
						iterator_pcs = pcs.iterator();
						}
					iterator_pms = pms.iterator();
					}
				iterator_rs = rs.iterator();
				}

	
			out.flush();
			out.close();				

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

}
class GA_multirun_Data{
	//experiment identification
	private GA_Configuration configuration;
	private int number_iterations;
	
	//fitness related data
	private double avg_max;
	private double avg_avg;
	private double avg_min;
	private double avg_stdDeviation;
	
	private double max_max;
	private double min_min;
	
	//alghorithm performance and execution time
	private double total_time_in_seconds;
	private double avg_time_in_seconds;
	private double min_time_in_seconds;
	private double max_time_in_seconds;
	private double stdDeviation_time_in_seconds;
	
	//alghorithm convergence
	private double total_number_of_convergence;
	private double average_number_of_convergence;
	private double average_convergence_percentage;
	
	private double min_iterations_to_converge;
	private double avg_iterations_to_converge;
	private double max_iterations_to_converge;
	
	public void calculate_data(ArrayList<GA_Data> data){

		this.configuration = data.get(0).getConfiguration();
		this.number_iterations = data.size();
		
		double MAXNUM = 9999999;
		
		//fitness statistics
		double avg_max=0;
		double avg_avg=0;
		double avg_min=0;
		double avg_stdDeviation=0;

		double max_max=0;
		double min_min=MAXNUM;

		//alghorithm performance and execution time
		double total_time_in_seconds=0;
		double avg_time_in_seconds=0;
		double min_time_in_seconds=MAXNUM;
		double max_time_in_seconds=0;
		double stdDeviation_time_in_seconds;
		ArrayList<Double> stdDeviation_time_in_seconds_aux = new ArrayList<Double>();

		//alghorithm convergence
		double total_number_of_convergence=0;
		double average_number_of_convergence=0;
		double average_convergence_percentage;

		double min_iterations_to_converge=MAXNUM;
		double avg_iterations_to_converge=0;
		double max_iterations_to_converge=0;


		Iterator<GA_Data> iterator = data.iterator();
		GA_Data elem;
		while(iterator.hasNext()){
			
			elem = iterator.next();
			
			//fitness statistics		
			avg_max += elem.getMax();
			avg_avg += elem.getAvg();
			avg_min += elem.getMin();
			avg_stdDeviation += elem.getStdDeviation();
			
			if(elem.getMax() > max_max)
				max_max = elem.getMax();
			if(elem.getMin() < min_min)
				min_min = elem.getMin();
			
			//alghorithm performance and execution time
			total_time_in_seconds += elem.getTime_in_seconds();
			
			if(elem.getTime_in_seconds() < min_time_in_seconds)
				min_time_in_seconds = elem.getTime_in_seconds();
			if(elem.getTime_in_seconds() > max_time_in_seconds)
				max_time_in_seconds = elem.getTime_in_seconds();
			
			stdDeviation_time_in_seconds_aux.add(elem.getTime_in_seconds());
			
			//alghorithm convergence
			if(elem.isConvergence()){
				total_number_of_convergence++;
				if(elem.getIterations_to_converge() < min_iterations_to_converge)
					min_iterations_to_converge = elem.getIterations_to_converge();
				if(elem.getIterations_to_converge() > max_iterations_to_converge)
					max_iterations_to_converge = elem.getIterations_to_converge();
				avg_iterations_to_converge += elem.getIterations_to_converge();
			}	
		}
		
		//fitness statistics
		avg_max /= data.size();
		avg_avg /= data.size();
		avg_min /= data.size();
		avg_stdDeviation /= data.size();
		
		//alghorithm performance and execution time
		avg_time_in_seconds = total_time_in_seconds / data.size();
		stdDeviation_time_in_seconds = Statistics.standardDeviation(stdDeviation_time_in_seconds_aux);
		
		//alghorithm convergence
		average_number_of_convergence = total_number_of_convergence/data.size();
		average_convergence_percentage = 100*total_number_of_convergence/data.size();
		avg_iterations_to_converge /= total_number_of_convergence;
		
		// FINAL ATTRIBUTION
		this.avg_max=avg_max;
		this.avg_avg=avg_avg;
		this.avg_min=avg_min;
		this.avg_stdDeviation=avg_stdDeviation;

		this.max_max=max_max;
		this.min_min=min_min;

		//alghorithm performance and execution time
		this.total_time_in_seconds=total_time_in_seconds;
		this.avg_time_in_seconds=avg_time_in_seconds;
		this.min_time_in_seconds=min_time_in_seconds;
		this.max_time_in_seconds=max_time_in_seconds;
		this.stdDeviation_time_in_seconds=stdDeviation_time_in_seconds;

		//alghorithm convergence
		this.total_number_of_convergence=total_number_of_convergence;
		this.average_number_of_convergence=average_number_of_convergence;
		this.average_convergence_percentage=average_convergence_percentage;

		this.min_iterations_to_converge=min_iterations_to_converge;
		this.avg_iterations_to_converge=avg_iterations_to_converge;
		this.max_iterations_to_converge=max_iterations_to_converge;
	}
	
	public void print(boolean showConfiguration){	
		
		if(showConfiguration){
			System.out.println();
			System.out.println("N  (population size): " + configuration.getN());
			System.out.println("l  (element length): " + configuration.getL());
			//System.out.print("s (tournament selection): " + configuration.getS());
			System.out.println("Pm (mutation probability): " + configuration.getPm());
			System.out.println("Pc (crossover probability): " + configuration.getPc());
			System.out.println("r  (replacement factor: " + configuration.getR());
			System.out.println("Selection method: " + configuration.getSEL_method_name());
			System.out.println("CrossOver method: " + configuration.getCO_method_name());
			System.out.println();
			
		}
		
		System.out.println("EXPERIMENT IDENTIFICATION");
		System.out.println("- configuration: ");
		System.out.println("- # iterations: " + number_iterations);
		System.out.println();
		
		System.out.println("FITNESS RELATED DATA");
		System.out.println("- avg_Max: " + getAvg_max());
		System.out.println("- avg_Avg: " + getAvg_avg());
		System.out.println("- avg_Min: " + getAvg_min());
		System.out.println("- avg_stdDeviation: " + getAvg_stdDeviation());
		System.out.println();
		System.out.println("- max_Max: " + getMax_max());
		System.out.println("- min_Min: " + getMin_min());
		System.out.println();
		
		System.out.println("ALGORITHM PERFORMANCE AND EXECUTION TIME");
		System.out.println("- total_time_in_seconds: " + this.getTotal_time_in_seconds());
		System.out.println("- avg_time_in_seconds: " + getAvg_time_in_seconds());
		System.out.println("- min_time_in_seconds: " + getMin_time_in_seconds());
		System.out.println("- max_time_in_seconds: " + getMax_time_in_seconds());
		System.out.println("- stdDeviation_time_in_seconds: " + getStdDeviation_time_in_seconds());
		System.out.println();
		
		System.out.println("ALGORITHM CONVERGENCE ");
		System.out.println("- total_number_of_convergence: " + getTotal_number_of_convergence());
		System.out.println("- average_number_of_convergence: " + getAverage_number_of_convergence());
		System.out.println("- average_convergence_percentage: " + getAverage_convergence_percentage());
		
		System.out.println("- min_iterations_to_converge: " + getMin_iterations_to_converge());
		System.out.println("- avg_iterations_to_converge: " + getAvg_iterations_to_converge());
		System.out.println("- max_iterations_to_converge: " + getMax_iterations_to_converge());
		
	}
	public GA_Configuration getConfiguration() {
		return configuration;
	}
	public void setConfiguration(GA_Configuration configuration) {
		this.configuration = configuration;
	}
	public int getNumber_iterations() {
		return number_iterations;
	}
	public void setNumber_iterations(int number_iterations) {
		this.number_iterations = number_iterations;
	}
	public double getAvg_max() {
		return avg_max;
	}
	public void setAvg_max(double avg_max) {
		this.avg_max = avg_max;
	}
	public double getAvg_avg() {
		return avg_avg;
	}
	public void setAvg_avg(double avg_avg) {
		this.avg_avg = avg_avg;
	}
	public double getAvg_min() {
		return avg_min;
	}
	public void setAvg_min(double avg_min) {
		this.avg_min = avg_min;
	}
	public double getAvg_stdDeviation() {
		return avg_stdDeviation;
	}
	public void setAvg_stdDeviation(double avg_stdDeviation) {
		this.avg_stdDeviation = avg_stdDeviation;
	}
	public double getMax_max() {
		return max_max;
	}
	public void setMax_max(double max_max) {
		this.max_max = max_max;
	}
	public double getMin_min() {
		return min_min;
	}
	public void setMin_min(double min_min) {
		this.min_min = min_min;
	}
	public double getTotal_time_in_seconds() {
		return total_time_in_seconds;
	}
	public void setTotal_time_in_seconds(double total_time_in_seconds) {
		this.total_time_in_seconds = total_time_in_seconds;
	}
	public double getAvg_time_in_seconds() {
		return avg_time_in_seconds;
	}
	public void setAvg_time_in_seconds(double avg_time_in_seconds) {
		this.avg_time_in_seconds = avg_time_in_seconds;
	}
	public double getMin_time_in_seconds() {
		return min_time_in_seconds;
	}
	public void setMin_time_in_seconds(double min_time_in_seconds) {
		this.min_time_in_seconds = min_time_in_seconds;
	}
	public double getMax_time_in_seconds() {
		return max_time_in_seconds;
	}
	public void setMax_time_in_seconds(double max_time_in_seconds) {
		this.max_time_in_seconds = max_time_in_seconds;
	}
	public double getStdDeviation_time_in_seconds() {
		return stdDeviation_time_in_seconds;
	}
	public void setStdDeviation_time_in_seconds(double stdDeviation_time_in_seconds) {
		this.stdDeviation_time_in_seconds = stdDeviation_time_in_seconds;
	}
	public double getTotal_number_of_convergence() {
		return total_number_of_convergence;
	}
	public void setTotal_number_of_convergence(double total_number_of_convergence) {
		this.total_number_of_convergence = total_number_of_convergence;
	}
	public double getAverage_number_of_convergence() {
		return average_number_of_convergence;
	}
	public void setAverage_number_of_convergence(
			double average_number_of_convergence) {
		this.average_number_of_convergence = average_number_of_convergence;
	}
	public double getAverage_convergence_percentage() {
		return average_convergence_percentage;
	}
	public void setAverage_convergence_percentage(
			double average_convergence_percentage) {
		this.average_convergence_percentage = average_convergence_percentage;
	}
	public double getMin_iterations_to_converge() {
		return min_iterations_to_converge;
	}
	public void setMin_iterations_to_converge(double min_iterations_to_converge) {
		this.min_iterations_to_converge = min_iterations_to_converge;
	}
	public double getAvg_iterations_to_converge() {
		return avg_iterations_to_converge;
	}
	public void setAvg_iterations_to_converge(double avg_iterations_to_converge) {
		this.avg_iterations_to_converge = avg_iterations_to_converge;
	}
	public double getMax_iterations_to_converge() {
		return max_iterations_to_converge;
	}
	public void setMax_iterations_to_converge(double max_iterations_to_converge) {
		this.max_iterations_to_converge = max_iterations_to_converge;
	}

	

	
	
}
class GA_Data{
	//Configuration used, each element represents one initialization parameter:
	private GA_Configuration configuration;
	
	//min, max and avg fitness of the population
	private double max;
	private double avg;
	private double min;
	private double stdDeviation;
	
	//time executing the iterations until convergeance or until iteration limit
	private double time_in_seconds;
	
	//Did the alghorithm converged before reaching iteration limit?
	private boolean convergence;
	
	//If the alghorithm converged, how many iterations did it took?
	private int iterations_to_converge;

	public double getMin() {
		return min;
	}
	public void setMin(double min) {
		this.min = min;
	}
	public double getMax() {
		return max;
	}
	public void setMax(double max) {
		this.max = max;
	}
	public double getAvg() {
		return avg;
	}
	public void setAvg(double avg) {
		this.avg = avg;
	}
	public double getTime_in_seconds() {
		return time_in_seconds;
	}
	public void setTime_in_seconds(double time_in_seconds) {
		this.time_in_seconds = time_in_seconds;
	}
	public boolean isConvergence() {
		return convergence;
	}
	public void setConvergence(boolean convergence) {
		this.convergence = convergence;
	}
	public int getIterations_to_converge() {
		return iterations_to_converge;
	}
	public void setIterations_to_converge(int iterations_to_converge) {
		this.iterations_to_converge = iterations_to_converge;
	}
	public GA_Configuration getConfiguration() {
		return configuration;
	}
	public void setConfiguration(GA_Configuration configuration) {
		this.configuration = configuration;
	}
	public double getStdDeviation() {
		return stdDeviation;
	}
	public void setStdDeviation(double stdDeviation) {
		this.stdDeviation = stdDeviation;
	}
	
	public void print(boolean showConfiguration){	
		System.out.println("Max: " + getMax());
		System.out.println("Avg: " + getAvg());
		System.out.println("Min: " + getMin());
		System.out.println("StdDeviation: " + getStdDeviation());
		System.out.println();
		System.out.println("Time: " + getTime_in_seconds());
		System.out.println("Convergence: " + (isConvergence()==true ? "Yes!" : "Nop... stopped at iteration: " + getConfiguration().getG()));
		if(isConvergence())
			System.out.println("Iterations to Converge:" + getIterations_to_converge());
	
		if(showConfiguration){
			System.out.println();
			System.out.println("N  (population size): " + configuration.getN());
			System.out.println("l  (element length): " + configuration.getL());
			//System.out.print("s (tournament selection): " + configuration.getS());
			System.out.println("Pm (mutation probability): " + configuration.getPm());
			System.out.println("Pc (crossover probability): " + configuration.getPc());
			System.out.println("r  (replacement factor: " + configuration.getR());
			System.out.println("Selection method: " + configuration.getSEL_method_name());
			System.out.println("CrossOver method: " + configuration.getCO_method_name());
			
		}
	}
	
}

class GA_Configuration{
	/**1st - N - number of elements (must be a multiple of s)
	 * 2nd - l - lenght of the string representation
	 * 3rd - s - tournament size
	 * 4th - Pm - probability of mutation
	 * 5th - Pc - probability of crossover
	 * 6th - r - fraction of the worst elements from the 
		   old population that is going to be replaced.
		   The actual number of individuals should be round up.
	 * 7th - g - number of generations
			
	 * 8th - SEL_method
	 * 9th - CO_method
	**/
	private int N;
	private int l;
	private int s;
	private double Pm;
	private double Pc;
	private double r;
	private int g;
	private int SEL_method;
	private int CO_method;
	
	public int getN() {
		return N;
	}
	public void setN(int n) {
		N = n;
	}
	public int getL() {
		return l;
	}
	public void setL(int l) {
		this.l = l;
	}
	public int getS() {
		return s;
	}
	public void setS(int s) {
		this.s = s;
	}
	public double getPm() {
		return Pm;
	}
	public void setPm(double pm) {
		Pm = pm;
	}
	public double getPc() {
		return Pc;
	}
	public void setPc(double pc) {
		Pc = pc;
	}
	public double getR() {
		return r;
	}
	public void setR(double r) {
		this.r = r;
	}
	public int getG() {
		return g;
	}
	public void setG(int g) {
		this.g = g;
	}
	public int getSEL_method() {
		return SEL_method;
	}
	public void setSEL_method(int SEL_method) {
		this.SEL_method = SEL_method;
	}
	public String getSEL_method_name(){
		// SEL_method
		/*
		 * 0 -> Roulette Wheel
		 * 1 -> SUS
		 * 2 -> Tournament with s = 2
		 * 3 -> Tournament with s = 4
		 * 4 -> Truncation with 10% cutoff
		 * 5 -> Truncation with 50% cutoff
		 */
		switch (getSEL_method()){
			case 0: return "Roulette Wheel";
			case 1: return "SUS";
			case 2: return "Tournament with s = 2";
			case 3: return "Tournament with s = 4";
			case 4: return "Truncation with 10% cutoff";
			case 5: return "Truncation with 50% cutoff";
			default: return "not specified";
		}
	}
	public int getCO_method() {
		return CO_method;
	}
	public void setCO_method(int CO_method) {
		this.CO_method = CO_method;
	}
	public String getCO_method_name(){
		// CO_method
		/*
		 * 0 -> Uniform Crossover
		 * 1 -> One-point Crossover
		 * 2 -> Two-point Crossover
		 */
		switch (getCO_method()){
			case 0: return "Uniform Crossover";
			case 1: return "One-point Crossover";
			case 2: return "Two-point Crossover";
			default: return "not specified";
		}
	}



}

class CustomGeneticAlghorithm{
	
	/**
	 *  Iteratively Simulates one generation of a genetic algorithm operating on the onemax problem.
	 *  
	 * 	@param N -  number of elements
	   	@param l -  lenght of the string representation
	   	@param s -  tournament size
		@param pm - probability of mutation
		@param pc - probability of crossover
		@param r -  fraction of the worst elements from the 
		    		old population that is going to be replaced.
					The actual number of individuals should be round up.
		@param M - the amount of random numbers given. Uniformly generated in [0..1)
		@param u - list with random numbers	
	 * @return 
	 * @return 
	*/
	public static GA_Data geneticONEMAX(int SEL_method, int CO_method, int N, int l, int s, double pm, double pc, double r, int g){
		
		//DecimalFormat df = new DecimalFormat("#0.00");
		
		boolean convergence = false;
		int iterations_to_converge = -1;
		
		//Stores the configuration
		GA_Configuration config = new GA_Configuration();
		config.setN(N);
		config.setL(l);
		config.setS(s);
		config.setPm(pm);
		config.setPc(pc);
		config.setPm(pm);
		config.setR(r);
		config.setG(g);
		config.setSEL_method(SEL_method);
		config.setCO_method(CO_method);
				
		long startTime = System.nanoTime();
		/** BEGINNING OF THE ALGHORITHM **/
		
		// INITIAL POPULATION GENERATION	
		// generates the population and calculates the fitness for each element
		ArrayList<String> population = BinaryString.generatePopulation(N, l);
		ArrayList<BinaryString> populationBS = BinaryString.onemaxPopulationFitness(population);

		// prints the status of the initial population
		double max = BinaryString.fitness_max(populationBS).getFitness();
		double avg;// = BinaryString.fitness_average(populationBS);
		double min;// = BinaryString.fitness_min(populationBS).getFitness();
		double stdDeviation;// = BinaryString.fitness_stdDeviation(populationBS);

		
		//System.out.println("0: " + df.format(max) + " " + df.format(avg) + " " + df.format(min));
		
		for(int i=0; i<g;i++){
			populationBS = geneticONEMAX_iteration(SEL_method,CO_method,populationBS, N, l, s, pm, pc, r);
			
			max = BinaryString.fitness_max(populationBS).getFitness();
			
			if(max == l){
				convergence = true;
				iterations_to_converge = i;
				break;
			}
			//System.out.println(i+1 + ": " + df.format(max) + " " + df.format(avg) + " " + df.format(min));
		}
		
		/** END OF THE ALGHORITHM **/
		
		//Stores the alghorithm iteration data
		GA_Data data = new GA_Data();
		
		//Stores configuration
		data.setConfiguration(config);
		
		//Stores running time
		long estimatedTime = System.nanoTime() - startTime;
		double elapsedTime = (double)estimatedTime/1000000000;
		//System.out.println("Time:" + elapsedTime);
		data.setTime_in_seconds(elapsedTime);
		
		//Store min, max and avg
		max = BinaryString.fitness_max(populationBS).getFitness();
		avg = BinaryString.fitness_average(populationBS);
		min = BinaryString.fitness_min(populationBS).getFitness();
		stdDeviation = BinaryString.fitness_stdDeviation(populationBS);
		
		data.setMin(min);
		data.setMax(max);
		data.setAvg(avg);
		data.setStdDeviation(stdDeviation);
		
		//Store convergence confirmation and the number of iterations necessary (stores -1 if convergence was not achieved)
		data.setConvergence(convergence);
		data.setIterations_to_converge(iterations_to_converge);
		
		return data;
	}
	
	private static ArrayList<BinaryString> geneticONEMAX_iteration(int SEL_method, int CO_method, ArrayList<BinaryString> population, int N, int l, int s, double pm, double pc, double r){
		
		//SELECTION	
		/*
		 * 0 -> Roulette Wheel
		 * 1 -> SUS
		 * 2 -> Tournament with s = 2
		 * 3 -> Tournament with s = 4
		 * 4 -> Truncation with 10% cutoff
		 * 5 -> Truncation with 50% cutoff
		 */
		ArrayList<BinaryString> selected;
		
		switch (SEL_method) {
        case 0:	selected = RouletteWheel.selectN(population, N);
        	//System.out.println("using: Roulette Wheel");
    		break;
        case 1: selected = SUS.selectN(population, N);
        	//System.out.println("using: SUS");
        	break;
        case 2: selected = Tournament.selectN(population, 2, N);
        	//System.out.println("using: Tournament with s = 2");
        	break;
        case 3: selected = Tournament.selectN(population, 4, N);
        	//System.out.println("using: Tournament with s = 4");
    		break;
        case 4: selected = Truncation.selectNfromFittestPercentage(population, N, 10);
        	//System.out.println("using: Truncation with 10% cutoff");
    		break;
        case 5: selected = Truncation.selectNfromFittestPercentage(population, N, 50);
        	//System.out.println("using: Truncation with 50% cutoff");
    		break;
        default: selected = Truncation.selectNfromFittestPercentage(population, N, 50);
        	//System.out.println("using: Truncation with 50% cutoff");
    		break;	
		}
		
		//CROSS-OVER
		// CO_method
		/*
		 * 0 -> Uniform Crossover
		 * 1 -> One-point Crossover
		 * 2 -> Two-point Crossover
		 */
		//apply One Point Crossover to every two elements of the selected list
		ArrayList<String> population_new = new ArrayList<String>();
		
		// for every 2 elements of the selected list, according to crossover probability
		// cross over and add to the new population
		for(int i=0; i<selected.size()/2; i++){
			if(RandomGenerator.generate() < pc){
				ArrayList<String> offspring = new ArrayList<String>();
				switch (CO_method) {
		        case 0:	offspring = CrossOver.uniformCrossOver(selected.get(i*2).getBinaryString(), selected.get(i*2+1).getBinaryString());
		    		break;
		        case 1: offspring = CrossOver.onePointCrossOver(selected.get(i*2).getBinaryString(), selected.get(i*2+1).getBinaryString());
		        	break;
		        case 2: offspring = CrossOver.twoPointCrossOver(selected.get(i*2).getBinaryString(), selected.get(i*2+1).getBinaryString());
		        	break;
		        default: offspring = CrossOver.uniformCrossOver(selected.get(i*2).getBinaryString(), selected.get(i*2+1).getBinaryString());
		    		break;	
				}
				
				population_new.addAll(offspring);
			}
			else{
				population_new.add(selected.get(i*2).getBinaryString());
				population_new.add(selected.get(i*2+1).getBinaryString());
			}
		}
		
		//MUTATION
		//apply mutation to every new element
		for(int i=0;i<population_new.size();i++)				
			population_new.set(i, Mutation.binaryStringMutation(population_new.get(i), pm));
		
		// REPLACE		
		ArrayList<BinaryString> populationBS_new = BinaryString.onemaxPopulationFitness(population_new);
		ArrayList<BinaryString> populationBS_final = Replace.replace(population, populationBS_new, r);
		
				
		return populationBS_final;
		
	}
	
}

class SimpleGeneticAlghorithm{
		
	/**
	 * Simulate one generation of a genetic algorithm operating on the onemax problem.
	 *  
	 * 	@param N -  number of elements
	   	@param l -  lenght of the string representation
	   	@param s -  tournament size
		@param pm - probability of mutation
		@param pc - probability of crossover
		@param r -  fraction of the worst elements from the 
		    		old population that is going to be replaced.
					The actual number of individuals should be round up.
		@param M - the amount of random numbers given. Uniformly generated in [0..1)
		@param u - list with random numbers	
	*/
	public void oneGenerationONEMAX(int N, int l, int s, double pm, double pc, double r){
		
		DecimalFormat df = new DecimalFormat("#0.00");

		
		// INITIAL POPULATION GENERATION	
		// generates the population and calculates the fitness for each element
		ArrayList<String> population = BinaryString.generatePopulation(N, l);
		ArrayList<BinaryString> populationBS = BinaryString.onemaxPopulationFitness(population);		
		
		// prints the status of the initial population
		double max = BinaryString.fitness_max(populationBS).getFitness();
		double avg = BinaryString.fitness_average(populationBS);
		double min = BinaryString.fitness_min(populationBS).getFitness();
		
		System.out.println("0: " + df.format(max) + " " + df.format(avg) + " " + df.format(min));
		
		//SELECTION		
		TournamentWithoutReplacement twr = new TournamentWithoutReplacement();
		ArrayList<BinaryString> selected =  twr.runTournaments(N, s, populationBS);
		
		//CROSS-OVER
		//apply One Point Crossover to every two elements of the selected list
		ArrayList<String> population_new = new ArrayList<String>();
		
		// for every 2 elements of the selected list, according to crossover probability
		// cross over and add to the new population
		for(int i=0; i<selected.size()/2; i++)
			if(RandomGenerator.generate() < pc)
				population_new.addAll(
						CrossOver.onePointCrossOver(
								selected.get(i*2).getBinaryString(),
								selected.get(i*2+1).getBinaryString()));
			else{
				population_new.add(selected.get(i*2).getBinaryString());
				population_new.add(selected.get(i*2+1).getBinaryString());
			}

		//MUTATION
		//apply mutation to every new element
		for(int i=0;i<population_new.size();i++){	
			population_new.set(i, Mutation.binaryStringMutation(population_new.get(i), pm));
		}		
		
		
		// REPLACE	
		ArrayList<BinaryString> populationBS_new = BinaryString.onemaxPopulationFitness(population_new);
		ArrayList<BinaryString> populationBS_final = Replace.replace(populationBS, populationBS_new, r);
				
		// prints the status of the new population
		max = BinaryString.fitness_max(populationBS_final).getFitness();
		avg = BinaryString.fitness_average(populationBS_final);
		min = BinaryString.fitness_min(populationBS_final).getFitness();

		System.out.println("1: " + df.format(max) + " " + df.format(avg) + " " + df.format(min));
		
	}
	
	/**
	 *  Iteratively Simulates one generation of a genetic algorithm operating on the onemax problem.
	 *  
	 * 	@param N -  number of elements
	   	@param l -  lenght of the string representation
	   	@param s -  tournament size
		@param pm - probability of mutation
		@param pc - probability of crossover
		@param r -  fraction of the worst elements from the 
		    		old population that is going to be replaced.
					The actual number of individuals should be round up.
		@param M - the amount of random numbers given. Uniformly generated in [0..1)
		@param u - list with random numbers	
	 * @return 
	*/
	public static void geneticONEMAX(int N, int l, int s, double pm, double pc, double r, double g){
		
		// INITIAL POPULATION GENERATION	
		// generates the population and calculates the fitness for each element
		ArrayList<String> population = BinaryString.generatePopulation(N, l);
		ArrayList<BinaryString> populationBS = BinaryString.onemaxPopulationFitness(population);

		// prints the status of the initial population
		double max = BinaryString.fitness_max(populationBS).getFitness();
		double avg = BinaryString.fitness_average(populationBS);
		double min = BinaryString.fitness_min(populationBS).getFitness();

		DecimalFormat df = new DecimalFormat("#0.00");
		System.out.println("0: " + df.format(max) + " " + df.format(avg) + " " + df.format(min));
		
		for(int i=0; i<g;i++){
			populationBS = geneticONEMAX_iteration(populationBS, N, l, s, pm, pc, r);
			
			// prints the status of the initial population
			max = BinaryString.fitness_max(populationBS).getFitness();
			avg = BinaryString.fitness_average(populationBS);
			min = BinaryString.fitness_min(populationBS).getFitness();

			System.out.println(i+1 + ": " + df.format(max) + " " + df.format(avg) + " " + df.format(min));
		}
	}
	
	private static ArrayList<BinaryString> geneticONEMAX_iteration(ArrayList<BinaryString> population, int N, int l, int s, double pm, double pc, double r){
		
		//SELECTION	
		TournamentWithoutReplacement twr = new TournamentWithoutReplacement();
		ArrayList<BinaryString> selected =  twr.runTournaments(N, s, population);
		
		//CROSS-OVER
		//apply One Point Crossover to every two elements of the selected list
		ArrayList<String> population_new = new ArrayList<String>();
		
		// for every 2 elements of the selected list, according to crossover probability
		// cross over and add to the new population
		for(int i=0; i<selected.size()/2; i++)
			if(RandomGenerator.generate() < pc)
				population_new.addAll(
						CrossOver.onePointCrossOver(
								selected.get(i*2).getBinaryString(),
								selected.get(i*2+1).getBinaryString()));
			else{
				population_new.add(selected.get(i*2).getBinaryString());
				population_new.add(selected.get(i*2+1).getBinaryString());
			}
		
		//MUTATION
		//apply mutation to every new element
		for(int i=0;i<population_new.size();i++)				
			population_new.set(i, Mutation.binaryStringMutation(population_new.get(i), pm));
		
		// REPLACE		
		ArrayList<BinaryString> populationBS_new = BinaryString.onemaxPopulationFitness(population_new);
		ArrayList<BinaryString> populationBS_final = Replace.replace(population, populationBS_new, r);
		
				
		return populationBS_final;
	}
	
}

class Replace{
	
	public static ArrayList<BinaryString> replace (ArrayList<BinaryString> population_old, ArrayList<BinaryString>population_new, double replacement_factor){
		
		ArrayList<BinaryString> population_final = new ArrayList<BinaryString>();
		
		int keepers = population_old.size()-(int)(replacement_factor*population_old.size()); 
		
		Collections.sort(population_old);
		//Collections.sort(population_new);
		
		population_final.addAll(population_old.subList(0, keepers));
		population_final.addAll(population_new.subList(0, population_old.size()-keepers));	
			
		return population_final;
	}

}
class CrossOver{
	// Given two BinaryStrings of size l and a random number, returns the results of onePointCrossover for those binaryStrings (which is another pair of binaryStrings)
	public static ArrayList<String> onePointCrossOver (String bsX, String bsY){
		
		// Generates pattern r
		String r = onePointPatternGenerator(bsX.length());
		
		return childs(bsX,bsY,r);
		
	}
	
	// Given two BinaryStrings of size l and a random number, returns the results of onePointCrossover for those binaryStrings (which is another pair of binaryStrings)
	public static ArrayList<String> twoPointCrossOver (String bsX, String bsY){

		// Generates pattern r
		String r = twoPointPatternGenerator(bsX.length());
		
		//System.out.println();
		//System.out.println(r);
		return childs(bsX,bsY,r);

	}

	// Given two BinaryStrings of size l and l random numbers, returns the results of uniformCrossOver for those binaryStrings (which is another pair of binaryStrings)
	public static ArrayList<String> uniformCrossOver (String bsX, String bsY){
		
		// String manipulator
		//BinaryStringManipulator bsm = new BinaryStringManipulator();
		
		int length = bsX.length();
		// Generates pattern r and the inverse r_
		String r = uniformPatternGenerator(length);
		/*String r_ 	= bsm.NOT(r);
		
		// child 1 from X and Y = r_.x ^ r.y   (. is the "AND" operation and ^ the "XOR")
		// child 2 from X and Y = r.x  ^ r_.y  (. is the "AND" operation and ^ the "XOR") 
		String bsXY1 = bsm.XOR(bsm.AND(bsX, r), bsm.AND(bsY, r_)); 
		String bsXY2 = bsm.XOR(bsm.AND(bsY, r), bsm.AND(bsX, r_));
		
		
		// Store and return the resulting children from the crossover
		ArrayList<String> childs = new ArrayList<String>();
		childs.add(bsm.missingZeros(bsXY1, r.length()));
		childs.add(bsm.missingZeros(bsXY2, r.length()));
		
		return childs;*/
		
		return childs(bsX,bsY,r);
	}

	// Return a random onePoint pattern with the given length
	private static String twoPointPatternGenerator(int l) {
		
		// calculates the pattern splitting point
		int point_a = RandomGenerator.generateInRange_int(0,l-1);
		int point_b = RandomGenerator.generateInRange_int(0,l-1);
		while(point_b == point_a)
			point_b = RandomGenerator.generateInRange_int(0,l-1);
		
		
		// bits before the smallest and after the largest of the two splitting points get a 1,
		// bits between the splitting points get a 0
		String pattern = "";
		for(int i=0; i<=Math.min(point_a,point_b);i++ )
			pattern += "1";
		for(int i=Math.min(point_a,point_b); i<Math.max(point_a,point_b);i++ )
			pattern += "0";
		for(int i=Math.max(point_a,point_b); i<l-1;i++ )
			pattern += "1";
		
		// return the pattern
		return pattern;
	}
	
	// Return a random onePoint pattern with the given length
	private static String onePointPatternGenerator(int l) {
		
		// calculates the pattern splitting point
		int onePoint = RandomGenerator.generateInRange_int(0,l-1);
		
		// bits before the splitting point get a 1, bits at the splitting point and after get a 0 
		String pattern = "";
		for(int i=0; i<=onePoint;i++ )
			pattern += "1";
		for(int i=onePoint; i<l-1;i++ )
			pattern += "0";
		
		// return the pattern
		return pattern;
	}

	// Given a list of random numbers returns the crossover application pattern
	private static String uniformPatternGenerator(int length){
		return BinaryString.generate(length);
	}
	
	private static ArrayList<String> childs (String bsX, String bsY, String r){
		// Initialize childs
		String bsX_="";
		String bsY_ = "";
		
		// In positions where pattern r has the value 0, keep bits the same as corresponding parent
		// 						where r has the value 1, swap the bit with the other parent
		for(int i=0; i<r.length(); i++)
			if(r.charAt(i)=='1'){
				//swap bits
				bsX_ += bsY.charAt(i);
				bsY_ += bsX.charAt(i);
			}else{
				//maintain bits
				bsX_ += bsX.charAt(i);
				bsY_ += bsY.charAt(i);
			}
		
		// Store and return the resulting children from the crossover
		ArrayList<String> childs = new ArrayList<String>();
		childs.add(bsY_);
		childs.add(bsX_);
		
		return childs;
	}
}


class Mutation{
	
	// Computes the mutation of a binaryString bit by bit given the probability of mutation for all the bits and a random value for each bit 
	public static String binaryStringMutation(String binaryString, double pMut){
						
		// The returning bynaryString
		String binaryString_mutated = binaryString;
		
		// for all the bits in the string, check if the corresponding random value is less than the probability of mutation, if so mutate the bit.
		for(int i=0; i<binaryString.length();i++)
			if(RandomGenerator.generate()<pMut)
				binaryString_mutated = BinaryString.FLIP(binaryString_mutated, i);
		
		return binaryString_mutated;
		
	}
}

/////// SELECTION /////

class Truncation{
	
	// Returns a list of N elements from the specified percentage of the fittest elements of the list
	public static ArrayList<BinaryString> selectNfromFittestPercentage(ArrayList<BinaryString> binaryString, int N, double percentage){
		
		ArrayList<BinaryString> fittest = selectFittestPercentage(binaryString, percentage);
		ArrayList<BinaryString> selected = selectNfromFittestPercentage(fittest, N);
		
		return selected;
	}
	
	private static ArrayList<BinaryString> selectNfromFittestPercentage(ArrayList<BinaryString> fittestBinaryStrings, int N){
		
		ArrayList<BinaryString> selected = new ArrayList<BinaryString>();
		
		//iterate through the selected elements
		Iterator<BinaryString> iterator = fittestBinaryStrings.iterator();
		
		for(int i=0; i<N; i++){
			//if there is no more selected elements, re-assign the iterator
			if(!iterator.hasNext())
				iterator = fittestBinaryStrings.iterator();
			
			selected.add(iterator.next());
		}
		
		assert selected.size() == N;
		
		return selected;
		
		
	}
	
	// Returns the specified percentage of the fittest elements of the list
	private static ArrayList<BinaryString> selectFittestPercentage(ArrayList<BinaryString> binaryStrings, double percentage){
		
		Collections.sort(binaryStrings);
		
		ArrayList<BinaryString> result = new ArrayList<BinaryString>();
		result.addAll(binaryStrings.subList(0, (int)Math.floor(binaryStrings.size()*percentage/100)));
		
		return result;
		
	}
	
}

//Stochastic Universal Sampling
class SUS{
		
	// u is a random number and numberElements is the number of elements to keep
	public static ArrayList<BinaryString> selectN(ArrayList<BinaryString> binaryStrings, int numberElements){
		
		// List containing the selected elements
		ArrayList<BinaryString> selected = new ArrayList<BinaryString>();
		
		// Wheel Distribution list
		ArrayList<Double> wheelDistribution = generateWheelDistribution(binaryStrings);
		
		// Calculates population total fitness
		double totalFitness = calculateTotalFitness(binaryStrings);
		
		// Calculates the distance between the pointers (equal distance between all points) 
		double pointerDistance = totalFitness/numberElements;
		
		// Chooses the starting point of the random points
		double startingPointer = RandomGenerator.generateInRange(0,pointerDistance);
		
		//scales the pointerDistance and startingPointer into 0 to 1 values
		double pointerDistanceScaled = pointerDistance/totalFitness;
		double startingPointerScaled = startingPointer/totalFitness;
		
		// Searches the wheel distribution for the element corresponding to each point
		Iterator<Double> iterator = wheelDistribution.iterator();
		
		int pointerNumber = 0;
		double currentPointer = startingPointerScaled;
		int index = 0;
		
		double currentWheel = iterator.next();
		
		while(pointerNumber < numberElements ){
			// add the element corresponding to the currrent slice if the random number is bigger than the current slice limit, but wasn't bigger than the previous
			if(currentPointer<=currentWheel){
				 selected.add(binaryStrings.get(index));
				 pointerNumber++;
				 currentPointer += pointerDistanceScaled; 
			}else{
				currentWheel = iterator.next();
				index++;
			}
		}
				
		return selected;
	}

	// Calculates the sum of the fitness of all elements (basic accumulation function)
	private static Double calculateTotalFitness (ArrayList<BinaryString> binaryStrings){
		double totalFitness = 0;
		Iterator<BinaryString> iterator = binaryStrings.iterator();
		while(iterator.hasNext()){
			totalFitness += iterator.next().getFitness();
		}

		return totalFitness;
	}

	// Given a list of BinaryStrings (with corresponding fitnesses), generates the ordered list with the corresponding wheel "slice" for each element. 
	// Offset taken into consideration so, each element corresponds to the maximum value between 0 and 1 that belongs to the indexed "slice"
	// Ex: if element 1 has p=0.15 and element 2 has p=0.20... The list will have at least 3 elements, being the first two 0.15 and 0.35. Where 0.35 = 0.15 + 0.20.
	private static ArrayList<Double> generateWheelDistribution (ArrayList<BinaryString> binaryStrings){

		// Stores each slice value (slice probability + all the previous slices probabilities)
		ArrayList<Double> wheelDistribution = new ArrayList<Double>();

		// Calculates Total Fitness
		double totalFitness = calculateTotalFitness(binaryStrings);

		// Offset accumulator
		double offset = 0;

		// Computes the wheelDistribution considering each element
		Iterator<BinaryString> iterator = binaryStrings.iterator();
		while(iterator.hasNext()){
			offset += iterator.next().getFitness() / totalFitness;
			wheelDistribution.add(offset);
		}

		// The offset value has to be exactly 1 uppon accumulating all the elements.
		assert(offset==1);

		return wheelDistribution;
	}

}

// Can be improved, generating the Wheel Distribution can be done only once for all the population,
// it's not necessary to generate it every time the wheel is executed
class RouletteWheel{
	
	// O(3*n)
	/*public static ArrayList<BinaryString> selectN(ArrayList<BinaryString> binaryStrings, int N){
					
			// Wheel Distribution list
			ArrayList<Double> wheelDistribution = generateWheelDistribution(binaryStrings);
			
			// List with the selected elements
			ArrayList<BinaryString> selected = new ArrayList<BinaryString>();
			
			for(int i=0; i<N; i++){
				// Searches for the slice corresponding to the given random number;
				int index = 0;
				Iterator<Double> iterator = wheelDistribution.iterator();
				
				double random = RandomGenerator.generate();
				while(iterator.hasNext())
					// return the element corresponding to the currrent slice if the random number is bigger than the current slice limit, but wasn't bigger than the previous
					if(random<iterator.next()){
						selected.add(binaryStrings.get(index));
						break;
					}
					else
						index++;
			}
			
				
			return selected;
			
		}
		*/
	
	public static ArrayList<BinaryString> selectN(ArrayList<BinaryString> binaryStrings, int N){
		ArrayList<BinaryString> selected = new ArrayList<BinaryString>();
		for(int i=0; i<N; i++)
			selected.add(selectOne(binaryStrings));
		
		return selected;
	}
	
	// O(3*n)
	public static BinaryString selectOne(ArrayList<BinaryString> binaryStrings){
				
		// Wheel Distribution list
		ArrayList<Double> wheelDistribution = generateWheelDistribution(binaryStrings);
		
		// Searches for the slice corresponding to the given random number;
		int index = 0;
		Iterator<Double> iterator = wheelDistribution.iterator();
		
		double random = RandomGenerator.generate();
		
		while(iterator.hasNext())
			// return the element corresponding to the currrent slice if the random number is bigger than the current slice limit, but wasn't bigger than the previous
			if(random<iterator.next())
				return binaryStrings.get(index);
			else
				index++;
		
		return null;
		
	}
	
	// Given a list of BinaryStrings (with corresponding fitnesses), generates the ordered list with the corresponding wheel "slice" for each element. 
	// Offset taken into consideration so, each element corresponds to the maximum value between 0 and 1 that belongs to the indexed "slice"
	// Ex: if element 1 has p=0.15 and element 2 has p=0.20... The list will have at least 3 elements, being the first two 0.15 and 0.35. Where 0.35 = 0.15 + 0.20.
	private static ArrayList<Double> generateWheelDistribution (ArrayList<BinaryString> binaryStrings){
		
		// Stores each slice value (slice probability + all the previous slices probabilities)
		ArrayList<Double> wheelDistribution = new ArrayList<Double>();
		
		// Calculates Total Fitness
		double totalFitness = calculateTotalFitness(binaryStrings);
		
		// Offset accumulator
		double offset = 0;
		
		// Computes the wheelDistribution considering each element
		Iterator<BinaryString> iterator = binaryStrings.iterator();
		while(iterator.hasNext()){
			offset += iterator.next().getFitness() / totalFitness;
			wheelDistribution.add(offset);
		}
		
		// The offset value has to be exactly 1 uppon accumulating all the elements.
		assert(offset==1);
		
		return wheelDistribution;
	}
	
	// Calculates the sum of the fitness of all elements (basic accumulation function)
	private static Double calculateTotalFitness (ArrayList<BinaryString> binaryStrings){
		double totalFitness = 0;
		Iterator<BinaryString> iterator = binaryStrings.iterator();
		while(iterator.hasNext()){
			totalFitness += iterator.next().getFitness();
		}
		
		return totalFitness;
	}
	
}

class TournamentWithoutReplacement {
	// N/s tournaments
	// s -> number of chosen elements
	// s * (N - 1) random numbers
	// Runs the tournament for the elements on binaryStrings list, randomly chooses s elements, and from those winner elements choose the fittest as the ultimate winner
	public ArrayList<BinaryString> runTournaments(int N, int s, ArrayList<BinaryString> binaryStrings){
		
		// List of winners
		ArrayList<BinaryString> winners = new ArrayList<BinaryString>();
		
		// You will	need to generate a random permutation of the population s times
		ArrayList<Permutation> permutations = generatePermutations(N,s);
		assert(permutations.size()==s); 
		
		// For each permutation
		for(int i=0; i<s;i++){
			
			// Perform a N/s group of selections, with N/s elements
			ArrayList<ArrayList<Integer>> selections = permutationTournamentSelections(N,s,permutations.get(i));
			
			for(int j = 0; j<N/s; j++){
				ArrayList<BinaryString> tournamentSelected = selectedBinaryStringFromIndexes(binaryStrings, selections.get(j));
				BinaryString winner = selectFittestBinaryString(tournamentSelected);
				winners.add(winner);
			}
		}
		
		
		return winners;
			
			
	}
	
	// Generate s random permutations of size N-1 
	private ArrayList<Permutation> generatePermutations (int N, int s){
		
		// Generate a random permutation of the population s times
		ArrayList<Permutation> permutations = new ArrayList<Permutation>();
		for(int i = 0; i<s; i++){
			
			// Generate the permutation and add it to the list
			Permutation permutation = Permutation.generateRandom(N);
			permutations.add(permutation);
		}
		return permutations;
	}

	// Given a permutation and the number of tournaments it represents, return an array with each group of selected elements
	private ArrayList<ArrayList<Integer>> permutationTournamentSelections(int N, int s, Permutation permutation){
		
		ArrayList<ArrayList<Integer>> tournamentSelections = new ArrayList<ArrayList<Integer>>();
		//int t = N/s;
		
		// selects a group of elements
		for(int j=0; j<N/s; j++){

			// select the ith group of N/s values from the permutation
			ArrayList<Integer> selected = new ArrayList<Integer>();
			int from = j*(s);
			int to = (j+1)*(s);
			selected.addAll(permutation.getPermutation().subList(from, to));
			tournamentSelections.add(selected);
		}

		return tournamentSelections;
	}

	// Given a list of binaryStrings and a list of indexes with the positions of selected binaryStrings, return a list of the indexed binaryStrings
	private ArrayList<BinaryString> selectedBinaryStringFromIndexes(ArrayList<BinaryString> binaryStrings, ArrayList<Integer> selected){
		
		//The list wich will store the elements to return
		ArrayList<BinaryString> tournamentSelected = new ArrayList<BinaryString>();
		
		//for each index element, select the corresponding binaryString
		for(int i = 0 ; i<selected.size();i++)
			tournamentSelected.add(binaryStrings.get(selected.get(i)));
		
		return tournamentSelected;
	}

	// Receives a non-empty list of BinaryString and returns the fittest element
	private BinaryString selectFittestBinaryString(ArrayList<BinaryString> binaryStrings){
		Iterator<BinaryString> iterator = binaryStrings.iterator();

		// If the list is not empty search for the fittest element
		if(!binaryStrings.isEmpty()){
			BinaryString fittest, current;
			fittest = current = iterator.next();

			// while there is elementes on the list, update the fittest
			while(iterator.hasNext()){
				current = iterator.next();
				if(current.getFitness() > fittest.getFitness())
					fittest = current;
			}

			return fittest;


		}else // if the list is empty return null
			return null;
	}

}

class Tournament{
	
	// Runs the tournament for the elements on binaryStrings list, randomly chooses s winner elements
	// From these s winner elements choose the fittest as the ultimate winner
	public static ArrayList<BinaryString> selectN(ArrayList<BinaryString> binaryStrings, int s, int N){
		ArrayList<BinaryString> selected = new ArrayList<BinaryString>();
		for(int i=0; i<N; i++)
			selected.add(selectOne(binaryStrings,s));
		
		return selected;
	}
	
	// Runs the tournament for the elements on binaryStrings list, randomly chooses s winner elements
	// From these s winner elements choose the fittest as the ultimate winner
	public static BinaryString selectOne(ArrayList<BinaryString> binaryStrings, int s){
		ArrayList<BinaryString> randomBinaryStrings = new ArrayList<BinaryString>();
		
		// Selects s random binaryStrings
		for(int i = 0; i<s; i++){
			randomBinaryStrings.add(selectRandomBinaryString(binaryStrings));
		}
		
		return selectFittestBinaryString(randomBinaryStrings);
		
	}
		
	// Receives a list of BinaryString returns a random element from the list
	private static BinaryString selectRandomBinaryString(ArrayList<BinaryString> binaryStrings){
		int index = RandomGenerator.generateInRange_int(1,binaryStrings.size());
		return binaryStrings.get(index-1);
	}
	
	// Receives a non-empty list of BinaryString and returns the fittest element
	private static BinaryString selectFittestBinaryString(ArrayList<BinaryString> binaryStrings){
		Iterator<BinaryString> iterator = binaryStrings.iterator();

		// If the list is not empty search for the fittest element
		if(!binaryStrings.isEmpty()){
			BinaryString fittest, current;
			fittest = current = iterator.next();
			
			// while there is elementes on the list, update the fittest
			while(iterator.hasNext()){
				current = iterator.next();
				if(current.getFitness() > fittest.getFitness())
					fittest = current;
			}
			
			return fittest;

			
		}else // if the list is empty return null
			return null;
	}

}

class Permutation {
	
	private ArrayList<Integer> permutation;
	
	/* CONSTRUCTORS */
	public Permutation(){
		this.permutation = new ArrayList<Integer>();
	}
	public Permutation(ArrayList<Integer> permutation){
		this.permutation = permutation;
	}

	/* GETTERS & SETTERS */
	public ArrayList<Integer> getPermutation() {
		return permutation;
	}
	public void setPermutation(ArrayList<Integer> permutation) {
		this.permutation = permutation;
	}
	
	/* METHODS */
	// prints one element of the permutation per line
	public void print_1perLine(){
		//
		for(int i=0;i<permutation.size();i++)
			System.out.println(permutation.get(i));
	}
	
	// swaps the elements of the referred positions in the permutation list
	public void swap(int a, int b){
		int aux = permutation.get(a);
		permutation.set(a, permutation.get(b));
		permutation.set(b,aux);
	}

	// generates a random permutation of size N using N-1 random numbers
	public static Permutation generateRandom(int N){

		//The permutation
		ArrayList<Integer> permutation_ = new ArrayList<Integer>();

		//Generate the permutation ordered 
		//for( i=0; i<N; i++)
		for(int i=0; i<N;i++)
			// v[i] = i;
			permutation_.add(i);

		// Instantiate the permutation object
		Permutation permutation = new Permutation(permutation_);

		// Swap the elements to generate a random permutation
		//for(i=0;i<N-1;i++)
		for(int i=0; i<N-1;i++){
			// r = obtain a uniformly generated random integer in [i..N-1];
			int mappedPosition = randomMap(i, N-1, RandomGenerator.generate());

			// exchange the contents of v[i] with the contents of v[r]
			permutation.swap(i, mappedPosition);
		}

		//returns the random permutation
		return permutation;
	}

	// Receives a random number and the smaller and larger limits. Maps the random number to an integer between those limits
	private static int randomMap(double a, double b, double u){
		return (int)(a + (Math.floor(u * (double)(b - a + 1))));
	}
}


class BinaryString implements Comparator<BinaryString>, Comparable<BinaryString>{ 
	
	private String binaryString;
	private double fitness;
	
	/* CONSTRUCTORS */
	public String getBinaryString() {
		return binaryString;
	}
	
	/* SETTERS & GETTERS*/
	public void setBinaryString(String binaryString) {
		this.binaryString = binaryString;
	}
	public double getFitness() {
		return fitness;
	}
	public void setFitness(double d) {
		this.fitness = d;
	}

	// COMPARATOR
	@Override
	public int compareTo(BinaryString bs) {
		if(this.getFitness()>bs.getFitness())
			return -1;
		else if(this.getFitness()<bs.getFitness())
			return 1;
		else
			return 0;
	}
	@Override
	public int compare(BinaryString bs0, BinaryString bs1) {
		if(bs0.getFitness()>bs1.getFitness())
			return -1;
		else if(bs0.getFitness()<bs1.getFitness())
			return 1;
		else
			return 0;
	}

	/* METHODS */
	// Generates a population of N random binaryStrings with the given length 
	public static ArrayList<String> generatePopulation (int N, int length){

		// List that will store the generated population of elements
		ArrayList<String> population = new ArrayList<String>();

		//Generate N elements
		for(int i=0; i<N; i++)
			population.add(generate(length));

		// Return the generated population
		return population;
	}

	// Generates a random binaryString with the given length
	public static String generate(int length){

		// the default binary string
		String binaryString = "";

		for(int i=0; i<length; i++)
			if(RandomGenerator.generate() < 0.5)
				binaryString = binaryString + "0";
			else
				binaryString = binaryString + "1";

		// return the binaryString
		return binaryString;		
	}

	// Given a list with binaryString representations, calculate each individual fitness and generate a
	// list with all the population as BinaryStrings
	public static ArrayList<BinaryString> onemaxPopulationFitness (ArrayList<String> population){

		// ArrayList that will store the population of binaryStrings
		ArrayList<BinaryString> populationBS = new ArrayList<BinaryString>();

		// For every representation
		Iterator<String> iterator = population.iterator();
		while(iterator.hasNext()){

			// get the string element
			String stringElement = iterator.next();

			//create the binaryString element and set its representation and fitness
			BinaryString bsElement = new BinaryString();
			bsElement.setBinaryString(stringElement);
			bsElement.setFitness(onemax(stringElement, stringElement.length()));	// fitness is given by onemax function

			// add the BinaryString element to the BinaryString population
			populationBS.add(bsElement);
		}

		// return the population of BinaryStrings
		return populationBS;
	}

	// Receives a binaryString and calculates the number of 0 bits
	public static int onemax(String binaryString, int length){
		int count=0;

		for(int i = 0; i<length;i++)
			if(binaryString.charAt(i) == '1')
				count++;

		return count;
	}

	// Converts a binaryString into a decimal number
	public static int binaryStringtoDecimal(String binaryString){
		int result = 0;
		int length = binaryString.length();
		for(int i =0; i<length; i++)
			result+=Math.pow(2, i) * Integer.parseInt(binaryString.charAt(length-i-1)+"");

		return result;
	}

	// Converts a decimal number into a binaryString
	public static String decimaltoBinaryString(int value){
		return Integer.toBinaryString(value);
	}

	// Calculates the fitness of a binaryString based on the function f(x)=x^2
	public static int fitSquareX(String binaryString, int length){
		return (int) Math.pow(binaryStringtoDecimal(binaryString),2);
	}

	// Print a list of BinaryStrings, one per line
	public static void printListOnePerLine(ArrayList<BinaryString> binaryStrings){

		Iterator<BinaryString> iterator = binaryStrings.iterator();

		while(iterator.hasNext())
			System.out.println(iterator.next().getBinaryString());
	}

	// bitwise operation "BIT FLIP" (not actually computed with bitwise operators)
	public static String FLIP(String binaryString, int index){

		// store the sub strings of bits preceding and proceeding the bit we want to flip 
		String pre = binaryString.substring(0, index);
		String post = binaryString.substring(index+1, binaryString.length());

		// return the concatenation of the preceding substring, the fliped bit and the proceeding substring of bits
		if(binaryString.charAt(index)=='0')
			return pre + "1" + post;
		else
			return pre + "0" + post;
	}

	// bitwise operation "AND"
	public static String AND(String binaryString1,String binaryString2){
		// Convert Strings to decimal
		int binaryString1_10 = binaryStringtoDecimal(binaryString1);
		int binaryString2_10 = binaryStringtoDecimal(binaryString2);

		// Realize bitwise operation
		int bsResult_10 = binaryString1_10 & binaryString2_10;

		// Return corresponding String
		return decimaltoBinaryString(bsResult_10);
	}

	// bitwise operation "OR"
	public static String OR(String binaryString1,String binaryString2){
		// Convert Strings to decimal
		int binaryString1_10 = binaryStringtoDecimal(binaryString1);
		int binaryString2_10 = binaryStringtoDecimal(binaryString2);

		// Realize bitwise operation
		int bsResult_10 = binaryString1_10 | binaryString2_10;

		// Return corresponding String
		return decimaltoBinaryString(bsResult_10);
	}

	// bitwise operation "exclusive-OR"
	public static String XOR(String binaryString1,String binaryString2){
		// Convert Strings to decimal
		int binaryString1_10 = binaryStringtoDecimal(binaryString1);
		int binaryString2_10 = binaryStringtoDecimal(binaryString2);

		// Realize bitwise operation
		int bsResult_10 = binaryString1_10 ^ binaryString2_10;

		// Return corresponding String
		return decimaltoBinaryString(bsResult_10);
	}

	// bitwise operation "NOT"
	public static String NOT(String binaryString){
		// Convert String to decimal
		int binaryString_10 = binaryStringtoDecimal(binaryString);

		// Realize bitwise operation
		int bsResult_10 = ~binaryString_10;

		// Return corresponding String
		return decimaltoBinaryString(bsResult_10);
	}

	// bitwise operation "SHIFT LEFT"
	public static String SHIFT_L(String binaryString, int n){
		// Convert String to decimal
		int binaryString_10 = binaryStringtoDecimal(binaryString);

		// Realize bitwise operation
		int bsResult_10 = binaryString_10 << n;

		// Return corresponding String
		return decimaltoBinaryString(bsResult_10);
	}

	// bitwise operation "SHIFT RIGHT"
	public static String SHIFT_R(String binaryString, int n){
			// Convert String to decimal
			int binaryString_10 = binaryStringtoDecimal(binaryString);

			// Realize bitwise operation
			int bsResult_10 = binaryString_10 >> n;

			// Return corresponding String
			return decimaltoBinaryString(bsResult_10);
		}

	// Search for the fittest element
	public static BinaryString fitness_max(ArrayList<BinaryString> binaryStrings){
		//the fittest so far
		BinaryString fittest = binaryStrings.get(0);
		BinaryString current;
		
		Iterator<BinaryString> iterator = binaryStrings.iterator();
		while(iterator.hasNext()){
			current = iterator.next();
			if(current.getFitness()>fittest.getFitness())
				fittest = current;
		}
		
		return fittest;
		
	}
	
	// Search for the fittest element
	public static BinaryString fitness_min(ArrayList<BinaryString> binaryStrings){
		//the most unsuited so far
		BinaryString unsuited = binaryStrings.get(0);
		BinaryString current;

		Iterator<BinaryString> iterator = binaryStrings.iterator();
		while(iterator.hasNext()){
			current = iterator.next();
			if(current.getFitness()<unsuited.getFitness())
				unsuited = current;
		}

		return unsuited;

	}
	
	// Calculate the average fitness for a list of binaryStrings
	public static Double fitness_average(ArrayList<BinaryString> binaryStrings){
		//the most unsuited so far
		double accumulator = 0;

		Iterator<BinaryString> iterator = binaryStrings.iterator();
		while(iterator.hasNext()){
			accumulator += iterator.next().getFitness();
		}

		return (double)accumulator/(double)binaryStrings.size();

	}
	
	// Calculate the standard deviation fitness for a list of binaryStrings
	public static Double fitness_stdDeviation(ArrayList<BinaryString> binaryStrings){
		double average = fitness_average(binaryStrings);
		
		//accumulator
		double variance_sum = 0;

		Iterator<BinaryString> iterator = binaryStrings.iterator();

		//calculate and sum the difference of each data point from the mean, and square the result of each:
		while(iterator.hasNext())
			variance_sum += Math.pow((double)(iterator.next().getFitness()-average),2);
		
		double variance_mean = (double)variance_sum/(double)binaryStrings.size();
		
		return Math.sqrt(variance_mean);
	}

}

class RandomGenerator{
	
	public static double generate(){
		return Math.random();
	}

	// Generates an Integer random number inside the given limits (using floor)
	public static int generateInRange_int(int min, int max){
		return (int)(min + (Math.floor(RandomGenerator.generate() * (double)(max - min + 1))));
	}
	
	// Generates a random number inside the given limits
	public static double generateInRange(double min, double max){
		return generate() * (double)(max - min);
	}
}

class Statistics{
	
	//Returns the average value of a list of values
	public static double average(ArrayList<Double> list){
		//accumulator
		double accumulator = 0;

		Iterator<Double> iterator = list.iterator();

		//calculate sum of all elements
		while(iterator.hasNext())
			accumulator += iterator.next();
		
		//return average
		return (double)accumulator/(double)list.size();
	}
	
	public static double standardDeviation(ArrayList<Double> list){
		double average = average(list);
		
		//accumulator
		double variance_sum = 0;

		Iterator<Double> iterator = list.iterator();

		//calculate and sum the difference of each data point from the mean, and square the result of each:
		while(iterator.hasNext())
			variance_sum += Math.pow((double)(iterator.next()-average),2);
		
		double variance_mean = (double)variance_sum/(double)list.size();
		
		return Math.sqrt(variance_mean);

	}
	
}


