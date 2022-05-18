package T3dD;

//package d3;

import java.text.DecimalFormat;
import java.util.*;

public class T3dD {

	public static void main(String[] args) {
		
		Scanner sc = new Scanner(System.in);
		
		// Read input constants
		int N = sc.nextInt();			//4 N - number of elements
		int l = sc.nextInt();			//3 l - lenght of the string representation
		int s = sc.nextInt();			//2 s - tournament size
		double pm = sc.nextDouble();	// Pm - probability of mutation
		double pc = sc.nextDouble();	// Pc - probability of crossover
		double r = 1;					// r - fraction of the worst elements from the 
										// old population that is going to be replaced.
										// The actual number of individuals should be round up.
		double g = sc.nextInt();		// g - number of generations
		double M = sc.nextInt();		// M - the amount of random numbers given. Uniformly generated in [0..1)
		
		int u1 = N*l; 		// 20 number of random elements for generating the random population
		int u2 = s*(N-1);	// 6 number of random elements for tournament without replacement
		
		// Read M random numbers
		ArrayList<Double> u = new ArrayList<Double>();
		for(int i=0;i<M;i++)
			u.add(sc.nextDouble());

		sc.close();
		
		SimpleGeneticAlghorithm.geneticONEMAX(N, l, s, pm, pc, r, g, M, u);

	}
		
}


class SimpleGeneticAlghorithm{
	
	private static ArrayList<Double> u;
	private static int u_offset; 
	
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
	public void oneGenerationONEMAX(int N, int l, int s, double pm, double pc, double r, double M, ArrayList<Double> u){
		
		DecimalFormat df = new DecimalFormat("#0.00");
		
		// INITIAL POPULATION GENERATION
		// select the random numbers destined to generating the population (N*l random numbers)
		ArrayList<Double> u_ = new ArrayList<Double>();
		
		u_.addAll(u.subList(0, u_offset));
		
		// generates the population and calculates the fitness for each element
		ArrayList<String> population = BinaryString.generatePopulation(N, l, u_);
		ArrayList<BinaryString> populationBS = BinaryString.onemaxPopulationFitness(population);		
		
		// prints the status of the initial population
		double max = BinaryString.fitness_max(populationBS).getFitness();
		double avg = BinaryString.fitness_average(populationBS);
		double min = BinaryString.fitness_min(populationBS).getFitness();
		
		System.out.println("0: " + 
				df.format(max) + " " +
				df.format(avg) + " " +
				df.format(min));
		
		//SELECTION
		//update offset to retrieve correct number of random numbers
		u_.clear();
		u_.addAll(u.subList(u_offset, u_offset+s*(N-1)));
		u_offset += s*(N-1); //26
		
		TournamentWithoutReplacement twr = new TournamentWithoutReplacement();
		ArrayList<BinaryString> selected =  twr.runTournaments(N, s, populationBS, u_);
		
		//CROSS-OVER
		//apply One Point Crossover to every two elements of the selected list
		CrossOver cO = new CrossOver();
		ArrayList<String> population_new = new ArrayList<String>();
		
		//for every 2 elements of the selected list
		for(int i=0; i<selected.size()/2; i++){
			
			//cross over and add to the new population
			if(u.get(u_offset++) < pc)
				population_new.addAll(
						cO.onePointCrossOver(
								selected.get(i*2).getBinaryString(),
								selected.get(i*2+1).getBinaryString(), 
								u.get(u_offset++)));
			else{
				population_new.add(selected.get(i*2).getBinaryString());
				population_new.add(selected.get(i*2+1).getBinaryString());
			}

		}

		
		//MUTATION
		//apply mutation to every new element
		Mutation mtt = new Mutation();
		for(int i=0;i<population_new.size();i++){
			//prepare random numbers
			u_.clear();
			u_.addAll(u.subList(u_offset,u_offset+l));
			
			//apply mutation
			population_new.set(i, mtt.binaryStringMutation(population_new.get(i), pm, u_));
			
			//update random numbers offset
			u_offset+=l;
		}		
		
		
		// REPLACE
		Replace rep = new Replace();
		
		ArrayList<BinaryString> populationBS_new = BinaryString.onemaxPopulationFitness(population_new);
		ArrayList<BinaryString> populationBS_final = rep.replace(populationBS, populationBS_new, r);
		
		
		
		// prints the status of the new population
		max = BinaryString.fitness_max(populationBS_final).getFitness();
		avg = BinaryString.fitness_average(populationBS_final);
		min = BinaryString.fitness_min(populationBS_final).getFitness();

		System.out.println("1: " + 
				df.format(max) + " " +
				df.format(avg) + " " +
				df.format(min));
				
		//System.out.println();
		
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
	public static void geneticONEMAX(int N, int l, int s, double pm, double pc, double r, double g, double M, ArrayList<Double> U){
		
		u = U;
		u_offset = 0;
		
		// INITIAL POPULATION GENERATION
		// select the random numbers destined to generating the population (N*l random numbers)
		ArrayList<Double> u_ = new ArrayList<Double>();
		u_offset = N*l; //20
		u_.addAll(u.subList(0, u_offset));
		
		// generates the population and calculates the fitness for each element
		ArrayList<String> population = BinaryString.generatePopulation(N, l, u_);
		ArrayList<BinaryString> populationBS = BinaryString.onemaxPopulationFitness(population);

		
		// prints the status of the initial population
		double max = BinaryString.fitness_max(populationBS).getFitness();
		double avg = BinaryString.fitness_average(populationBS);
		double min = BinaryString.fitness_min(populationBS).getFitness();

		DecimalFormat df = new DecimalFormat("#0.00");
		System.out.println("0: " + 
				df.format(max) + " " +
				df.format(avg) + " " +
				df.format(min));
		
		for(int i=0; i<g;i++){
			populationBS = geneticONEMAX_iteration(populationBS, N, l, s, pm, pc, r);
			
			// prints the status of the initial population
			max = BinaryString.fitness_max(populationBS).getFitness();
			avg = BinaryString.fitness_average(populationBS);
			min = BinaryString.fitness_min(populationBS).getFitness();

			System.out.println(i+1 + ": " + 
					df.format(max) + " " +
					df.format(avg) + " " +
					df.format(min));
		}
	}
	
	private static ArrayList<BinaryString> geneticONEMAX_iteration(ArrayList<BinaryString> population, int N, int l, int s, double pm, double pc, double r){
		
		ArrayList<Double> u_ = new ArrayList<Double>();
		
		//SELECTION
		//update offset to retrieve correct number of random numbers
		u_.addAll(u.subList(u_offset, u_offset+s*(N-1)));
		u_offset += s*(N-1);
		
		TournamentWithoutReplacement twr = new TournamentWithoutReplacement();
		ArrayList<BinaryString> selected =  twr.runTournaments(N, s, population, u_);
		
		//CROSS-OVER
		//apply One Point Crossover to every two elements of the selected list
		CrossOver cO = new CrossOver();
		ArrayList<String> population_new = new ArrayList<String>();
		
		//for every 2 elements of the selected list
		for(int i=0; i<selected.size()/2; i++){
			
			//cross over and add to the new population
			if(u.get(u_offset++) < pc)
				population_new.addAll(
						cO.onePointCrossOver(
								selected.get(i*2).getBinaryString(),
								selected.get(i*2+1).getBinaryString(), 
								u.get(u_offset++)));
			else{
				population_new.add(selected.get(i*2).getBinaryString());
				population_new.add(selected.get(i*2+1).getBinaryString());
			}

		}
		
		//MUTATION
		//apply mutation to every new element
		Mutation mtt = new Mutation();
		for(int i=0;i<population_new.size();i++){
			//prepare random numbers
			u_.clear();
			u_.addAll(u.subList(u_offset,u_offset+l));
			
			//apply mutation
			population_new.set(i, mtt.binaryStringMutation(population_new.get(i), pm, u_));
			
			//update random numbers offset
			u_offset+=l;
		}		
		
		
		// REPLACE
		Replace rep = new Replace();
		
		ArrayList<BinaryString> populationBS_new = BinaryString.onemaxPopulationFitness(population_new);
		ArrayList<BinaryString> populationBS_final = rep.replace(population, populationBS_new, r);
		
				
		return populationBS_final;
	}
	
}

class Replace{
	
	public ArrayList<BinaryString> replace (ArrayList<BinaryString> population_old, ArrayList<BinaryString>population_new, double replacement_factor){
		
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
	public ArrayList<String> onePointCrossOver (String bsX, String bsY, double u){
		
		// Generates pattern r
		String r = onePointPatternGenerator(u, bsX.length());
		
		return childs(bsX,bsY,r);
		
	}
	
	// Given two BinaryStrings of size l and l random numbers, returns the results of uniformCrossOver for those binaryStrings (which is another pair of binaryStrings)
	public ArrayList<String> uniformCrossOver (String bsX, String bsY, ArrayList<Double> u){
		
		// String manipulator
		//BinaryStringManipulator bsm = new BinaryStringManipulator();
		
		// Generates pattern r and the inverse r_
		String r = uniformPatternGenerator(u);
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
	
	// Given a random number and the desired length, return a onePoint pattern
	private String onePointPatternGenerator(double u, double l) {
		
		// calculates the pattern splitting point
		int onePoint = randomMap(0,l-1,u);
		
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
	private String uniformPatternGenerator(ArrayList<Double> u){
		return BinaryString.generate(u);
	}
	
	// Recieves a random number and the smaller and larger limits. Maps the random number to an integer between those limits
	private int randomMap(double a, double b, double u){
		return (int)Math.floor(u * (double)(b - a));
	}

	private ArrayList<String> childs (String bsX, String bsY, String r){
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
	public String binaryStringMutation(String binaryString, double pMut, ArrayList<Double> u){
		
		// The number of random numbers provided must be equal to the number of bits on the binaryString
		assert(binaryString.length() == u.size());
				
		// The returning bynaryString
		String binaryString_mutated = binaryString;
		
		// for all the bits in the string, check if the corresponding random value is less than the probability of mutation, if so mutate the bit.
		for(int i=0; i<binaryString.length();i++)
			if(u.get(i)<pMut)
				binaryString_mutated = BinaryString.FLIP(binaryString_mutated, i);
		
		return binaryString_mutated;
		
	}
}


//Stochastic Universal Sampling
class SUS{
	// u is a random number and numberElements is the number of elements to keep
	public ArrayList<BinaryString> run(ArrayList<BinaryString> binaryStrings, double u, int numberElements){
		
		// List containing the selected elements
		ArrayList<BinaryString> selected = new ArrayList<BinaryString>();
		
		// Wheel Distribution list
		ArrayList<Double> wheelDistribution = generateWheelDistribution(binaryStrings);
		
		// Calculates population total fitness
		double totalFitness = calculateTotalFitness(binaryStrings);
		
		// Calculates the distance between the pointers (equal distance between all points) 
		double pointerDistance = totalFitness/numberElements;
		
		// Chooses the starting point of the random points
		double startingPointer = randomMap(0,pointerDistance,u);
		
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
	private Double calculateTotalFitness (ArrayList<BinaryString> binaryStrings){
		double totalFitness = 0;
		Iterator<BinaryString> iterator = binaryStrings.iterator();
		while(iterator.hasNext()){
			totalFitness += iterator.next().getFitness();
		}

		return totalFitness;
	}

	// Recieves a random number and the smallest and biggest index numbers, and returns a selected one
	private double randomMap(double a, double b, double u){
		return u * (double)(b - a);
	}

	// Given a list of BinaryStrings (with corresponding fitnesses), generates the ordered list with the corresponding wheel "slice" for each element. 
	// Offset taken into consideration so, each element corresponds to the maximum value between 0 and 1 that belongs to the indexed "slice"
	// Ex: if element 1 has p=0.15 and element 2 has p=0.20... The list will have at least 3 elements, being the first two 0.15 and 0.35. Where 0.35 = 0.15 + 0.20.
	private ArrayList<Double> generateWheelDistribution (ArrayList<BinaryString> binaryStrings){

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

//Can be improved, generating the Wheel Distribution can be done only once for all the population,
//it's not necessary to generate it every time the wheel is executed
class RouletteWheel{
	
	// O(3*n)
	public BinaryString run(ArrayList<BinaryString> binaryStrings, double u){
		
		// assert u is a random number
		assert(u>=0);
		assert(u<1);
		
		// Wheel Distribution list
		ArrayList<Double> wheelDistribution = generateWheelDistribution(binaryStrings);
		
		// Searches for the slice corresponding to the given random number;
		int index = 0;
		Iterator<Double> iterator = wheelDistribution.iterator();
		while(iterator.hasNext())
			// return the element corresponding to the currrent slice if the random number is bigger than the current slice limit, but wasn't bigger than the previous
			if(u<iterator.next())
				return binaryStrings.get(index);
			else
				index++;
		
		return null;
		
	}
	
	// Given a list of BinaryStrings (with corresponding fitnesses), generates the ordered list with the corresponding wheel "slice" for each element. 
	// Offset taken into consideration so, each element corresponds to the maximum value between 0 and 1 that belongs to the indexed "slice"
	// Ex: if element 1 has p=0.15 and element 2 has p=0.20... The list will have at least 3 elements, being the first two 0.15 and 0.35. Where 0.35 = 0.15 + 0.20.
	private ArrayList<Double> generateWheelDistribution (ArrayList<BinaryString> binaryStrings){
		
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
	private Double calculateTotalFitness (ArrayList<BinaryString> binaryStrings){
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
	// Runs the tournament for the elements on binaryStrings list, randomly chooses u.size elements, and from those winner elements choose the fittest as the ultimate winner
	public ArrayList<BinaryString> runTournaments(int N, int s, ArrayList<BinaryString> binaryStrings, ArrayList<Double> u){
		
		// List of winners
		ArrayList<BinaryString> winners = new ArrayList<BinaryString>();
		
		// You will	need to generate a random permutation of the population s times
		ArrayList<Permutation> permutations = generatePermutations(N,s,u);
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
	
	// Generate s permutations for the s groups of N-1 provided random values. # random values = 
	private ArrayList<Permutation> generatePermutations (int N, int s, ArrayList<Double> u){
		
		// You will	need to generate a random permutation of the population s times
		ArrayList<Permutation> permutations = new ArrayList<Permutation>();
		for(int i = 0; i<s;i++){
			
			// select the necessary random values (N-1) values
			ArrayList<Double> u_ = new ArrayList<Double>();
			int from = i*(N-1);
			int to = (i+1)*(N-1);
			u_.addAll(u.subList(from,to));
			
			//generate the permutation and add it to the list
			Permutation permutation = Permutation.generateRandom(N, u_);
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
	
	// Runs the tournament for the elements on binaryStrings list, randomly chooses u.size elements, and from those winner elements choose the fittest as the ultimate winner
	public BinaryString runTournament(ArrayList<BinaryString> binaryStrings, ArrayList<Double> u){
		ArrayList<BinaryString> randomBinaryStrings = new ArrayList<BinaryString>();
		
		// Selects u.size random binaryStrings
		for(int i = 0; i<u.size(); i++){
			randomBinaryStrings.add(this.selectRandomBinaryString(binaryStrings, u.get(i)));
		}
		
		return this.selectFittestBinaryString(randomBinaryStrings);
		
	}
	
	// Receives a random number and the smallest and biggest index numbers, and returns a selected one
	public int selectRandom(int a, int b, double u){
		return (int)(a + (Math.floor(u * (double)(b - a + 1))));
	}
	
	// Receives a list of BinaryString and a random number and returns the selected BinaryString
	public BinaryString selectRandomBinaryString(ArrayList<BinaryString> binaryStrings, double u){
		int index = this.selectRandom(1,binaryStrings.size(), u);
		return binaryStrings.get(index-1);
	}
	
	// Receives a non-empty list of BinaryString and returns the fittest element
	public BinaryString selectFittestBinaryString(ArrayList<BinaryString> binaryStrings){
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

	// generates a random permutation of size N, given a list of random numbers of size N-1
	public static Permutation generateRandom(int N, ArrayList<Double> u){

		assert(u.size() == N-1);

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
			int mappedPosition = randomMap(i
					, N-1, u.get(i));

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
	
	// Generates a population of random binaryStrings (based of function generate) when given
	// the number of elements N, the length l, and a list with the necessary number of random elements
	public static ArrayList<String> generatePopulation (int N, int l, ArrayList<Double> u){

		// Make sure the list contains the necessary number of random elements
		assert(u.size() == N*l);

		// List that will contain the population
		ArrayList<String> population = new ArrayList<String>();

		//Generate N elements
		for(int i=0; i<N; i++){

			// define the portion of the random elements list to retrieve
			// based on the String length of l
			int from = i*l;
			int to = (i+1)*l;
			ArrayList<Double> u_ = new ArrayList<Double>();
			u_.addAll(u.subList(from,to));

			// generate an element from the random list portion with size l
			population.add(generate(u_));
		}

		// Return the generated population
		return population;
	}

	// Receives a list with l random numbers and generates a binary string from those numbers
	// Random number u represents a 0 if u < 0.5 or a 1 if u>=0.5
	public static String generate(List<Double> u){

		// the default binary string
		String binaryString = "";

		// while there is random numbers, iterate through them (there should be l of them)
		Iterator<Double> iterator = u.iterator();
		double current;
		while(iterator.hasNext()){

			// gets the next random number
			current = iterator.next();

			// if the random number is smaller than 0.5 add 0 to the binary string, else add 1
			if(current < 0.5)
				binaryString = binaryString + "0";
			else
				binaryString = binaryString + "1";
		}

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

	// Receives a binaryString and calculates the number of 1 bits
	public static int onemax(String binaryString, int length){
		int count=0;

		for(int i = 0; i<length;i++)
			if(binaryString.charAt(i) == '0')
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

	

}
