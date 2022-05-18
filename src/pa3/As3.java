package pa3;

import java.util.*;

public class As3 {
	public static void main(String[] args) {

		Scanner sc = new Scanner(System.in);

		// Read N
		int N = sc.nextInt();

		// Read N-1 random numbers
		ArrayList<Double> u = new ArrayList<Double>();
		for(int i=0;i<N-1;i++)
			u.add(sc.nextDouble());

		// Randomly generates the permutation
		PermutationManipulator permM = new PermutationManipulator();
		Permutation permutation = permM.generateRandom(N,u);

		// Print the permutation
		permutation.print_1perLine();

	}	

}

class Permutation {
	
	private ArrayList<Integer> permutation;
	private int range;
	
	/* CONSTRUCTORS */
	public Permutation(){
		this.permutation = new ArrayList<Integer>();
	}
	public Permutation(ArrayList<Integer> permutation){
		this.permutation = permutation;
	}
	public Permutation(ArrayList<Integer> permutation, int range){
		this.permutation = permutation;
		this.range = range;
	}
	
	/* GETTERS & SETTERS */
	public ArrayList<Integer> getPermutation() {
		return permutation;
	}
	public void setPermutation(ArrayList<Integer> permutation) {
		this.permutation = permutation;
	}
	
	/* FUNCTIONS */
	// prints one element of the permutation per line
	public void print_1perLine(){
		//
		for(int i=0;i<permutation.size();i++)
			System.out.println(permutation.get(i));
	}
	
	// swaps the elements of the refered positions in the permutation list
	public void swap(int a, int b){
		int aux = permutation.get(a);
		permutation.set(a, permutation.get(b));
		permutation.set(b,aux);
	}
}

class PermutationManipulator {
	
	public Permutation generateRandom(int N, ArrayList<Double> u){
		
		//assert(N == u.size()+1);
		
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
			int mappedPosition = this.randomMap(i, N-1, u.get(i));
			
			// exchange the contents of v[i] with the contents of v[r]
			permutation.swap(i, mappedPosition);
		}
			
		//returns the random permutation
		return permutation;
	}
	
	// Recieves a random number and the smaller and larger limits. Maps the random number to an integer between those limits
	private int randomMap(double a, double b, double u){
		return (int)(a + (Math.floor(u * (double)(b - a + 1))));
	}
	
}
