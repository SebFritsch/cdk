//TODO: CDK copyright info goes here
package org.openscience.cdk.tools;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Queue;

import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.GraphUtil.EdgeToBondMap;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

//TODO class info goes here
public class FunctionalGroupsFinder {
	/**
	 * Defines the working mode.
	 */
	public static enum Mode{
		/**
		 * Default mode including generalization step
		 */
		DEFAULT,
		/**
		 * Skips the generalization step. Functional groups will keep their full "environment".
		 */
		NO_GENERALIZATION;
	}
	
	/**
	 * Flag to mark environmental C's in not-yet-generalized functional groups.
	 */
	private final static String ENVIRONMENTAL_C_FLAG = "ENVIRONMENTAL_C"; 
    
    /**
     * Default initial capacity of the list that takes the functional groups
     */
    private final static int DEFAULT_INITIAL_OUTPUT_CAPACITY = 8;
    
    /**
     * Initial size of the collections containing each functional group. 
     */
    private final static int FUNCTIONAL_GROUP_INITIAL_CAPACITY = 15;
    
    // Predicates
    // TODO ?
    
    private final Mode 		mode; 
    private final int		initialOutputCapacity;
    
    private EdgeToBondMap 		bondMap;
    private int[][] 			adjList;
    private HashSet<Integer>	markedAtoms;
    private List<Integer>		unmarkedCAtoms;
    private List<Integer> 		hAtoms;
    
    public FunctionalGroupsFinder() {
    	this(Mode.DEFAULT, DEFAULT_INITIAL_OUTPUT_CAPACITY);
    }
    
    public FunctionalGroupsFinder(Mode mode, int initialOutputCapacity) {
    	this.mode = mode;
    	this.initialOutputCapacity = initialOutputCapacity;
    }
    
    public  List<IAtomContainer> find(IAtomContainer molecule){
    	// init GraphUtil+EdgeToBondMap
    	bondMap = EdgeToBondMap.withSpaceFor(molecule);
    	adjList = GraphUtil.toAdjList(molecule, bondMap);
    	
    	markAtoms(molecule);
    	
    	// List<IAtomContainer> groupsWithEnvironment = partitionGroups(molecule);
    	
    	// if(mode) List<IAtomContainer> groupsGeneralized = generalizeGroups(groupsWithEnvironment);
    	
    	//FIXME
    	return new ArrayList<IAtomContainer>(initialOutputCapacity);
    }

    // TODO:	*	check if it is worth it to mark connected Atoms (Hetero, C=C, C≡C, 3-rings) right away and check if 
    // 				atom has already been marked at the beginning of the loop (HashSet<Integer idx>.contains)
    //			*	is is better to mark atoms by setting a property instead of using an external HashSet?
    private void markAtoms(IAtomContainer molecule) {
    	// TODO: FOR DEBUGGING ONLY!
    	// stores which condition in the algorithm got the atom(index) marked
    	Map<Integer, String> debugAtomConditionMap = new HashMap<>(molecule.getAtomCount());
    	
    	// classify marked, unmarked C's & rest (H's)
    	markedAtoms = new HashSet<>(molecule.getAtomCount());	// TODO: initial capacity!
    	unmarkedCAtoms = new ArrayList<>();						// TODO: initial capacity!	TODO: remove?!
    	hAtoms = new ArrayList<>();								// TODO: initial capacity!	TODO: remove?!
    	
    	for(int idx = 0; idx < molecule.getAtomCount(); idx++) {
    		IAtom cAtom = molecule.getAtom(idx);
    		int atomicNr = cAtom.getAtomicNumber();
    		
    		// NOTE: order assuming #Cs > #Hs > #HeteroAtoms
    		// if C...
    		if(atomicNr == 6) {
    			// TODO: see condition for marking C's
    			boolean isMarked = false;		// to detect if foor loop ran with or without marking the C atom
    			int oNSCounter = 0;				// count for the number of connected O, N & S atoms
    			for(int connectedIdx : adjList[idx]) {
    				//TODO: scratch connectedAtom, replace by connectedBond.getOther? 
    				IAtom connectedAtom = molecule.getAtom(connectedIdx); 
    				IBond connectedBond = bondMap.get(idx, connectedIdx);
    				
    				// if connected to Heteroatom or C in aliphatic double or triple bond... [CONDITIONS 2.1 & 2.2]
    				if(connectedAtom.getAtomicNumber() != 1 
    						&& ((connectedBond.getOrder() == Order.DOUBLE || connectedBond.getOrder() == Order.TRIPLE)
    						&& connectedBond.isAromatic() == false)) {
    					// set as marked and break out of connected atoms
    					debugAtomConditionMap.put(idx, "Condition 2.1/2.2");		// FIXME debug only!
    					isMarked = true;
    					break;
    				}
    				// if connected to O/N/S in single bond... (FIXME: & aliphatic?!)
    				else if((connectedAtom.getAtomicNumber() == 7 
    						|| connectedAtom.getAtomicNumber() == 8
    						|| connectedAtom.getAtomicNumber() == 16)
    						&& connectedBond.getOrder() == Order.SINGLE){
    					// if "acetal C" (2+ O/N/S in single bonds connected tp sp3-C)... [CONDITION 2.3]
    					oNSCounter++;
    					if(oNSCounter > 1 && adjList[idx].length == 4) {
    						// set as marked and break out of connected atoms
    						debugAtomConditionMap.put(idx, "Condition 2.3");		// FIXME debug only!
    						isMarked = true;
    						break;
    					}
    					
    					// if part of oxirane, aziridine and thiirane ring... [CONDITION 2.4]
    					// FIXME not the best way to check every C like this!
    					for(int connectedInSphere2Idx : adjList[connectedIdx]) {
    						IAtom connectedInSphere2Atom = molecule.getAtom(connectedInSphere2Idx);
    						if(connectedInSphere2Atom.getAtomicNumber() == 6) {
    							for(int connectedInSphere3Idx : adjList[connectedInSphere2Idx]) {
    								IAtom connectedInSphere3Atom = molecule.getAtom(connectedInSphere3Idx);
    								if(connectedInSphere3Atom.equals(cAtom)) {
    									// set as marked and break out of connected atoms
    									debugAtomConditionMap.put(idx, "Condition 2.4");		// FIXME debug only!
    									isMarked = true;
    									break;
    								}
    							}
    						}
    					}
    				}
    			}
    			if(isMarked) {
    				markedAtoms.add(idx);
    				continue;
    			}
    			// if none of the conditions 2.X apply, store as unmarked C for now
    			else {
    				unmarkedCAtoms.add(idx);
    			}
    			
    		}
    		// if H...
    		else if (atomicNr == 1){
    			// TODO: is collecting H's any useful?
    			hAtoms.add(idx);
    			continue;
    		}
    		// if Hereroatom... (CONDITION 1)
    		else {
    			debugAtomConditionMap.put(idx, "Condition 1");		// FIXME debug only!
    			markedAtoms.add(idx);
    			continue;
    		}
    	}
    	
    	//DEBUG CODE:
    	printMarkedAtomsDebugInfo(debugAtomConditionMap, molecule);
    }
    
    // TODO:	*	AtomContainerManipulator.extractSubstructure(IAtomContainer atomContainer, int... atomIndices)
    //				vs. cloning mol., deleting bonds/atoms & using ConnectivityChecker.partitionIntoMolecules(disconnectedContainer) ?
    //			* 	DFS should be more memory efficient than BFS, but does not work with a boolean lastMarked. 
    public List<IAtomContainer> extractGroups(IAtomContainer molecule) throws CloneNotSupportedException{
    	List<IAtomContainer> functionalGroups = new ArrayList<>(initialOutputCapacity);
    	
    	while(!markedAtoms.isEmpty()) { 				// TODO: = for with iterator?
    		// get next markedAtom as the starting node for the search 
    		int beginIdx = markedAtoms.iterator().next();
    		List<Integer> fGroupIndices = new ArrayList<>(FUNCTIONAL_GROUP_INITIAL_CAPACITY); //TODO: better to use #atoms-size int[] ?
    		
    		System.out.println("###### SEARCHING NEW FUNCTIONAL GROUP FROM: "+ molecule.getAtom(beginIdx).getSymbol() + beginIdx);
    		
    		// do a BFS from there
    		boolean[] visited = new boolean[molecule.getAtomCount()];
    		Queue<Integer> queue = new ArrayDeque<>();
    		queue.add(beginIdx);
    		visited[beginIdx] = true;
    		boolean lastMarked = true;
    		
    		while(!queue.isEmpty()) {
    			Integer idx = queue.poll();
    			// if we find a marked atom ...
    			if(markedAtoms.contains(idx)) {
    				
    				System.out.println("VISITED MARKED ATOM: "+ molecule.getAtom(idx).getSymbol() + idx);
    				
    				// add its index to the functional group
    				fGroupIndices.add(idx);
    				// also scratch the index from markedAtoms
    				markedAtoms.remove(idx);
    				// make note that we had a marked atom
    				if(!lastMarked) lastMarked = true;
    				
    				System.out.println("SET LASTMARKED: : "+ lastMarked);
    				
    				// and take look at the connected atoms
    				for(int connectedIdx : adjList[idx]) {
    					if(!visited[connectedIdx]) {
    						queue.add(connectedIdx);
    						visited[connectedIdx] = true;
    					}
    				}
    			}
    			// if we find a C and the atom we came from was marked (= environmental C)...
    			else if(lastMarked && molecule.getAtom(idx).getAtomicNumber() == 6) {
    				
    				System.out.println("VISITED ENVIRONMENTAL C ATOM: "+ molecule.getAtom(idx).getSymbol() + idx);
    				
    				// add its index to the functional group
    				fGroupIndices.add(idx);
    				// if generalization is wanted later on, flag as environmental C
    				if(mode != Mode.NO_GENERALIZATION) {
    					molecule.getAtom(idx).setProperty(ENVIRONMENTAL_C_FLAG, true);
    				}
    				// make note that this was not a marked atom
    				if(lastMarked) lastMarked = false;
    				// end here and do not look at connected atoms
    				
    				System.out.println("SET LASTMARKED: : "+ lastMarked);
    			}
    			else {
    				System.out.println("VISITED NON_MARKED ATOM: "+ molecule.getAtom(idx).getSymbol() + idx);
    			}
    		}
    		System.out.println("###### END OF FUNCTIONAL GROUP");
    		
    		// extract functional group from the collected indices
    		// FIXME: not a nice workaround, see todo above! 
    		int[] fGroupIndicesArray = fGroupIndices.stream().mapToInt(i->i).toArray();
    		IAtomContainer fGroup = AtomContainerManipulator.extractSubstructure(molecule, fGroupIndicesArray);
    		
    		// create ID
    		IDCreator.createIDs(fGroup);
    		
    		// add functional group to list
    		functionalGroups.add(fGroup);
    	}
    	
    	return functionalGroups;
    }
    
    // TODO:	*	ATM: environmental C's are identified by an attached String property-flag
    //				option 1: store all of them (across all FGs) in some data structure (HashMap?)
    //				option 2: introduce inner class functionalGroup that stores an AC plus a HashSet with the env.C's idices.
    //				option 3: do not remove marked atoms from map and search for C's not in the map.
    //			*	is it faster to do pattern matching? see class Pattern & VentoFoggia
    //			*	is it faster to do (unique) SMILES-String (regex) matching?
    public void generalizeEnvironments(List<IAtomContainer> fGroups) {
    	for(IAtomContainer fGroup : fGroups) {
    		int atomCount = fGroup.getAtomCount();
    		
    		// STEP 1: 		Remove environments on marked C's
    		// EXCEPTION: 	Substisubstituents on carbonyl (to distinguish aldehydes & ketones)
    		boolean isAldehydeOrKetone = false;
    		if(atomCount == 3 || atomCount == 4) {
    			for(IBond bond : fGroup.bonds()) {
    				if(bond.getOrder() == Order.DOUBLE 
    						// FIXME: is there a simpler way?
    						&& (bond.getBegin().getAtomicNumber() + bond.getEnd().getAtomicNumber() == 14 
    						&& (bond.getBegin().getAtomicNumber() == 8 || bond.getBegin().getAtomicNumber() == 6))) {
    					isAldehydeOrKetone = true;
    				}
    			}
    		}
    		if(!isAldehydeOrKetone) {
    			for(IAtom atom : fGroup.atoms()) {
    				if(atom.getProperty(ENVIRONMENTAL_C_FLAG) != null) {
    					fGroup.removeAtom(atom);
    				}
    			}
    		}
    		
    		// STEP 2: 		Fill valences on heteroatoms by R-Atoms (representing C or H)
    		// EXCEPTION 1:	H's on -OH
    		// EXCEPTION 2:	H's on simple amines and thiols (to distinguish secondary & tertiary amines and thiols & sulfides)

    		// TODO!!!

    	}
    	
    	// FIXME: remove when finished
    	throw new UnsupportedOperationException("Method is WIP!");
    }

    private static final boolean isHeteroatom(IAtom atom) {
    	// TODO does this make any difference to calling getAtomicNumber() each time?
    	int atomicNr = atom.getAtomicNumber();
    	return atomicNr != 1 && atomicNr != 6;
    }
    
    // TODO: for testing only, remove!
    public static IAtom[] getConnectedAtoms(int idx,IAtomContainer molecule, int[][] adjList) {
    	IAtom[] connectedAtoms = new IAtom[adjList[idx].length];
    	int counter = 0;
    	for(int connectedIdx : adjList[idx]) {
    		connectedAtoms[counter] = molecule.getAtom(connectedIdx);
    		counter++;
    	}
    	return connectedAtoms;
    }
    
    // FOR DEBUGGING ONLY!
    public static void printMarkedAtomsDebugInfo(Map<Integer, String> debugAtomConditionMap, IAtomContainer mol) {
    	System.out.println("################# ATOM CONDITION DATA ###################");
    	System.out.println("Size check... " + (mol.getAtomCount() == debugAtomConditionMap.size() ? "DONE" : "NUMBER OF ATOMS DOES NOT MATCH!"));
    	for(Map.Entry<Integer, String> e : debugAtomConditionMap.entrySet()) {
    		System.out.printf("#%d:%s   %s%n", e.getKey(), mol.getAtom(e.getKey()).getSymbol(), e.getValue());
    	}
    	System.out.println("##########################################################");
    }
    

}
