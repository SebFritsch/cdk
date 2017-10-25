//TODO: CDK copyright info goes here
package org.openscience.cdk.tools;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;

import org.openscience.cdk.Atom;
import org.openscience.cdk.PseudoAtom;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.GraphUtil.EdgeToBondMap;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

//TODO class info goes here
public class FunctionalGroupsFinder {
	private static ILoggingTool log = LoggingToolFactory.createLoggingTool(FunctionalGroupsFinder.class);
	
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
	 * Defines whether an environmental C is connected to a heteroatom or another C.
	 */
	private static enum EnvironmentalCType{ C_ON_HETEROATOM, C_ON_C };
	
	/**
	 * Flag to mark environmental C's in not-yet-generalized functional groups.
	 */
	private final static String ENVIRONMENTAL_C_FLAG = "ENVIRONMENTAL_C"; 
    
    /**
     * Initial size of the collections containing each functional group. 
     */
    private final static int FUNCTIONAL_GROUP_INITIAL_CAPACITY = 15;
    
    private final Mode 		mode; 
    
    private EdgeToBondMap 		bondMap;
    private int[][] 			adjList;
    private HashSet<Integer>	markedAtoms;
    
    /**
     * Default constructor for FunctionalGroupsFinder.
     */
    public FunctionalGroupsFinder() {
    	this(Mode.DEFAULT);
    }
    
    /**
     * Constructor for FunctionalGroupsFinder.
     * 
     * @param mode working mode (see {@code FunctionalGroupsFinder.Mode}).
     */
    public FunctionalGroupsFinder(Mode mode) {
    	this.mode = mode;
    }
    
    /**
     * Find all functional groups contained in a molecule.
     * 
     * @param molecule the molecule which contains the functional groups
     * @return a list with all functional groups found in the molecule.
     * @throws CloneNotSupportedException an exception thrown if the molecule cannot be cloned
     */
    public  List<IAtomContainer> find(IAtomContainer molecule) throws CloneNotSupportedException{
    	// init GraphUtil+EdgeToBondMap
    	bondMap = EdgeToBondMap.withSpaceFor(molecule);
    	adjList = GraphUtil.toAdjList(molecule, bondMap);
    	
    	// TODO!
    	
//    	markAtoms(molecule);
    	
//    	List<IAtomContainer> fGroups = extractGroups(molecule);
    	
//    	if generalization is wanted do it
    	
    	//FIXME replace!
    	return new LinkedList<IAtomContainer>();
    }

    /**
     * TODO / NOTES
     * 
     * TODO:	*	change public -> private
     *			*	check if it is worth it to mark connected Atoms (Hetero, C=C, C≡C, 3-rings) right away and check if 
     * 				atom has already been marked at the beginning of the loop (HashSet<Integer idx>.contains)
     *			*	is is better to mark atoms by setting a property instead of using an external HashSet?
     */
    /**
     * Mark all atoms according to the conditions below and store them in a set for further processing.
     * 
     * 	1. mark all heteroatoms in a molecule, including halogens
     * 	2. mark also the following carbon atoms:
     * 		• atoms connected by non-aromatic double or triple bond to any heteroatom
     * 		• atoms in nonaromatic carbon–carbon double or triple bonds
     * 		• acetal carbons, i.e. sp3 carbons connected to two or more oxygens, nitrogens or sulfurs; these O, N or S atoms must have only single bonds
     *		• all atoms in oxirane, aziridine and thiirane rings (such rings are traditionally considered to be functional groups due to their high reactivity)
     * TODO: add citation!
     * 
     * @param molecule molecule with atoms to mark
     */
    public void markAtoms(IAtomContainer molecule) {
    	log.debug("########## Starting search for atoms to mark ... ##########");
    	
    	// store marked atoms
    	markedAtoms = new HashSet<>(molecule.getAtomCount());	// TODO: initial capacity!
    	
    	for(int idx = 0; idx < molecule.getAtomCount(); idx++) {
    		IAtom cAtom = molecule.getAtom(idx);
    		int atomicNr = cAtom.getAtomicNumber();
    		
    		// NOTE: order assuming #Cs > #Hs > #HeteroAtoms
    		// if C...
    		if(atomicNr == 6) {
    			boolean isMarked = false;		// to detect if foor loop ran with or without marking the C atom
    			int oNSCounter = 0;				// count for the number of connected O, N & S atoms
    			for(int connectedIdx : adjList[idx]) {
    				//TODO: scratch connectedAtom, replace by connectedBond.getOther? 
    				IAtom connectedAtom = molecule.getAtom(connectedIdx); 
    				IBond connectedBond = bondMap.get(idx, connectedIdx);
    				
    				// if connected to Heteroatom or C in aliphatic double or triple bond... [CONDITIONS 2.1 & 2.2]
    				if(connectedAtom.getAtomicNumber() != 1 
    						&& ((connectedBond.getOrder() == Order.DOUBLE || connectedBond.getOrder() == Order.TRIPLE)
    						&& !connectedBond.isAromatic())) {
    					// set as marked and break out of connected atoms
    					log.debug(String.format("Marking Atom #%d (%s) - Met condition 2.1/2.2", idx, cAtom.getSymbol())); 
    					isMarked = true;
    					break;
    				}
    				// if connected to O/N/S in single bond...
    				else if((connectedAtom.getAtomicNumber() == 7 
    						|| connectedAtom.getAtomicNumber() == 8
    						|| connectedAtom.getAtomicNumber() == 16)
    						&& connectedBond.getOrder() == Order.SINGLE){
    					// if "acetal C" (2+ O/N/S in single bonds connected tp sp3-C)... [CONDITION 2.3]
    					oNSCounter++;
    					if(oNSCounter > 1 && adjList[idx].length == 4) {
    						// set as marked and break out of connected atoms
    						log.debug(String.format("Marking Atom #%d (%s) - Met condition 2.3", idx, cAtom.getSymbol())); 
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
    									log.debug(String.format("Marking Atom #%d (%s) - Met condition 2.4", idx, cAtom.getSymbol())); 
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
    			// if none of the conditions 2.X apply, we have an unmarked C (not relevant here)
    		}
    		// if H...
    		else if (atomicNr == 1){
    			// convert to implicit H
    			IAtom connectedAtom = molecule.getAtom(adjList[idx][0]);
    			
    			if(connectedAtom.getImplicitHydrogenCount() == null) {
    				connectedAtom.setImplicitHydrogenCount(1);
    			}
    			else {
    				connectedAtom.setImplicitHydrogenCount(connectedAtom.getImplicitHydrogenCount() + 1);
    			}
    			continue;
    		}
    		// if heteroatom... (CONDITION 1)
    		else {
    			log.debug(String.format("Marking Atom #%d (%s) - Met condition 1", idx, cAtom.getSymbol())); 
    			markedAtoms.add(idx);
    			continue;
    		}
    	}
    	log.debug(String.format("########## End of search. Marked %d/%d atoms. ##########", markedAtoms.size(), molecule.getAtomCount()));
    }
    
    /**
     * TODO / NOTES
     * 
     * TODO:	*	change public -> private
     *			*	Environmental C's & type could better be stored in HashMap instead of atom's property!
     *			*	AtomContainerManipulator.extractSubstructure(IAtomContainer atomContainer, int... atomIndices)
     *				vs. cloning mol., deleting bonds/atoms & using ConnectivityChecker.partitionIntoMolecules(disconnectedContainer) ?
     *			* 	DFS should be more memory efficient than BFS
     * NOTE:	*	returns groups with implicit hydrogens!
     */
    /**
     * Searches the molecule for groups of connected marked atoms and extracts each as a new functional group.
     * The extraction includes connected unmarked carbon atoms that are not part of the functional group
     * but form its "environment".
     * 
     * @param molecule the molecule which contains the functional groups
     * @return a list of all functional groups (including "environments") extracted from the molecule
     * @throws CloneNotSupportedException an exception thrown if the molecule cannot be cloned
     */
    public List<IAtomContainer> extractGroups(IAtomContainer molecule) throws CloneNotSupportedException{
    	log.debug("########## Starting identification & extraction of functional groups... ##########");
    	
    	List<IAtomContainer> functionalGroups = new LinkedList<>();
    	while(!markedAtoms.isEmpty()) {
    		// get next markedAtom as the starting node for the search 
    		int beginIdx = markedAtoms.iterator().next();
    		List<Integer> fGroupIndices = new ArrayList<>(FUNCTIONAL_GROUP_INITIAL_CAPACITY);
    		
    		log.debug(String.format("Searching new functional group from atom #%d (%s)...", beginIdx,  molecule.getAtom(beginIdx).getSymbol()));
    		
    		// do a BFS from there
    		boolean[] visited = new boolean[molecule.getAtomCount()];
    		Queue<Integer> queue = new ArrayDeque<>();
    		queue.add(beginIdx);
    		visited[beginIdx] = true;
    		boolean isParentC = false; // value does not matter, will be reset anyway
    		int resetParentCount = 0;
    		
    		while(!queue.isEmpty()) {
    			Integer idx = queue.poll();
    			IAtom currentAtom = molecule.getAtom(idx);
    			
//    			if(mode != Mode.NO_GENERALIZATION) {
//    				// handle identification of whether or not the current parent is a C
//    				if(resetParentCount < 2) {
//    					isParentC = currentAtom.getAtomicNumber() == 6;
//    					resetParentCount = adjList[idx].length;
//    					log.debug("		resetParentCount elapsed. Checking atom for C: " + currentAtom.getSymbol() + idx + " -> " + isParentC + ". New resetParentCount: " + resetParentCount); // FIXME debug only
//    				}
//    				else {
//    					resetParentCount--;
//    					log.debug("		resetParentCount: " + resetParentCount + "(isParentC: " + isParentC + ")"); // FIXME debug only
//    				}
//    			}
    			
    			// if we find a marked atom ...
    			if(markedAtoms.contains(idx)) {
    				log.debug(String.format("	visited marked atom: #%d (%s)", idx, currentAtom.getSymbol()));
    				
    				// add its index to the functional group
    				fGroupIndices.add(idx);
    				// also scratch the index from markedAtoms
    				markedAtoms.remove(idx);
    				
    				// and take look at the connected atoms
    				for(int connectedIdx : adjList[idx]) {
    					if(!visited[connectedIdx]) {
    						queue.add(connectedIdx);
    						visited[connectedIdx] = true;
    					}
    				}
    			}
    			// if we find a C...
    			else if(currentAtom.getAtomicNumber() == 6) {
    				// add its index to the functional group
    				fGroupIndices.add(idx);
    				// if generalization is wanted later on, flag as environmental C (either connected to another C or a heteroatom)
    				if(mode != Mode.NO_GENERALIZATION) {
    					EnvironmentalCType type = isParentC ? EnvironmentalCType.C_ON_C : EnvironmentalCType.C_ON_HETEROATOM;
    					currentAtom.setProperty(ENVIRONMENTAL_C_FLAG, type);
    				}
    				// end here and do not look at connected atoms
    				
    				log.debug(String.format("	visited env. C atom: #%d (%s), type: %s", idx, currentAtom.getSymbol(), isParentC ? "C_ON_C" : "C_ON_HETEROATOM"));
    			}
    			else {
    				log.debug(String.format("	visited non-marked atom: #%d (%s)", idx, currentAtom.getSymbol()));
    			}
    			
    			if(mode != Mode.NO_GENERALIZATION) {
    				// handle identification of whether or not the current parent is a C
    				if(resetParentCount < 2) {
    					isParentC = currentAtom.getAtomicNumber() == 6;
    					resetParentCount = adjList[idx].length;
//    					log.debug("		resetParentCount elapsed. Checking atom for C: " + currentAtom.getSymbol() + idx + " -> " + isParentC + ". New resetParentCount: " + resetParentCount); // FIXME debug only
    				}
    				else {
    					resetParentCount--;
//    					log.debug("		resetParentCount: " + resetParentCount + "(isParentC: " + isParentC + ")"); // FIXME debug only
    				}
    			}
    		}
    		log.debug("	search completed.");
    		
    		// extract functional group from the collected indices
    		// FIXME: workaround, see todo above! 
    		int[] fGroupIndicesArray = fGroupIndices.stream().mapToInt(i->i).toArray();
    		IAtomContainer fGroup = AtomContainerManipulator.extractSubstructure(molecule, fGroupIndicesArray);
    		
    		log.debug("Extracting functional group by atom indices: ", Arrays.toString(fGroupIndicesArray));
    		
    		// create IDs for the group
    		IDCreator.createIDs(fGroup);
    		
    		// add functional group to list
    		functionalGroups.add(fGroup);
    	}
    	log.debug(String.format("########## Found & extracted %d functional groups. ##########", functionalGroups.size()));
    	
    	return functionalGroups;
    }
    
    /**
     * TODO / NOTES
     * 
     * TODO:	*	change to private void!
     *			*	ATM: environmental C's are identified by an attached String property-flag
     *				option 1: store all of them (across all FGs) in a list of Hashmaps (idx - type)
     *				option 2: introduce inner class functionalGroup that stores an AC plus the HashSet.
     *			*	is it faster to do pattern matching? see class Pattern & VentoFoggia
     */
    /**
     * Generalizes the full environments of functional groups, providing a good balance between preserving 
     * meaningful detail and generalization.
     * 
     * @param fGroups the list of functional groups including full "environments"
     * @return a list of the same functional groups with generalized "environments"
     */
    public List<IAtomContainer> generalizeEnvironments(List<IAtomContainer> fGroups){
    	log.debug("########## Starting generalization of functional groups... ##########");
    	
    	fGroupLoop:
    	for(IAtomContainer fGroup : fGroups) {
    		int atomCount = fGroup.getAtomCount();
    		
    		log.debug(String.format("Generalizing functional group (%d atoms)...", atomCount));
    		
    		// PRECHECKING FOR EXCEPTION CASES
    		boolean isAldehydeOrKetone = false, isSingleO = false, isSecAmineOrSimpleThiol = false, isSingleN = false;
    		
    		if(atomCount == 2) {
    			IBond bond = fGroup.bonds().iterator().next();
    			if(bond.getOrder() == Order.SINGLE) {
    				int atomicNumberSum = bond.getBegin().getAtomicNumber() + bond.getEnd().getAtomicNumber();
    				switch(atomicNumberSum) {
    				case 13: 	isSingleN = true; // C-N
    				break;
    				case 14:	isSingleO = true; // C-O
    				break;
    				case 22:	isSecAmineOrSimpleThiol = true; // C-S
    				break;
    				}
    			}
    		}
    		else if(atomCount == 3 || atomCount == 4) {
    			byte bondOrderSum = 0;
    			for(IBond bond : fGroup.bonds()) {
    				// check for C=O double bond
    				if(bond.getOrder() == Order.DOUBLE && bond.getBegin().getAtomicNumber() + bond.getEnd().getAtomicNumber() == 14) {
    					isAldehydeOrKetone = true; // aldehyde or ketone
    					continue;
    				}
    				bondOrderSum += bond.getOrder().numeric();
    			}
    			// if its neither aldehyde or ketone , all two bonds (-> 3 atoms) are single bonds...
    			if(!isAldehydeOrKetone && bondOrderSum == 2) {
    				int atomicNumberSum = 0;
    				for(IAtom atom : fGroup.atoms()) {
    					atomicNumberSum += atom.getAtomicNumber();
    				}
    				if(atomicNumberSum == 19) {
    					isSecAmineOrSimpleThiol = true; // sec. amide 
    				}
    			}
    		}
    		log.debug(String.format("	precheck result: 	isSingleO: %b, isSingleN: %b, isSecAmineOrSimpleThiol: %b, isAldehydeOrKetone: %b",
    				isSingleO, isSingleN, isSecAmineOrSimpleThiol, isAldehydeOrKetone));

    		
    		// PROCESSING
    		if(isSingleN || isSingleO || isSecAmineOrSimpleThiol) {
    			for(IAtom atom : fGroup.atoms()) {
    				if(atom.getAtomicNumber() != 6) {
    					restoreExplicitHydrogens(fGroup, atom);
    					if(!isSecAmineOrSimpleThiol) {
    						continue fGroupLoop; // continue because we want to keep the C's in C-OH & C-NH2
    					}
    				}
    			}
    		}
    		for(IAtom atom : fGroup.atoms()) {
    			if(atom.getAtomicNumber() == 6) {
    				// STEP 1: delete environments on carbons... see exceptions!
    				EnvironmentalCType type = atom.getProperty(ENVIRONMENTAL_C_FLAG);
    				if(type == EnvironmentalCType.C_ON_C && !isAldehydeOrKetone) {
    					log.debug("	removing env. C (type C_ON_C)...");
    					fGroup.removeAtom(atom);
    					continue;
    				}
    				// STEP 3: replace other environmental C's by R-Atoms (represent H or C)... see exceptions!
    				// TODO: double check exception for carbonyls!
    				else if(type == EnvironmentalCType.C_ON_HETEROATOM && !isSingleN && !isSingleO) {
    					IPseudoAtom rAtom = new PseudoAtom("R");
    					rAtom.setAttachPointNum(1);
    					AtomContainerManipulator.replaceAtomByAtom(fGroup, atom, rAtom);
    					log.debug("	replacing env. C (type C_ON_HETEROATOM) by R-atom...");
    				}
    			}
    			// STEP 2: fill free valences on heteroatoms by R-Atoms... see exceptions!
    			else if(isHeteroatom(atom)){
    				if(atom.getImplicitHydrogenCount() != null) {
    					int rAtomsToAdd = atom.getImplicitHydrogenCount();
    					for(int i = 0; i < rAtomsToAdd; i++) {
    						IPseudoAtom rAtom = new PseudoAtom("R");
    						rAtom.setAttachPointNum(1);
    						IBond bond = atom.getBuilder().newBond();
    						bond.setAtoms(new IAtom[] {rAtom, atom});
    						bond.setOrder(Order.SINGLE); // TODO: necessary?
    						fGroup.addAtom(rAtom);
    						fGroup.addBond(bond);
    						log.debug("	filling valence on heteroatom with R-atom...");
    					}
    				}
    			}
    		}
    	}
    	log.debug("########## Generalization of functional groups completed. ##########");
    
    	return fGroups;
    }
    
    /**
     * Returns true if an atom is a heteroatom (not hydrogen or carbon).
     * 
     * @param atom atom to check
     * @return true if it is a heteroatom
     */
    private static final boolean isHeteroatom(IAtom atom) {
    	int atomicNr = atom.getAtomicNumber();
    	return atomicNr != 1 && atomicNr != 6;
    }
    
    // TODO: optimize?
    /**
     * Converts implicit hydrogens to explicit for a single atom in a molecule.
     * 
     * @param molecule molecule in which the conversion should happen
     * @param atom atom with implicit hydrogens that should be converted
     */
    private static void restoreExplicitHydrogens(IAtomContainer molecule, IAtom atom) {
    	int formerAtomCount = molecule.getAtomCount();
    	for(int i = 0; i < atom.getImplicitHydrogenCount(); i++) {
    		IAtom hydrogen = new Atom(1);
    		IBond bond = atom.getBuilder().newBond();
    		bond.setAtoms(new IAtom[]{hydrogen, atom});
    		bond.setOrder(Order.SINGLE); // TODO: necessary?
    		molecule.addAtom(hydrogen);
    		molecule.addBond(bond);
    	}
    	log.debug("	restored "+ (molecule.getAtomCount()-formerAtomCount) +" explicit hydrogens at "+atom.getSymbol()+molecule.indexOf(atom)); // FIXME debug only
    	atom.setImplicitHydrogenCount(0);
    }
}
