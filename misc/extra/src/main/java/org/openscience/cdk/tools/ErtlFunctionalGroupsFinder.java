//TODO: CDK copyright info goes here
package org.openscience.cdk.tools;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.PseudoAtom;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.GraphUtil.EdgeToBondMap;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import com.google.common.primitives.Ints;

//TODO class info goes here
public class ErtlFunctionalGroupsFinder {
	private static ILoggingTool log = LoggingToolFactory.createLoggingTool(ErtlFunctionalGroupsFinder.class);
	
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
	//private static enum EnvironmentalCType{ C_ON_HETEROATOM, C_ON_C, CARBONYL_C};
	private static enum EnvironmentalCType{ C_ON_HETEROATOM, C_ON_C};
	
	private static enum searchParentType{HETEROATOM, C, CARBONYL_C};
	
	/**
	 * Flag to mark environmental C's in not-yet-generalized functional groups.
	 */
	private final static String ENVIRONMENTAL_C_FLAG = "ENVIRONMENTAL_C"; 
	
	/**
	 * Flag to mark carbons in carbonyl groups.
	 */
	private final static String CARBONYL_C_FLAG = "CARBONYL_C"; 
    
    /**
     * Initial size of the collections containing each functional group. 
     */
    private final static int FUNCTIONAL_GROUP_INITIAL_CAPACITY = 15;
    
    private final Mode 		mode; 
    
    private EdgeToBondMap 				bondMap;
    private int[][] 					adjList;
    private HashSet<Integer>			markedAtoms;
    private HashMap<Integer, Boolean>	aromaticHeteroAtoms; // <atomIdx, isInGroup>
    
    /**
     * Default constructor for FunctionalGroupsFinder.
     */
    public ErtlFunctionalGroupsFinder() {
    	this(Mode.DEFAULT);
    }
    
    /**
     * Constructor for FunctionalGroupsFinder.
     * 
     * @param mode working mode (see {@code FunctionalGroupsFinder.Mode}).
     */
    public ErtlFunctionalGroupsFinder(Mode mode) {
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
    	markedAtoms = new HashSet<>(Math.max((int) (molecule.getAtomCount()/.75f) + 1, 16));
    	// store aromatic heteroatoms
    	aromaticHeteroAtoms = new HashMap<>();
    	
    	for(int idx = 0; idx < molecule.getAtomCount(); idx++) {
    		// skip atoms that already got marked in a previous iteration
    		if(markedAtoms.contains(idx)) {
    			continue;
    		}
    		IAtom cAtom = molecule.getAtom(idx);
    		// skip aromatic atoms and add aromatic heteroatoms to map
    		if(cAtom.isAromatic()) {
    			if(isHeteroatom(cAtom)) {
    				aromaticHeteroAtoms.put(idx, false);
    			}
    			continue;
    		}
    		int atomicNr = cAtom.getAtomicNumber();
    		
    		// NOTE: order assuming #Cs > #Hs > #HeteroAtoms
    		// if C...
    		if(atomicNr == 6) {
    			boolean isMarked = false;		// to detect if foor loop ran with or without marking the C atom
    			int oNSCounter = 0;				// count for the number of connected O, N & S atoms
    			for(int connectedIdx : adjList[idx]) {
    				IAtom connectedAtom = molecule.getAtom(connectedIdx); 
    				IBond connectedBond = bondMap.get(idx, connectedIdx);
    				
    				// if connected to Heteroatom or C in aliphatic double or triple bond... [CONDITIONS 2.1 & 2.2]
    				if(connectedAtom.getAtomicNumber() != 1 
    						&& ((connectedBond.getOrder() == Order.DOUBLE || connectedBond.getOrder() == Order.TRIPLE)
    						&& !connectedBond.isAromatic())) {
    					// set the connected atom as marked
    					String connectedAtomCondition = connectedAtom.getAtomicNumber() == 6 ? "2.1/2.2" : "1";
    					log.debug(String.format("Marking Atom #%d (%s) - Met condition %s", connectedIdx, connectedAtom.getSymbol(), connectedAtomCondition));
    					markedAtoms.add(connectedIdx);
    					// set the current atom as marked and break out of connected atoms
    					log.debug(String.format("Marking Atom #%d (%s) - Met condition 2.1/2.2", idx, cAtom.getSymbol())); 
    					isMarked = true;
    					
    					
    					if(connectedAtom.getAtomicNumber() == 8 && connectedBond.getOrder() == Order.DOUBLE && adjList[idx].length == 3) {
    						log.debug("                     - was flagged as Carbonly-C");
    						cAtom.setProperty(CARBONYL_C_FLAG, true);
    					}
    					
    					
    					break;
    				}
    				// if connected to O/N/S in single bond...
    				else if((connectedAtom.getAtomicNumber() == 7 
    						|| connectedAtom.getAtomicNumber() == 8
    						|| connectedAtom.getAtomicNumber() == 16)
    						&& connectedBond.getOrder() == Order.SINGLE){
    					// if connected O/N/S is not aromatic...
    					if(!connectedAtom.isAromatic()) {
    						// set the connected O/N/S atom as marked
    						log.debug(String.format("Marking Atom #%d (%s) - Met condition 1", connectedIdx, connectedAtom.getSymbol()));
    						markedAtoms.add(connectedIdx);
    						
    						// if "acetal C" (2+ O/N/S in single bonds connected to sp3-C)... [CONDITION 2.3]
    						boolean isAllSingleBonds = true;
    						for(int connectedInSphere2Idx : adjList[connectedIdx]) {
    							IBond sphere2Bond = bondMap.get(connectedIdx, connectedInSphere2Idx);
    							if(sphere2Bond.getOrder() != Order.SINGLE) {
    								isAllSingleBonds = false;
    								break;
    							}
    						}
    						if(isAllSingleBonds) {
    							oNSCounter++;
    							// FIXME: getImplicitHydrogenCount can be null / CDKConstants.UNSET for unset Atoms -> check?!
    							if(oNSCounter > 1 && adjList[idx].length + cAtom.getImplicitHydrogenCount() == 4) {
    								// set as marked and break out of connected atoms
    								log.debug(String.format("Marking Atom #%d (%s) - Met condition 2.3", idx, cAtom.getSymbol())); 
    								isMarked = true;
    								break;
    							}
    						}
    					}
    					// if part of oxirane, aziridine and thiirane ring... [CONDITION 2.4]
    					for(int connectedInSphere2Idx : adjList[connectedIdx]) {
    						IAtom connectedInSphere2Atom = molecule.getAtom(connectedInSphere2Idx);
    						if(connectedInSphere2Atom.getAtomicNumber() == 6) {
    							for(int connectedInSphere3Idx : adjList[connectedInSphere2Idx]) {
    								IAtom connectedInSphere3Atom = molecule.getAtom(connectedInSphere3Idx);
    								if(connectedInSphere3Atom.equals(cAtom)) {
    									// set connected atoms as marked
    									log.debug(String.format("Marking Atom #%d (%s) - Met condition 2.4", connectedInSphere2Idx, connectedInSphere2Atom.getSymbol())); 
    									log.debug(String.format("Marking Atom #%d (%s) - Met condition 2.4", connectedInSphere3Idx, connectedInSphere3Atom.getSymbol()));
    									markedAtoms.add(connectedInSphere2Idx);
    									markedAtoms.add(connectedInSphere3Idx);
    									// set current atom as marked and break out of connected atoms
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
     * 			* 	find a better solution for parent-C-check?
     *			*	Environmental C's & type could better be stored in HashMap instead of atom's property!
     *			*	AtomContainerManipulator.extractSubstructure(IAtomContainer atomContainer, int... atomIndices)
     *				vs. cloning mol., deleting bonds/atoms & using ConnectivityChecker.partitionIntoMolecules(disconnectedContainer) ?
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
    		
    		// stores if the "parent" atom (atom we are coming from) is a carbon or not
    		Queue<Boolean> isParentCQueue = new ArrayDeque<>();
    		isParentCQueue.add(molecule.getAtom(beginIdx).getAtomicNumber() == 6);
    		
    		while(!queue.isEmpty()) {
    			Integer idx = queue.poll();
    			IAtom currentAtom = molecule.getAtom(idx);
    			
    			boolean isParentC = false;
    			if(mode != Mode.NO_GENERALIZATION) {
    				isParentC = isParentCQueue.poll();
    			}
    			
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
    						if(mode != Mode.NO_GENERALIZATION) {
    							//isParentCQueue.add(currentAtom.getAtomicNumber() == 6);
    							isParentCQueue.add(currentAtom.getAtomicNumber() == 6 && currentAtom.getProperty(CARBONYL_C_FLAG) != null);
    						}
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
    			// if we find an aromatic heteroatom
    			//TODO isHeteroatom should not be necessary?!
    			else if(isHeteroatom(currentAtom) && currentAtom.isAromatic()) {
    				// add its index to the functional group
    				fGroupIndices.add(idx);
    				// note that this aromatic heteroatom has been added to a group
    				aromaticHeteroAtoms.put(idx, true);
    			}
    			else {
    				log.debug(String.format("	visited non-marked atom: #%d (%s)", idx, currentAtom.getSymbol()));
    			}
    		}
    		log.debug("	search completed.");
    		
    		// extract functional group from the collected indices
    		int[] fGroupIndicesArray = Ints.toArray(fGroupIndices);
    		IAtomContainer fGroup = extractGroupByIndices(molecule, fGroupIndicesArray);
    		
    		log.debug("Extracting functional group by atom indices: ", Arrays.toString(fGroupIndicesArray));
    		
    		// create IDs for the group
    		IDCreator.createIDs(fGroup);
    		
    		// add functional group to list
    		functionalGroups.add(fGroup);
    	}
    	// extract aromatic heteroatoms that are not connected to an FG as single atoms
    	for(int idx : aromaticHeteroAtoms.keySet()) {
    		if(!aromaticHeteroAtoms.get(idx)) {
    			log.debug("Extracting single aromatic heteroatom by atom index: ", idx);
    			IAtomContainer fGroup = extractGroupByIndices(molecule, new int[] {idx});
        		IDCreator.createIDs(fGroup);
        		functionalGroups.add(fGroup);
    		}
    	}
    	
    	
    	log.debug(String.format("########## Found & extracted %d functional groups. ##########", functionalGroups.size()));
    	
    	return functionalGroups;
    }
    
    /**
     * TODO / NOTES
     * 
     * TODO:	*	change to private void!
     * 			*	FIXME: overhaul detection of env.C's on carbonyl-C's! works but its a quick & dirty fix!
     *			*	ATM: environmental C's are identified by an attached String property-flag
     *				option 1: store all of them (across all FGs) in a list of Hashmaps (idx - type)
     *				option 2: introduce inner class functionalGroup that stores an AC plus the HashSet.
     *			*	is it faster to do pattern matching for carbonyls? see class Pattern & VentoFoggia
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
    		boolean isSingleO = false, isSecAmineOrSimpleThiol = false, isSingleN = false;
    		
    		if(atomCount == 2) {
    			IBond bond = fGroup.bonds().iterator().next();
    			int atomicNumberSum = bond.getBegin().getAtomicNumber() + bond.getEnd().getAtomicNumber();
    			switch(atomicNumberSum) {
    			case 13: 	isSingleN = true; // C-N
    			break;
    			case 14:	isSingleO = true; // C-O
    			break;
    			case 22:	if(bond.getOrder() == Order.SINGLE) {
    							isSecAmineOrSimpleThiol = true; // C-S
    						}
    			break;
    			}
    		}
    		else if(atomCount == 3) {
    			// checking for two carbons and one nitrogen should be enough 
    			// (bond order has to be single or it should be more than three atoms anyway)
    			byte cCount = 0, nCount = 0;
    			for(IAtom atom : fGroup.atoms()) {
    				if(atom.getAtomicNumber() == 6) {
    					cCount++;
    					continue;
    				}
    				else if(atom.getAtomicNumber() == 7) {
    					nCount++;
    					continue;
    				}
    				break;
    			}
    			isSecAmineOrSimpleThiol = cCount == 2 && nCount == 1; 
    		}
    		log.debug(String.format("	precheck result: 	isSingleO: %b, isSingleN: %b, isSecAmineOrSimpleThiol: %b",
    				isSingleO, isSingleN, isSecAmineOrSimpleThiol));

    		
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
    		List<IAtom> atomsToRemove = new LinkedList<>();
    		for(IAtom atom : fGroup.atoms()) {
    			if(atom.getAtomicNumber() == 6) {
    				// STEP 1: delete environments on carbons... see exceptions!
    				EnvironmentalCType type = atom.getProperty(ENVIRONMENTAL_C_FLAG);
    				if(type == EnvironmentalCType.C_ON_C) {
    					//FIXME quick & dirty fix
//    					IAtom connectedC = atom.bonds().iterator().next().getOther(atom);
    					IAtom connectedC = fGroup.getConnectedBondsList(atom).get(0).getOther(atom);
    					boolean isCarbonylC = false;
//    					for(IBond bond : fGroup.getConnectedBondsList(connectedC)) {
//    						if(bond.getOrder() == Order.DOUBLE && bond.getBegin().getAtomicNumber() + bond.getEnd().getAtomicNumber() == 14) {
//    							log.debug("	found env. C (type CARBONLY_C)...");
//    							isCarbonylC = true;
//    							type = EnvironmentalCType.CARBONYL_C;
//    							atom.setProperty(ENVIRONMENTAL_C_FLAG, type);
//    							break;
//    						}
//    					}
    					if(!isCarbonylC) {
    						log.debug("	removing env. C (type C_ON_C)...");
    						atomsToRemove.add(atom);
    						
    					}
    					else {
    						IPseudoAtom rAtom = createRAtom();
        					AtomContainerManipulator.replaceAtomByAtom(fGroup, atom, rAtom);
        					log.debug("	replacing env. C (type C_ON_C, carbonyl-C) by R-atom...");
    					}
    					continue;
    				}
    				// STEP 3: replace other environmental C's by R-Atoms (represent H or C)... see exceptions!
    				// TODO: double check exception for carbonyls!
    				else if(type == EnvironmentalCType.C_ON_HETEROATOM && !isSingleN && !isSingleO) {
    					IPseudoAtom rAtom = createRAtom();
    					AtomContainerManipulator.replaceAtomByAtom(fGroup, atom, rAtom);
    					log.debug("	replacing env. C (type C_ON_HETEROATOM) by R-atom...");
    				}
    			}
    			// STEP 2: fill free valences on heteroatoms by R-Atoms... see exceptions!
    			else if(isHeteroatom(atom)){
    				if(atom.getImplicitHydrogenCount() != null) {
    					int rAtomsToAdd = atom.getImplicitHydrogenCount();
    					if(atom.getAtomicNumber() == 8 || isSecAmineOrSimpleThiol) {
    						rAtomsToAdd = 0;
    						restoreExplicitHydrogens(fGroup, atom);
    					}
    					else {
    						atom.setImplicitHydrogenCount(0);
    					}
    					for(int i = 0; i < rAtomsToAdd; i++) {
    						IPseudoAtom rAtom = createRAtom();
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
    		for(IAtom atomToRemove : atomsToRemove) {
    			fGroup.removeAtom(atomToRemove);
    		}
    		checkAndResolveRCycles(fGroup);
    	}
    	
    	log.debug("########## Generalization of functional groups completed. ##########");
    
    	return fGroups;
    }
    
    /**
     * Extract a functional group from an atom container, in the form of a new
     * cloned atom container with only the atoms with indices in atomIndices and
     * bonds that connect these atoms (except for the bonds between two 
     * environmental carbons).
     * 
     * @param molecule the source container to extract from
     * @param atomIndices the atom indices of the group
     * @return a cloned atom container with the group of the source
     * @throws CloneNotSupportedException if the source container cannot be cloned
     */
    public IAtomContainer extractGroupByIndices(IAtomContainer molecule, int... atomIndices) throws CloneNotSupportedException {    	
    	 IAtomContainer extract = AtomContainerManipulator.extractSubstructure(molecule, atomIndices);
         // remove all bonds between atoms marked as environmental C's
    	 for(Iterator<IBond> bondIter = extract.bonds().iterator(); bondIter.hasNext();) {
    		 IBond bond = bondIter.next();
    		 if(bond.getBegin().getProperty(ENVIRONMENTAL_C_FLAG) != null
        			 && bond.getEnd().getProperty(ENVIRONMENTAL_C_FLAG) != null) {
        		 extract.removeBond(bond);
        	 }
    	 }
         return extract;
    }
    
    public IAtomContainer extractGroupByIndicesV2(IAtomContainer molecule, int... atomIndices) throws CloneNotSupportedException {
    	// copy-paste from AtomContainerManipulator
    	//int
    	IAtomContainer group = (IAtomContainer) molecule.clone();
        int numberOfAtoms = group.getAtomCount();
        IAtom[] atoms = new IAtom[numberOfAtoms];
        for (int atomIndex = 0; atomIndex < numberOfAtoms; atomIndex++) {
            atoms[atomIndex] = group.getAtom(atomIndex);
        }
        
        Arrays.sort(atomIndices);
        for (int index = 0; index < numberOfAtoms; index++) {
            if (Arrays.binarySearch(atomIndices, index) < 0) {
                IAtom atom = atoms[index];
                group.removeAtom(atom);
            }
        }
        
        // map alt->neu
        
        
        
        return group;
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
    	if(atom.getImplicitHydrogenCount() == null) {
    		System.out.println("WARNING: atom with unset implicit hydrogen count");
    		return;
    	}
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
    
    /**
     * Creates a new R-Atom (representing C or H)
     * @return R-Atom
     */
    private static IPseudoAtom createRAtom() {
    	IPseudoAtom rAtom = new PseudoAtom("R");
		rAtom.setAttachPointNum(1);
		rAtom.setImplicitHydrogenCount(0);
		return rAtom;
    }
    
    /**
     * Checks generalised functional groups for ring-over-R-Atoms cases and resolves them (e.g. X1-R-X2 -> X1-R R-X2)
     * 
     * TODO: Optimise ?!
     * 
     * @param fGroup functional group
     */
    private static void checkAndResolveRCycles(IAtomContainer fGroup) {
    	Set<IAtom> rAtomsInBondSet = new HashSet<>();
    	for(Iterator<IBond> bondIter = fGroup.bonds().iterator(); bondIter.hasNext(); ) {
    		IBond bond = bondIter.next();
    		for(int atomInBondIdx = 0; atomInBondIdx < bond.getAtomCount(); atomInBondIdx++) {
    			IAtom atomInBond = bond.getAtom(atomInBondIdx);
    			// if atom in bond is an R-Atom
    			if(atomInBond instanceof PseudoAtom) {
    				if(rAtomsInBondSet.contains(atomInBond)) {
    					// if R-Atom is already in another bond we need to interfere
    					IAtom otherAtomInBond = bond.getOther(atomInBond);
    					if(!(otherAtomInBond instanceof PseudoAtom)) {
    						if(!(isHeteroatom(atomInBond) || atomInBond.getProperty(CARBONYL_C_FLAG) != null))
    						{
    							// if it is connected to a regular atom, swap the R-Atom in the bond for a newly created additional R-Atom
    							IPseudoAtom additionalRAtom = createRAtom();
    							fGroup.addAtom(additionalRAtom);
    							bond.setAtom(additionalRAtom, atomInBondIdx);
    						}
    						else {
    							fGroup.removeBond(atomInBond, otherAtomInBond);
    						}
    						
    						
    					}
    					else {
    						// if the R-Atom is connected to another R-Atom, remove the bond
    						fGroup.removeBond(atomInBond, otherAtomInBond);
    					}
    					continue;
    				}
    				else {
    					rAtomsInBondSet.add(atomInBond);
    				}
    			}
    		}
    	}
    }
}
