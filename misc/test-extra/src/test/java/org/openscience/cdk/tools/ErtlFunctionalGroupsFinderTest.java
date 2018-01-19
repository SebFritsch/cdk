package org.openscience.cdk.tools;


import static org.hamcrest.Matchers.is;

import java.beans.ExceptionListener;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import javax.lang.model.util.Elements;

import org.apache.commons.math.stat.Frequency;
import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.CDKTestCase;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.PseudoAtom;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.io.IChemObjectReader.Mode;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.isomorphism.Ullmann;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.isomorphism.VentoFoggia;
import org.openscience.cdk.silent.Atom;
import org.openscience.cdk.silent.Bond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
//import org.openscience.cdk.tools.ErtlFunctionalGroupsFinder.RAtom;
import org.openscience.cdk.tools.diff.ElementDiff;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multimaps;
import com.google.common.collect.Multiset;
import com.google.common.collect.Multisets;

import net.sf.jniinchi.INCHI_RET;


/**
 * @cdk.module test-standard
 * @author Sebastian Fritsch
 */
public class ErtlFunctionalGroupsFinderTest extends CDKTestCase {

    public ErtlFunctionalGroupsFinderTest() {
        super();
    }
    
    @Test
    public void testFind1() throws Exception {
    	String moleculeSmiles = "Cc1cc(C)nc(NS(=O)(=O)c2ccc(N)cc2)n1";
    	String[] expectedFGs = new String[] {"R-:N=:R", "[R]N([R])S(=O)(=O)[R]", "C_ar-NH2", "R-:N=:R"};
    	testFind(moleculeSmiles, expectedFGs);
    }
    
    @Test
    public void testFind2() throws Exception{
    	String moleculeSmiles = "NC(=N)c1ccc(\\\\C=C\\\\c2ccc(cc2O)C(=N)N)cc1";
    	String[] expectedFGs = new String[] {"NR=C-NR2", "C=C", "C_ar-OH", "NR=C-NR2"};
    	testFind(moleculeSmiles, expectedFGs);
    }
    
	@Test
	public void testFind3() throws Exception {
		String moleculeSmiles = "CC(=O)Nc1nnc(s1)S(=O)(=O)N";
    	String[] expectedFGs = new String[] {"[R]N([R])C(=O)[R]", "R-:N-:N-:R", "R-:S-:R", "[R]S(=O)(=O)N([R])[R]"};
    	testFind(moleculeSmiles, expectedFGs);
	}

	@Test
	public void testFind4() throws Exception {
		String moleculeSmiles = "NS(=O)(=O)c1cc2c(NCNS2(=O)=O)cc1Cl";
    	String[] expectedFGs = new String[] {"[R]S(=O)(=O)N([R])[R]", "[R]S(=O)(=O)N([R])[C]N([R])[R]", "Cl-R"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind5() throws Exception {
		String moleculeSmiles = "CNC1=Nc2ccc(Cl)cc2C(=N(=O)C1)c3ccccc3";
    	String[] expectedFGs = new String[] {"[R]N([R])[C]=N[R]", "Cl-R", "[R]N(=O)=[C]"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	//FIXME:	Group 2: env.C in 4-ring is BOTH carbonyl-c and C_ON_C! -> gets removed because of the latter
	@Test
	public void testFind6() throws Exception {
		String moleculeSmiles = "Cc1onc(c2ccccc2)c1C(=O)N[C@H]3[C@H]4SC(C)(C)[C@@H](N4C3=O)C(=O)O";
    	String[] expectedFGs = new String[] {"R-:O-:N=:R", "O=C([R])N([R])[R]", "O=C([R])N([R])[C]S[R]", "O=C([R])O"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind7() throws Exception {
		String moleculeSmiles = "Clc1ccccc1C2=NCC(=O)Nc3ccc(cc23)N(=O)=O";
    	String[] expectedFGs = new String[] {"Cl-R", "[R]N=[C]", "[R]C(=O)N([R])[R]", "O=N([R])=O"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind8() throws Exception {
		String moleculeSmiles = "COc1cc(cc(C(=O)NCC2CCCN2CC=C)c1OC)S(=O)(=O)N";
    	String[] expectedFGs = new String[] {"[R]O[R]", "[R]N([R])C(=O)[R]", "N([R])([R])[R]", "[C]=[C]", "[R]O[R]", "[R]S(=O)(=O)N([R])[R]"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	//FIXME:	Group 2: looks like isAldehydeOrKetone prevent the replacement of the env.C's by R-Atoms!
	@Test
	public void testFind9() throws Exception {
		String moleculeSmiles = "Cc1ccc(Cl)c(Nc2ccccc2C(=O)O)c1Cl";
    	String[] expectedFGs = new String[] {"Cl-R", "[R]N(H)[R]", "O=C(O[R])[R]", "Cl-R"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	//TODO: C in C=N not marked in paper! -> 4 instead of 5 groups
	@Test
	public void testFind10() throws Exception {
		String moleculeSmiles = "Clc1ccc2Oc3ccccc3N=C(N4CCNCC4)c2c1";
    	String[] expectedFGs = new String[] {"[R]Cl", "[R]O[R]", "[R]N([R])[C]=N[R]", "[R]N([H])[R]"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind11() throws Exception {
		String moleculeSmiles = "FC(F)(F)CN1C(=O)CN=C(c2ccccc2)c3cc(Cl)ccc13";
    	String[] expectedFGs = new String[] {"[R]F", "[R]F", "[R]F", "O=C([R])N([R])[R]", "[R]N=[C]", "[R]Cl"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	 
	@Test
	public void testFind12() throws Exception {
		String moleculeSmiles = "OC[C@H]1O[C@H](C[C@@H]1O)n2cnc3[C@H](O)CNC=Nc23";
    	String[] expectedFGs = new String[] {"[C]O[H]", "[R]O[R]", "[C]OH", "R-:N(-R)-:R" , "R-:N=:R", "[C]OH", "[R]N=CN([R])[R]"}; // FIXME 7 according to paper
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind13() throws Exception {
		String moleculeSmiles = "CCN[C@H]1C[C@H](C)S(=O)(=O)c2sc(cc12)S(=O)(=O)N";
    	String[] expectedFGs = new String[] {"[R]N([R])H", "O=S(=O)([R])[R]", "R-:S-:R", "[R]S(=O)(=O)N([R])[R]"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind14() throws Exception {
		String moleculeSmiles = "C[C@@H](O)[C@@H]1[C@H]2[C@@H](C)C(=C(N2C1=O)C(=O)O)S[C@@H]3CN[C@@H](C3)C(=O)N(C)C";
    	String[] expectedFGs = new String[] {"[C]O[H]", "O=C([R])N([R])C(C(=O)(O[R]))=[C]S[R]", "[R]N(H)[R]", "[R]N([R])C([R])=O"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	//TODO carbon in C=O & carbons in C=C not marked in paper! -> 6 instead of 8 groups
	@Test
	public void testFind15() throws Exception {
		String moleculeSmiles = "C[C@@H]1CN(C[C@H](C)N1)c2c(F)c(N)c3C(=O)C(=CN(C4CC4)c3c2F)C(=O)O";
    	String[] expectedFGs = new String[] {"[R]N([R])[R]", "[R]N([H])[R]", "[R]F", "C_ar-NH2", "[R]OC(=O)\\C(=C/N([R])[R])\\C(=O)[R]", "[R]F"}; //FIXME: 8 according to paper
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	//FIXME should the R-atom yield a ring?
	@Test
	public void testFind16() throws Exception {
		String moleculeSmiles = "CC(=CCC1C(=O)N(N(C1=O)c2ccccc2)c3ccccc3)C";
    	String[] expectedFGs = new String[] {"C=C", "[R]C(=O)N([R])N([R])C(=O)[R]"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	//FIXME should the R-atom yield a ring?
	@Test
	public void testFind17() throws Exception {
		String moleculeSmiles = "Clc1ccc2N=C3NC(=O)CN3Cc2c1Cl";
    	String[] expectedFGs = new String[] {"Cl[R]", "[R]N=C(N([R])[R])N([R])C(=O)[R]", "Cl[R]"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind18() throws Exception {
		String moleculeSmiles = "CC(=O)N[C@@H]1[C@@H](NC(=N)N)C=C(O[C@H]1[C@H](O)[C@H](O)CO)C(=O)O";
    	String[] expectedFGs = new String[] {"[R]N([R])C(=O)[R]", "[R]N([R])C(=N[R])N([R])[R]", "O=C(O[R])C(=[C])O[R]" , "[C]OH", "[C]OH", "[C]OH"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind19() throws Exception {
		String moleculeSmiles = "C[C@H](O)[C@H](O)[C@H]1CNc2nc(N)nc(O)c2N1";
    	String[] expectedFGs = new String[] {"[C]OH", "[C]OH", "[R]N(H)[R]" , "R-:N=:R", "C_ar-NH2", "R-:N=:R", "C_ar-OH", "[R]N(H)[R]"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	//TODO carbons in C=C & carbon in C=O not marked in paper! -> 5 instead of 7 groups
	@Test
	public void testFind20() throws Exception {
		String moleculeSmiles = "N[C@@H]1CCCCN(C1)c2c(Cl)cc3C(=O)C(=CN(C4CC4)c3c2Cl)C(=O)O";
    	String[] expectedFGs = new String[] {"[C]N([H])[H]", "[R]N([R])[R]", "[R]Cl" , "[R]OC(=O)C(=[C]N([R])[R])C(=O)[R]", "[R]Cl"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFindExtraS() throws Exception {
		String moleculeSmiles = "SCCSCC=S";
    	String[] expectedFGs = new String[] {"HS[R]", "[R]S[R]", "[C]=S"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFindExtraRRings1() throws Exception {
		String moleculeSmiles = "C1N=CC=N1"; // avoid ring with X1-R-X2
    	String[] expectedFGs = new String[] {"[R]\\N=C\\C=N\\[R]"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFindExtraRRings2() throws Exception {
		String moleculeSmiles = "C1CN=CNO1"; // avoid ring with X1-R-R-X2
    	String[] expectedFGs = new String[] {"[R]ON([R])\\C=N\\[R]"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testAtomAromaticity() throws Exception {
		String moleculeSmiles = "Oc1ccccc1";
		
		SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer mol = smilesParser.parseSmiles(moleculeSmiles);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		Aromaticity.cdkLegacy().apply(mol);
		
		for(IAtom a : mol.atoms()) {
			if(a.isAromatic()) {
				System.out.println("aromatic " + a.getSymbol());
			}
		}
		System.out.println("---------------");
		mol = AtomContainerManipulator.extractSubstructure(mol, new int[]{0,1,2,4});
		for(IAtom a : mol.atoms()) {
			if(a.isAromatic()) {
				System.out.println("aromatic " + a.getSymbol());
			}
		}
	}
	
	@Test
	public void testFindAromaticityCases() throws Exception {
		// my expected results for...
		
    	// testFind10 = Case #1
		// OUUTCOME: all models give the same result
    	String moleculeSmiles10 = "Clc1ccc2Oc3ccccc3N=C(N4CCNCC4)c2c1";
    	String[] expectedFGs10 = new String[] {"[R]Cl", "[R]O[R]", "[R]N([R])[C]=N[R]", "[R]N([H])[R]"};
    	//testFind15 = Case #2
    	// OUTCOME: different result with dalight() TODO doc
    	String moleculeSmiles15 = "C[C@@H]1CN(C[C@H](C)N1)c2c(F)c(N)c3C(=O)C(=CN(C4CC4)c3c2F)C(=O)O";
    	String[] expectedFGs15 = new String[] {"[R]N([R])[R]", "[R]N([H])[R]", "[R]F", "C_ar-NH2", "[R]OC(=O)\\C(=C/N([R])[R])\\C(=O)[R]", "[R]F"};
    	// testFind20 = Case #4
		// OUTCOME: different result with daylight()! TODO doc
		String moleculeSmiles20 = "N[C@@H]1CCCCN(C1)c2c(Cl)cc3C(=O)C(=CN(C4CC4)c3c2Cl)C(=O)O";
    	String[] expectedFGs20 = new String[] {"[C]N([H])[H]", "[R]N([R])[R]", "[R]Cl" , "[R]OC(=O)C(=[C]N([R])[R])C(=O)[R]", "[R]Cl"};
    	
    	testFind(moleculeSmiles15, expectedFGs15, new Aromaticity(ElectronDonation.daylight(), Cycles.all()));
	}
	
	private void testFind(String moleculeSmiles, String[] fGStrings) throws Exception {
		testFind(moleculeSmiles, fGStrings, Aromaticity.cdkLegacy());
	}
	
	private void testFind(String moleculeSmiles, String[] fGStrings, Aromaticity aromaticity) throws Exception {
		SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer mol = smilesParser.parseSmiles(moleculeSmiles);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		aromaticity.apply(mol);
		
		this.printToConsoleWithIndices(mol);
		
		ErtlFunctionalGroupsFinder fgFinder = new ErtlFunctionalGroupsFinder();
        fgFinder.find(mol);
        fgFinder.markAtoms(mol);
        List<IAtomContainer> fGs = fgFinder.extractGroups(mol);
        
        Assert.assertThat("Number of (non-generalized) functional groups did not match", fGs.size(), is(fGStrings.length));
        
        System.out.println(fGs.size() + " FUNCTIONAL GROUPS:");
        this.printMolFiles(fGs);
        
        for(IAtomContainer fG : fGs) {
			System.out.println(String.format("FG: %d atoms and %d bonds", fG.getAtomCount(), fG.getBondCount()));
			this.assertNoNonseseBonds(fG);
		}
        
        fGs = fgFinder.generalizeEnvironments(fGs);
        
		for(IAtomContainer fG : fGs) {
			System.out.println(String.format("FG: %d atoms and %d bonds", fG.getAtomCount(), fG.getBondCount()));
			this.assertNoNonseseBonds(fG);
		}
        
        System.out.println("GENERALIZED FUNCTIONAL GROUPS:");
		this.printMolFiles(fGs);
        
		Assert.assertThat("Number of (generalized) functional groups did not match", fGs.size(), is(fGStrings.length));
        
        // EXPECTED functional groups:
        List<IAtomContainer>	expectedFGs = new LinkedList<>();
        for(String fGString : fGStrings) {
        	expectedFGs.add(buildFunctionalGroup(fGString));
        }
        
        this.assertIsomorphism(expectedFGs, fGs);
	}
	
	// TODO: tidy up!
	@Test
	public void testPerformance1() throws Exception {
		List<String> moleculeSmilesList = new LinkedList<>();
		List<String[]> expectedFGsList = new LinkedList<>();
		
		// example 1
		moleculeSmilesList.add("Cc1cc(C)nc(NS(=O)(=O)c2ccc(N)cc2)n1");
		expectedFGsList.add(new String[] {"R-:N=:R", "[R]N([R])S(=O)(=O)[R]", "C_ar-NH2", "R-:N=:R"});
		// example 2
		moleculeSmilesList.add("NC(=N)c1ccc(\\\\\\\\C=C\\\\\\\\c2ccc(cc2O)C(=N)N)cc1");
		expectedFGsList.add(new String[] {"NR=C-NR2", "C=C", "C_ar-OH", "NR=C-NR2"});
		// example 3
		moleculeSmilesList.add("CC(=O)Nc1nnc(s1)S(=O)(=O)N");
		expectedFGsList.add(new String[] {"[R]N([R])C(=O)[R]", "R-:N-:N-:R", "R-:S-:R", "[R]S(=O)(=O)N([R])[R]"});
		// example 4
		moleculeSmilesList.add("NS(=O)(=O)c1cc2c(NCNS2(=O)=O)cc1Cl");
		expectedFGsList.add(new String[] {"[R]S(=O)(=O)N([R])[R]", "[R]S(=O)(=O)N([R])[C]N([R])[R]", "Cl-R"});
		// example 5
		moleculeSmilesList.add("CNC1=Nc2ccc(Cl)cc2C(=N(=O)C1)c3ccccc3");
		expectedFGsList.add(new String[] {"[R]N([R])[C]=N[R]", "Cl-R", "[R]N(=O)=[C]"});
		
		
		
		// parse smiles
		List<IAtomContainer> molecules = new LinkedList<>();
		List<List<IAtomContainer>> expectedFGs = new LinkedList<>();
		
		SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		for(String smiles : moleculeSmilesList) {
			IAtomContainer mol = smilesParser.parseSmiles(smiles);
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
			Aromaticity.cdkLegacy().apply(mol);
			molecules.add(mol);
		}
		
		for(String[] a : expectedFGsList) {
			List<IAtomContainer> expFGs = new LinkedList<>();
			for(String s : a) {
				expFGs.add(buildFunctionalGroup(s));
			}
			expectedFGs.add(expFGs);
		}
		
		// ensure that all examples run
		ErtlFunctionalGroupsFinder fgFinder = new ErtlFunctionalGroupsFinder();
		for(int i = 0; i < molecules.size(); i++) {
			IAtomContainer mol = molecules.get(i);
			fgFinder.find(mol);
			fgFinder.markAtoms(mol);
			List<IAtomContainer> fGs = fgFinder.extractGroups(mol);
			fGs = fgFinder.generalizeEnvironments(fGs);
			
			List<IAtomContainer> expFGs = expectedFGs.get(i);
			Assert.assertThat("Number of functional groups did not match", fGs.size(), is(expFGs.size()));
			assertIsomorphism(expFGs, fGs);
		}
		
		// perf. test
		int totalIterations = 20;
		int singleIterations = 2000;
		
		long startTime, endTime;
		for(int i = 0; i < totalIterations; i++) {
			startTime = System.currentTimeMillis();
			for(int j = 0; j < singleIterations; j++) {
				// run all examples
				for(IAtomContainer mol : molecules) {
					fgFinder.find(mol);
					fgFinder.markAtoms(mol);
					List<IAtomContainer> fGs = fgFinder.extractGroups(mol);
					fGs = fgFinder.generalizeEnvironments(fGs);
				}
			}
			endTime = System.currentTimeMillis();
			System.out.println(molecules.size()*singleIterations + " Iterations (" + molecules.size() + " cases): " + (endTime-startTime) +" ms");
		}	
	}
	
	@Test
	public void testChEBI_lite_3star() throws Exception {
		List<IAtomContainer> molecules = new LinkedList<>();
		//String filename = "data/mdl/ChEBI/ChEBI_lite_3star.sdf";
		String filename = "C:/Users/ETPCO/Lokal Documents/ChEBI/ChEBI_lite_3star_minicrop.sdf";
		//String filename = "data/mdl/BremserPredictionTest.mol";
        //InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
		InputStream in = new FileInputStream(filename);
        MDLV2000Reader reader = new MDLV2000Reader(in, Mode.STRICT);
        
        int readerStopsCount = 0;
        int basicFilteredCount = 0;
        int oranometallicCount = 0; // counts structures with 1 or more metals/metaloids
        int chargeCount = 0; // counts structures with 1 or more charged atoms
        int unconnectedCount = 0; // counts the atom containers that are not connected
        int invalidCount = 0;
//        for(IAtomContainer mol = reader.read(new AtomContainer()); mol != null;){
//        	n++;
//        	System.out.println(n);
//        }
        IAtomContainer mol = new AtomContainer();
        HashSet<String> nameSet = new HashSet<>();
        HashMap<String, Integer> exceptionCountMap = new HashMap<>();
        HashMap<String, Integer> metalsCountMap = new HashMap<>();
        while(true){
        	boolean isChargeCounted = false;
        	boolean isMetalsCounted = false;
        	boolean isUnconnected = false;
        	boolean isInvalid = false;
        	readerStopsCount++;
        	try {
        		mol = reader.read(new AtomContainer());
        	}
        	catch(Exception e) {
        		exceptionCountMap.put(e.getMessage(), exceptionCountMap.get(e.getMessage()) == null ? 1 : exceptionCountMap.get(e.getMessage()) + 1);
        		System.out.println("Reading error in enty " + readerStopsCount);
        		continue;
        	}
        	if(mol == null) {
        		break;
        	}
        	if(mol.getAtomCount() == 0 || mol.getBondCount() == 0) {
        		Exception e= new Exception("Atom and/or Bond count is zero!");
        		exceptionCountMap.put(e.getMessage(), exceptionCountMap.get(e.getMessage()) == null ? 1 : exceptionCountMap.get(e.getMessage()) + 1);
        		System.out.println("Skipping enty " + readerStopsCount);
        		continue;
        	}
        	basicFilteredCount++;
        	String name = mol.getProperty("ChEBI Name");
        	//nameSet.add(name);
        	//use mol
        	System.out.println(readerStopsCount + " ("+name+")");
        	
        	
        	
        	AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
        	HashSet<Integer> nonmetallicAtomicNumbers = new HashSet<>(Arrays.asList( //also non metalloid
        			1, 2, 6, 7, 8, 9, 10, 15, 16, 17, 18, 34, 35, 36, 53, 54, 86));
        	
        	for(IAtom atom : mol.atoms()) {
        		if(!nonmetallicAtomicNumbers.contains(atom.getAtomicNumber())) {
        			if(!isMetalsCounted) {
        				oranometallicCount++;
        				metalsCountMap.put(atom.getSymbol(), metalsCountMap.get(atom.getSymbol()) == null ? 1 : metalsCountMap.get(atom.getSymbol()) + 1);
        				isMetalsCounted = true;
        			}
        		}
        		if(atom.getFormalCharge() != 0) {
        			if(!isChargeCounted) {
        				chargeCount++;
        				isChargeCounted = true;
        			}
        		}
        	}
        	if(!ConnectivityChecker.isConnected(mol)) {
        		unconnectedCount++;
        		isUnconnected = true;
        	}
        	
        	if(isChargeCounted || isUnconnected) {
        		try {
        			mol = fragmentAndNeutraliseCharges(mol);
        		}
        		catch (CDKException e) {
        			invalidCount++;
        			isInvalid = true;
        		}
        	}
        	
        	if(!isMetalsCounted && !isInvalid) {
        		Aromaticity.cdkLegacy().apply(mol);
        		molecules.add(mol);
        	}
        }
        
        

		//IteratingSDFReader reader = new IteratingSDFReader(in, DefaultChemObjectBuilder.getInstance(), true);
		
//		List<IAtomContainer> atomContainerList = new LinkedList<>();
//		while(reader.hasNext()) {
//			IAtomContainer c = reader.next();
//			if(c == null)
//				System.out.println("Could not read entry.");
//			else
//				atomContainerList.add(reader.next());
//		}
		
        //IChemFile chemFile = reader.read(new ChemFile());
        reader.close();
        System.out.println("Found " + molecules.size() + " valid structures.");
        System.out.println("Read " + basicFilteredCount + " structures."); //Name set size: " + nameSet.size());
        System.out.println("with " + oranometallicCount + " organometallic structures & " + chargeCount + " structures with at least 1 charged atom & " + unconnectedCount + " unconnected structures");
        System.out.println("Number of invalid structures: " + invalidCount);
//        System.out.println("Metals and Metaloids log:");
//        for(String symbol : metalsCountMap.keySet()) {
//        	System.out.println(metalsCountMap.get(symbol) + "x: " + symbol);
//        }
        //System.out.println("Exception log:");
        //for(String message : exceptionCountMap.keySet()) {
        //	System.out.println(exceptionCountMap.get(message) + "x: " + message);
        //}
        System.out.println("######## start processing ... ##########");
        List<IAtomContainer> allFGs = new LinkedList<>();
        ErtlFunctionalGroupsFinder fgFinder = new ErtlFunctionalGroupsFinder();
        for(IAtomContainer aC : molecules) {
        	
            fgFinder.find(aC);
            fgFinder.markAtoms(aC);
            List<IAtomContainer> fGs = fgFinder.extractGroups(aC);
            fGs = fgFinder.generalizeEnvironments(fGs);
            allFGs.addAll(fGs);
        }
        System.out.println("######## processing done ##########");
        
        printMostCommonFunctionalGroups(allFGs);
	}
    
    @Test
    public void uitTest() throws Exception {
    	UniversalIsomorphismTester uiTester = new UniversalIsomorphismTester();
    	
    	SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer mol = smilesParser.parseSmiles("[R]O[R]");
        
        
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		Aromaticity.cdkLegacy().apply(mol);
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
		for(IAtom a : mol.atoms()) {
        	if(a instanceof PseudoAtom) {
        		System.out.println(a.getSymbol());
        		a.setSymbol("R");
        	}
        }
		
		System.out.println(mol.getAtomCount());
		
		IAtomContainer c = new AtomContainer();
		IAtom a1 = new Atom("O");
		IPseudoAtom a2 = new PseudoAtom();
		System.out.println(a2.getSymbol());
		a2.setSymbol("R");
		IPseudoAtom a3 = new PseudoAtom();
		a3.setSymbol("R");
		IBond b1 = new Bond(new IAtom[] {a1,  a2}, Order.SINGLE);
		IBond b2 = new Bond(new IAtom[] {a1,  a3}, Order.SINGLE);
		c.addAtom(a1);
		c.addAtom(a2);
		c.addAtom(a3);
		c.addBond(b1);
		c.addBond(b2);
		
		//PATTERN
		Pattern ventoFoggia = VentoFoggia.findIdentical(mol);
		Pattern ullmann = Ullmann.findIdentical(mol);
		
		long startTime, endTime;
		
		for(int k = 0; k < 10; k++) {
			startTime = System.currentTimeMillis();
			for(int i= 0; i < 100000; i++) {
				ullmann.matches(c);
			};
			endTime = System.currentTimeMillis();
			System.out.println("Ullmann: " + (endTime - startTime) + " ms");	
			
		startTime = System.currentTimeMillis();
		for(int i= 0; i < 100000; i++) {
			ventoFoggia.matches(c);
		};
		endTime = System.currentTimeMillis();
		System.out.println("VentoFoggia: " + (endTime - startTime) + " ms");
		
		startTime = System.currentTimeMillis();
		for(int i= 0; i < 100000; i++) {
			uiTester.isIsomorph(mol, c);
		};
		endTime = System.currentTimeMillis();
		System.out.println("UITester: " + (endTime - startTime) + " ms");
		}
    }
    
    @Test
    public void testFragmentAndNeutraliseCharges() throws Exception{
    	String smiles = "CC[O-].C";
    	
    	SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer mol = smilesParser.parseSmiles(smiles);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		Aromaticity.cdkLegacy().apply(mol);
		
		mol = fragmentAndNeutraliseCharges(mol);
		
		SmilesGenerator sg = SmilesGenerator.unique();
		Assert.assertEquals("OCC", sg.create(mol));
    }
    
    private IAtomContainer fragmentAndNeutraliseCharges(IAtomContainer mol) throws CDKException {
    	if(!ConnectivityChecker.isConnected(mol)) {
    		IAtomContainerSet fragmentSet = ConnectivityChecker.partitionIntoMolecules(mol);
    		IAtomContainer biggestFragment = null;
    		for(IAtomContainer fragment : fragmentSet.atomContainers()) {
    			 if(biggestFragment == null || biggestFragment.getAtomCount() < fragment.getAtomCount()) {
    				 biggestFragment = fragment;
    			 }
    		}
    		mol = biggestFragment;
    	}
    	
    	CDKHydrogenAdder hydrogenAdder = CDKHydrogenAdder.getInstance(mol.getBuilder());
    	CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(mol.getBuilder());
    	for(IAtom atom : mol.atoms()) {
    		int formalCharge = atom.getFormalCharge();
    		if(formalCharge != 0) {
    			atom.setFormalCharge(0);
    			IAtomType matched = matcher.findMatchingAtomType(mol, atom);
    			if (matched != null) AtomTypeManipulator.configure(atom, matched);
				hydrogenAdder.addImplicitHydrogens(mol, atom);
    		}
    	}
    	
    	return mol;
    }
    
    private void printMostCommonFunctionalGroups(List<IAtomContainer> molecules) throws CDKException {
    	Multiset<String> occurences = HashMultiset.create(); // counts occurences of inchiKeys
    	HashMap<String,String> inchiKeyToSmilesMap = new HashMap<>(); // maps inchiKeys to unique SMILES
    	
    	SmilesGenerator smilesGen = SmilesGenerator.unique();
    	InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
    	
    	int skippedCount = 0;
    	
    	molLoop:
    	for(IAtomContainer mol : molecules) {
//    		//replace R-Atoms by U
    		for(IAtom atom : mol.atoms()) {
    			if(atom instanceof IPseudoAtom) {
    				atom.setImplicitHydrogenCount(0);
//    				IAtom replacementAtom = new Atom("U");
//    				replacementAtom.setImplicitHydrogenCount(0);
//    				AtomContainerManipulator.replaceAtomByAtom(mol, atom, replacementAtom);
    			}
    			//replace aromatic C's by Au
    			else if (atom.getAtomicNumber() == 6 && atom.isAromatic()) {
    				IAtom replacementAtom = new Atom("Au");
    				replacementAtom.setImplicitHydrogenCount(0);
    				AtomContainerManipulator.replaceAtomByAtom(mol, atom, replacementAtom);
    				//System.out.println("###### replaced aromatic" + atom);
    			}
    		}
//    		for(IAtom atom : mol.atoms()) {
//    			if(atom.getImplicitHydrogenCount() == null) {
//    				skippedCount ++;
//    				System.out.println("UNSET IMPL. HYDROGEN COUNT ON " + atom.getSymbol());
//    				continue molLoop; //FIXME!!!
//    			}
//    		}
    		// strip aromaticity info from bonds
    		for(IBond bond : mol.bonds()) {
    			bond.setIsAromatic(false);
    		}
    		// strip aromaticity info from atoms & remove implicit hydrogens
    		for(IAtom atom : mol.atoms()) {
    			atom.setImplicitHydrogenCount(0);
    			atom.setIsAromatic(false);
    		}
    		// create unique smiles
    		String smiles = smilesGen.create(mol);
    		// create inchiKey
    		InChIGenerator gen = factory.getInChIGenerator(mol, "FixedH");
    		String inchiKey = gen.getInchiKey();
    		INCHI_RET returnType = gen.getReturnStatus();
    		if(returnType == INCHI_RET.WARNING) {
    			// warning but still works
    			System.out.println("### Inchi WARNING: " + gen.getMessage());
    		}
    		else if(returnType != INCHI_RET.OKAY) {
    			// skip fG on error
    			System.out.println("### Inchi ERROR !!! (->Skipping): " + gen.getMessage());
    			skippedCount++;
    			continue molLoop;
    		}
    		// count 
    		occurences.add(inchiKey);
    		// map inchiKey to SMILES if not done before
    		if(inchiKeyToSmilesMap.containsKey(inchiKey)) {
    			// if it already contains the inchiKey, double check if the SMILES for the inchiKey are equal
    			String smilesInMap = inchiKeyToSmilesMap.get(inchiKey);
    			if(!smiles.equals(smilesInMap)) {
    				System.out.println("### Inchi/Smiles WARNING: InchiKey " + inchiKey + " mapped to SMILES " + smilesInMap + ". Other proposes SMILES " + smiles);
    			}
    		}
    		else {
    			inchiKeyToSmilesMap.put(inchiKey, smiles);
    		}
    		
    		// print to console
    		int frequencyThreshold = 10; // threshold for printout
    		int totalFGCount = occurences.size();
    		System.out.println("-----------------------------------------------------------------");
    		System.out.println("Counted " + occurences.size() + " unique FGs.");
        	System.out.println("Skipped " + skippedCount);
    		String format = "%-20s%-20s%s%s";
    		for(String key : Multisets.copyHighestCountFirst(occurences)) {
    			String smi = inchiKeyToSmilesMap.get(key);
    			int frequency = occurences.count(key);
    			double percentage = frequency/totalFGCount * 100;
    			
    			if(frequency >= frequencyThreshold)
    				System.out.printf(format, frequency + "x", percentage + "%", smi);
    		}
    	}
    }
    
    private void printToConsoleWithIndices(IAtomContainer mol) {
    	IAtomContainer c = new AtomContainer();
    	try {
			c = mol.clone();
		} catch (CloneNotSupportedException e1) {
			System.out.println("#ERROR# Could not clone molecule!");
		}
    	
    	for(int i = 0; i < c.getAtomCount(); i++) {
    		IAtom a = c.getAtom(i);
    		a.setMassNumber(i);
    	}
		
		SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.AtomicMass);
		String smiles = new String();
		try {
			smiles = smiGen.create(c);
		} catch (CDKException e) {
			System.out.println("#ERROR# Could not generate SMILES!");
		}
		System.out.println("------------------------------------------------------");
		System.out.println(smiles);
		System.out.println("------------------------------------------------------");
    }
    
    private void printMolFiles(List<IAtomContainer> fGs) throws IOException {
		StringWriter stringWriter = new StringWriter();
		MDLV2000Writer writer = new MDLV2000Writer(stringWriter);
		for(IAtomContainer fG : fGs) {
				try {
					writer.writeMolecule(fG);
				} catch (Exception e) {
					System.out.println(e.getMessage());
					Assert.assertFalse(true);
				}
		}
		writer.close();
		String molFile = stringWriter.toString();
		System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		System.out.println(molFile);
		System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    }
    
    private void assertNoNonseseBonds(IAtomContainer mol) {
    	for(IBond b : mol.bonds()) {
    		IAtom begin = b.getBegin();
    		IAtom end = b.getEnd();
    		Assert.assertNotNull(begin);
    		Assert.assertNotNull(end);
    		Assert.assertTrue(mol.contains(begin));
    		Assert.assertTrue(mol.contains(end));
    	}
    }
    
    //FIXME: input needs same order!
    //TODO: use VentoFoggia.findIdentical instead of UIT? -> much faster
    /**
     * Uses asserts to check the following (actual vs. expected):
     * 	-	Number of functional groups
     * 	-	Atom counts for each group
     * 	-	Bond counts for each group
     * 	-	Isomorphism of each group pair (as interpreted by UniversalIsomorphismTester)
     * 	-	Same aromaticities on carbons
     * 	-	Same aromaticities for bonds
     * 
     * @param expectedFGs 	list of expected functional groups
     * @param fGs			list of actual functional groups
     * @throws Exception	if anything does not work as planned
     */
    private void assertIsomorphism(List<IAtomContainer> expectedFGs, List<IAtomContainer> fGs) throws Exception {
    	Assert.assertThat("Number of functional groups does not match the expected number of groups", fGs.size(), is(expectedFGs.size()));
    	UniversalIsomorphismTester uiTester = new UniversalIsomorphismTester();
    	for(int i = 0; i < expectedFGs.size(); i++) {
    		IAtomContainer cExp = expectedFGs.get(i);
    		IAtomContainer cAct = fGs.get(i);
    		
    		Assert.assertThat("Groups #" + i + ": different atom count", cAct.getAtomCount(), is(cExp.getAtomCount()));
    		Assert.assertThat("Groups #" + i + ": different bond count", cAct.getBondCount(), is(cExp.getBondCount()));

    		Assert.assertThat("Groups #" + i + ": not isomorph", uiTester.isIsomorph(cExp, cAct), is(true));

    		for(int j = 0; j < cExp.getAtomCount(); j++) {
    			if(cExp.getAtom(j).getAtomicNumber() == 6) {
    				Assert.assertThat("Groups #" + i + ": Carbon aromaticity does not match", cAct.getAtom(j).isAromatic(), is(cExp.getAtom(j).isAromatic()));
    			}
    		}
    		for(int j = 0; j < cExp.getBondCount(); j++) {
    			Assert.assertThat("Groups #" + i + ": Bond aromaticity does not match", cAct.getBond(j).isAromatic(), is(cExp.getBond(j).isAromatic()));
    		}
    	}
    }
    
    private IAtomContainer buildFunctionalGroup(String string) {
        IAtom 					a1, a2, a3, a4, a5, a6, a7, a8, a9;
        IBond 					b1, b2, b3, b4, b5, b6, b7, b8, b9;
        IChemObjectBuilder 		builder = DefaultChemObjectBuilder.getInstance();
        IAtomContainer 			container;
        
        switch(string) {
        case "R-:N=:R":
        	a1 = builder.newInstance(IPseudoAtom.class, "R");
        	a2 = builder.newInstance(IPseudoAtom.class, "R");
        	a3 = builder.newInstance(IAtom.class, "N");

        	b1 = builder.newInstance(IBond.class, a1, a3, Order.SINGLE);
        	b1.setIsAromatic(true);
        	b2 = builder.newInstance(IBond.class, a2, a3, Order.DOUBLE);
        	b2.setIsAromatic(true);

        	container = new AtomContainer();
        	container.setAtoms(new IAtom[] {a1, a2, a3});
        	container.setBonds(new IBond[] {b1, b2});
        	return container;
        	
        case "R-:N(-R)-:R":
        	a1 = builder.newInstance(IPseudoAtom.class, "R");
        	a2 = builder.newInstance(IPseudoAtom.class, "R");
        	a3 = builder.newInstance(IPseudoAtom.class, "R");
        	a4 = builder.newInstance(IAtom.class, "N");

        	b1 = builder.newInstance(IBond.class, a2, a4, Order.SINGLE);
        	b1.setIsAromatic(true);
        	b2 = builder.newInstance(IBond.class, a3, a4, Order.SINGLE);
        	b2.setIsAromatic(true);
        	b3 = builder.newInstance(IBond.class, a1, a4, Order.SINGLE);
        	//b3.setIsAromatic(true);
        	
        	container = new AtomContainer();
        	container.setAtoms(new IAtom[] {a1, a2, a3, a4});
        	container.setBonds(new IBond[] {b3, b2, b1});
        	return container;
        
        case "C_ar-NH2":
            a1 = builder.newInstance(IAtom.class, "C");
            a1.setIsAromatic(true);
            a2 = builder.newInstance(IAtom.class, "N");
            a3 = builder.newInstance(IAtom.class, "H");
            a4 = builder.newInstance(IAtom.class, "H");
            
            b1 = builder.newInstance(IBond.class, a1, a2, Order.SINGLE);
            b2 = builder.newInstance(IBond.class, a2, a3, Order.SINGLE);
            b3 = builder.newInstance(IBond.class, a2, a4, Order.SINGLE);
            
            container = new AtomContainer();
            container.setAtoms(new IAtom[] {a1, a2, a3, a4});
            container.setBonds(new IBond[] {b1, b2, b3});
            return container;
        
        // TODO remove
        case "NR=C-NR2":
            a1 = builder.newInstance(IPseudoAtom.class, "R");
            a2 = builder.newInstance(IPseudoAtom.class, "R");
            a3 = builder.newInstance(IPseudoAtom.class, "R");
            a4 = builder.newInstance(IAtom.class, "N");
            a5 = builder.newInstance(IAtom.class, "N");
            a6 = builder.newInstance(IAtom.class, "C");
            
            b1 = builder.newInstance(IBond.class, a1, a4, Order.SINGLE);
            b2 = builder.newInstance(IBond.class, a2, a4, Order.SINGLE);
            b3 = builder.newInstance(IBond.class, a3, a5, Order.SINGLE);
            b4 = builder.newInstance(IBond.class, a4, a6, Order.SINGLE);
            b5 = builder.newInstance(IBond.class, a5, a6, Order.DOUBLE);
            
            container = new AtomContainer();
            container.setAtoms(new IAtom[] {a1, a2, a3, a4, a5, a6});
            container.setBonds(new IBond[] {b1, b2, b3, b4, b5});
            return container;
            
        // TODO remove & replace by [C]=[C]    
        case "C=C":
            a1 = builder.newInstance(IAtom.class, "C");
            a2 = builder.newInstance(IAtom.class, "C");
            
            b1 = builder.newInstance(IBond.class, a1, a2, Order.DOUBLE);
            
            container = new AtomContainer();
            container.setAtoms(new IAtom[] {a1, a2});
            container.setBonds(new IBond[] {b1});
            return container;
        	
        case "C_ar-OH":
            a1 = builder.newInstance(IAtom.class, "C");
            a1.setIsAromatic(true);
            a2 = builder.newInstance(IAtom.class, "O");
            a3 = builder.newInstance(IAtom.class, "H");
            
            b1 = builder.newInstance(IBond.class, a1, a2, Order.SINGLE);
            b2 = builder.newInstance(IBond.class, a2, a3, Order.SINGLE);
            
            container = new AtomContainer();
            container.setAtoms(new IAtom[] {a1, a2, a3});
            container.setBonds(new IBond[] {b1, b2});
            return container;
            
        case "R-:N-:N-:R":
            a1 = builder.newInstance(IPseudoAtom.class, "R");
            a2 = builder.newInstance(IPseudoAtom.class, "R");
            a3 = builder.newInstance(IAtom.class, "N");
            a4 = builder.newInstance(IAtom.class, "N");
            
            b1 = builder.newInstance(IBond.class, a1, a3, Order.SINGLE);
            b1.setIsAromatic(true);
            b2 = builder.newInstance(IBond.class, a2, a4, Order.SINGLE);
            b2.setIsAromatic(true);
            b3 = builder.newInstance(IBond.class, a3, a4, Order.SINGLE);
            b3.setIsAromatic(true);
            
            container = new AtomContainer();
            container.setAtoms(new IAtom[] {a1, a2, a3, a4});
            container.setBonds(new IBond[] {b1, b2, b3});
            return container;
            
        case "R-:S-:R":
            a1 = builder.newInstance(IPseudoAtom.class, "R");
            a2 = builder.newInstance(IPseudoAtom.class, "R");
            a3 = builder.newInstance(IAtom.class, "S");
            
            b1 = builder.newInstance(IBond.class, a1, a3, Order.SINGLE);
            b1.setIsAromatic(true);
            b2 = builder.newInstance(IBond.class, a2, a3, Order.SINGLE);
            b2.setIsAromatic(true);
            
            container = new AtomContainer();
            container.setAtoms(new IAtom[] {a1, a2, a3});
            container.setBonds(new IBond[] {b1, b2});
            return container;
            
        case "R-:O-:N=:R":
        	a1 = builder.newInstance(IPseudoAtom.class, "R");
            a2 = builder.newInstance(IPseudoAtom.class, "R");
            a3 = builder.newInstance(IAtom.class, "O");
            a4 = builder.newInstance(IAtom.class, "N");
            
            b1 = builder.newInstance(IBond.class, a1, a3, Order.SINGLE);
            b1.setIsAromatic(true);
            b2 = builder.newInstance(IBond.class, a3, a4, Order.SINGLE);
            b2.setIsAromatic(true);
            b3 = builder.newInstance(IBond.class, a4, a2, Order.DOUBLE);
            b3.setIsAromatic(true);
            
            container = new AtomContainer();
            container.setAtoms(new IAtom[] {a1, a2, a3, a4});
            container.setBonds(new IBond[] {b1, b2, b3});
            return container;
        
        // TODO remove & replace by Cl[R]
        case "Cl-R":
            a1 = builder.newInstance(IPseudoAtom.class, "R");
            a2 = builder.newInstance(IAtom.class, "Cl");
            
            b1 = builder.newInstance(IBond.class, a1, a2, Order.SINGLE);
                    
            container = new AtomContainer();
            container.setAtoms(new IAtom[] {a1, a2});
            container.setBonds(new IBond[] {b1});
            return container;
            
        default:
        	try {
        		SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        		container = smilesParser.parseSmiles(string);
                for(IAtom a : container.atoms()) {
                	if(a instanceof PseudoAtom) {
                		a.setSymbol("R");
                	}
                }
                return container;
        	}
        	catch(InvalidSmilesException e) {
        		throw new IllegalArgumentException("Input string '" + string + " could not be found as a template and is not a valid SMILES string.");
        	}
        }
    }
}