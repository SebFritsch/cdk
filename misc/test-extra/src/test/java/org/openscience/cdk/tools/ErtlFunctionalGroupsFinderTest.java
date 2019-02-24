package org.openscience.cdk.tools;


import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.CDKTestCase;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.PseudoAtom;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.isomorphism.VentoFoggia;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import static org.hamcrest.Matchers.is;


/**
 * Test for ErtlFunctionalGroupsFinder.
 *
 * @author Sebastian Fritsch
 */
public class ErtlFunctionalGroupsFinderTest extends CDKTestCase {

    public ErtlFunctionalGroupsFinderTest() {
        super();
    }
    
    @Test
    public void testFind1() throws Exception {
    	String moleculeSmiles = "Cc1cc(C)nc(NS(=O)(=O)c2ccc(N)cc2)n1";
    	String[] expectedFGs = new String[] {"[R]N([R])S(=O)(=O)[R]", "[c]N(H)H", "NarR3", "NarR3"};
    	testFind(moleculeSmiles, expectedFGs);
    }
    
    @Test
    public void testFind2() throws Exception{
    	String moleculeSmiles = "NC(=N)c1ccc(\\\\C=C\\\\c2ccc(cc2O)C(=N)N)cc1";
    	String[] expectedFGs = new String[] {"[R]N=C-N([R])[R]", "[C]=[C]", "[c]OH", "[R]N=C-N([R])[R]"};
    	testFind(moleculeSmiles, expectedFGs);
    }
    
	@Test
	public void testFind3() throws Exception {
		String moleculeSmiles = "CC(=O)Nc1nnc(s1)S(=O)(=O)N";
    	String[] expectedFGs = new String[] {"[R]N([R])C(=O)[R]", "[R]S(=O)(=O)N([R])[R]", "NarR3", "NarR3", "SarR2"};
    	testFind(moleculeSmiles, expectedFGs);
	}

	@Test
	public void testFind4() throws Exception {
		String moleculeSmiles = "NS(=O)(=O)c1cc2c(NCNS2(=O)=O)cc1Cl";
    	String[] expectedFGs = new String[] {"[R]S(=O)(=O)N([R])[R]", "[R]S(=O)(=O)N([R])[C]N([R])[R]", "[R]Cl"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind5() throws Exception {
		String moleculeSmiles = "CNC1=Nc2ccc(Cl)cc2C(=N(=O)C1)c3ccccc3";
    	String[] expectedFGs = new String[] {"[R]N([R])[C]=N[R]", "[R]Cl", "[R]N(=O)=[C]"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind6() throws Exception {
		String moleculeSmiles = "Cc1onc(c2ccccc2)c1C(=O)N[C@H]3[C@H]4SC(C)(C)[C@@H](N4C3=O)C(=O)O";
    	String[] expectedFGs = new String[] {"O=C([R])N([R])[R]",  "O=C([R])N([R])[C]S[R]", "O=C([R])OH", "OarR2", "NarR3"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind7() throws Exception {
		String moleculeSmiles = "Clc1ccccc1C2=NCC(=O)Nc3ccc(cc23)N(=O)=O";
		String[] expectedFGs = new String[] {"[R]Cl", "[R]N=[C]", "[R]C(=O)N([R])[R]", "O=N([R])=O"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind8() throws Exception {
		String moleculeSmiles = "COc1cc(cc(C(=O)NCC2CCCN2CC=C)c1OC)S(=O)(=O)N";
    	String[] expectedFGs = new String[] {"[R]O[R]", "[R]N([R])C(=O)[R]", "N([R])([R])[R]", "[C]=[C]", "[R]O[R]", "[R]S(=O)(=O)N([R])[R]"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind9() throws Exception {
		String moleculeSmiles = "Cc1ccc(Cl)c(Nc2ccccc2C(=O)O)c1Cl";
    	String[] expectedFGs = new String[] {"[R]Cl", "[R]N(H)[R]", "O=C(OH)[R]", "[R]Cl"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
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
		String moleculeSmiles = "OC[C@H]1O[C@H](C[C@@H]1O)n2cnc3[C@H](O)CNC=Nc23";;
    	String[] expectedFGs = new String[] {"[C]O[H]", "[R]O[R]", "[C]OH", "[C]OH", "[R]N=CN([R])[R]", "NarR3", "NarR3"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind13() throws Exception {
		String moleculeSmiles = "CCN[C@H]1C[C@H](C)S(=O)(=O)c2sc(cc12)S(=O)(=O)N";
    	String[] expectedFGs = new String[] {"[R]N([R])H", "O=S(=O)([R])[R]", "[R]S(=O)(=O)N([R])[R]", "SarR2"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind14() throws Exception {
		String moleculeSmiles = "C[C@@H](O)[C@@H]1[C@H]2[C@@H](C)C(=C(N2C1=O)C(=O)O)S[C@@H]3CN[C@@H](C3)C(=O)N(C)C";
    	String[] expectedFGs = new String[] {"[C]O[H]", "O=C([R])N([R])C(C(=O)(OH))=[C]S[R]", "[R]N(H)[R]", "[R]N([R])C([R])=O"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind15() throws Exception {
		String moleculeSmiles = "C[C@@H]1CN(C[C@H](C)N1)c2c(F)c(N)c3C(=O)C(=CN(C4CC4)c3c2F)C(=O)O";
    	String[] expectedFGs = new String[] {"[R]N([R])[R]", "[R]N([H])[R]", "[R]F", "[c]N(H)H", "[c]=O", "[R]F", "[R]C(=O)OH", "NarR3"};
		testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind16() throws Exception {
		String moleculeSmiles = "CC(=CCC1C(=O)N(N(C1=O)c2ccccc2)c3ccccc3)C";
    	String[] expectedFGs = new String[] {"[C]=[C]", "[R]C(=O)N([R])N([R])C(=O)[R]"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind17() throws Exception {
		String moleculeSmiles = "Clc1ccc2N=C3NC(=O)CN3Cc2c1Cl";
    	String[] expectedFGs = new String[] {"Cl[R]", "[R]N=C(N([R])[R])N([R])C(=O)[R]", "Cl[R]"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind18() throws Exception {
		String moleculeSmiles = "CC(=O)N[C@@H]1[C@@H](NC(=N)N)C=C(O[C@H]1[C@H](O)[C@H](O)CO)C(=O)O";
		String[] expectedFGs = new String[] {"[R]N([R])C(=O)[R]", "[R]N([R])C(=N[R])N([R])[R]", "O=C(OH)C(=[C])O[R]" , "[C]OH", "[C]OH", "[C]OH"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind19() throws Exception {
		String moleculeSmiles = "C[C@H](O)[C@H](O)[C@H]1CNc2nc(N)nc(O)c2N1";
    	String[] expectedFGs = new String[] {"[C]OH", "[C]OH", "[R]N(H)[R]" , "[c]N(H)H",  "[c]OH", "[R]N(H)[R]", "NarR3", "NarR3"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFind20() throws Exception {
		String moleculeSmiles = "N[C@@H]1CCCCN(C1)c2c(Cl)cc3C(=O)C(=CN(C4CC4)c3c2Cl)C(=O)O";
    	String[] expectedFGs = new String[] {"[C]N([H])[H]", "[R]N([R])[R]", "[R]Cl" , "[c]=O", "[R]Cl", "[R]C(=O)OH", "NarR3"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFindExtraS() throws Exception { 
		String moleculeSmiles = "SCCSCC=S";
    	String[] expectedFGs = new String[] {"HS[R]", "[R]S[R]", "[C]=S"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFindExtraRNR() throws Exception {
		String moleculeSmiles = "CCCCN=CCCC";
    	String[] expectedFGs = new String[] {"[C]=N[R]"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFindExtraPyridine() throws Exception {
		String moleculeSmiles = "c1ccncc1";
    	String[] expectedFGs = new String[] {"NarR3"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFindExtraRRings1() throws Exception {
		String moleculeSmiles = "C1N=CC=N1"; // avoid ring with X1-R-X2
    	String[] expectedFGs = new String[] {"[R]N=CC=N[R]"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFindExtraRRings2() throws Exception {
		String moleculeSmiles = "C1CN=CNO1"; // avoid ring with X1-R-R-X2
    	String[] expectedFGs = new String[] {"[R]ON([R])C=N[R]"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFindExtraAromaticConnectedToAliphatic() throws Exception {
		String moleculeSmiles = "OC(=O)n1cccc1";
    	String[] expectedFGs = new String[] {"NarR2COOH"};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFindExtraEpoxide() throws Exception {
		String moleculeSmiles = "ClCC1CO1"; // Epichlorohydrin
    	String[] expectedFGs = new String[] {"[R]Cl", "[C]1[C]O1"};
    	testFind(moleculeSmiles, expectedFGs);
	}

	@Test
	public void testFindExtraExplicitHs() throws Exception {
		String moleculeSmiles = "[H]N([H])[H]";
		String[] expectedFGs = new String[] {"[R]N([R])[R]"};
		testFind(moleculeSmiles, expectedFGs);
	}
	
	@Test
	public void testFindExtraSingleH() throws Exception {
		String moleculeSmiles = "[H]";
    	String[] expectedFGs = new String[] {};
    	testFind(moleculeSmiles, expectedFGs);
	}
	
	private void testFind(String moleculeSmiles, String[] fGStrings) throws Exception {
		testFind(moleculeSmiles, fGStrings, new Aromaticity(ElectronDonation.daylight(), Cycles.all()));
	}
	
	private void testFind(String moleculeSmiles, String[] fGStrings, Aromaticity aromaticity) throws Exception {
    	// prepare input
		SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer mol = smilesParser.parseSmiles(moleculeSmiles);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		aromaticity.apply(mol);

		// find functional groups
		ErtlFunctionalGroupsFinder fgFinder = new ErtlFunctionalGroupsFinder();
		List<IAtomContainer> fGs = fgFinder.find(mol);

		// get expected groups
        List<IAtomContainer>	expectedFGs = new LinkedList<>();
        for(String fGString : fGStrings) {
        	expectedFGs.add(buildFunctionalGroup(fGString));
        }

        // compare
        this.assertIsomorphism(expectedFGs, fGs);
	}


    /**
     * NOTE: actual and expected functional groups must be in the same order!
	 *
     * @param expectedFGs 	list of expected functional groups
	 * @param actualFGs			list of actual functional groups
     * @throws Exception	if anything does not work as planned
     */
    private void assertIsomorphism(List<IAtomContainer> expectedFGs, List<IAtomContainer> actualFGs) {
    	Assert.assertThat("Number of functional groups does not match the expected number of groups",
				actualFGs.size(), is(expectedFGs.size()));

		for(int i = 0; i < expectedFGs.size(); i++) {
    		IAtomContainer cExp = expectedFGs.get(i);
    		IAtomContainer cAct = actualFGs.get(i);
    		
    		Assert.assertThat("Groups #" + i + ": different atom count",
					cAct.getAtomCount(), is(cExp.getAtomCount()));
    		Assert.assertThat("Groups #" + i + ": different bond count",
					cAct.getBondCount(), is(cExp.getBondCount()));

			Pattern pattern = VentoFoggia.findIdentical(cExp);

			Assert.assertThat("Groups #" + i + ": not isomorph", pattern.matches(cAct), is(true));
    		
    		Mappings mappings = pattern.matchAll(cAct);

    		Map<IAtom, IAtom> atomMap = mappings.toAtomMap().iterator().next();
    		for (Map.Entry<IAtom,IAtom> e : atomMap.entrySet()) {
    	         IAtom atomExp  = e.getKey();
    	         IAtom atomAct = e.getValue();
    	         Assert.assertThat("Groups #" + i + ": Atom aromaticity does not match" + atomAct.getSymbol() + atomAct.isAromatic() + atomExp.getSymbol() + atomExp.isAromatic(), atomAct.isAromatic(), is(atomExp.isAromatic()));
    	     }

    		Map<IBond, IBond> bondMap = mappings.toBondMap().iterator().next();
    		for (Map.Entry<IBond,IBond> e : bondMap.entrySet()) {
    	         IBond bondExp  = e.getKey();
    	         IBond bondAct = e.getValue();
    	         Assert.assertThat("Groups #" + i + ": Bond aromaticity does not match", bondAct.isAromatic(), is(bondExp.isAromatic()));
    	     }
    	}
    }
    
    private IAtomContainer buildFunctionalGroup(String string) {
        IAtom 					a1, a2, a3, a4, a5, a6, a7, a8, a9;
        IBond 					b1, b2, b3, b4, b5, b6, b7, b8, b9;
        IChemObjectBuilder 		builder = DefaultChemObjectBuilder.getInstance();
        IAtomContainer 			container;

        // custom templates
        switch(string) {
        case "NarR3":
        	a1 = builder.newInstance(IPseudoAtom.class, "R");
        	a2 = builder.newInstance(IPseudoAtom.class, "R");
        	a3 = builder.newInstance(IPseudoAtom.class, "R");
            a4 = builder.newInstance(IAtom.class, "N");
            a4.setIsAromatic(true);
            
            b1 = builder.newInstance(IBond.class, a1, a4, Order.SINGLE);
            b2 = builder.newInstance(IBond.class, a2, a4, Order.SINGLE);
            b3 = builder.newInstance(IBond.class, a3, a4, Order.SINGLE);
                    
            container = new AtomContainer();
            container.setAtoms(new IAtom[] {a1, a2, a3, a4});
            container.setBonds(new IBond[] {b1, b2, b3});
            return container;
            
        case "SarR2":
        	a1 = builder.newInstance(IPseudoAtom.class, "R");
        	a2 = builder.newInstance(IPseudoAtom.class, "R");
            a3 = builder.newInstance(IAtom.class, "S");
            a3.setIsAromatic(true);
            
            b1 = builder.newInstance(IBond.class, a1, a3, Order.SINGLE);
            b2 = builder.newInstance(IBond.class, a2, a3, Order.SINGLE);
                    
            container = new AtomContainer();
            container.setAtoms(new IAtom[] {a1, a2, a3});
            container.setBonds(new IBond[] {b1, b2});
            return container;
            
        case "OarR2":
        	a1 = builder.newInstance(IPseudoAtom.class, "R");
        	a2 = builder.newInstance(IPseudoAtom.class, "R");
            a3 = builder.newInstance(IAtom.class, "O");
            a3.setIsAromatic(true);
            
            b1 = builder.newInstance(IBond.class, a1, a3, Order.SINGLE);
            b2 = builder.newInstance(IBond.class, a2, a3, Order.SINGLE);
                    
            container = new AtomContainer();
            container.setAtoms(new IAtom[] {a1, a2, a3});
            container.setBonds(new IBond[] {b1, b2});
            return container;
            
        case "NarR2COOH":
        	a1 = builder.newInstance(IPseudoAtom.class, "R");
        	a2 = builder.newInstance(IPseudoAtom.class, "R");
            a3 = builder.newInstance(IAtom.class, "N");
            a3.setIsAromatic(true);
            a4 = builder.newInstance(IAtom.class, "C");
            a5 = builder.newInstance(IAtom.class, "O");
            a6 = builder.newInstance(IAtom.class, "O");
            a7 = builder.newInstance(IAtom.class, "H");
            
            b1 = builder.newInstance(IBond.class, a1, a3, Order.SINGLE);
            b2 = builder.newInstance(IBond.class, a2, a3, Order.SINGLE);
            b3 = builder.newInstance(IBond.class, a3, a4, Order.SINGLE);
            b4 = builder.newInstance(IBond.class, a4, a5, Order.DOUBLE);
            b5 = builder.newInstance(IBond.class, a4, a6, Order.SINGLE);
            b6 = builder.newInstance(IBond.class, a6, a7, Order.SINGLE);
                    
            container = new AtomContainer();
            container.setAtoms(new IAtom[] {a1, a2, a3, a4, a5, a6, a7});
            container.setBonds(new IBond[] {b1, b2, b3, b4, b5, b6});
            return container;

            // smiles
        default:
        	try {
        		SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        		try {
        			if(string.equals("[c]=O"))
        				smilesParser.kekulise(false);
        			container = smilesParser.parseSmiles(string);
        		}
        		catch(InvalidSmilesException e) {
        			smilesParser.kekulise(false);
        			container = smilesParser.parseSmiles(string);
        		}
        		
                for(IAtom a : container.atoms()) {
                	if(a instanceof PseudoAtom) {
                		a.setSymbol("R");
                	}
                }
                return container;
        	}
        	catch(InvalidSmilesException e) {
        		throw new IllegalArgumentException("Input string '" + string + " could not be found as a template " +
						"and is not a valid SMILES string.");
        	}
        }
    }
}