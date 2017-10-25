/* Copyright (C) 1997-2007  The Chemistry Development Kit (CDK) project
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.tools;

import java.io.StringWriter;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.CDKTestCase;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.PseudoAtom;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 * @cdk.module test-standard
 * @author Sebastian Fritsch
 */
public class FunctionalGroupsFinderTest extends CDKTestCase {

    public FunctionalGroupsFinderTest() {
        super();
    }
    
    @Test
    public void testExtractGroups() {
    	SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer mol = new AtomContainer();
		try {
			mol = smiPar.parseSmiles("NC(O)C\\C=C\\CC(=O)C1CO1");
			//mol = smiPar.parseSmiles("CC(=O)N[C@@H]1[C@@H](NC(=N)N)C=C(O[C@H]1[C@H](O)[C@H](O)CO)C(=O)O");
			
			addExplicitHydrogens(mol);
		} catch (Exception e) {
			Assert.assertFalse(true);
		}
		try {
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
			Aromaticity.cdkLegacy().apply(mol);
		} catch (Exception e) {
			Assert.assertFalse(true);
		}
		
		FunctionalGroupsFinder gF = new FunctionalGroupsFinder();
		try {
			gF.find(mol);
		} catch (CloneNotSupportedException e1) {
			Assert.assertFalse(true);
		}
		gF.markAtoms(mol);
		
		List<IAtomContainer> fGs = new ArrayList<>();
		try {
			fGs = gF.extractGroups(mol);
		} catch (CloneNotSupportedException e) {
			Assert.assertFalse(true);
		}
		
		printToConsoleWithIndices(mol);
		
		System.out.println(fGs.size() + " FUNCTIONAL GROUPS:");
		this.printMolFiles(fGs);

		
		fGs = gF.generalizeEnvironments(fGs);
		
		System.out.println("GENERALIZED FUNCTIONAL GROUPS:");
		this.printMolFiles(fGs);
    }
    
    @Test
    public void testFind1() throws Exception {
    	SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer mol = smilesParser.parseSmiles("Cc1cc(C)nc(NS(=O)(=O)c2ccc(N)cc2)n1");
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		Aromaticity.cdkLegacy().apply(mol);
        
        this.printToConsoleWithIndices(mol);
        
        FunctionalGroupsFinder gF = new FunctionalGroupsFinder();
        gF.find(mol);
        gF.markAtoms(mol);
        List<IAtomContainer> fGs = gF.extractGroups(mol);
        
        Assert.assertEquals(4, fGs.size());
        
        System.out.println(fGs.size() + " FUNCTIONAL GROUPS:");
        this.printMolFiles(fGs);
        
        for(IAtomContainer fG : fGs) {
			System.out.println(String.format("FG: %d atoms and %d bonds", fG.getAtomCount(), fG.getBondCount()));
			this.assertNoNonseseBonds(fG);
		}
        
        fGs = gF.generalizeEnvironments(fGs);
        
		for(IAtomContainer fG : fGs) {
			System.out.println(String.format("FG: %d atoms and %d bonds", fG.getAtomCount(), fG.getBondCount()));
			this.assertNoNonseseBonds(fG);
		}
        
        System.out.println("GENERALIZED FUNCTIONAL GROUPS:");
		this.printMolFiles(fGs);
        
        Assert.assertEquals(4, fGs.size());
        
        // EXPECTED functional groups:
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        SmilesParser sp = new SmilesParser(builder);
        List<IAtomContainer> expectedFGs = new LinkedList<>();
        
        IAtom a1, a2, a3, a4;
        IBond b1, b2, b3;
        IAtomContainer container;
        
        // R-:N=:R
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
        expectedFGs.add(container);
        
        // [R]N([R])S(=O)(=O)[R]
        container = sp.parseSmiles("[R]N([R])S(=O)(=O)[R]");
        for(IAtom a : container.atoms()) {
        	if(a instanceof PseudoAtom) {
        		a.setSymbol("R");
        	}
        }
        expectedFGs.add(container);
        
        // C_ar-NH2
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
        expectedFGs.add(container);
        
        // R-:N=:R
        container = expectedFGs.get(0).clone();
        expectedFGs.add(container);
        
        
        
        UniversalIsomorphismTester uiTester = new UniversalIsomorphismTester();
        ArrayList<IAtomContainer> expFGsCopy = new ArrayList<>(expectedFGs);
        for(IAtomContainer fG : fGs) {
        	for(int i = 0; i < expFGsCopy.size(); i++) {
        		if(uiTester.isIsomorph(fG, expFGsCopy.get(i))){
        			expFGsCopy.remove(i);
        			break;
        		}
        	}
        }
        Assert.assertEquals(0, expFGsCopy.size());
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
    
    private void printMolFiles(List<IAtomContainer> fGs) {
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
}
