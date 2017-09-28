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

import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.CDKTestCase;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.GraphUtil.EdgeToBondMap;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.io.MDLV3000Writer;
import org.openscience.cdk.sgroup.Sgroup;
import org.openscience.cdk.sgroup.SgroupType;
import org.openscience.cdk.silent.Element;
import org.openscience.cdk.silent.PseudoAtom;
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
    public void testGetConnectedAtoms() {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer mol = new AtomContainer();
		try {
			mol = smiPar.parseSmiles("CO");
			addExplicitHydrogens(mol);
		} catch (Exception e) {
			Assert.assertFalse(true);
		}
		
        
        int cIdx = -1;
        for(int idx = 0; idx < mol.getAtomCount(); idx++) {
        	if(mol.getAtom(idx).getAtomicNumber() == 6) {
        		cIdx = idx;
        		break;
        	}
        }
        
        EdgeToBondMap bondMap = GraphUtil.EdgeToBondMap.withSpaceFor(mol);
        int[][] adjList = GraphUtil.toAdjList(mol, bondMap);
        
        List<IAtom> actual = new ArrayList<>(Arrays.asList(FunctionalGroupsFinder.getConnectedAtoms(cIdx, mol, adjList)));
        List<IAtom> expected = mol.getConnectedAtomsList(mol.getAtom(cIdx));
        
        System.out.println(actual.stream().map(a -> a.getSymbol()).collect(Collectors.toList()));
        System.out.println(expected.stream().map(a -> a.getSymbol()).collect(Collectors.toList()));
        
        Assert.assertEquals(actual, expected);
    }
    
    @Test
    public void testDebugConditionMap() {
    	SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer mol = new AtomContainer();
		try {
			mol = smiPar.parseSmiles("NC(O)C\\C=C\\CC(=O)C1CO1");
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
		gF.find(mol);
		
		printToConsoleWithIndices(mol);
    }
    
    @Test
    public void testExtractGroups() {
    	SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer mol = new AtomContainer();
		try {
			mol = smiPar.parseSmiles("NC(O)C\\C=C\\CC(=O)C1CO1");
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
		gF.find(mol);
		List<IAtomContainer> fGs = new ArrayList<>();
		try {
			fGs = gF.extractGroups(mol);
		} catch (CloneNotSupportedException e) {
			Assert.assertFalse(true);
		}
		
		printToConsoleWithIndices(mol);
		
		System.out.println("FUNCTIONAL GROUPS:");
		SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.AtomicMass);
		for(IAtomContainer aC : fGs) {
			String smiles = new String();
			try {
				smiles = smiGen.create(aC);
			} catch (CDKException e) {
				System.out.println("#ERROR# Could not generate SMILES!");
			}
			System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
			System.out.println(smiles);
			System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
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
}
