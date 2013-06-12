/*
 * Copyright (c) 2013 European Bioinformatics Institute (EMBL-EBI)
 *                    John May <jwmay@users.sf.net>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version. All we ask is that proper credit is given
 * for our work, which includes - but is not limited to - adding the above
 * copyright notice to the beginning of your source code files, and to any
 * copyright notice that you may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 U
 */
package org.openscience.cdk.graph;

import org.junit.Test;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.*;
import static org.openscience.cdk.graph.InitialCyclesTest.anthracene;
import static org.openscience.cdk.graph.InitialCyclesTest.bicyclo;
import static org.openscience.cdk.graph.InitialCyclesTest.cyclophane_even;
import static org.openscience.cdk.graph.InitialCyclesTest.napthalene;

/**
 * @author John May
 * @cdk.module test-core
 */
public class MinimumCycleBasisTest {

    @Test(expected = NullPointerException.class)
    public void noGraph() {
        new RelevantCycles((int[][]) null);
    }

    @Test(expected = NullPointerException.class)
    public void noInitialCycles() {
        new RelevantCycles((InitialCycles) null);
    }

    @Test public void paths_bicyclo() {
        int[][] bicyclo = bicyclo();
        MinimumCycleBasis mcb = new MinimumCycleBasis(bicyclo);
        int[][] paths = mcb.paths();
        assertThat(paths.length, is(2));
        int[][] expected = new int[][]{{5, 0, 1, 2, 3, 4},
                                       {5, 0, 1, 2, 7, 6}};
        assertThat(paths, is(expected));
    }

    @Test public void paths_napthalene() {
        int[][] napthalene = napthalene();
        MinimumCycleBasis mcb = new MinimumCycleBasis(napthalene);
        int[][] paths = mcb.paths();
        assertThat(paths.length, is(2));
        int[][] expected = new int[][]{{5, 0, 1, 2, 3, 4},
                                       {5, 4, 7, 8, 9, 6}};
        assertThat(paths, is(expected));
    }

    @Test public void paths_anthracene() {
        int[][] anthracene = anthracene();
        MinimumCycleBasis mcb = new MinimumCycleBasis(anthracene);
        int[][] paths = mcb.paths();
        assertThat(paths.length, is(3));
        int[][] expected = new int[][]{{5, 0, 1, 2, 3, 4},
                                       {9, 6, 5, 4, 7, 8},
                                       {9, 8, 10, 11, 12, 13}};
        assertThat(paths, is(expected));
    }

    @Test public void paths_cyclophane_even() {
        int[][] cyclophane_even = cyclophane_even();
        MinimumCycleBasis mcb = new MinimumCycleBasis(cyclophane_even);
        int[][] paths = mcb.paths();
        assertThat(paths.length, is(2));
        int[][] expected = new int[][]{{3, 2, 1, 0, 5, 4},
                                       {3, 6, 7, 8, 9, 10, 11, 0, 1, 2}};
        assertThat(paths, is(expected));
    }

    @Test public void paths_cyclophane_odd() {
        int[][] cyclophane_even = cyclophane_even();
        MinimumCycleBasis mcb = new MinimumCycleBasis(cyclophane_even);
        int[][] paths = mcb.paths();
        assertThat(paths.length, is(2));
        int[][] expected = new int[][]{{3, 2, 1, 0, 5, 4},
                                       {3, 6, 7, 8, 9, 10, 11, 0, 1, 2}};
        assertThat(paths, is(expected));
    }

    @Test public void size_bicyclo() {
        int[][] bicyclo = bicyclo();
        MinimumCycleBasis mcb = new MinimumCycleBasis(bicyclo);
        assertThat(mcb.size(), is(2));
    }

    @Test public void size_napthalene() {
        int[][] napthalene = napthalene();
        MinimumCycleBasis mcb = new MinimumCycleBasis(napthalene);
        assertThat(mcb.size(), is(2));
    }

    @Test public void size_anthracene() {
        int[][] anthracene = anthracene();
        MinimumCycleBasis mcb = new MinimumCycleBasis(anthracene);
        assertThat(mcb.size(), is(3));
    }

    @Test public void size_cyclophane_even() {
        int[][] cyclophane_even = cyclophane_even();
        MinimumCycleBasis relevant = new MinimumCycleBasis(cyclophane_even);
        assertThat(relevant.size(), is(2));
    }

    @Test public void size_cyclophane_odd() {
        int[][] cyclophane_even = cyclophane_even();
        MinimumCycleBasis mcb = new MinimumCycleBasis(cyclophane_even);
        assertThat(mcb.size(), is(2));
    }
}
