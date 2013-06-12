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

import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;

import static com.google.common.base.Preconditions.checkNotNull;
import static org.openscience.cdk.graph.InitialCycles.Cycle;

/**
 * Compute the minimum cycle basis (MCB) of a graph. A cycle basis is a set of
 * cycles that can be combined to generate all the cycles of a graph. The MCB is
 * the basis of minimum weight (length). As an example, there are three cycles
 * in <a href="http://en.wikipedia.org/wiki/Naphthalene"><i>naphthalene</i></a>.
 * To generate the three cycles (cycle space) we actually only need two of the
 * three cycles. The third cycle can be generated by combining the other two.
 * Each combination of the two cycles is called a cycle basis. There is one
 * basis with the two rings of size six and two bases with either ring of size
 * six and the perimeter ring of size 10. The weights of each basis are 12
 * (6+6), 16 (6+10), 16 (6+10). The basis of the two six member rings has
 * minimum weight (12) and is thus the MCB of <i>naphthalene</i>.<p/>
 *
 * The Smallest Set of Smallest Rings (SSSR) is normally used interchangeably
 * with MCB although traditionally their meaning was different {@cdk.cite
 * Berger04}. Although the MCB can be used to generate all other cycles of a
 * graph it may not be unique. A compound with a bridged ring system, such as <a
 * href="http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:115239">
 * <i>3-quinuclidinol</i></a>, has three equally valid minimum cycle bases. From
 * any two six member rings the third ring can be generated by &oplus;-sum
 * (xoring) their edges. As there may be more than one MCB it should not be used
 * as a <i>unique</i> descriptor. The smallest (but potentially exponential)
 * <i>unique</i> set of short cycles which generates the cycle space is the
 * union of all minimum cycle bases. This set is known as the {@link
 * RelevantCycles}. The intersect of all minimum cycle bases ({@link
 * EssentialCycles}) is also unique. Although the number of {@link
 * EssentialCycles} is polynomial it can not always be used to generate the
 * cycle space.
 *
 * @author John May
 * @cdk.module core
 * @cdk.keyword sssr
 * @cdk.keyword smallest set of smallest rings
 * @cdk.keyword mcb
 * @cdk.keyword minimum cycle basis
 * @cdk.keyword cycle
 * @cdk.keyword ring
 * @cdk.keyword ring perception
 * @see RelevantCycles
 * @see org.openscience.cdk.ringsearch.SSSRFinder#findSSSR()
 * @see <a href="http://en.wikipedia.org/wiki/Cycle_space">Cycle Basis</a>
 */
@TestClass("org.openscience.cdk.graph.MinimumCycleBasisTest")
public final class MinimumCycleBasis {

    /** The minimum cycle basis. */
    private final GreedyBasis basis;

    /**
     * Generate the minimum cycle basis for a graph.
     *
     * @param graph undirected adjacency list
     * @see org.openscience.cdk.ringsearch.RingSearch#fused()
     * @see GraphUtil#subgraph(int[][], int[])
     */
    @TestMethod("noGraph")
    public MinimumCycleBasis(final int[][] graph) {
        this(new InitialCycles(graph));
    }

    /**
     * Generate the minimum cycle basis from a precomputed set of initial
     * cycles.
     *
     * @param initial set of initial cycles.
     * @throws NullPointerException null InitialCycles provided
     */
    @TestMethod("noInitialCycles")
    MinimumCycleBasis(final InitialCycles initial) {

        checkNotNull(initial, "No InitialCycles provided");

        this.basis = new GreedyBasis(initial.numberOfCycles(),
                                     initial.numberOfEdges());

        // processing by size add cycles which are independent of smaller cycles
        for (final Cycle cycle : initial.cycles()) {
            if (basis.isIndependent(cycle))
                basis.add(cycle);
        }
    }

    /**
     * The paths of all cycles in the minimum cycle basis.
     *
     * @return array of vertex paths
     */
    @TestMethod("paths_bicyclo,paths_napthalene,paths_anthracene," +
                        "paths_cyclophane_odd,paths_cyclophane_even")
    public int[][] paths() {
        final int[][] paths = new int[size()][0];
        int i = 0;
        for (final Cycle c : basis.members())
            paths[i++] = c.path();
        return paths;
    }

    /**
     * The number of the cycles in the minimum cycle basis.
     *
     * @return size of minimum cycle set
     */
    @TestMethod("size_bicyclo,size_napthalene,size_anthracene," +
                        "size_cyclophane_odd,size_cyclophane_even")
    public int size() {
        return basis.size();
    }

}
