/*
 * Copyright (c) 2013 John May <jwmay@users.sf.net>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 U
 */
package org.openscience.cdk.hash;

import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;

/**
 * A fast pseudorandom number generator based on feedback shift registers.
 *
 * @author John May
 * @see <a href="http://en.wikipedia.org/wiki/Xorshift">Xorshift</a>
 * @see <a href="http://www.javamex.com/tutorials/random_numbers/xorshift.shtml">Xorshift
 *      random number generators</a>
 * @cdk.githash
 * @cdk.module hash
 */
@TestClass("org.openscience.cdk.hash.XorshiftTest")
final class Xorshift extends Pseudorandom {

    /**
     * Generate the next pseudorandom number for the provided <i>seed</i>.
     *
     * @param seed random number seed
     * @return the next pseudorandom number
     */
    @TestMethod("testNext,testDistribution,demonstrateZeroLimitation")
    @Override
    long next(long seed) {
        seed = seed ^ seed << 21;
        seed = seed ^ seed >>> 35;
        return seed ^ seed << 4;
    }
}