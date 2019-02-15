/* Copyright (C) 1997-2018  The Chemistry Development Kit (CDK) project
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
package org.openscience.cdk.tools.fgfinderperformance;

/**
 * Main class for starting application.
 * 
 * @author Jonas Schaub
 */
public class Main {
    
    /**
     * Starts the application. Command line arguments must be the name of an SD-file to read (must be located in the 
     * same directory as the application's .jar file) and the number of different threads to use for calculation.
     * 
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            ErtlFunctionalGroupsFinderPerformanceTest tmpApplication = new ErtlFunctionalGroupsFinderPerformanceTest(args);
        } catch (Exception anException) {
            anException.printStackTrace(System.err);
            System.exit(1);
        }
    }
    
}


