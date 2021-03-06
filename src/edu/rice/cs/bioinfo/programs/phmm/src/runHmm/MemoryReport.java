/**
 * This file is part of PhyloNet-HMM.
 *
 * Copyright © 2013-2014 Kevin Liu, Jingxuan Dai, Kathy Truong, 
 * Ying Song, Michael H. Kohn, and Luay Nakhleh. <http://bioinfo.cs.rice.edu/>
 * 
 * PhyloNet-HMM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * PhyloNet-HMM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phmm.src.runHmm;

public class MemoryReport {
    private static final long MiB = 1024 * 1024;

    public static void report() {
        Runtime runtime = Runtime.getRuntime();

        System.err.println();

        long total = runtime.totalMemory() / MiB;
        long free = runtime.freeMemory() / MiB;
        long used = total - free;
        System.err.println("Total: " + total + " MiB; Free: " + free + " MiB; --> Used: " + used + " MiB.");

        long max = runtime.maxMemory() / MiB;
        System.err.println("Max: " + max + " MiB.");
        System.err.println();
    }
}
