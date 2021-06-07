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

package be.ac.ulg.montefiore.run.jahmm.test;

import java.io.IOException;

import junit.framework.TestCase;
import be.ac.ulg.montefiore.run.jahmm.*;
import be.ac.ulg.montefiore.run.jahmm.draw.GenericHmmDrawerDot;


public class GenerateTest
extends TestCase
{
    public final static String outputDir = "";

    private Hmm<ObservationInteger> hmm;


    protected void setUp()
    {
        hmm = new Hmm<ObservationInteger>(4, new OpdfIntegerFactory(2));
    }


    public void testDotGenerator()
    {
        GenericHmmDrawerDot hmmDrawer = new GenericHmmDrawerDot();

        try {
            hmmDrawer.write(hmm, outputDir + "hmm-generate.dot");
        }
        catch (IOException e) {
            assertTrue("Writing file triggered an exception: " + e, false);
        }
    }
}
