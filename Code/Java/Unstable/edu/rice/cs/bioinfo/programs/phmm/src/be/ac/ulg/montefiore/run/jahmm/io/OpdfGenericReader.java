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

package be.ac.ulg.montefiore.run.jahmm.io;

import java.io.IOException;
import java.io.StreamTokenizer;

import be.ac.ulg.montefiore.run.jahmm.Opdf;

public class OpdfGenericReader
extends OpdfReader<Opdf<?>>
{
    String keyword()
    {
        throw new AssertionError("Cannot call method");
    }


    public Opdf<?> read(StreamTokenizer st)
    throws IOException, FileFormatException
    {
        if (st.nextToken() != StreamTokenizer.TT_WORD)
            throw new FileFormatException("Keyword expected");

        for (OpdfReader r : new OpdfReader[] {
                new OpdfIntegerReader(),
                new OpdfGaussianReader(),
                new OpdfGaussianMixtureReader(),
                new OpdfMultiGaussianReader() })
            if (r.keyword().equals(st.sval)) {
                st.pushBack();
                return r.read(st);
            }

        throw new FileFormatException("Unknown distribution");
    }
}
