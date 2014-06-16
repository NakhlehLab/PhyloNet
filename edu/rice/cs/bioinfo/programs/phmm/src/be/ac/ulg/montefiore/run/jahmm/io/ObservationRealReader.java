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

/* jahmm package - v0.6.1 */

/*
  *  Copyright (c) 2004-2006, Jean-Marc Francois.
 *
 *  This file is part of Jahmm.
 *  Jahmm is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Jahmm is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Jahmm; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

 */

package be.ac.ulg.montefiore.run.jahmm.io;

import java.io.IOException;
import java.io.StreamTokenizer;

import be.ac.ulg.montefiore.run.jahmm.ObservationReal;


/**
 * Reads an {@link be.ac.ulg.montefiore.run.jahmm.ObservationReal
 * ObservationReal} up to (and including) a semi-colon.
 * <p>
 * The format of this observation is a simple real number
 * (following the format [+-]?[0123456789]+[.]?[0123456789]*, <i>e.g.</i>
 * 12.3).
 * <p>
 * For example, reading
 * <pre>76.3;</pre>
 * creates an observation such as the one generated by
 * <code>new ObservationReal(76.3);</code>
 */
public class ObservationRealReader
extends ObservationReader<ObservationReal>
{
    /**
     * Reads an {@link be.ac.ulg.montefiore.run.jahmm.ObservationReal
     * ObservationReal} reader, as explained in
     * {@link ObservationReader ObservationReader}.
     *
     * @param st A stream tokenizer.
     * @return The {@link be.ac.ulg.montefiore.run.jahmm.ObservationReal
     *         ObservationReal} read.
     */
    public ObservationReal read(StreamTokenizer st)
    throws IOException, FileFormatException
    {
        ObservationReal or;

        switch (st.nextToken()) {
        case StreamTokenizer.TT_EOL:
        case StreamTokenizer.TT_EOF:
        case StreamTokenizer.TT_WORD:
            throw new FileFormatException("Real value expected");

        case StreamTokenizer.TT_NUMBER:
            or = new ObservationReal(st.nval);
            break;

        default:
            throw new FileFormatException("Real value expected");
        }

        switch (st.nextToken()) {
        case ';':
            break;

        default:
            if (st.ttype != ';')
                throw new FileFormatException("';' expected");
        }

        return or;
    }
}