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
import java.io.Writer;

import be.ac.ulg.montefiore.run.jahmm.ObservationInteger;
import be.ac.ulg.montefiore.run.jahmm.OpdfInteger;

/**
 * This class can write a textual description of an {@link OpdfInteger}.
 * It is compatible with {@link OpdfIntegerReader}.
 */
public class OpdfIntegerWriter
extends OpdfWriter<OpdfInteger>
{
	public void write(Writer writer, OpdfInteger opdf)
	throws IOException
	{
		String s = "IntegerOPDF [";
		
		for (int i = 0; i < opdf.nbEntries(); i++)
			s += opdf.probability(new ObservationInteger(i)) + " ";
			
		writer.write(s + "]\n");
	}
}
