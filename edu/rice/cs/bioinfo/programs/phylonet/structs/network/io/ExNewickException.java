/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phylonet.structs.network.io;

/**
 * This class is used in parsing eNewick-format network input.
 */
public class ExNewickException extends RuntimeException {
	// Constructors
	/**
	 * Instantiates a new <code>ExNewickException</code> with the given
	 * message.
	 * 
	 * @param msg	the specific reason for this exception's existence
	 */
	public ExNewickException(String msg)
	{
		super(msg);
	}
	
	// Data members
	final static long serialVersionUID = 1;
}
