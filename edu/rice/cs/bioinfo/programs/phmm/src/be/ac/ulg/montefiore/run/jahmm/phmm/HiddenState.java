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

/**
 * Container for data associated with a hidden state.
 * Neater to do it this way.
 * Modularity good - may want to extend it later.
 */

package be.ac.ulg.montefiore.run.jahmm.phmm;

import phylogeny.*;

public class HiddenState
{
    protected EvoTree parentalTree;
    protected EvoTree geneGenealogy;

    public HiddenState (EvoTree inParentalTree, EvoTree inGeneGenealogy) {
	setParentalTree(inParentalTree);
	setGeneGenealogy(inGeneGenealogy);
    }

    public EvoTree getParentalTree () {
	return (parentalTree);
    }

    public EvoTree getGeneGenealogy () {
	return (geneGenealogy);
    }

    public void setParentalTree (EvoTree inParentalTree) {
	this.parentalTree = inParentalTree;
    }

    public void setGeneGenealogy (EvoTree inGeneGenealogy) {
	this.geneGenealogy = inGeneGenealogy;
    }

    public String toString () {
	return ("Parental tree:\n" +
		parentalTree.toString() + "\n" +
		"Gene genealogy:\n" +
		geneGenealogy.toString() + "\n");
    }
}
