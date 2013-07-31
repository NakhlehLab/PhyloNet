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

import java.util.Map;

import phylogeny.EvoTree;

public class HiddenState
{
    protected EvoTree parentalTree;
    protected EvoTree geneGenealogy;

    // maintain mapping between taxa in parentalTree and geneGenealogy
    // only store reference to a shared object
    protected Map<String,String> alleleToSpeciesMapping;

    public HiddenState (EvoTree inParentalTree, EvoTree inGeneGenealogy, Map<String,String> map) {
	setParentalTree(inParentalTree);
	setGeneGenealogy(inGeneGenealogy);
	setAlleleToSpeciesMapping(map);
    }
    
    // For test
    public HiddenState (EvoTree inGeneGenealogy) {
    	setGeneGenealogy(inGeneGenealogy);
    }

    public EvoTree getParentalTree () {
    	return (parentalTree);
    }

    public EvoTree getGeneGenealogy () {
    	return (geneGenealogy);
    }

    public Map<String,String> getAlleleToSpeciesMapping () {
    	return (alleleToSpeciesMapping);
    }

    public void setParentalTree (EvoTree inParentalTree) {
    	this.parentalTree = inParentalTree;
    }

    public void setGeneGenealogy (EvoTree inGeneGenealogy) {
    	this.geneGenealogy = inGeneGenealogy;
    }

    public void setAlleleToSpeciesMapping (Map<String,String> map) {
    	this.alleleToSpeciesMapping = map;
    }

    /**
     * To String Method --> optional arguments
     * @param displayBranchLengths
     * @param displayInternalNodeNames
     * @return
     */
    public String toString (boolean displayBranchLengths, boolean displayInternalNodeNames ) {
		return ("Parental tree:\n" +
			parentalTree.toNewickString(displayBranchLengths, displayInternalNodeNames) + "\n" +
			"Gene genealogy:\n" +
			geneGenealogy.toNewickString(displayBranchLengths, displayInternalNodeNames) + "\n");
    }
    
    /**
     * Default to String Method
     * Default will set display Branch lengths to true
     * and display Internal Node names to false
     */
    public String toString () {
		return this.toString(true,false);
    }
    
    
}
