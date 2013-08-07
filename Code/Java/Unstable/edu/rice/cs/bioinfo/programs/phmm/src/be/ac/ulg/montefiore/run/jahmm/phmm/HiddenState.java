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
import java.util.Vector;
import java.util.List;
import java.io.StringReader;
import java.io.IOException;
import phylogeny.EvoTree;

// kliu - Phylonet support libraries
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.ExNewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;


public class HiddenState
{
    protected EvoTree parentalTree;
    protected EvoTree geneGenealogy;
    protected int parentalTreeID;

    // maintain mapping between taxa in parentalTree and geneGenealogy
    // only store reference to a shared object
    protected Map<String,String> alleleToSpeciesMapping;

    /**
     * For coalescent model calculations.
     */
    protected GeneTreeProbability gtp;

    public HiddenState (EvoTree inParentalTree, int inParentalTreeID, EvoTree inGeneGenealogy, Map<String,String> map) {
    setParentalTree(inParentalTree);
    setGeneGenealogy(inGeneGenealogy);
    setAlleleToSpeciesMapping(map);
    gtp = new GeneTreeProbability();
    }

    // For test
    public HiddenState (EvoTree inGeneGenealogy) {
        setGeneGenealogy(inGeneGenealogy);
    }

    public EvoTree getParentalTree () {
        return (parentalTree);
    }

    public int getParentalTreeID() {
        return (parentalTreeID);
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

    public void setParentalTreeID(int inParentalTreeID) {
        this.parentalTreeID = inParentalTreeID;
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

    /**
     * Expose this method to other classes, since transition probability calculation will need to access this.
     * Default to no debug messages.
     */
    public double calculateProbabilityOfGeneGenealogyInParentalTree () {
    return (calculateProbabilityOfGeneGenealogyInParentalTree(false));
    }

    /**
     * Perform standard coalescent model calculation to obtain
     * probability P[g(s_i) | T(s_i), c_{T(s_i)}] of observing a gene genealogy given a parental tree.
     *
     * See writeup for details.
     *
     * Might move this later to HiddenState.
     *
     * Use code from ComputeGTProb.
     * WARNING - returns likelihood, *NOT* log likelihood!
     *
     */
    protected double calculateProbabilityOfGeneGenealogyInParentalTree (boolean debugFlag) {
    Network<Double> parentalTree = convertPHMMTreeToPhyloNetNetwork(getParentalTree());
    Map<String,String> alleleToSpeciesMapping = getAlleleToSpeciesMapping();
    Tree geneGenealogy = convertPHMMTreeToPhyloNetTree(getGeneGenealogy());
    // A list with one element. Inefficient - consider doing a one-shot approach later.
    Vector<Tree> geneGenealogies = new Vector<Tree>();
    geneGenealogies.add(geneGenealogy);

    gtp.emptyState();

    // look like the calculation is proceeding OK
    //
    // calculation under model from Yu et al. 2012
    // this method requires Network<Double>
    // Yun uses Double to store hybridization probabilities during calculation
        List<Double> probList = gtp.calculateGTDistribution(parentalTree, geneGenealogies, alleleToSpeciesMapping, debugFlag);

    // should only be a single entry
    if (probList.size() != 1) {
        System.err.println ("ERROR: GeneTreeProbability.calculateGTDistribution(...) didn't return exactly one probability. Returning -1 to signal error.");
        return (-1.0);
    }

        return (probList.get(0));
    }

    /**
     * Convert an EvoTree into a PhyloNet Tree object.
     */
    protected Tree convertPHMMTreeToPhyloNetTree (EvoTree etree) {
    NewickReader nr = new NewickReader(new StringReader(etree.toNewickString(true, true)));
    // kliu - change this over to a String data member
    // don't really use it anyways
    STITree<Double> newtr = new STITree<Double>(true);
    try {
        nr.readTree(newtr);
        return (newtr);
    }
    catch(Exception e) {
        System.err.println(e);
        e.printStackTrace();
        return (null);
    }
    }

    /**
     * Convert an EvoTree into a PhyloNet Network object
     */
    protected Network<Double> convertPHMMTreeToPhyloNetNetwork (EvoTree etree) {
    ExNewickReader<Double> enr = new ExNewickReader<Double>(new StringReader(etree.toNewickString(true, true)));
    // kliu - change this over to a String data member
    // don't really use it anyways
    try {
        Network<Double> network = enr.readNetwork();
        return (network);
    }
    catch(IOException ioe) {
        System.err.println(ioe);
        ioe.printStackTrace();
        return (null);
    }
    }

}
