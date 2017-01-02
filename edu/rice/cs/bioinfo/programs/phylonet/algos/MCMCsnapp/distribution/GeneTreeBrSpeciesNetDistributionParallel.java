package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.util.List;
import java.util.Map;

/**
 * Created by wendingqiao on 3/3/16.
 */
public class GeneTreeBrSpeciesNetDistributionParallel extends GeneTreeBrSpeciesNetDistribution {

    private List<UltrametricTree> _geneTrees;
    private int _totalTree;
    private int _currentTreeID = 0;

    public GeneTreeBrSpeciesNetDistributionParallel(Network network, List<UltrametricTree> gts,
                                                    Map<String, List<String>> species2alleles) {
        super(network, species2alleles);
        this._geneTrees = gts;
        this._totalTree = gts.size();
    }

    public void calculateGTDistribution(double[] likelihoodArray) {
        int treeID = getNextTreeID();
        while(treeID < _totalTree){
            UltrametricTree gt = _geneTrees.get(treeID);
            likelihoodArray[treeID] = super.calculateGTDistribution(gt);
            treeID = getNextTreeID();
        }
    }

    public synchronized int getNextTreeID(){
        return _currentTreeID++;
    }
}
