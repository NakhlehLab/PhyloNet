package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.distribution;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.TreeEmbedding;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.util.List;
import java.util.Map;

/**
 * Created by wendingqiao on 3/3/16.
 */
public class GeneTreeBrSpeciesNetDistributionParallel extends GeneTreeBrSpeciesNetDistribution {

    private List<UltrametricTree> _geneTrees;
    private List<TreeEmbedding> _embeddings;
    private int _totalTree;
    private int _currentTreeID = 0;

    public GeneTreeBrSpeciesNetDistributionParallel(Network network, List<UltrametricTree> gts,
                                                    List<TreeEmbedding> embeddings,
                                                    Map<String, List<String>> species2alleles) {
        super(network, species2alleles);
        this._geneTrees = gts;
        this._embeddings = embeddings;
        this._totalTree = gts.size();
    }

    public void calculateGTDistribution(double[] likelihoodArray) {
        int treeID = getNextTreeID();
        while(treeID < _totalTree){
            UltrametricTree gt = _geneTrees.get(treeID);
            TreeEmbedding embedding = _embeddings == null ? null : _embeddings.get(treeID);
            likelihoodArray[treeID] = super.calculateGTDistribution(gt, embedding);
            treeID = getNextTreeID();
        }
    }

    public synchronized int getNextTreeID(){
        return _currentTreeID++;
    }
}
