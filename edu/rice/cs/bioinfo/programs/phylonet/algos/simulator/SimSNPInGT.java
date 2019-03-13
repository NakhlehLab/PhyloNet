package edu.rice.cs.bioinfo.programs.phylonet.algos.simulator;

import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.DiscreteGammaDistribution;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 11/27/16
 * Time: 9:31 AM
 * To change this template use File | Settings | File Templates.
 */
public class SimSNPInGT {
    private BiAllelicGTR _BAGTRModel;
    private double _theta;
    private double _pi0;
    private double _pi1;
    private double _u;
    private double _v;
    private Random _random;
    private Long _seed = null;
    private DiscreteGammaDistribution _gamma = null;
    private boolean _rateVarAcrossMarkers = false;
    private boolean _rateVarAcrossLineages = false;
    private Double _invariantSitesProb = null;
    private double _gammaSampled = 0;


    public SimSNPInGT(BiAllelicGTR BAGTRModel) {
        _BAGTRModel = BAGTRModel;
        //_theta = theta;
        _pi0 = _BAGTRModel.getEquilibriumVector().get(0, 0);
        _pi1 = _BAGTRModel.getEquilibriumVector().get(1, 0);
        _u = 1.0 / (2.0 * _pi0);
        _v = 1.0 / (2.0 * _pi1);
    }

    /**
     * This function is to set seed to control the randomness
     */
    public void setSeed(Long seed){
        _seed = seed;
        if(_seed != null)
            _random = new Random(_seed);
        else
            _random = new Random();
    }

    public void setDiscreteGamma(DiscreteGammaDistribution gamma) {
        _gamma = gamma;
    }

    public void setRateVarAcrossMarkers(boolean setting) {
        _rateVarAcrossMarkers = setting;
    }

    public void setRateVarAcrossLineages(boolean setting) {
        _rateVarAcrossLineages = setting;
    }

    public void setInvariantSitesProb(double prob) {
        _invariantSitesProb = prob;
    }

    public Map<String, String> generateSingeSite(Tree tree, Map<TNode, Character> markers) {
        if(markers == null) {
            markers = new HashMap<>();
        }
        markers.clear();

        if(_seed == null){
            _random = new Random();
        }

        Map<String, String> res = new HashMap<>();
        res = new HashMap<>();
        TNode root = tree.getRoot();
        // Need to burn a random number, otherwise the bahavior is weird.
        _random.nextDouble();
        markers.put(root, _random.nextDouble() < _pi0 ? '0' : '1');

        if(_rateVarAcrossMarkers) {
            _gammaSampled = _gamma.sample();
        }

        for (TNode child : root.getChildren()) {
            boolean invariant = (_invariantSitesProb != null) && (_random.nextDouble() < _invariantSitesProb);
            processNode(child, markers, invariant);
        }

        for (String leafName : tree.getLeaves()) {
            res.put(leafName, markers.get(tree.getNode(leafName)).toString());
        }
        return res;
    }

    public Map<String, String> generateSingeSite(Tree tree) {
        return generateSingeSite(tree, null);
    }

    private void processNode(TNode node, Map<TNode, Character> markers, boolean invariant) {
        double t = node.getParentDistance();

        if(_rateVarAcrossLineages) {
            t *= _gamma.sample();
        } else if(_rateVarAcrossMarkers) {
            t *= _gammaSampled;
        }

        double x = 1.0 - Math.exp(-(_u + _v) * t);
        char parentState = markers.get(node.getParent());
        char curState = parentState == '0' ?
                (_random.nextDouble() < (1.0 - _pi1 * x) ? '0' : '1')
                :(_random.nextDouble() < (_pi0 * x) ? '0' : '1');

        if(invariant) {
            curState = parentState;
        }
        markers.put(node, curState);
        for(TNode child : node.getChildren()) {
            processNode(child, markers, invariant);
        }
    }

    public static void main(String [] args) {

    }
}