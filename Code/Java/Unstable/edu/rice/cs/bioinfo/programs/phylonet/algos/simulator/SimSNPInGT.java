package edu.rice.cs.bioinfo.programs.phylonet.algos.simulator;

import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
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


    public Map<String, String> generateSingeSite(Tree tree) {
        if(_seed == null){
            _random = new Random();
        }

        Map<String, String> res = new HashMap<>();
        res = new HashMap<>();
        Map<TNode, Character> markers = new HashMap<>();
        TNode root = tree.getRoot();
        markers.put(root, _random.nextDouble() < _pi0 ? '0' : '1');
        for (TNode child : root.getChildren()) {
            processNode(child, markers);
        }

        for (String leafName : tree.getLeaves()) {
            res.put(leafName, markers.get(tree.getNode(leafName)).toString());
        }
        return res;
    }

    private void processNode(TNode node, Map<TNode, Character> markers) {
        double t = node.getParentDistance();
        double x = 1.0 - Math.exp(-(_u + _v) * t);
        char parentState = markers.get(node.getParent());
        char curState = parentState == '0' ?
                (_random.nextDouble() < (1.0 - _pi1 * x) ? '0' : '1')
                :(_random.nextDouble() < (_pi0 * x) ? '0' : '1');
        markers.put(node, curState);
        for(TNode child : node.getChildren()) {
            processNode(child, markers);
        }
    }

    public static void main(String [] args) {

    }
}