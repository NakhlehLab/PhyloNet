package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.all;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.List;

/**
 * Created by wendingqiao on 4/29/16.
 * Multiply times of trees and network, population size of network by a random factor
 * TODO by dw20: doesn't change population sizes due to correlations (cite: MCMCcoal)
 */
public class ScaleAll extends AllOperator {

    private double _upperLimit = 1.0 - 1e-8;
    private double _lowerLimit;
    private double _scaleFactor;

    private Double scale;
    private List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> _edges;

    public ScaleAll(List<UltrametricTree> trees, UltrametricNetwork net) {
        super(trees, net);
        _lowerLimit = Math.min(_upperLimit,
                1.0 / Math.exp(20.0 / (trees.size()+1.0) / trees.get(0).getInternalNodes().size())); // TODO: 50
        _scaleFactor = _lowerLimit;
    }

    @Override
    public double propose() {
        scale = getScaler();
        int dimension = 0;
//        int dimension = 1;
//        if(Utils._CONST_POP_SIZE) {
//            _edges = null;
//        } else {
//            _edges = Networks.getAllEdges(_speciesNet.getNetwork());
//            dimension += _edges.size();
//        }
//        scalePopSize (scale);

        dimension += _speciesNet.getInternalNodeCount();
        scaleHeight(scale, false);

        for (UltrametricTree ut : _allTrees) {
            dimension += ut.getInternalNodes().size();
            ut.scale(scale, false);
        }

        return Math.log(scale) * (dimension - 2);
    }

    @Override
    public void undo() {
//        scalePopSize(1.0 / scale);
        scaleHeight(1.0 / scale, true);
        for (UltrametricTree ut : _allTrees) {
            ut.scale(1.0 / scale, true);
        }
    }

    @Override
    public String getName() {
        return "Scale-All";
    }

    @Override
    public void optimize(double logAlpha) {
        double delta = Utils.calcDelta(this, logAlpha);
        delta += Math.log(1.0 / _scaleFactor - 1.0);
        setCoercableParameterValue(1.0 / (Math.exp(delta) + 1.0));
    }

    private void setCoercableParameterValue(final double value) {
        _scaleFactor = Math.max(Math.min(value, _upperLimit), _lowerLimit);
    }

    private double getScaler() {
        return _scaleFactor + Randomizer.getRandomDouble() * (1.0 / _scaleFactor - _scaleFactor);
    }

//    private void scalePopSize(double scale) {
//        NetNode root = _speciesNet.getNetwork().getRoot();
//        root.setRootPopSize(root.getRootPopSize() * scale);
//        if(Utils._CONST_POP_SIZE) {
//            return;
//        }
//        if(_edges == null) {
//            throw new RuntimeException("not CONST population size!!!!");
//        }
//        for(Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge : _edges) {
//            edge.Item1.setParentSupport(edge.Item2, edge.Item1.getParentSupport(edge.Item2) * scale);
//        }
//    }

    private void scaleHeight(double scale, boolean undo) {
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(_speciesNet.getNetwork())) {
            if(node.isLeaf()) continue;
            if(undo) {
                node.getData().recoverHeight();
            } else {
                node.getData().storeHeight();
                double newH = node.getData().getHeight() * scale;
                node.getData().setHeight(newH);
            }
            for(NetNode<NetNodeInfo> child : node.getChildren()) {
                child.setParentDistance(node, node.getData().getHeight() - child.getData().getHeight());
            }
        }
    }
}
