package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move.network.topo;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.distribution.SpeciesNetPriorDistribution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move.network.NetworkOperator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by wendingqiao on 3/7/16.
 * select a random internal node, select a random child of node
 * select a random height between [child, node]
 * select a random edge cross with the height where edge.child != child and edge.parent != node
 * swap parents
 * the same as SMARTIE
 * http://sysbio.oxfordjournals.org/cgi/pmidlookup?view=long&pmid=20525618
 */
public class SwapNodes extends NetworkOperator {

    private double _logHR;
    private Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> _edge1, _edge2;

    public SwapNodes(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        _edge1 = _edge2 = null; //reset

        List<NetNode<NetNodeInfo>> internalNodes = Networks.getInternalNodes(_network.getNetwork());
        NetNode<NetNodeInfo> parent, child;
        do {
            parent = internalNodes.get(Randomizer.getRandomInt(internalNodes.size()));
        } while (parent.isRoot());

        List<NetNode<NetNodeInfo>> children = IterableHelp.toList(parent.getChildren());
        child = children.get(Randomizer.getRandomInt(children.size()));
        _edge1 = new Tuple<>(child, parent);

        List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> edges = getEdgesGivenEdgeHeight(_edge1);
        if(edges.size() == 0) {
            _logHR = Utils.INVALID_MOVE;
            _violate = false;
        } else {
            _edge2 = edges.get(Randomizer.getRandomInt(edges.size()));

            if(_edge1.Item1.hasParent(_edge2.Item2) || _edge2.Item1.hasParent(_edge1.Item2)) {
                _logHR = Utils.INVALID_MOVE;
                _violate = false;
            } else {
                swap(_edge1.Item1, _edge1.Item2, _edge2.Item1, _edge2.Item2);

                _logHR = 0;
                _violate = true;
            }
        }
        return _logHR;
    }

    @Override
    public void undo() {
        if(_logHR == Utils.INVALID_MOVE) return;
        swap(_edge1.Item1, _edge2.Item2, _edge2.Item1, _edge1.Item2);
    }

    @Override
    public String getName() {
        return "Swap-Nodes";
    }

    private void swap(NetNode<NetNodeInfo> c1, NetNode<NetNodeInfo> p1, NetNode<NetNodeInfo> c2, NetNode<NetNodeInfo> p2) {

        double rootPopSize = _network.getNetwork().getRoot().getRootPopSize();
        if(rootPopSize == Double.NaN) throw new RuntimeException("Invalid root pop size in SwapKernel!!!");

        double[] paramP1C1 = getParameters(p1, c1);
        double[] paramP2C2 = getParameters(p2, c2);

        p1.removeChild(c1);
        p2.removeChild(c2);

        adopt(p1, c2, paramP2C2);
        adopt(p2, c1, paramP1C1);
    }

    private List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>>
    getEdgesGivenEdgeHeight(Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge) {

        double height = edge.Item2.getData().getHeight();
        List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> edges = new ArrayList<>();
        for(NetNode<NetNodeInfo> child : _network.getNetwork().dfs()) {
            if(child.getData().getHeight() >= height || child.hasParent(edge.Item2)) continue;
            for(NetNode<NetNodeInfo> par : child.getParents()) {
                if(par.getData().getHeight() < height) continue;
                edges.add(new Tuple<>(child, par));
            }
        }
        return edges;
    }

}
