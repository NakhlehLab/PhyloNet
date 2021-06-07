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
 * Created by wendingqiao on 3/16/16.
 * Move the tail of a random selected edge to a random source
 */
public class MoveTail extends NetworkOperator {

    private double _logHR;
    private NetNode<NetNodeInfo> _v1, _v2, _v3, _v4, _v5, _v6;
    private double _oldHeight;

    public MoveTail(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        _v1 = _v2 = _v3 = _v4 = _v5 = _v6 = null; // reset

        List<NetNode<NetNodeInfo>> internalTreeNodes = Networks.getInternalTreeNodes(_network.getNetwork());
        do {
            _v1 = internalTreeNodes.get(Randomizer.getRandomInt(internalTreeNodes.size()));
        } while(_v1.isRoot());

        List<NetNode<NetNodeInfo>> children = IterableHelp.toList(_v1.getChildren());
        _v2 = children.get(Randomizer.getRandomInt(children.size()));

        _v3 = _v1.getParents().iterator().next();
        _v4 = Networks.getOtherChild(_v1, _v2);

        if(_v4.hasParent(_v3)) {
            _violate = false;
            _logHR = Utils.INVALID_MOVE;
        } else {
            List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> edges = getEdgesAboveNode(_v2, _v1);
            if(edges.size() == 0) {
                _violate = false;
                _logHR = Utils.INVALID_MOVE;
            } else {
                _oldHeight = _v1.getData().getHeight();

                Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge = edges.get(Randomizer.getRandomInt(edges.size()));
                _v5 = edge.Item2;
                _v6 = edge.Item1;

                double t5 = _v5.getData().getHeight();
                double tl = Math.max(_v2.getData().getHeight(), _v6.getData().getHeight());
                double newHeight = tl + Randomizer.getRandomDouble() * (t5 - tl);

                moveTail(newHeight, _v3, _v4, _v5, _v6);

                double l2 = t5 - tl;
                double l1 = _v3.getData().getHeight() - Math.max(_v2.getData().getHeight(), _v4.getData().getHeight());

                _violate = true;
                _logHR = Math.log(l2 / l1);
            }
        }
        return _logHR;
    }

    @Override
    public void undo() {
        if(_logHR == Utils.INVALID_MOVE) return;
        moveTail(_oldHeight, _v5, _v6, _v3, _v4);
    }

    @Override
    public String getName() {
        return "Move-Tail";
    }

    private void moveTail(double height,
                          NetNode<NetNodeInfo> v3, NetNode<NetNodeInfo> v4,
                          NetNode<NetNodeInfo> v5, NetNode<NetNodeInfo> v6) {

        _v1.getData().setHeight(height);
        _v2.setParentDistance(_v1, height - _v2.getData().getHeight());

        double[] paramV3V1 = getParameters(v3, _v1);
        double[] paramV1V4 = getParameters(_v1, v4);
        double[] paramV5V6 = getParameters(v5, v6);

        v3.removeChild(_v1);
        _v1.removeChild(v4);
        v5.removeChild(v6);

        adopt(v3, v4, paramV1V4);
        adopt(v5, _v1, paramV3V1);
        adopt(_v1, v6, paramV5V6);
    }

    private List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> getEdgesAboveNode(NetNode<NetNodeInfo> v2
            , NetNode<NetNodeInfo> v1 ) {

        List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> list = new ArrayList<>();
        for(NetNode<NetNodeInfo> node : _network.getNetwork().bfs()) {
            for(NetNode<NetNodeInfo> par: node.getParents()) {
                if(par.getData().getHeight() > v2.getData().getHeight()
                        && par != v1 && node != v1 && node != v2) {
                    list.add(new Tuple<>(node, par));
                }
            }
        }
        return list;
    }

}
