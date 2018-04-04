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
 * Move the head of a random selected edge from to a random selected destination.
 */
public class MoveHead extends NetworkOperator {

    private NetNode<NetNodeInfo> _v1, _v2, _v3, _v4, _v5, _v6;
    private double _oldHeight;
    private double _logHR;

    public MoveHead(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        _v1 = _v2 = _v3 = _v4 = _v5 = _v6 = null; // reset

        List<NetNode<NetNodeInfo>> netNodes = IterableHelp.toList(_network.getNetwork().getNetworkNodes());
        _v2 = netNodes.get(Randomizer.getRandomInt(netNodes.size()));
        List<NetNode<NetNodeInfo>> parents = IterableHelp.toList(_v2.getParents());
        _v1 = parents.get(Randomizer.getRandomInt(parents.size()));

        _v3 = Networks.getOtherParent(_v2, _v1);
        _v4 = _v2.getChildren().iterator().next();

        if(_v4.hasParent(_v3)) {
            _violate = false;
            _logHR = Utils.INVALID_MOVE;
        } else {
            List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> edges = getEdgesUnderNode(_v1, _v2);
            Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge = edges.get(Randomizer.getRandomInt(edges.size()));
            _v5 = edge.Item2;
            _v6 = edge.Item1;

            double t6 = _v6.getData().getHeight();
            double th = Math.min(_v1.getData().getHeight(), _v5.getData().getHeight());
            double newHeight = t6 + Randomizer.getRandomDouble() * (th - t6);
            _oldHeight = _v2.getData().getHeight();

            moveHead(newHeight, _v3, _v4, _v5, _v6);

            double t4 = _v4.getData().getHeight();
            double thPrime = Math.min(_v1.getData().getHeight(), _v3.getData().getHeight());

            _violate = true;
            _logHR = Math.log(th - t6) - Math.log(thPrime - t4);
        }
        return _logHR;
    }

    @Override
    public void undo() {
        if(_logHR == Utils.INVALID_MOVE) return;
        moveHead(_oldHeight, _v5, _v6, _v3, _v4);
    }

    @Override
    public String getName() {
        return "Move-Head";
    }

    private void moveHead(double height,
                          NetNode<NetNodeInfo> v3, NetNode<NetNodeInfo> v4,
                          NetNode<NetNodeInfo> v5, NetNode<NetNodeInfo> v6) {

        _v2.getData().setHeight(height);
        _v2.setParentDistance(_v1, _v1.getData().getHeight() - height);

        double[] paramV2V4 = getParameters(_v2, v4);
        double[] paramV3V2 = getParameters(v3, _v2);
        double[] paramV5V6 = getParameters(v5, v6);

        v3.removeChild(_v2);
        _v2.removeChild(v4);
        v5.removeChild(v6);

        adopt(v3, v4, paramV2V4);
        adopt(v5, _v2, paramV3V2);
        adopt(_v2, v6, paramV5V6);
    }

    private List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> getEdgesUnderNode(NetNode<NetNodeInfo> v1
            , NetNode<NetNodeInfo> v2 ) {

        List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> list = new ArrayList<>();
        for(NetNode<NetNodeInfo> node : _network.getNetwork().bfs()) {
            for(NetNode<NetNodeInfo> par: node.getParents()) {
                if(node.getData().getHeight() < v1.getData().getHeight()
                        && node != v2 && par != v2 && par != v1) {
                    list.add(new Tuple<>(node, par));
                }
            }
        }
        return list;
    }


}
