package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.network.dimension;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.List;

/**
 * Created by wendingqiao on 3/6/16.
 * Add reticulation
 */
public class AddReticulationToBackbone extends AddReticulation {

    Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> retiEdge = null;

    public AddReticulationToBackbone(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        double logHasting = super.propose();

        if(logHasting != Utils.INVALID_MOVE) {
            logHasting += Math.log(2.0);
            List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> retiEdges = _network.getRetiEdges();
            if (_v1.getData().getHeight() > _v2.getData().getHeight()) {
                retiEdge = new Tuple<>(_v2, _v1);
            } else {
                retiEdge = new Tuple<>(_v1, _v2);
            }
            retiEdges.add(retiEdge);
        }

        return logHasting;
    }

    @Override
    public void undo() {
        super.undo();

        List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> retiEdges = _network.getRetiEdges();
        if(retiEdge != null) {
            retiEdges.remove(retiEdge);
            retiEdge = null;
        }
    }

    @Override
    public String getName() {
        return "Add-Reticulation-To-Backbone";
    }

}
