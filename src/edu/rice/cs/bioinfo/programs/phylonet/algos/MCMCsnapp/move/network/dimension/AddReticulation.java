package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.dimension;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import org.apache.commons.math3.distribution.BetaDistribution;

import java.util.List;

/**
 * Created by wendingqiao on 3/6/16.
 * Add reticulation
 */
public class AddReticulation extends DimensionChange {

    protected NetNode<NetNodeInfo> _v1, _v2, _v3, _v4, _v5, _v6;

    public AddReticulation(UltrametricNetwork net) {
        super(net);
    }

    @Override
    public double propose() {
        _v1 = _v2 = _v3 = _v4 = _v5 = _v6 = null; // reset

        Networks.autoLabelNodes(_network.getNetwork());

        List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> edges = Networks.getAllEdges(_network.getNetwork());
        edges.add(new Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>(_network.getNetwork().getRoot(), null));

        for(Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> tuple : edges) {
            if(!tuple.Item1.isRoot() && !Double.isNaN(tuple.Item1.getRootPopSize())) {
                System.err.println("Not root but has root pop size!");
            }
        }

        int numEdges = edges.size();
        int numRetiNodes = _network.getNetwork().getReticulationCount();

        Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge1, edge2;
        edge1 = edges.get(Randomizer.getRandomInt(numEdges));
        do {
            edge2 = edges.get(Randomizer.getRandomInt(numEdges));
        } while (edge1 == edge2);

        //System.out.println(edge1.Item1.getName() + "-" + (edge1.Item2!=null?edge1.Item2.getName():"") + " " + edge2.Item1.getName() + "-" + (edge2.Item2!=null?edge2.Item2.getName():""));

        _v3 = edge1.Item2;
        _v4 = edge1.Item1;
        _v5 = edge2.Item2;
        _v6 = edge2.Item1;

        double t3 = _v3 != null ? _v3.getData().getHeight() : Utils.ROOT_TIME_UPPER_BOUND;
        double t4 = _v4.getData().getHeight();
        double l1 = t3 - t4;
        double t1 = t4 + Randomizer.getRandomDouble() * l1;

        double t5 = _v5 != null ? _v5.getData().getHeight() : Utils.ROOT_TIME_UPPER_BOUND;
        double t6 = _v6.getData().getHeight();
        double l2 = t5 - t6;
        double t2 = t6 + Randomizer.getRandomDouble() * l2;

        double gamma = 0.0;
        if(Utils._MCMC) {
            gamma = Randomizer.getRandomDouble();
        } else {
            gamma = Randomizer.getRandomDouble();
            //gamma = new BetaDistribution(2.0, 2.0).inverseCumulativeProbability(Randomizer.getRandomDouble());
        }

        _v1 = new BniNetNode<>();
        _v1.setData(new NetNodeInfo(t1));
        _v2 = new BniNetNode<>();
        _v2.setData(new NetNodeInfo(t2));

        double logPopSize = (t1 > t2) ? addReticulation(_v1, _v2, _v3, _v4, _v5, _v6, gamma) // v3 v4 t1 v5 v6 t2
                : addReticulation(_v2, _v1, _v5, _v6, _v3, _v4, gamma);  // v5 v6 t2 v3 v4 t1

        double pda = numRetiNodes == 0 ? 0.5 : 1.0;
        double numRetiEdges = 2 * (numRetiNodes + 1);
        double hr = pda * l1 * l2 * numEdges * (numEdges-1.0) / numRetiEdges;

        //if(Networks.hasTheSameTopology(_network.getNetwork(), Networks.readNetwork("[0.02](((((C:0.01)#H1:0.02::0.6,B:0.03):0.01)#H2:0.02::0.6,A:0.06):0.03,((D:0.02,#H1:0.01::0.4):0.03,#H2:0.01::0.4):0.04);"))) {
        //    System.err.println("Got it!");
        //}

        _violate = false;
        return Math.log(hr) + logPopSize;
    }

    @Override
    public void undo() {
        if(_v1.getData().getHeight() > _v2.getData().getHeight()) {
            removeReticulation(_v1, _v2, _v3, _v4, _v5, _v6);
        } else {
            removeReticulation(_v2, _v1, _v5, _v6, _v3, _v4);
        }
    }

    @Override
    public String getName() {
        return "Add-Reticulation";
    }

}
