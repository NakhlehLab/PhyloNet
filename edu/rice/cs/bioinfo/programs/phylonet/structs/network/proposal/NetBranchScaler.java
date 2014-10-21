package edu.rice.cs.bioinfo.programs.phylonet.structs.network.proposal;

import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.List;
import java.util.Random;

/**
 * Created by Dingqiao Wen on 10/18/14.
 * Re-scale the distance between several tree nodes to their parents
 * Randomly choose a number from uniform distribution of tree nodes size.
 */
public class NetBranchScaler<T> extends NetworkBranchScaler<T> {

    public NetBranchScaler(Network net) {
        super(net);
    }

    @Override
    protected void setLength(NetNode<T> node, NetNode<T> parent, double dist) {
        node.setParentDistance(parent, dist);
        parent.adoptChild(node, dist);
    }

    @Override
    protected void selectNodes() {
        List<NetNode<T>> nodes = IterableHelp.toList(_net.getTreeNodes());
        nodes.remove(_net.getRoot());
        Random rand = new Random();
        int num = rand.nextInt(nodes.size());
        System.out.println(num);
        while(num > 0) {
            num--;
            selectedNodes.add(nodes.get(rand.nextInt(nodes.size())));
        }
    }

    @Override
    protected void scaleLengths() {
        NormalDistribution normal = new NormalDistribution(0.0, 0.25);
        for(NetNode<T> node : selectedNodes) {
            // tree node - only one parent
            NetNode<T> par = node.getParents().iterator().next();
            double dist = node.getParentDistance(par);
            originals.add(dist);
            double scale = dist * (1.0 + normal.sample());
            while(scale <= 0.0) scale = dist * (1.0 + normal.sample());
            setLength(node, par, scale);
        }
    }
}
