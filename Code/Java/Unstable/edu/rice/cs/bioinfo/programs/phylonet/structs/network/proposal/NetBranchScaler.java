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

    private double windowSize = 0.1;

    public NetBranchScaler(Network net, Random rand) {
        super(net, rand);
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
        int num = random.nextInt(nodes.size());
        while(num > 0) {
            num--;
            selectedNodes.add(nodes.get(random.nextInt(nodes.size())));
        }
    }

    @Override
    protected void scaleLengths() {
        for(NetNode<T> node : selectedNodes) {
            // tree node - only one parent
            NetNode<T> par = node.getParents().iterator().next();
            double dist = node.getParentDistance(par);
            originals.add(dist);
            double scale = windowSize * 2.0 * (random.nextDouble() - 0.50) + dist;
            while(scale <= 0.0) scale = windowSize * 2.0 * (random.nextDouble() - 0.50) + dist;
            setLength(node, par, scale);
        }
    }
}
