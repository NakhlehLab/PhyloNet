package edu.rice.cs.bioinfo.programs.phylonet.structs.network.proposal;

import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.Iterator;
import java.util.List;
import java.util.Random;

/**
 * Created by Dingqiao Wen on 10/18/14.
 * Re-scale the interitance probability of several network nodes to their parents
 */
public class NetInheritScaler<T> extends NetworkBranchScaler<T> {

    private double windowSize = 0.1;

    public NetInheritScaler(Network net, Random rand) {
        super(net, rand);
    }

    @Override
    protected void setLength(NetNode<T> node, NetNode<T> parent, double prob) {
        node.setParentProbability(parent, prob);
    }

    @Override
    protected void selectNodes() {
        List<NetNode<T>> nodes = IterableHelp.toList(_net.getNetworkNodes());
        int num = nodes.size();
        while(num < 1 && nodes.size() > 0) {
            num = random.nextInt(nodes.size());
        }
        while(num > 0) {
            num--;
            selectedNodes.add(nodes.get(random.nextInt(nodes.size())));
        }
    }

    @Override
    protected void scaleLengths() {
        for(NetNode<T> node : selectedNodes) {
            // network node - two parents
            Iterator<NetNode<T>> iterator = node.getParents().iterator();
            NetNode<T> par1 = iterator.next(), par2 = iterator.next();
            double p1 = node.getParentProbability(par1);
            double p2 = node.getParentProbability(par2);
            originals.add(p1);
            originals.add(p2);
            double p3 =  windowSize * 2.0 * (random.nextDouble() - 0.50) + p1;
            while(p3 <= 0.0 || p3 >= 1.0) {
                p3 =  windowSize * 2.0 * (random.nextDouble() - 0.50) + p1;
            }
            setLength(node, par1, p3);
            setLength(node, par2, 1.0-p3);
        }
    }
}
