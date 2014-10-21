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

    public NetInheritScaler(Network net) {
        super(net);
    }

    @Override
    protected void setLength(NetNode<T> node, NetNode<T> parent, double prob) {
        node.setParentProbability(parent, prob);
    }

    @Override
    protected void selectNodes() {
        List<NetNode<T>> nodes = IterableHelp.toList(_net.getNetworkNodes());
        Random rand = new Random();
        int num = nodes.size();
        while(num < 1 && nodes.size() > 0) {
            num = rand.nextInt(nodes.size());
        }
        while(num > 0) {
            num--;
            selectedNodes.add(nodes.get(rand.nextInt(nodes.size())));
        }
    }

    @Override
    protected void scaleLengths() {
        NormalDistribution normal = new NormalDistribution(0.0, 0.25);
        for(NetNode<T> node : selectedNodes) {
            // network node - two parents
            Iterator<NetNode<T>> iterator = node.getParents().iterator();
            NetNode<T> par1 = iterator.next(), par2 = iterator.next();
            double p1 = node.getParentProbability(par1);
            double p2 = node.getParentProbability(par2);
            originals.add(p1);
            originals.add(p2);
            if((new Random()).nextDouble() < 0.5) {
                double newp = p1 * (1.0 + normal.sample());
                while(newp >= 1.0 || newp <= 0.0){
                    newp = p1 * (1.0 + normal.sample());
                }
                setLength(node, par1, newp);
                setLength(node, par2, 1.0-newp);
            } else {
                double newp = p2 * (1.0 + normal.sample());
                while(newp >= 1.0 || newp <= 0.0){
                    newp = p2 * (1.0 + normal.sample());
                }
                setLength(node, par1, 1.0-newp);
                setLength(node, par2, newp);
            }

        }
    }
}
