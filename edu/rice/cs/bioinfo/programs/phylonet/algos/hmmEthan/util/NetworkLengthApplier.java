package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.util;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

public class NetworkLengthApplier
{
    final private Network net;
    final double[] networkLengths;

    private int netLengthIndex = 0;

    public NetworkLengthApplier(Network net, double[] networkLengths)
    {
        this.net =net;
        this.networkLengths = networkLengths;
    }

    @SuppressWarnings("unchecked")
    private static <T> Network<T> uglyCast(Network net)
    {
        return net;
    }

    private void applyNetNodeLengths(NetNode<Double> net, double[] lengths)
    {
        if (net.getData() != null)
            return;

        if (net.isLeaf())
        {
            net.setData(0.0);
        }
        else
        {
            double maxChildHeight = 0;
            for (NetNode<Double> child : net.getChildren())
            {
                if (child.getData() == null)
                    applyNetNodeLengths(child, lengths);

                double childHeight = child.getData();
                maxChildHeight = Math.max(maxChildHeight, childHeight);
            }

            double myHeight = maxChildHeight + lengths[netLengthIndex++];
            net.setData(myHeight);

            for (NetNode<Double> child : net.getChildren())
            {
                child.setParentDistance(net, myHeight - child.getData());
            }

        }

    }


    public void apply()
    {
        netLengthIndex = 0;
        Network<Double> foo = uglyCast(net);

        for (NetNode<Double> node:  foo.bfs())
        {
            node.setData(null);
        }

        for (NetNode<Double> node : foo.bfs())
        {
            applyNetNodeLengths(node, networkLengths);
        }

        if (netLengthIndex != networkLengths.length)
            throw new RuntimeException("FAIL adding netlengths");

    }
}
