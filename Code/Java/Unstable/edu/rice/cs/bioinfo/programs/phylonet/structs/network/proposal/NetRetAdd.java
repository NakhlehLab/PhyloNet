package edu.rice.cs.bioinfo.programs.phylonet.structs.network.proposal;

import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;

import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Random;

/**
 * Created by Dingqiao Wen on 10/18/14.
 */
public class NetRetAdd<T> implements NetworkOperation {

    // try up to [count] times
    private int count = 10;

    private Network _net;

    private Random random;

    private double logQ = 1.0;

    /**
     * one of the parents of network node
     */
    private NetNode<T> source;
    private NetNode<T> srcPar;
    private double srcDist;
    private NetNode<T> parentNode;

    /**
     * Network node
     */
    private NetNode<T> destination;
    private NetNode<T> dstPar;
    private double dstDist;
    private NetNode<T> netNode;

    public NetRetAdd(Network net, Random rand) {
        this._net = net;
        this.random = rand;
    }

    @Override
    public double operate() {
        List<NetNode<T>> treeNodes = IterableHelp.toList(_net.getTreeNodes());
        treeNodes.remove(_net.getRoot());
        int r1 = random.nextInt(treeNodes.size()), r2 = random.nextInt(treeNodes.size());
        // r1 != r2
        if(r1 == r2) {
            count--;
            if(count == 0) {
                logQ = Double.MIN_VALUE;
                return logQ;
            }
            return operate();
        }
        source = treeNodes.get(r1);
        destination = treeNodes.get(r2);
        // destination node cannot be a parent, grand parent... of add node
        // avoid cycle
        Queue<NetNode<T>> queue = new LinkedList<NetNode<T>>();
        queue.add(source);
        while(!queue.isEmpty()) {
            NetNode<T> cur = queue.poll();
            if(cur.equals(destination)) { // end this loop
                count--;
                if(count == 0) {
                    logQ = Double.MIN_VALUE;
                    return logQ;
                }
                return operate();
            }
            for(NetNode<T> p : cur.getParents()) queue.add(p);
        }
        // add parent node between source node and its parent
        srcPar = source.getParents().iterator().next();
        srcDist = source.getParentDistance(srcPar);
        // add parent node
        parentNode = new BniNetNode<T>("", null);
        double d = srcDist * random.nextDouble();
        // source parent node -> parent node
        srcPar.adoptChild(parentNode, d);
        parentNode.setParentDistance(srcPar, d);
        // parent node -> source
        parentNode.adoptChild(source, srcDist - d);
        source.setParentDistance(parentNode, srcDist - d);
        // [remove] source parent node -> source [remove]
        srcPar.removeChild(source);

        // add network node between destination node and its parent
        dstPar = destination.getParents().iterator().next();
        dstDist = destination.getParentDistance(dstPar);
        // add network node
        netNode = new BniNetNode<T>("", null);
        d = dstDist * random.nextDouble();
        // destination parent node -> network node
        dstPar.adoptChild(netNode, d);
        netNode.setParentDistance(dstPar, d);
        netNode.setParentProbability(dstPar, 0.5); // set probability
        // network node -> destination node
        netNode.adoptChild(destination, dstDist - d);
        destination.setParentDistance(netNode, dstDist - d);
        // remove destination parent node - > destination node
        dstPar.removeChild(destination);

        // parent node -> network node
        parentNode.adoptChild(netNode, 1.0);
        netNode.setParentDistance(parentNode, 1.0);
        netNode.setParentProbability(parentNode, 0.5);

        return logQ;
    }

    @Override
    public void undo() {
        // no-op
        if(logQ == Double.MIN_VALUE) return;
        // remove parent node -> network node
        parentNode.removeChild(netNode);
        // remove network node
        dstPar.adoptChild(destination, dstDist);
        destination.setParentDistance(dstPar, dstDist);
        netNode.removeChild(destination);
        dstPar.removeChild(netNode);
        // remove parent node
        srcPar.adoptChild(source, srcDist);
        source.setParentDistance(srcPar, srcDist);
        parentNode.removeChild(source);
        srcPar.removeChild(parentNode);
    }
}
