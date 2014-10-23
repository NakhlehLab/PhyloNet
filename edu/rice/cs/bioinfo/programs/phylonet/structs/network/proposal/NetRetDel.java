package edu.rice.cs.bioinfo.programs.phylonet.structs.network.proposal;

import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Created by Dingqiao Wen on 10/18/14.
 */
public class NetRetDel<T> implements NetworkOperation {

    // try up to [count] times
    private int count = 10;

    private Network _net;

    private Random random;

    private double logQ = 1.0;

    /**
     * store deleted reticulation nodes
     */
    private NetNode<T> deleted;
    private List<NetNode<T>> children;
    private List<Double> childrenProbs = new ArrayList<Double>();

    /**
     * store one of its two parents
     */
    private NetNode<T> deletedParent;
    private double prob1;
    private double dist1;
    private NetNode<T> dpp;
    private double distp;
    private List<NetNode<T>> dPchildren;
    private List<Double> dPCProbs = new ArrayList<Double>();

    /**
     * One of its two parents and the only parent after operation
     */
    private NetNode<T> parent;
    private double prob2;
    private double dist2;


    public NetRetDel(Network net, Random rand) {
        this._net = net;
        this.random = rand;
    }

    @Override
    public double operate() {
        List<NetNode<T>> netNodes = IterableHelp.toList(_net.getNetworkNodes());
        if(netNodes.size() == 0) {
            logQ = Double.MIN_VALUE;
            return logQ;
        }
        deleted = netNodes.get(random.nextInt(netNodes.size()));

        List<NetNode<T>> parents = IterableHelp.toList(deleted.getParents());
        // delete parent node cannot be network node
        // cannot delete two network nodes at once
        if(parents.get(0).isNetworkNode() && parents.get(1).isNetworkNode()){
            count--;
            if(count == 0) {
                logQ = Double.MIN_VALUE;
                return logQ;
            }
            return operate();
        } else if(parents.get(0).isNetworkNode()) {
            deletedParent = parents.get(1);
            parent = parents.get(0);
        } else if(parents.get(1).isNetworkNode()) {
            deletedParent = parents.get(0);
            parent = parents.get(1);
        } else {
            if(random.nextDouble() < 0.50) {
                deletedParent = parents.get(0);
                parent = parents.get(1);
            } else {
                deletedParent = parents.get(1);
                parent = parents.get(0);
            }
        }

        // constraint for delete parent: its parent should not be a parent of its children
        if(!deletedParent.getParents().iterator().hasNext()) {
            count--;
            if(count == 0) {
                logQ = Double.MIN_VALUE;
                return logQ;
            }
            return operate();
        }
        dpp = deletedParent.getParents().iterator().next();
        for(NetNode<T> pc : deletedParent.getChildren()) {
            if(pc.hasParent(dpp)){
                count--;
                if(count == 0) {
                    logQ = Double.MIN_VALUE;
                    return logQ;
                }
                return operate();
            }
        }
        // 1:: store values between deleted parent - deleted node
        prob1 = deleted.getParentProbability(deletedParent);
        dist1 = deleted.getParentDistance(deletedParent);
        deletedParent.removeChild(deleted);
        // store values between another parent - deleted node
        prob2 = deleted.getParentProbability(parent);
        dist2 = deleted.getParentDistance(parent);
        parent.removeChild(deleted);

        // 2:: remove deleted node from the network
        children = IterableHelp.toList(deleted.getChildren());
        for(NetNode<T> child : children) {
            // 2.1:: add all its children under parent
            double d = child.getParentDistance(deleted) + dist2;
            parent.adoptChild(child, d);
            child.setParentDistance(parent, d);
            if(child.isNetworkNode()) {
                double prob = child.getParentProbability(deleted);
                childrenProbs.add(prob);
                child.setParentProbability(parent, prob);
            }
            // 2.2:: remove all its children from deleted
            deleted.removeChild(child);
        }

        // 3:: remove deleted parent from the network
        dPchildren = IterableHelp.toList(deletedParent.getChildren());
        distp = deletedParent.getParentDistance(dpp);
        for(NetNode<T> child : dPchildren) {
            double d = child.getParentDistance(deletedParent) + distp;
            dpp.adoptChild(child, d);
            child.setParentDistance(dpp, d);
            if(child.isNetworkNode()) {
                double prob = child.getParentProbability(deletedParent);
                dPCProbs.add(prob);
                child.setParentProbability(dpp, prob);
            }
            deletedParent.removeChild(child);
        }
        dpp.removeChild(deletedParent);
        return logQ;
    }

    @Override
    public void undo() {
        // no-op
        if(logQ == Double.MIN_VALUE) return;
        // 1:: restore deleted parent
        dpp.adoptChild(deletedParent, distp);
        deletedParent.setParentDistance(dpp, distp);
        int idx = 0;
        for(NetNode<T> child : dPchildren) {
            double d = child.getParentDistance(dpp) - distp;
            deletedParent.adoptChild(child, d);
            child.setParentDistance(deletedParent, d);
            if(child.isNetworkNode()) {
                child.setParentProbability(deletedParent, dPCProbs.get(idx++));
            }
            dpp.removeChild(child);
        }
        assert (idx == dPCProbs.size());
        if(idx != dPCProbs.size())
            throw new IllegalArgumentException("Network node size != probability size!!!");

        // 2:: restore deleted node
        deletedParent.adoptChild(deleted, dist1);
        deleted.setParentDistance(deletedParent, dist1);
        deleted.setParentProbability(deletedParent, prob1);

        parent.adoptChild(deleted, dist2);
        deleted.setParentDistance(parent, dist2);
        deleted.setParentProbability(parent, prob2);
        // 3:: restore deleted node's children
        idx = 0;
        for(NetNode<T> child : children) {
            double d = child.getParentDistance(parent) - dist2;
            deleted.adoptChild(child, d);
            child.setParentDistance(deleted, d);
            if(child.isNetworkNode()) {
                child.setParentProbability(deleted, childrenProbs.get(idx++));
            }
            parent.removeChild(child);
        }
        assert (idx == childrenProbs.size());
        if(idx != childrenProbs.size())
            throw new IllegalArgumentException("Network node size != probability size!!!");
    }
}
