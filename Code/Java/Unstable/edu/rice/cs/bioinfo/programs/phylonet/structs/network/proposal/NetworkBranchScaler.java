package edu.rice.cs.bioinfo.programs.phylonet.structs.network.proposal;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import java.util.*;

/**
 * Created by Dingqiao Wen on 10/20/14.
 */
public abstract class NetworkBranchScaler<T> implements NetworkOperation {

    protected Network _net;

    protected Random random;

    protected double logQ = 0.0;

    protected List<NetNode<T>> selectedNodes = new ArrayList<NetNode<T>>();

    protected List<Double> originals = new ArrayList<Double>();

    public NetworkBranchScaler(Network net, Random random) {
        this._net = net;
        this.random = random;
    }

    public double operate() {
        selectNodes();
        scaleLengths();
        if(selectedNodes.size() == 0) logQ = Double.MIN_VALUE;
        return logQ;
    }

    public void undo() {
        if(logQ == Double.MIN_VALUE) return;
        int size = selectedNodes.size();
        int idx = originals.size();
        for(int i = size-1; i >= 0; i--) {
            Iterable<NetNode<T>> parents = selectedNodes.get(i).getParents();
            for(NetNode<T> par : parents) {
                setLength(selectedNodes.get(i), par, originals.get(--idx));
            }
        }
        assert(idx == 0);
        //if(idx != 0) System.out.println("wrong");
    }

    protected abstract void setLength(NetNode<T> node, NetNode<T> parent, double dist);

    /**
     * Selects nodes in netNodes to reset branch lengths
     * Save original lengths in originals for undo().
     */
    protected abstract void selectNodes();

    /**
     * Re-scale branch lengths and store them in lengths
     */
    protected abstract void scaleLengths();
}
