package edu.rice.cs.bioinfo.programs.phylonet.structs.network.proposal;

import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Random;

/**
 * Created by Dingqiao Wen on 10/18/14.
 * Swap two sub-net
 */
public class NetworkSwap<T> implements NetworkOperation {

    private Network net;

    private Random random;

    private double logQ = 1.0;

    private NetNode<T> n1;

    private NetNode<T> n2;

    public NetworkSwap(Network net, Random rand) {
        this.net = net;
        this.random = rand;
    }

    /**
     * Swap two tree nodes with no direct relationship (not siblings, not child/ grandchild)
     * @return  hastings ratio
     */
    public double operate() {
        op();
        return logQ;
    }

    private void op() {
        List<NetNode<T>> nodes = IterableHelp.toList(net.getTreeNodes());
        nodes.remove(net.getRoot());
        int count = 10; // try up to 10 times
        while(count > 0) {
            count--;
            int r1 = random.nextInt(nodes.size());
            int r2 = random.nextInt(nodes.size());
            // n1 != n2
            if(r1 == r2) continue;
            n1 = nodes.get(r1);
            n2 = nodes.get(r2);
            //System.out.println(r1 + "  " + r2 + "  " + n1.getName() + "  " + n2.getName());
            NetNode<T> n1p = n1.getParents().iterator().next(),
                    n2p = n2.getParents().iterator().next();
            // n1p != n2p
            if(n1p.equals(n2p)) continue;
            // n1 != n2..p
            boolean valid = true;
            Queue<NetNode<T>> queue = new LinkedList<NetNode<T>>();
            queue.add(n2p);
            while(!queue.isEmpty()) {
                NetNode<T> cur = queue.poll();
                if(n1.equals(cur)) {
                    valid = false;
                    break;
                }
                for(NetNode<T> n : cur.getParents()){
                    queue.add(n);
                }
            }
            if(!valid) continue;
            // n2 != n1..p
            queue.add(n1p);
            while(!queue.isEmpty()) {
                NetNode<T> cur = queue.poll();
                if(n2.equals(cur)) {
                    valid = false;
                    break;
                }
                for(NetNode<T> n : cur.getParents()){
                    queue.add(n);
                }
            }
            if(!valid) continue;
            swap();
            return;
        }
        logQ = Double.MIN_VALUE;
    }

    public void undo() {
        // no-op
        if(logQ == Double.MIN_VALUE) return;
        swap();
    }

    private void swap() {
        NetNode<T> n1p = n1.getParents().iterator().next(),
                n2p = n2.getParents().iterator().next();
        double n1dist = n1.getParentDistance(n1p);
        double n2dist = n2.getParentDistance(n2p);
        n1p.removeChild(n1);
        n2p.removeChild(n2);
        n1p.adoptChild(n2, n2dist);
        n2p.adoptChild(n1, n1dist);
    }
}
