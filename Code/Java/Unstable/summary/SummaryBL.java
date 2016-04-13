package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.summary;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.*;
import java.util.*;

/**
 * Created by wendingqiao on 2/17/15.
 */
public class SummaryBL {

    class Branch{

        private NetNode parent;
        private NetNode child;
        private double minBL, maxBL, meanBL, stdBL;
        private double minProb, maxProb, meanProb, stdProb;
        private int size;

        public Branch(NetNode p, NetNode c, double dist, double prob) {
            this.parent = p;
            this.child = c;
            this.meanBL = this.maxBL = 0;
            this.minBL = Double.MAX_VALUE;
            this.stdBL = 0;
            this.size = 0;
            if(prob == NetNode.NO_PROBABILITY) {
                this.meanProb = NetNode.NO_PROBABILITY;
                this.stdProb = NetNode.NO_PROBABILITY;
            } else {
                this.meanProb = this.maxProb = 0;
                this.minProb = Double.MAX_VALUE;
                this.stdProb = 0;
            }
        }

        public void addInfo(double dist, double prob) {
            this.meanBL += dist;
            this.stdBL += dist * dist;
            this.minBL = Math.min(this.minBL, dist);
            this.maxBL = Math.max(this.maxBL, dist);

            if(this.meanProb == this.child.NO_PROBABILITY) return;

            this.meanProb += prob;
            this.stdProb += prob * prob;
            this.minProb = Math.min(this.minProb, prob);
            this.maxProb = Math.max(this.maxProb, prob);
            this.size++;
        }

        public int getSize() {return this.size;}

        public String report() {
            //System.out.println(this.size);
            this.meanBL /= (double)size;
            this.stdBL = Math.sqrt(this.stdBL / (double)size - this.meanBL * this.meanBL);
            String s = this.parent.getName() + " -> " + this.child.getName() + ": " + this.meanBL + ", " +
                    this.minBL + ", " + this.maxBL + ", " + this.stdBL;
            if(this.meanProb == this.child.NO_PROBABILITY) return s;
            this.meanProb /= (double) size;
            this.stdProb = Math.sqrt(this.stdProb / (double)size - this.meanProb * this.meanProb);
            return s + ": " + this.meanProb + ", " +
                    this.minProb + ", " + this.maxProb + ", " + this.stdProb;
        }
    }


    private Network _net;
    private Map<String, Branch> _branches;
    private int _size;

    public SummaryBL(String s) {
        this._net = Networks.readNetwork(s);
        this._branches = getBranches(this._net);
        this._size = 0;
    }

    public void addFile(String file) {
        try{
            BufferedReader in = new BufferedReader(new FileReader(file));
            String s;
            String[] ss;
            boolean start = false;
            while((s = in.readLine()) != null) {
                ss = s.split("\\s+");
                if(ss.length == 7 && ss[0].charAt(0) != 'I') {
                    s = in.readLine();
                    if(ss[2].charAt(0) == '1' || start) {
                        start = true;
                        Network net = Networks.readNetwork(s);
                        if(net == null) continue;
                        //System.out.println(net.toString());
                        if(Networks.hasTheSameTopology(this._net, net)) {
                            Map<NetNode, NetNode> map = Networks.mapTwoNetworks(this._net, net);
                            this._size++;
                            addInformation(this._net, map);
                        }
                    }
                }
            }
            in.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void report() {
        System.out.println(_size);
        System.out.println("edge: branch length - mean min max std : inheritance probability - mean min max std");
        for(String key : this._branches.keySet()) {
            System.out.println(this._branches.get(key).report());
        }
    }

    private Map<String, Branch> getBranches(Network net) {
        Map<String, Branch> res = new HashMap<String, Branch>();
        for(Object n1 : Networks.postTraversal(net)) {
            NetNode child = (NetNode) n1;
            for(Object n2 : child.getParents()) {
                NetNode parent = (NetNode) n2;
                String id = parent.getName() + "-" + child.getName();
                if(res.containsKey(id)) continue;
                res.put(id, new Branch(parent, child, child.getParentDistance(parent), child.getParentProbability(parent)));
            }
        }
        return res;
    }

    private void addInformation(Network net, Map<NetNode, NetNode> map) {
        for(Object n1 : Networks.postTraversal(net)) {
            NetNode child = (NetNode) n1;
            for(Object n2 : child.getParents()) {
                NetNode parent = (NetNode) n2;
                String id = parent.getName() + "-" + child.getName();
                Branch br = this._branches.get(id);
                if(br.getSize() == this._size) continue;
                NetNode ch = map.get(child), par = map.get(parent);
                br.addInfo(ch.getParentDistance(par), ch.getParentProbability(par));
            }
        }
    }
}