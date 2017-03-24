package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.summary;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by wendingqiao on 2/17/15.
 */
public class SummaryBranch {

    private class Info{
        private int count;
        private double avgPost;
        private double stdPost;

        public Info(double post) {
            count = 1;
            avgPost = post;
            stdPost = post * post;
        }

        public void add(double post) {
            count++;
            avgPost += post;
            stdPost += post * post;
        }

        public String toString() {
            double avg = avgPost / count;
            double std = Math.sqrt(stdPost / count - avg * avg);
            return count + "    " + avg + "  " + std;
        }
    }

    private class Branch{

        private NetNode parent;
        private NetNode child;
        private double minBL, maxBL, meanBL, stdBL;
        private double minProb, maxProb, meanProb, stdProb;
        private double minSp, maxSp, meanSp, stdSp;
        private int size;

        public Branch(NetNode p, NetNode c, double dist, double prob, double support) {
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
                this.meanProb = this.maxProb = this.stdProb = 0;
                this.minProb = Double.MAX_VALUE;
            }
            if(support == NetNode.NO_SUPPORT) {
                this.meanSp = NetNode.NO_SUPPORT;
                this.stdSp = NetNode.NO_SUPPORT;
            } else {
                this.meanSp = this.maxSp = this.stdSp = 0;
                this.minSp = Double.MAX_VALUE;
            }
        }

        public void addInfo(double dist, double prob, double support) {
            this.size++;
            this.meanBL += dist;
            this.stdBL += dist * dist;
            this.minBL = Math.min(this.minBL, dist);
            this.maxBL = Math.max(this.maxBL, dist);
            if(this.meanProb != this.child.NO_PROBABILITY) {
                this.meanProb += prob;
                this.stdProb += prob * prob;
                this.minProb = Math.min(this.minProb, prob);
                this.maxProb = Math.max(this.maxProb, prob);
            }
            if(this.meanSp != this.child.NO_SUPPORT) {
                this.meanSp += support;
                this.stdSp += support * support;
                this.minSp = Math.min(this.minSp, support);
                this.maxSp = Math.max(this.maxSp, support);
            }
        }

        public int getSize() {return this.size;}

        public String report() {
            return this.report(1.0);
        }

        public String report(double scale) {
            return this.report(scale, 1.0);
        }

        public String report(double scale, double popScale) {
            this.meanBL /= (double)size;
            this.stdBL = Math.sqrt(this.stdBL / (double)size - this.meanBL * this.meanBL);
            String s = this.parent.getName() + " -> " + this.child.getName() + ": " + (this.meanBL/scale) + ", "
                    + (this.stdBL/scale);
            if(this.meanProb == this.child.NO_PROBABILITY) {
                s += " : NO_PROBABILITY";
            } else {
                this.meanProb /= (double) size;
                this.stdProb = Math.sqrt(this.stdProb / (double)size - this.meanProb * this.meanProb);
                s += " : " + this.meanProb + ", "
                        + this.stdProb;
            }
            if(this.meanSp == this.child.NO_SUPPORT) {
                s += " : NO_POP_SIZE";
            } else {
                this.meanSp /= (double) size;
                this.stdSp = Math.sqrt(this.stdSp / (double)size - this.meanSp * this.meanSp);
                s += " : " + (this.meanSp / popScale) + ", "
                        + (this.stdSp / popScale);
            }
            return s;

        }
    }

    private Network _net;
    private Map<String, Branch> _branches;
    private int _size;
    private double _rootPopSizeAvg = 0;
    private double _rootPopSizeStd = 0;
    private int _rootPopSizeSize = 0;
    private Map<Network, Info> _topologyCount = new HashMap<>();

    public SummaryBranch(String s) {
        this._net = Networks.readNetwork(s);
        Networks.autoLabelNodes(this._net);
        this._branches = getBranches(this._net);
        this._size = 0;
    }

    public void addFile(String file) {
        this.addFile(file, -1, Integer.MAX_VALUE);
    }

    public void addFile(String file, int startIter, int endIter) {
        try{
            BufferedReader in = new BufferedReader(new FileReader(file));
            String s;
            String[] ss;
            boolean start = false;
            while((s = in.readLine()) != null) {
                ss = s.trim().split("\\s+");
                if(ss.length != 7 || ss[0].startsWith("I") || !ss[0].endsWith(";")) continue;
                s = in.readLine();
                int iter = Integer.parseInt(ss[0].substring(0, ss[0].length() - 1));
                if(iter < startIter) {
                    continue;
                }
                if(iter > endIter) {
                    break;
                }
                double posterior = Double.parseDouble(ss[3].substring(0, ss[1].length() - 1));
                double ess = Double.parseDouble(ss[2].substring(0, ss[2].length() - 1));
                if(ess <= 0 && !start) continue;
                start = true;
                Network net = Networks.readNetwork(s);
                if(net == null) {
                    continue;
                }
                if(Networks.hasTheSameTopology(this._net, net)) {
                    Map<NetNode, NetNode> map = Networks.mapTwoNetworks(this._net, net);
                    this._size++;
                    addInformation(this._net, map);
                    if(s.startsWith("[")) {
                        double popSize = Double.parseDouble(s.substring(1, s.indexOf("]")));
                        _rootPopSizeAvg += popSize;
                        _rootPopSizeStd += popSize * popSize;
                        s = s.substring(s.indexOf("]") + 1);
                        _rootPopSizeSize++;
                    }
                } else {
                    boolean add = false;
                    for(Network key : _topologyCount.keySet()) {
                        if(Networks.hasTheSameTopology(key, net)) {
                            _topologyCount.get(key).add(posterior);
                            add = true;
                            break;
                        }
                    }
                    if(!add) _topologyCount.put(net, new Info(posterior));
                }
            }
            in.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void report(double scale, double popScale) {
        System.out.println(_size);
        _rootPopSizeAvg /= _rootPopSizeSize;
        _rootPopSizeStd = Math.sqrt(_rootPopSizeStd / _rootPopSizeSize - _rootPopSizeAvg * _rootPopSizeAvg);
        System.out.printf("sample_size = %d, root_pop_size_mean = %f, root_pop_size_std = %f\n",
                _rootPopSizeSize, (_rootPopSizeAvg / popScale), (_rootPopSizeStd / popScale));
        System.out.println("edge: tau_mean, tau_std, gamma_mean, gamma_std");
        for(String key : this._branches.keySet()) {
            Branch br = this._branches.get(key);
            System.out.println(br.report(scale, popScale));
            br.child.setParentDistance(br.parent, br.meanBL);
            br.child.setParentSupport(br.parent, br.child.NO_SUPPORT);
        }
        System.out.println(this._net.toString());
    }

    private Map<String, Branch> getBranches(Network net) {
        Map<String, Branch> res = new HashMap<String, Branch>();
        for(Object n1 : Networks.postTraversal(net)) {
            NetNode child = (NetNode) n1;
            for(Object n2 : child.getParents()) {
                NetNode parent = (NetNode) n2;
                String id = parent.getName() + "-" + child.getName();
                if(res.containsKey(id)) continue;
                res.put(id, new Branch(parent, child, child.getParentDistance(parent),
                        child.getParentProbability(parent), child.getParentSupport(parent)));
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
                br.addInfo(
                        ch.getParentDistance(par), ch.getParentProbability(par), ch.getParentSupport(par)
                );
            }
        }
    }

}