package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.summary;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 1/19/17
 * Time: 6:13 PM
 * To change this template use File | Settings | File Templates.
 */
import edu.rice.cs.bioinfo.programs.phylonet.algos.clustering.ParentalTreeOperation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class SummaryBL {

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
        private List<Double> listBL;
        private List<Double> listProb;
        private List<Double> listSp;

        public Branch(NetNode p, NetNode c, double dist, double prob, double support) {
            this.listBL = new ArrayList<>();
            this.listBL.add(dist);
            this.listProb = new ArrayList<>();
            this.listSp = new ArrayList<>();
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
                this.listProb.add(prob);
            }
            if(support == NetNode.NO_SUPPORT) {
                this.meanSp = NetNode.NO_SUPPORT;
                this.stdSp = NetNode.NO_SUPPORT;
            } else {
                this.meanSp = this.maxSp = this.stdSp = 0;
                this.minSp = Double.MAX_VALUE;
                this.listSp.add(support);
            }
        }

        public void addInfo(double dist, double prob, double support) {
            this.size++;
            this.listBL.add(dist);
            this.meanBL += dist;
            this.stdBL += dist * dist;
            this.minBL = Math.min(this.minBL, dist);
            this.maxBL = Math.max(this.maxBL, dist);
            if(this.meanProb != this.child.NO_PROBABILITY) {
                this.listProb.add(prob);
                this.meanProb += prob;
                this.stdProb += prob * prob;
                this.minProb = Math.min(this.minProb, prob);
                this.maxProb = Math.max(this.maxProb, prob);
            }
            if(this.meanSp != this.child.NO_SUPPORT) {
                this.listSp.add(support);
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
//                   + this.minBL + ", " + this.maxBL + ", "
                    + (this.stdBL/scale);
            if(this.meanProb == this.child.NO_PROBABILITY) {
                s += " : NO_PROBABILITY";
            } else {
                this.meanProb /= (double) size;
                this.stdProb = Math.sqrt(this.stdProb / (double)size - this.meanProb * this.meanProb);
                s += " : " + this.meanProb + ", "
//                        + this.minProb + ", " + this.maxProb + ", "
                        + this.stdProb;
            }
            if(this.meanSp == this.child.NO_SUPPORT) {
                s += " : NO_POP_SIZE";
            } else {
                this.meanSp /= (double) size;
                this.stdSp = Math.sqrt(this.stdSp / (double)size - this.meanSp * this.meanSp);
                s += " : " + (this.meanSp / popScale) + ", "
//                        + this.minSp + ", " + this.maxSp + ", "
                        + (this.stdSp / popScale);
            }
            return s;

        }

        public void reportForDensityPlot(PrintWriter out, double scale, double popScale) {
            out.println(this.parent.getName() + " -> " + this.child.getName() + " : branch length");
            for(Double bl : listBL) {
                out.println(bl / scale);
            }
            if(this.listProb.size() > 0) {
                out.println(this.parent.getName() + " -> " + this.child.getName() + " : inheritance prob");
                for(Double prob : listProb ) {
                    out.println(prob);
                }
            }
            if(this.listSp.size() > 0) {
                out.println(this.parent.getName() + " -> " + this.child.getName() + " : pop size");
                for(Double sp : listSp ) {
                    out.println(sp / popScale);
                }
            }
        }
    }

    private Network _net;
    private Map<String, Branch> _branches;
    private int _size;
    private double _rootPopSizeAvg = 0;
    private double _rootPopSizeStd = 0;
    private int _rootPopSizeSize = 0;
    private Map<Network, Info> _topologyCount = new HashMap<>();

    public SummaryBL(String s) {
        this._net = Networks.readNetwork(s);
        Networks.autoLabelNodes(this._net);
        this._branches = getBranches(this._net);
        this._size = 0;
    }

    public void addFile(String file, boolean notBeast) {
        this.addFile(file, notBeast, -1, Integer.MAX_VALUE);
    }

    public void addFile(String file, boolean notBeast, int startIter, int endIter) {
        try{
            BufferedReader in = new BufferedReader(new FileReader(file));
            String s;
            String[] ss;
            boolean start = false;
            String currentChain = "";
            boolean loggerStart = false;
            while((s = in.readLine()) != null) {
                ss = s.split("\\s+");
                if (notBeast) {
                    if(s.contains("Logger")) {
                        loggerStart = true;
                        continue;
                    } else if(s.contains("Summarization"))  {
                        loggerStart = false;
                        continue;
                    }
                    else if(!loggerStart)
                        continue;

                    if(s.startsWith("Temp")) {
                        currentChain = s;
                        continue;
                    } else if(!currentChain.contains("main")){
                        continue;
                    }
                    if(ss.length == 7 && ss[0].charAt(0) != 'I') {
                        s = in.readLine();
                        int iter = Integer.parseInt(ss[0].substring(0, ss[0].length() - 1));
                        if(iter < startIter) {
                            continue;
                        }
                        if(iter > endIter) {
                            break;
                        }
                        double posterior = Double.parseDouble(ss[3].substring(0, ss[1].length() - 1));
                        double ess = 1; //Double.parseDouble(ss[2].substring(0, ss[2].length() - 1));
                        if(ess > 0 || start) {
                            double popSize = Double.NaN;
                            start = true;
//                            if(s.startsWith("[")) {
//                                popSize = Double.parseDouble(s.substring(1, s.indexOf("]")));
//                                _rootPopSizeAvg += popSize;
//                                _rootPopSizeStd += popSize * popSize;
//                                s = s.substring(s.indexOf("]") + 1);
//                                _rootPopSizeSize++;
//                            }
                            Network net = Networks.readNetwork(s);
                            if(net == null) {
                                continue;
                            }
                            if(Networks.hasTheSameTopology(this._net, net)) {
                                Map<NetNode, NetNode> map = Networks.mapTwoNetworks(this._net, net);
                                this._size++;
                                addInformation(this._net, map);
                                if(s.startsWith("[")) {
                                    popSize = Double.parseDouble(s.substring(1, s.indexOf("]")));
                                    _rootPopSizeAvg += popSize;
                                    _rootPopSizeStd += popSize * popSize;
                                    s = s.substring(s.indexOf("]") + 1);
                                    _rootPopSizeSize++;
                                }
                            } else {
                                boolean add = false;
                                for(Network key : _topologyCount.keySet()) {
                                    if(Networks.hasTheSameTopology(key, net)) {
//                                        _topologyCount.get(key).add(popSize);
                                        _topologyCount.get(key).add(posterior);
                                        add = true;
                                        break;
                                    }
                                }
//                                if(!add) _topologyCount.put(net, new Info(popSize));
                                if(!add) _topologyCount.put(net, new Info(posterior));
                            }
                        }
                    }
                } else {
                    if (ss.length == 4 && ss[0].compareTo("tree") == 0 && Integer.parseInt(ss[1].split("_")[1]) > -1) {

                        String REGEX = "\\[&dmv=[0-9]+.[0-9]+\\]";
                        String INPUT = ss[3];
                        String REPLACE = "";
                        Pattern p = Pattern.compile(REGEX);
                        Matcher m = p.matcher(INPUT);
                        StringBuffer sb = new StringBuffer();
                        while(m.find()){
                            m.appendReplacement(sb, REPLACE);
                        }
                        m.appendTail(sb);
//                        System.out.println(sb.toString());
//                        String tree = sb.toString().replace("1:", "Ash:").replace("2:", "Asp:").replace("3:", "At:")
//                                .replace("4:", "TaA:").replace("5:", "TaB:").replace("6:", "TaD:")
//                                .replace("7:", "Tm:").replace("8:", "Tu:");
                        String tree = sb.toString()
                                .replace("1:", "A:").replace("2:", "B:").replace("3:", "C:");
//                                .replace("1:", "A:").replace("2:", "C:").replace("3:", "G:").replace("4:", "L:").replace("5:", "Q:").replace("6:", "R:");
//                                .replace("1:", "Sbay:").replace("2:", "Scas:").replace("3:", "Scer:").replace("4:", "Sklu:").replace("5:", "Skud:").replace("6:", "Smik:").replace("7:", "Spar:");

//                        System.out.println(tree);
                        Network net = Networks.readNetwork(tree);
                        if(net == null) continue;
                        if(Networks.hasTheSameTopology(this._net, net)) {
                            Map<NetNode, NetNode> map = Networks.mapTwoNetworks(this._net, net);
                            this._size++;
                            addInformation(this._net, map);
                        } else {
                            boolean add = false;
                            for(Network key : _topologyCount.keySet()) {
                                if(Networks.hasTheSameTopology(key, net)) {
                                    _topologyCount.get(key).add(1);
                                    add = true;
                                    break;
                                }
                            }
                            if(!add) _topologyCount.put(net, new Info(1));
                        }
                    }
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
        System.out.println((_rootPopSizeAvg / popScale) + " " + (_rootPopSizeStd / popScale) + " " + _rootPopSizeSize);
        System.out.println("edge: branch length - mean min max std : inheritance probability - mean min max std");
        for(String key : this._branches.keySet()) {
            Branch br = this._branches.get(key);
            System.out.println(br.report(scale, popScale));
            br.child.setParentDistance(br.parent, br.meanBL);
            br.child.setParentSupport(br.parent, br.child.NO_SUPPORT);
        }
        System.out.println(this._net.toString());

        for(Network key : _topologyCount.keySet()) {
            if(_topologyCount.get(key).count < _rootPopSizeSize * 0.005) {
                continue;
            }
            System.out.println(_topologyCount.get(key).toString() + "  " + key.toString());
        }
    }

    public void reportForDensityPlot(String filename, double scale, double popScale) {
        System.out.println(_size);
        try {
            PrintWriter out = new PrintWriter(filename);
            for(String key : this._branches.keySet()) {
                Branch br = this._branches.get(key);
                br.reportForDensityPlot(out, scale, popScale);
            }
            out.println("->");
            out.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

        for(Network key : _topologyCount.keySet()) {
            if(_topologyCount.get(key).count < _rootPopSizeSize * 0.005) {
                continue;
            }
            System.out.println(_topologyCount.get(key).toString() + "  " + key.toString());
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

    public static void main(String[] args) {
        String netA = "[0.036](((C:0.025199999999999997,(A:0.012599999999999998)I6#H1:0.012599999999999998::0.2)I5:0.010799999999999999,G:0.036)I1:0.05399999999999999,(R:0.072,(L:0.05399999999999999,(Q:0.036,I6#H1:0.0234::0.8)I4:0.018)I3:0.018)I2:0.018)I0;";
        String netB = "[0.036](((C:0.036,G:0.036)I1:0.036,((L:0.05399999999999999,(A:0.036,Q:0.036)I4:0.018)I3:0.009)I8#H1:0.009::0.3)I7:0.018,(R:0.072,I8#H1:0.009::0.7)I2:0.018)I0;";
        String netC = "[0.036]((C:0.018,G:0.018)I1:0.072,((R:0.026999999999999996,(Q:0.009)I8#H1:0.018::0.3)I7:0.045,(L:0.036,(A:0.018,I8#H1:0.009::0.7)I4:0.018)I3:0.036)I2:0.018)I0;";
        String netD = "[0.036]((G:0.036,(C:0.018,(A:0.009)I6#H1:0.009::0.2)I5:0.018)I1:0.05399999999999999,((R:0.026999999999999996,(Q:0.009)I8#H2:0.018::0.3)I7:0.045,(L:0.036,(I6#H1:0.009::0.8,I8#H2:0.009::0.7)I4:0.018)I3:0.036)I2:0.018)I0;";
        String net = "((R:0.07418337143189109,(L:0.06326041807321096,(A:0.04879879316357197,Q:0.04879879316357197):0.014461624909638995):0.010922953358680126):0.015047126703711472,(G:0.03786213356113884,C:0.03786213356113884):0.05136836457446372);";
        int start = 1000, end = 4001;
        SummaryBL sbl = new SummaryBL(netD);
        String file = "/scratch/jz55/usePolyMono/run/slurm-3213054_14.out";
        sbl.addFile(file, true, start, end);
        //sbl.report( 0.036 / 2, 1);
        sbl.reportForDensityPlot(file + ".DP", 0.036 / 2, 1);

//        SummaryBL sbl = new SummaryBL(net);
//        for(int i = 0; i < 10; i++) {
////            System.out.println("id = " + i);
////            SummaryBL sbl = new SummaryBL(net);
//            String file = "/Users/dw20/Downloads/phases/sim-mosq/scale/out/0.5/" + i + "/run.txt";
//            sbl.addFile(file, true, start, end);
////            sbl.report(
////                    1, 1
////                    //8E-7, 8E-8
////            );
//        }
//        sbl.report(
//                1, 1
//                //8E-7, 8E-8
//        );
    }

}