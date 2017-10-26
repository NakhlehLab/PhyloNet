package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.summary;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 1/19/17
 * Time: 6:13 PM
 * To change this template use File | Settings | File Templates.
 */
import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.clustering.ParentalTreeOperation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetworkMetricNakhleh;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

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
            if(this.parent == null)
                return "";
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
            if(this.parent == null) {
                out.println("[KEY]-> root : pop size");
                for(Double sp : listSp ) {
                    out.println(sp / popScale);
                }
                return;
            }
            out.println("[KEY]" + this.parent.getName() + " -> " + this.child.getName() + " : branch length");
            for(Double bl : listBL) {
                out.println(bl / scale);
            }
            if(this.listProb.size() > 0) {
                out.println("[KEY]" + this.parent.getName() + " -> " + this.child.getName() + " : inheritance prob");
                for(Double prob : listProb ) {
                    out.println(prob);
                }
            }
            if(this.listSp.size() > 0) {
                out.println("[KEY]" + this.parent.getName() + " -> " + this.child.getName() + " : pop size");
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
    private List<Network> _samples = new ArrayList<>();

    public SummaryBL(String s) {
        this._net = Networks.readNetwork(s);
        if(s.startsWith("[")) {
            double popSize = Double.parseDouble(s.substring(1, s.indexOf("]")));
            this._net.getRoot().setRootPopSize(popSize);
        }
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
                    if(s.contains("MC BEGINS")) {
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
                            _samples.add(Networks.readNetwork(s));
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
                                    _branches.get("root").addInfo(NetNode.NO_DISTANCE, NetNode.NO_PROBABILITY, popSize);
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
            if(br.parent == null) continue;
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

    class BoundingBox{
        public double x1 = 0;
        public double y1 = 0;
        public double x2 = 0;
        public double y2 = 0;
    }

    class TextBox {
        BoundingBox rect;
        public String info;
        public boolean frame;
    }

    BoundingBox getBoundingBoxOf(List<TextBox> boxes) {
        BoundingBox ret = new BoundingBox();
        for(TextBox box : boxes) {
            ret.x1 = Math.min(ret.x1, box.rect.x1);
            ret.y1 = Math.min(ret.y1, box.rect.y1);
            ret.x2 = Math.max(ret.x2, box.rect.x2);
            ret.y2 = Math.max(ret.y2, box.rect.y2);
        }
        return ret;
    }

    List<TextBox> drawNode(Branch br) {
        BoundingBox selfBox = new BoundingBox();
        TextBox textBox = new TextBox();

        selfBox.x1 = selfBox.y1 = 0;
        selfBox.x2 = br.meanSp;
        selfBox.y2 = br.meanBL;

        if(br.child.isLeaf()) {


        } else if(br.child.isTreeNode()) {
            List<BoundingBox> childrenBoundingBox = new ArrayList<>();
            List<TextBox> childDrawing1 = null;
            List<TextBox> childDrawing2 = null;

            for(String key : this._branches.keySet()) {
                Branch nextbr = this._branches.get(key);
                if(nextbr.parent != br.child) continue;
                List<TextBox> childDrawing = drawNode(nextbr);
                if(childDrawing1 == null)
                    childDrawing1 = childDrawing;
                else
                    childDrawing2 = childDrawing;
            }

            BoundingBox childBoundingBox1 = getBoundingBoxOf(childDrawing1);
        }



        return null;
    }

    public void reportForNetworkDrawing(double scale, double popScale) {
        double totalWidth = 0.0;
        double totalHeight = 0.0;

        double rootLength = 0.1;

        List<TextBox> elements = new ArrayList<>();

        _rootPopSizeAvg /= _rootPopSizeSize;
        _rootPopSizeStd = Math.sqrt(_rootPopSizeStd / _rootPopSizeSize - _rootPopSizeAvg * _rootPopSizeAvg);

        totalWidth += _rootPopSizeAvg / popScale;
        totalHeight += 0.1;
        Branch reticulationEdge1 = null;
        Branch reticulationEdge2 = null;

        for(String key : this._branches.keySet()) {
            Branch br = this._branches.get(key);
            if(br.parent == null) continue;
            br.child.setParentDistance(br.parent, br.meanBL);
            br.child.setParentSupport(br.parent, br.child.NO_SUPPORT);
            if(br.child.isNetworkNode()) {
                if(reticulationEdge1 == null)
                    reticulationEdge1 = br;
                else
                    reticulationEdge2 = br;
            }
        }


    }

    public void reportForDensityPlot(String filename, double scale, double popScale) {
        System.out.println("Samples: " + _size);
        try {
            PrintWriter out = new PrintWriter(filename);
            for(String key : this._branches.keySet()) {
                Branch br = this._branches.get(key);
                br.reportForDensityPlot(out, scale, popScale);
            }
            out.println("[END]");
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

    public void reportForTopology(String filename) {
        System.out.println("Samples: " + _samples.size());
        System.out.println("Correct topologies: " + _size);

        Tree backbone = getBackbone(_net);

        List<Double> treeDistances = new ArrayList<>();
        List<Double> networkDistances = new ArrayList<>();
        for(Network key : _samples) {
            if(key.getReticulationCount() == 0) {
                double minDist = 1e99;
                String s = Networks.getTopologyString(key);
                Tree t = Trees.readTree(s);
                minDist = Math.min(minDist, getRootedRobinsonFouldsDistance(backbone, t));
                treeDistances.add(minDist);
            } else {
                NetworkMetricNakhleh metric = new NetworkMetricNakhleh();
                double dist = metric.computeDistanceBetweenTwoNetworks(key, _net);
                networkDistances.add(dist);
            }
        }

        try {
            PrintWriter out = new PrintWriter(filename);
            out.println("[KEY]Tree");
            out.println(0.0);
            for(double t : treeDistances) {
                out.println(t);
            }
            out.println("[KEY]Network");
            out.println(0.0);
            for(double t : networkDistances) {
                out.println(t);
            }
            out.println("[END]");
            out.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    public void reportForHeight(String filename) {
        System.out.println("Samples: " + _samples.size());
        System.out.println("Correct topologies: " + _size);

        try {
            PrintWriter out = new PrintWriter(filename);
            out.println("[KEY]Height");
            out.println(getNetworkHeight(_net));

            for(Network key : _samples) {
                out.println(getNetworkHeight(key));
            }

            out.println("[END]");
            out.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    public void reportForReticulation(String filename) {
        System.out.println("Samples: " + _samples.size());
        System.out.println("Correct topologies: " + _size);

        try {
            PrintWriter out = new PrintWriter(filename);
            out.println("[KEY]Reticulation");
            out.println(_net.getReticulationCount());

            for(Network key : _samples) {
                out.println(key.getReticulationCount());
            }

            out.println("[END]");
            out.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
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
        res.put("root", new Branch(null, net.getRoot(), NetNode.NO_DISTANCE, NetNode.NO_PROBABILITY, net.getRoot().getRootPopSize()));
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

    public double getRootedRobinsonFouldsDistance(Tree tree1, Tree tree2) {
        SymmetricDifference symmetricDifference = new SymmetricDifference();
        symmetricDifference.computeDifference(tree1, tree2, true);
        double diff = symmetricDifference.getWeightedAverage();

        return diff;
    }

    public double getNetworkHeight(Network net) {
        NetNode node = (NetNode) net.getLeaves().iterator().next();
        double result = 0.0;
        while(!node.isRoot()) {
            NetNode parent = (NetNode) node.getParents().iterator().next();
            result += node.getParentDistance(parent);
            node = parent;
        }
        return result;
    }

    public Tree getBackbone(Network net) {
        Network curNet = Networks.readNetwork(net.toString());
        curNet.getRoot().setRootPopSize(net.getRoot().getRootPopSize());
        List<NetNode> networkNodes = new ArrayList<>();
        for(Object nodeObj : curNet.getNetworkNodes()) {
            networkNodes.add((NetNode) nodeObj);
        }
        for(NetNode node : networkNodes) {
            Iterator it = node.getParents().iterator();
            NetNode p1 = (NetNode) it.next();
            NetNode p2 = (NetNode) it.next();
            if(node.getParentProbability(p1) > node.getParentProbability(p2)) {
                p1 = p2;
            }
            p1.removeChild(node);
        }
        String s = Networks.getTopologyString(curNet);
        MutableTree t = (MutableTree) Trees.readTree(s);
        Trees.removeBinaryNodes(t);
        return t;
    }

    public static void main(String[] args) {
        String netA = "[0.036](((C:0.025199999999999997:0.036,(A:0.012599999999999998:0.036)I6#H1:0.012599999999999998:0.036:0.2)I5:0.010799999999999999:0.036,G:0.036:0.036)I1:0.05399999999999999:0.036,(R:0.072:0.036,(L:0.05399999999999999:0.036,(Q:0.036:0.036,I6#H1:0.0234:0.036:0.8)I4:0.018:0.036)I3:0.018:0.036)I2:0.018:0.036)I0;";
        String netB = "[0.036](((C:0.036,G:0.036)I1:0.036,((L:0.05399999999999999,(A:0.036,Q:0.036)I4:0.018)I3:0.009)I8#H1:0.009::0.3)I7:0.018,(R:0.072,I8#H1:0.009::0.7)I2:0.018)I0;";
        String netC = "[0.036]((C:0.018:0.036,G:0.018:0.036)I1:0.072:0.036,((R:0.026999999999999996:0.036,(Q:0.009:0.036)I8#H1:0.018:0.036:0.3)I7:0.045:0.036,(L:0.036:0.036,(A:0.018:0.036,I8#H1:0.009:0.036:0.7)I4:0.018:0.036)I3:0.036:0.036)I2:0.018:0.036)I0;";
        String netD = "[0.036]((G:0.036:0.036,(C:0.018:0.036,(A:0.009:0.036)I6#H1:0.009:0.036:0.2)I5:0.018:0.036)I1:0.05399999999999999:0.036,((R:0.026999999999999996:0.036,(Q:0.009:0.036)I8#H2:0.018:0.036:0.3)I7:0.045:0.036,(L:0.036:0.036,(I6#H1:0.009:0.036:0.8,I8#H2:0.009:0.036:0.7)I4:0.018:0.036)I3:0.036:0.036)I2:0.018:0.036)I0;";
        String netS = "[0.005]((((C:0.005:0.005)I1#H1:0.006:0.005:0.8,D:0.011:0.005):0.009:0.005,(B:0.014:0.005,I1#H1:0.009:0.005:0.2):0.006:0.005):0.005:0.005,A:0.025:0.005);";
        String netO = "[0.21545722449921065](((((HYB:0.003914413112003172:0.4176146908155346)#H1:0.017795921488291882:0.43115993308468703:0.27243274959227726,mcpmcp:0.021710334600295055:0.39534479563469904):0.05006671672056976:0.33564197409296903,mcplac:0.07177705132086482:0.2106105837692384):0.022926534912054042:0.19532354605930963,((#H1:0.024245525688431534:0.06424656769681:0.7275672504077227,vul:0.028159938800434707:0.1060232965423852):0.04314872010473817:0.47969536332364177,mcccal:0.07130865890517288:0.22924196754984502):0.02339492732774598:0.08797760028321075):0.011203153707864497:0.24137432580356583,cro:0.10590673994078335:0.28202492943752444);";
        //String netO="[0.12681791544104312](((mcccal:0.10324554898151703:0.07498756905531118,((HYB:0.018668688439637513:0.03429908136988121)#H1:0.052148254039290104:0.045084565479595835:0.573826907820398,vul:0.07081694247892761:0.03812038618556495):0.03242860650258941:0.03261690207000356):0.019549498068229462:0.037603755462335915,((#H1:0.029604547063352134:0.032054637614948284:0.426173092179602,mcpmcp:0.04827323550298965:0.021273866063836175):0.03578119089377865:0.033957981610605474,mcplac:0.0840544263967683:0.07805001108990976):0.03874062065297819:0.09196942861228893):0.00658368681952827:0.03596084334798809,cro:0.12937873386927476:0.11626211262353897);";
        //String net = "((R:0.07418337143189109,(L:0.06326041807321096,(A:0.04879879316357197,Q:0.04879879316357197):0.014461624909638995):0.010922953358680126):0.015047126703711472,(G:0.03786213356113884,C:0.03786213356113884):0.05136836457446372);";
        String netO8 = "[0.03678804256592915]((mccmcc:0.14889368885248497:0.0401585123579523,(sesspl:0.07484846132667702:0.05334010792803318,((HYB:0.029738544838728925:0.045259539795100646)#H1:0.004500135083003898:0.05739433597670084:0.5786813009538128,sesses:0.03423867992173282:0.03912830945064944):0.0406097814049442:0.021553843442927427):0.07404522752580794:0.11405841692866706):0.058943734376352225:0.04611915359113329,(#H1:0.010502859957408392:0.06686347677981164:0.42131869904618724,cae:0.04024140479613732:0.059065877126376044):0.1675960184326999:0.07103450968740839);";
        String netO4 = "[0.031269291383741966]((cae:0.023498520577909784:0.11453515690711917,(HYB:0.009167547496130975:0.03805835880867663)#H1:0.01433097308177881:0.03423044524104861:0.5101383020949305):0.17382587887270695:0.024235535337852533,(mccmcc:0.18439345369566149:0.06925787527744415,(mcplac:0.1465258072175068:0.08235404912736036,(#H1:0.014564359287868745:0.03939682755966605:0.4898616979050695,mcccal:0.02373190678399972:0.04882348145145685):0.12279390043350707:0.04534154430390398):0.03786764647815469:0.05099477836461758):0.012930945754955264:0.06687609363382643);";
        //String net = "[0.011246994934593928](wal:0.16856764465578844:0.024643805160552118,(((the57:0.005532380717127541:0.0328604191625876)#H1:0.008861256779885193:0.0349981605141176:0.5697578143531001,c513:0.014393637497012734:0.06774556777311104):0.02033133765658179:0.05536861556843345,(m523:0.03221433499579413:0.05265524780165408,(agla569:0.010982981528325176:0.04211202761569184,(#H1:4.459081325119832E-4:0.044983217465569256:0.43024218564689987,amar48:0.0059782888496395245:0.06404852435647355):0.005004692678685652:0.05142074487290094):0.021231353467468954:0.05110140055230508):0.0025106401578003923:0.05027578096009688):0.13384266950219392:0.02824200898517967);";
        String netR1 = "[0.006](((((Q:0.004:0.006)I5#H1:0.002:0.005:0.7,A:0.006:0.006)I3:0.006:0.005,L:0.012:0.006)I2:0.012:0.005,(I5#H1:0.003:0.005:0.3,R:0.007:0.006)I4:0.017:0.005)I1:0.016:0.005,C:0.04:0.006);";
        int start = 400, end = 3000;
        SummaryBL sbl = new SummaryBL(netO4);
        String file = "/Users/zhujiafan/Documents/PhyloDataResults/Ourisia/test4_2.txt";
        sbl.addFile(file, true, start, end);
        //sbl.report( 0.036 / 2, 1);
        sbl.report( 1, 1);
        //sbl.reportForDensityPlot(file + ".DP", 1, 1);

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