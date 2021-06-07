package edu.rice.cs.bioinfo.programs.phylonet.algos.simulator;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution.SNAPPLikelihoodSampling;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.Splitting;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.LineageConfiguration;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.R;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.RPattern;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 3/6/19
 * Time: 11:29 AM
 * To change this template use File | Settings | File Templates.
 */
public class SimLCInNetwork {
    private BiAllelicGTR _model;
    private SimGTInNetworkWithTheta _gtsim;
    private SimSNPInGT _snpsim;
    private Long _seed = null;
    private Random _random;

    public SimLCInNetwork(BiAllelicGTR model, Long seed) {
        _model = model;
        _seed = seed;
        _random = new Random(_seed);
        _gtsim = new SimGTInNetworkWithTheta();
        _gtsim.setSeed(_seed);
        _snpsim = new SimSNPInGT(model);
        _snpsim.setSeed(_seed);
    }

    public LineageConfiguration generateLC(Network network, Map<String, List<String>> species2alleles) {
        if(species2alleles == null) {
            species2alleles = new HashMap<>();
            for(Object nodeObj : network.getLeaves()) {
                NetNode leaf = (NetNode) nodeObj;
                species2alleles.put(leaf.getName(), new ArrayList<>());
                species2alleles.get(leaf.getName()).add(leaf.getName());
            }
        }

        BiAllelicGTR model = _model;
        double pi0 = model.getEquilibriumVector().get(0, 0);
        double pi1 = model.getEquilibriumVector().get(1, 0);
        double u = 1.0 / (2.0 * pi0);
        double v = 1.0 / (2.0 * pi1);
        double mu = 2.0 * u * v / (u + v);

        // <parent, child>
        Map<Tuple<NetNode,NetNode>, List<STINode>> netEdgeTop2geneLineages = new HashMap<>();
        Map<Tuple<NetNode,NetNode>, List<STINode>> netEdgeBottom2geneLineages = new HashMap<>();
        // Data of each node is the height of the node.
        STITree<Double> gt = _gtsim.generateOneGT(network, species2alleles, mu, netEdgeTop2geneLineages, netEdgeBottom2geneLineages);

        // Init heights of nodes.
        for(Object treenodeObj : gt.postTraverse()) {
            STINode<Double> treenode = (STINode<Double>) treenodeObj;
            if(treenode.isLeaf()) {
                treenode.setData(0.0);
            } else {
                STINode<Double> child = treenode.getChildren().iterator().next();
                treenode.setData(child.getData() + child.getParentDistance());
            }
        }

        Map<STINode<Double>, List<Tuple<NetNode,NetNode>>> geneLineage2netedges = new HashMap<>();

        // Get inverse mapping from gt node to net edges.
        for(Tuple<NetNode,NetNode> edge : netEdgeBottom2geneLineages.keySet()) {
            for(STINode treenode : netEdgeBottom2geneLineages.get(edge)) {
                if(!geneLineage2netedges.containsKey(treenode)) {
                    geneLineage2netedges.put(treenode, new ArrayList<>());
                }
                geneLineage2netedges.get(treenode).add(edge);
            }
        }

        for(Tuple<NetNode,NetNode> edge : netEdgeTop2geneLineages.keySet()) {
            for(STINode treenode : netEdgeTop2geneLineages.get(edge)) {
                if(!geneLineage2netedges.containsKey(treenode)) {
                    geneLineage2netedges.put(treenode, new ArrayList<>());
                }
                if(!geneLineage2netedges.get(treenode).contains(edge))
                    geneLineage2netedges.get(treenode).add(edge);
            }
        }

        // Init height of nodes.
        Map<NetNode, Double> netnode2height = new HashMap<>();
        for(Object nodeObj : Networks.postTraversal(network)) {
            NetNode node = (NetNode) nodeObj;
            if(node.isLeaf()) {
                netnode2height.put(node, 0.0);
            } else {
                NetNode child = (NetNode) node.getChildren().iterator().next();
                netnode2height.put(node, netnode2height.get(child) + child.getParentDistance(node));
            }
        }

        Map<NetNode, Map<NetNode, List<STINode<Double>>>> netnodeBottom2treenode = new HashMap<>();
        Map<NetNode, Map<NetNode, List<STINode<Double>>>> netnodeTop2treenode = new HashMap<>();
        for(Object nodeObj : Networks.postTraversal(network)) {
            NetNode node = (NetNode) nodeObj;
            netnodeBottom2treenode.put(node, new HashMap<>());
            netnodeTop2treenode.put(node, new HashMap<>());

            for(Object parentObj : node.getParents()) {
                NetNode parent = (NetNode) parentObj;
                netnodeBottom2treenode.get(node).put(parent, new ArrayList<>());
                netnodeTop2treenode.get(node).put(parent, new ArrayList<>());
                if(node.isLeaf()) {
                    for (String allele : species2alleles.get(node.getName())) {
                        //netnodeBottom2treenode.get(node).get(parent).add(gt.getNode(allele));
                    }
                }
            }


        }

        for(STINode<Double> treenode : geneLineage2netedges.keySet()) {
            Collections.sort(geneLineage2netedges.get(treenode), (Tuple<NetNode,NetNode> a, Tuple<NetNode,NetNode> b)->Double.compare(netnode2height.get(a.Item2), netnode2height.get(b.Item2)));
            STINode<Double> cur = treenode;
            // parent, child
            for(Tuple<NetNode,NetNode> edge : geneLineage2netedges.get(treenode)) {
                NetNode parent = edge.Item1;
                NetNode node = edge.Item2;
                if(node.isRoot()) break;
                if(cur.getParent().getData() < netnode2height.get(parent)) {
                    netnodeBottom2treenode.get(node).get(parent).add(cur);
                    break;
                }

                STINode<Double> newtreenode = cur.getParent().createChild();
                newtreenode.adoptChild(cur);


                newtreenode.setData(netnode2height.get(parent));
                newtreenode.setParentDistance(newtreenode.getParent().getData() - newtreenode.getData());
                cur.setParentDistance(newtreenode.getData() - cur.getData());

                if(!node.isRoot()) {
                    if(cur != treenode || cur.isLeaf())
                        netnodeBottom2treenode.get(node).get(parent).add(cur);
                    netnodeTop2treenode.get(node).get(parent).add(newtreenode);
                } else {
                    netnodeBottom2treenode.get(node).put(node, new ArrayList<>());
                    netnodeBottom2treenode.get(node).get(node).add(cur);
                }

                cur = newtreenode;
            }
        }

        Map<TNode, Character> markers = new HashMap<>();
        _snpsim.generateSingeSite(gt, markers);
        LineageConfiguration lc = new LineageConfiguration();

        for(Object nodeObj : Networks.postTraversal(network)) {
            NetNode node = (NetNode) nodeObj;
            for(Object parentObj : node.getParents()) {
                NetNode parent = (NetNode) parentObj;
                int c0 = 0, c1 = 0;
                for(STINode<Double> treenode : netnodeBottom2treenode.get(node).get(parent)) {

                    if(markers.get(treenode).equals('0')) {
                        c0++;
                    } else if(markers.get(treenode).equals('1')) {
                        c1++;
                    } else {
                        throw new RuntimeException("!!!!!!!");
                    }
                }
                lc.setBottom(node, parent, new R(c0 + c1, new int[]{c0}));

                c0 = 0;
                c1 = 0;
                for(STINode<Double> treenode : netnodeTop2treenode.get(node).get(parent)) {

                    if(markers.get(treenode).equals('0')) {
                        c0++;
                    } else if(markers.get(treenode).equals('1')) {
                        c1++;
                    } else {
                        throw new RuntimeException("!!!!!!!");
                    }
                }
                lc.setTop(node, parent, new R(c0 + c1, new int[]{c0}));
            }


        }

        int c0 = 0, c1 = 0;
        for(Object childObj : network.getRoot().getChildren()) {
            for (STINode<Double> treenode : netnodeTop2treenode.get(childObj).get(network.getRoot())) {

                if (markers.get(treenode).equals('0')) {
                    c0++;
                } else if (markers.get(treenode).equals('1')) {
                    c1++;
                } else {
                    throw new RuntimeException("!!!!!!!");
                }
            }
        }
        lc.setBottom(network.getRoot(), network.getRoot(), new R(c0 + c1, new int[]{c0}));

        return lc;
    }

    public static void main(String args[]) {
        double pi0 = 0.5;
        double pi1 = 1- pi0;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {pi0, pi1}, new double[] {1.0/ (2.0 * pi0)});


        double alpha1 = 0.5;
        double alpha2 = 0.4;
        double alpha3 = 0.6;
        double alpha4 = 0.8;
        double a = 0;
        double b = 0;
        double c = 0;
        double theta = 0.1;
        double gamma = 0.55;

        Map<Tuple<String, String>, Map<R, Integer>> countBottom = new HashMap<>();
        Map<Tuple<String, String>, Map<R, Integer>> countTop = new HashMap<>();

        String netstring = String.format("[%f](B:%f,A:%f)I0;", theta, alpha2, alpha2);
        //String netstring = String.format("[%f](((B:%f)I3#H1:%f::%f,A:%f)I1:%f,(C:%f,I3#H1:%f::%f)I2:%f)I0;", theta, alpha2, alpha1-alpha2, gamma, alpha1, alpha4-alpha1, alpha3, alpha3-alpha2, 1.0-gamma, alpha4-alpha3);

        Map<String, List<String>> species2alleles = new HashMap<>();
        species2alleles.put("A", new ArrayList<>());
        species2alleles.get("A").add("A_0");
        //species2alleles.get("A").add("A_1");
        species2alleles.put("B", new ArrayList<>());
        species2alleles.get("B").add("B_0");
        //species2alleles.get("B").add("B_1");
        species2alleles.put("C", new ArrayList<>());
        species2alleles.get("C").add("C_0");
        //species2alleles.get("C").add("C_1");

        SimLCInNetwork simulator = new SimLCInNetwork(BAGTRModel, new Random().nextLong());
        Map<Tuple<R, R>, Integer> count = new HashMap<>();
        Map<RPattern, Integer> countRPatterns = new HashMap<>();

        for(int i = 0 ; i < 100000 ; i++) {

            Network trueNetwork = Networks.readNetworkWithRootPop(netstring);
            LineageConfiguration lc = simulator.generateLC(trueNetwork, species2alleles);

            //for (NetNode node : lc.keySet()) {
            //    node.setName(lc.get(node).getNum(0) + "_" + lc.get(node).getNum(1));
            //}
            for(Object nodeObj : Networks.postTraversal(trueNetwork)) {
                NetNode node = (NetNode) nodeObj;
                if(node.isRoot()) {
                    R r = lc.getBottom(node, node);
                    if(!countBottom.containsKey(new Tuple<>(node.getName(), node.getName()))) {
                        countBottom.put(new Tuple<>(node.getName(), node.getName()), new HashMap<>());
                    }
                    if(!countBottom.get(new Tuple<>(node.getName(), node.getName())).containsKey(r)) {
                        countBottom.get(new Tuple<>(node.getName(), node.getName())).put(r, 0);
                    }

                    countBottom.get(new Tuple<>(node.getName(), node.getName())).put(r, countBottom.get(new Tuple<>(node.getName(), node.getName())).get(r) + 1);
                } else {
                    for(Object parentObj : node.getParents()) {
                        NetNode parent = (NetNode) parentObj;
                        R r = lc.getBottom(node, parent);
                        if(!countBottom.containsKey(new Tuple<>(node.getName(), parent.getName()))) {
                            countBottom.put(new Tuple<>(node.getName(), parent.getName()), new HashMap<>());
                        }
                        if(!countBottom.get(new Tuple<>(node.getName(), parent.getName())).containsKey(r)) {
                            countBottom.get(new Tuple<>(node.getName(), parent.getName())).put(r, 0);
                        }

                        countBottom.get(new Tuple<>(node.getName(), parent.getName())).put(r, countBottom.get(new Tuple<>(node.getName(), parent.getName())).get(r) + 1);
                    }
                }
            }

            for(Object nodeObj : Networks.postTraversal(trueNetwork)) {
                NetNode node = (NetNode) nodeObj;
                for(Object parentObj : node.getParents()) {
                    NetNode parent = (NetNode) parentObj;
                    R r = lc.getTop(node, parent);
                    if(!countTop.containsKey(new Tuple<>(node.getName(), parent.getName()))) {
                        countTop.put(new Tuple<>(node.getName(), parent.getName()), new HashMap<>());
                    }
                    if(!countTop.get(new Tuple<>(node.getName(), parent.getName())).containsKey(r)) {
                        countTop.get(new Tuple<>(node.getName(), parent.getName())).put(r, 0);
                    }

                    countTop.get(new Tuple<>(node.getName(), parent.getName())).put(r, countTop.get(new Tuple<>(node.getName(), parent.getName())).get(r) + 1);
                }
            }

            Map<String, R> leafR = new HashMap<>();
            for(Object nodeObj : Networks.postTraversal(trueNetwork)) {
                NetNode node = (NetNode) nodeObj;
                if(node.isLeaf()) {
                    NetNode parent = (NetNode) node.getParents().iterator().next();
                    leafR.put(node.getName(), lc.getBottom(node, parent));
                }
            }
            RPattern patttern = new RPattern(leafR);
            Map<RPattern, double[]> patterns = new HashMap<>();
            patterns.put(patttern, new double[]{1.0, 0.0});

            Splitting splitting = new Splitting(null, species2alleles, patttern, BAGTRModel);
            splitting.lineageConfigurations = lc;
            List<Splitting> splittings = new ArrayList<>();
            splittings.add(splitting);

            //double logL = SNAPPLikelihoodSampling.computeSNAPPLikelihoodST(trueNetwork, splittings, patterns, BAGTRModel);
            //System.out.println(Math.exp(logL));

            R ATop = lc.getTop(trueNetwork.findNode("A"), trueNetwork.findNode("I0"));
            R ABottom = lc.getBottom(trueNetwork.findNode("A"), trueNetwork.findNode("I0"));
            Tuple<R, R> ATuple = new Tuple<>(ATop, ABottom);
            if(!count.containsKey(ATuple)) count.put(ATuple, 0);
            count.put(ATuple, count.get(ATuple) + 1);

            if(!countRPatterns.containsKey(patttern)) {
                countRPatterns.put(patttern, 0);
            }
            countRPatterns.put(patttern, countRPatterns.get(patttern) + 1);
        }

        for(Tuple<String, String> tuple : countBottom.keySet()) {
            for(R r : countBottom.get(tuple).keySet()) {
                System.out.println(tuple.Item1 + " " + tuple.Item2 + " Bottom " + r + " " + countBottom.get(tuple).get(r));
            }
        }
        System.out.println();
        for(Tuple<String, String> tuple : countTop.keySet()) {
            for(R r : countTop.get(tuple).keySet()) {
                System.out.println(tuple.Item1 + " " + tuple.Item2 + " Top " + r + " " + countTop.get(tuple).get(r));
            }
        }
        System.out.println();

        for(Tuple<R, R> ATuple : count.keySet()) {
            System.out.println(ATuple.Item1 + " " + ATuple.Item2 + " " + count.get(ATuple));
        }

        System.out.println();

        for(RPattern rPattern : countRPatterns.keySet()) {
            System.out.println(rPattern + " " + countRPatterns.get(rPattern));
        }

        // Weird bahavior! If not burning the initial random double.
        int s = 0;
        SimLCInNetwork simulator0 = new SimLCInNetwork(BAGTRModel, new Random().nextLong());
        for(int i = 0 ; i < 10000 ; i++) {
            Network trueNetwork = Networks.readNetworkWithRootPop(netstring);

            double _u = 1.0;
            double _v = 1.0;
            double t = 0.1;
            double _pi1 = 0.5;
            double _pi0 = 0.5;
            double mu = 2.0 * _u * _v / (_u + _v);

            STITree tree = simulator0._gtsim.generateOneGT(trueNetwork, null, mu, null, null);
            t = (double )tree.getNode("A").getParentDistance() - 0.4;
            double x = 1.0 - Math.exp(-(_u + _v) * t);

            char parentState =  simulator0._random.nextDouble() < _pi0 ? '0' : '1';
            char curState = parentState == '0' ?
                    (simulator0._random.nextDouble() < (1.0 - _pi1 * x) ? '0' : '1')
                    :(simulator0._random.nextDouble() < (_pi0 * x) ? '0' : '1');
            if(curState == '1') s++;
        }

        System.out.println(s);
    }
}
