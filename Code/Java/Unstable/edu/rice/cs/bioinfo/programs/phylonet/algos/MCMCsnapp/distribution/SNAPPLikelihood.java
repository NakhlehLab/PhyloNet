package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.R;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.SNAPPAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.OneNucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.DoubleAdder;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 11/23/16
 * Time: 3:17 PM
 * To change this template use File | Settings | File Templates.
 */
public class SNAPPLikelihood {
    public static int ALGORITHM = 0;

    static public double verify(Network network, Map<String, String> alleles2species, List<Alignment> alignments, BiAllelicGTR BAGTRModel) {
        SNAPPLikelihood.ALGORITHM = 0;
        double t0 = SNAPPLikelihood.computeSNAPPLikelihood(network, null, alignments, BAGTRModel);
        SNAPPLikelihood.ALGORITHM = 1;
        double t1 = SNAPPLikelihood.computeSNAPPLikelihood(network, null, alignments, BAGTRModel);
        SNAPPLikelihood.ALGORITHM = 2;
        double t2 = SNAPPLikelihood.computeSNAPPLikelihood(network, null, alignments, BAGTRModel);
        if(Math.abs(t0 - t1) > 1e-6) System.out.println("ALGO0 != ALGO1");
        if(Math.abs(t0 - t2) > 1e-6) System.out.println("ALGO0 != ALGO2");
        if(Math.abs(t1 - t2) > 1e-6) System.out.println("ALGO1 != ALGO2");

        return t0;
    }

    static public double computeSNAPPLikelihood(Network network, Map<String, String> alleles2species, List<Alignment> alignments, BiAllelicGTR BAGTRModel) {

        /*double [] pi =  new double[2];
        double [] rate = new double[1];

        int totalSites = 0;
        for(Alignment alg : alignments) {
            int count = 0;
            for(String s : alg.getAlignment().values()) {
                for(int i = 0 ; i < s.length() ; i++) {
                    count += s.charAt(i) == '0' ? 0 : 2;
                }
            }
            totalSites += alg.getSiteCount() * alg.getTaxonSize();
            pi[1] += 1.0 * count;
        }
        pi[1] /= 2.0 * totalSites;
        pi[0] = 1 - pi[1];

        rate[0] = 1 / (2*pi[1]);

        R.dims = 1;
        _BAGTRModel = new BiAllelicGTR(pi, rate);*/
        double sum = 0;
        DoubleAdder adder = new DoubleAdder();

        int nameCount = 0;
        Network net = Networks.readNetwork(network.toString());
        for(Object node : net.dfs()) {
            NetNode mynode = (NetNode) node;
            if(mynode.getName().equals("")) {
                mynode.setName("II" + nameCount);
                nameCount++;
            }
        }
        String netstring = net.toString();



        for(Alignment alg : alignments) {
            ExecutorService executor = Executors.newFixedThreadPool(Utils._NUM_THREADS);

            List<String> names = alg.getTaxaNames();

            if(ALGORITHM == 0) {
                for (Integer represent : alg.getCache().keySet()) {
                    executor.execute(new Runnable() {
                        public void run() {
                            int cur = represent;
                            Integer count = alg.getCache().get(cur);
                            Map<String, Character> colorMap = new HashMap<String, Character>();
                            for (int i = names.size() - 1; i >= 0; i--) {
                                Character site = cur % 2 == 1 ? '1' : '0';
                                colorMap.put(names.get(i), site);
                                cur /= 2;
                            }
                            OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);
                            Network cloneNetwork = Networks.readNetwork(netstring);
                            cloneNetwork.getRoot().setRootPopSize(network.getRoot().getRootPopSize());
                            SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, BAGTRModel, (Utils._CONST_POP_SIZE ? network.getRoot().getRootPopSize() : null));
                            double likelihood = 0;
                            try {
                                likelihood = run.getProbability(converter);
                            } catch (Exception e) {
                                e.printStackTrace();
                                System.out.println("Exceptional network: " + netstring);
                            }
                            double logL = Math.log(likelihood);
                            adder.add(logL * count);
                        }
                    });
                }

                try {
                    executor.shutdown();
                    executor.awaitTermination(1000, TimeUnit.SECONDS);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            } else if(ALGORITHM == 1){
                int minDeep = Integer.MAX_VALUE;
                String bestName = null;
                Network cloneNetwork = Networks.readNetwork(netstring);
                cloneNetwork.getRoot().setRootPopSize(network.getRoot().getRootPopSize());
                int nameIndex = names.size() - 1;
                int bestNameIndex = 0;
                Map<NetNode,Boolean> node2update = new HashMap<>();
                for(String name : names) {
                    NetNode node = cloneNetwork.findNode(name);
                    Map<NetNode,Boolean> curnode2update = new HashMap<>();
                    int deep = 0;
                    while(!node.isRoot()) {
                        deep++;
                        curnode2update.put(node, true);
                        if(node.isNetworkNode()) { //should find longest path
                            deep = Integer.MAX_VALUE;
                            break;
                        }
                        node = (NetNode) node.getParents().iterator().next();
                    }
                    if(deep < minDeep) {
                        minDeep = deep;
                        bestName = name;
                        bestNameIndex = nameIndex;
                        node2update = curnode2update;
                    }
                    nameIndex--;
                }
                node2update.put(cloneNetwork.getRoot(), true);
                SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, BAGTRModel, (Utils._CONST_POP_SIZE ? network.getRoot().getRootPopSize() : null));

                for(int represent = 0 ; represent < (1 << names.size()) ; represent++) {
                    if(bestName != null && (represent & (1 << bestNameIndex)) > 0) continue;
                    int cur = represent;
                    Integer count = alg.getCache().get(cur);
                    Map<String, Character> colorMap = new HashMap<String, Character>();
                    for (int i = names.size() - 1; i >= 0; i--) {
                        Character site = cur % 2 == 1 ? '1' : '0';
                        colorMap.put(names.get(i), site);
                        cur /= 2;
                    }
                    OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);
                    double likelihood = 0;

                    if(count != null) {
                        try {
                            likelihood = run.getProbability(converter);
                        } catch (Exception e) {
                            e.printStackTrace();
                            System.out.println("Exceptional network: " + netstring);
                        }
                        double logL = Math.log(likelihood);
                        adder.add(logL * count);
                        //System.out.println(represent + " " + likelihood);
                    }
                    if(bestName == null) continue;

                    cur = represent | ((1 << bestNameIndex));
                    count = alg.getCache().get(cur);
                    if(count == null) continue;

                    colorMap.put(bestName, '1');
                    converter = new OneNucleotideObservation(colorMap);
                    if(alg.getCache().get(represent) == null) {
                        try {
                            likelihood = run.getProbability(converter);
                        } catch (Exception e) {
                            e.printStackTrace();
                            System.out.println("Exceptional network: " + netstring);
                        }
                    } else
                        likelihood = run.updateProbability(0, converter, null, node2update);
                    double logL = Math.log(likelihood);
                    adder.add(logL * count);
                    //System.out.println(cur + " " + likelihood);
                }
            } else {
                int minDeep = Integer.MAX_VALUE;
                String bestName = null;
                Network cloneNetwork = Networks.readNetwork(netstring);
                cloneNetwork.getRoot().setRootPopSize(network.getRoot().getRootPopSize());
                int nameIndex = names.size() - 1;
                int bestNameIndex = 0;
                Map<NetNode,Boolean> node2update = new HashMap<>();
                for(String name : names) {
                    NetNode node = cloneNetwork.findNode(name);
                    Map<NetNode,Boolean> curnode2update = new HashMap<>();
                    int deep = 0;
                    while(!node.isRoot()) {
                        deep++;
                        curnode2update.put(node, true);
                        if(node.isNetworkNode()) { //should find longest path
                            deep = Integer.MAX_VALUE;
                            break;
                        }
                        node = (NetNode) node.getParents().iterator().next();
                    }
                    if(deep < minDeep) {
                        minDeep = deep;
                        bestName = name;
                        bestNameIndex = nameIndex;
                        node2update = curnode2update;
                    }
                    nameIndex--;
                }
                node2update.put(cloneNetwork.getRoot(), true);

                if(bestName == null) {
                    for (Integer represent : alg.getCache().keySet()) {
                        executor.execute(new Runnable() {
                            public void run() {
                                int cur = represent;
                                Integer count = alg.getCache().get(cur);
                                Map<String, Character> colorMap = new HashMap<String, Character>();
                                for (int i = names.size() - 1; i >= 0; i--) {
                                    Character site = cur % 2 == 1 ? '1' : '0';
                                    colorMap.put(names.get(i), site);
                                    cur /= 2;
                                }
                                OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);
                                Network cloneNetwork = Networks.readNetwork(netstring);
                                cloneNetwork.getRoot().setRootPopSize(network.getRoot().getRootPopSize());
                                SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, BAGTRModel, (Utils._CONST_POP_SIZE ? network.getRoot().getRootPopSize() : null));
                                double likelihood = 0;
                                try {
                                    likelihood = run.getProbability(converter);
                                } catch (Exception e) {
                                    e.printStackTrace();
                                    System.out.println("Exceptional network: " + netstring);
                                }
                                double logL = Math.log(likelihood);
                                adder.add(logL * count);
                            }
                        });
                    }
                } else {
                    List<Integer> representations = new ArrayList<>();
                    for (Integer represent : alg.getCache().keySet()) {
                        if( (represent & (1 << bestNameIndex)) == 0) {
                            representations.add(represent);
                        } else {
                            if(!alg.getCache().containsKey(represent ^ (1 << bestNameIndex)))
                                representations.add(represent);
                        }
                    }
                    final int bestNameIndex1 = bestNameIndex;
                    final String bestName1 = bestName;
                    for (Integer represent : representations) {
                        final int represent1 = represent;
                        executor.execute(new Runnable() {
                            public void run() {
                                int cur = represent1;

                                Integer count = alg.getCache().get(cur);
                                Map<String, Character> colorMap = new HashMap<String, Character>();
                                for (int i = names.size() - 1; i >= 0; i--) {
                                    Character site = cur % 2 == 1 ? '1' : '0';
                                    colorMap.put(names.get(i), site);
                                    cur /= 2;
                                }
                                OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);
                                Network cloneNetwork = Networks.readNetwork(netstring);
                                cloneNetwork.getRoot().setRootPopSize(network.getRoot().getRootPopSize());
                                Map<NetNode,Boolean> node2update = new HashMap<NetNode, Boolean>();
                                NetNode node = cloneNetwork.findNode(bestName1);
                                while(!node.isRoot()) {
                                    node2update.put(node, true);
                                    node = (NetNode) node.getParents().iterator().next();
                                }
                                node2update.put(node, true);
                                SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, BAGTRModel, network.getRoot().getRootPopSize());

                                double likelihood = 0;

                                if(count != null) {
                                    try {
                                        likelihood = run.getProbability(converter);
                                    } catch (Exception e) {
                                        e.printStackTrace();
                                        System.out.println("Exceptional network: " + netstring);
                                    }
                                    double logL = Math.log(likelihood);
                                    adder.add(logL * count);
                                    //System.out.println(represent1 + " " + likelihood);
                                }
                                cur = represent1 ^ (1 << bestNameIndex1);
                                count = alg.getCache().get(cur);
                                colorMap.put(bestName1, '1');
                                converter = new OneNucleotideObservation(colorMap);
                                if( (represent1 & (1 << bestNameIndex1)) == 0 && count != null) {
                                    likelihood = run.updateProbability(0, converter, null, node2update);
                                    double logL = Math.log(likelihood);
                                    adder.add(logL * count);
                                    //System.out.println(cur + " " + likelihood);
                                }
                            }
                        });
                    }
                }

                try {
                    executor.shutdown();
                    executor.awaitTermination(1000, TimeUnit.SECONDS);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }


            sum += adder.sumThenReset();
        }

        return sum;

    }
}
