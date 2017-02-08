package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.R;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.SNAPPAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimSNPInNetwork;
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
    public static boolean useOnlyPolymorphic = false;
    public static boolean debugMode = false ;

    static public double verify(Network network, Map<String, String> alleles2species, List<Alignment> alignments, BiAllelicGTR BAGTRModel) {
        SNAPPLikelihood.ALGORITHM = 0;
        double t0 = SNAPPLikelihood.computeSNAPPLikelihood1(network, null, alignments, BAGTRModel);
        SNAPPLikelihood.ALGORITHM = 1;
        double t1 = SNAPPLikelihood.computeSNAPPLikelihood1(network, null, alignments, BAGTRModel);
        SNAPPLikelihood.ALGORITHM = 2;
        double t2 = SNAPPLikelihood.computeSNAPPLikelihood1(network, null, alignments, BAGTRModel);
        if(Math.abs(t0 - t1) > 1e-6) System.out.println("ALGO0 != ALGO1");
        if(Math.abs(t0 - t2) > 1e-6) System.out.println("ALGO0 != ALGO2");
        if(Math.abs(t1 - t2) > 1e-6) System.out.println("ALGO1 != ALGO2");

        return t0;
    }
    static public double computeSNAPPLikelihood(Network network, Map<String, String> alleles2species, List<Alignment> alignments, BiAllelicGTR BAGTRModel) {
        if(debugMode) {
            double result = verify(network, alleles2species, alignments, BAGTRModel);
            int numSites = 100000;

            Network trueNetwork = network.clone();
            trueNetwork.getRoot().setRootPopSize(network.getRoot().getRootPopSize());

            for(Object nodeObject : Networks.postTraversal(trueNetwork)) {
                NetNode node = (NetNode) nodeObject;
                for(Object parentObject :  node.getParents()) {
                    NetNode parent = (NetNode) parentObject;
                    node.setParentSupport(parent, trueNetwork.getRoot().getRootPopSize());
                }
            }
            int nameCount = 0;
            for(Object node : trueNetwork.dfs()) {
                NetNode mynode = (NetNode) node;
                if(mynode.getName().equals("")) {
                    mynode.setName("II" + nameCount);
                    nameCount++;
                }
            }

            SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, null);
            Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, numSites, !useOnlyPolymorphic);



            List<Map<String, String>> snpdata = new ArrayList<>();
            snpdata.add(onesnp);
            List<Alignment> alns = new ArrayList<>();
            for(Map<String, String> input : snpdata) {
                Alignment aln = new Alignment(input);

                Map<Integer, Integer> cache = new HashMap<>();

                for(int i = 0 ; i < aln.getSiteCount() ; i++) {
                    Map<String, Character> colorMap = new TreeMap<String, Character>();
                    for(String taxon : aln.getAlignment().keySet()) {
                        colorMap.put(taxon, aln.getAlignment().get(taxon).charAt(i));
                    }
                    Integer represent = 0;
                    for(String s : colorMap.keySet()) {
                        represent = (represent << 1) + (colorMap.get(s) == '0' ? 0 : 1);
                    }
                    if(!cache.containsKey(represent)) {
                        cache.put(represent, 0);
                    }
                    cache.put(represent, cache.get(represent) + 1);
                }

                if(useOnlyPolymorphic) {
                    cache.put(0, 0);
                    cache.put((1 << input.keySet().size()) - 1, 0);
                }

                aln.setCache(cache);
                alns.add(aln);



                for(Alignment alg : alns) {
                    double P0 = 0;
                    double P1 = 0;
                    Map<Integer, Double> likelihoods = new TreeMap<>();
                    List<String> names = alg.getTaxaNames();
                    double sum = 0;

                    for(Integer represent : alg.getCache().keySet()) {
                        int cur = represent;
                        Integer count = alg.getCache().get(cur);
                        Map<String, Character> colorMap = new HashMap<String, Character>();
                        for(int i = names.size() - 1 ; i >= 0 ; i--) {
                            Character site = cur % 2 == 1 ? '1' : '0';
                            colorMap.put(names.get(i), site);
                            cur /= 2;
                        }
                        OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);
                        Network cloneNetwork = Networks.readNetwork(trueNetwork.toString());
                        cloneNetwork.getRoot().setRootPopSize(trueNetwork.getRoot().getRootPopSize());
                        SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, BAGTRModel, null);
                        double likelihood = 0;
                        try {
                            likelihood = run.getProbability(converter);
                            if(represent == 0) P0 = likelihood;
                            if(represent == ((1 << names.size()) - 1)) P1 = likelihood;
                            sum += likelihood;
                            likelihoods.put(represent, likelihood);
                            //System.out.println(represent + " " + likelihood  + " " + count );
                        } catch(Exception e) {
                            e.printStackTrace();
                            System.out.println("Exceptional network");
                        }


                    }

                    for(Integer represent : likelihoods.keySet()) {
                        Integer count = alg.getCache().get(represent);
                        double likelihood = likelihoods.get(represent);
                        if(Math.abs( 1.0 * count / numSites - likelihood) > 1e-2) {
                            throw new RuntimeException("Likelihood error: " + trueNetwork.toString());
                        }

                    }


                }
            }
            return result;
        } else
            return computeSNAPPLikelihood1(network, alleles2species, alignments, BAGTRModel);
    }

    static public double computeSNAPPLikelihood1(Network network, Map<String, String> alleles2species, List<Alignment> alignments, BiAllelicGTR BAGTRModel) {

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
        DoubleAdder sitecounter = new DoubleAdder();
        DoubleAdder P0 = new DoubleAdder();
        DoubleAdder P1 = new DoubleAdder();
        DoubleAdder probAdder = new DoubleAdder();

        int nameCount = 0;
        Network net = Networks.readNetwork(network.toString());
        for(Object node : net.dfs()) {
            NetNode mynode = (NetNode) node;
            if(mynode.getName().equals("")) {
                mynode.setName("II" + nameCount);
                nameCount++;
            }
        }

        if(!Utils._ESTIMATE_POP_SIZE) {
            network.getRoot().setRootPopSize(Utils._POP_SIZE_MEAN);
        }
        String netstring = net.toString();


        for(Alignment alg : alignments) {
            ExecutorService executor = Executors.newFixedThreadPool(Utils._NUM_THREADS);

            List<String> names = alg.getTaxaNames();


            if(useOnlyPolymorphic) {
                alg.getCache().put(0, 0);
                alg.getCache().put((1 << names.size()) - 1, 0);

            }

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
                                if(debugMode) probAdder.add(likelihood);
                                if(useOnlyPolymorphic && represent == 0) P0.add(likelihood);
                                if(useOnlyPolymorphic && represent == ((1 << names.size()) - 1)) P1.add(likelihood);
                            } catch (Exception e) {
                                e.printStackTrace();
                                System.out.println("Exceptional network: " + netstring);
                            }
                            double logL = Math.log(likelihood);
                            adder.add(logL * count);
                            sitecounter.add(count);
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
                            if(debugMode) probAdder.add(likelihood);
                            if(useOnlyPolymorphic)
                                throw new Exception("not supported");
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
                    if(debugMode) probAdder.add(likelihood);
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
                                    if(debugMode) probAdder.add(likelihood);
                                    if(useOnlyPolymorphic && represent == 0) P0.add(likelihood);
                                    if(useOnlyPolymorphic && represent == ((1 << names.size()) - 1)) P1.add(likelihood);
                                } catch (Exception e) {
                                    e.printStackTrace();
                                    System.out.println("Exceptional network: " + netstring);
                                }
                                double logL = Math.log(likelihood);
                                adder.add(logL * count);
                                sitecounter.add(count);
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
                                        if(debugMode) probAdder.add(likelihood);
                                        if(useOnlyPolymorphic && represent1 == 0) P0.add(likelihood);
                                        if(useOnlyPolymorphic && represent1 == ((1 << names.size()) - 1)) P1.add(likelihood);
                                    } catch (Exception e) {
                                        e.printStackTrace();
                                        System.out.println("Exceptional network: " + netstring);
                                    }
                                    double logL = Math.log(likelihood);
                                    adder.add(logL * count);
                                    sitecounter.add(count);
                                    //System.out.println(represent1 + " " + likelihood);
                                }
                                cur = represent1 ^ (1 << bestNameIndex1);
                                count = alg.getCache().get(cur);
                                colorMap.put(bestName1, '1');
                                converter = new OneNucleotideObservation(colorMap);
                                if( (represent1 & (1 << bestNameIndex1)) == 0 && count != null) {
                                    likelihood = run.updateProbability(0, converter, null, node2update);
                                    double logL = Math.log(likelihood);
                                    if(debugMode) probAdder.add(likelihood);
                                    if(useOnlyPolymorphic && cur == 0) P0.add(likelihood);
                                    if(useOnlyPolymorphic && cur == ((1 << names.size()) - 1)) P1.add(likelihood);
                                    adder.add(logL * count);
                                    sitecounter.add(count);
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
            if(useOnlyPolymorphic) {
                sum -= sitecounter.sumThenReset() * Math.log(1.0 - P0.sumThenReset() - P1.sumThenReset());
            }
            if(debugMode) {
                double check = probAdder.sumThenReset();
                if(Math.abs(check - 1.0) > 1e-3) {
                    throw new RuntimeException("Likelihood error!");
                }
            }

        }
        //System.out.println(sum + " " + netstring);
        return sum;

    }
}
