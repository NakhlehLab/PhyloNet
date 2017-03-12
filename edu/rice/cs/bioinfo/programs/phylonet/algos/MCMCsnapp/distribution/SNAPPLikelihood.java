package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.R;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.RPattern;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.SNAPPAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimSNPInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.OneNucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import org.apache.commons.math3.util.ArithmeticUtils;

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

    static public BiAllelicGTR getModel(List<Alignment> alignments) {
        double [] pi =  new double[2];
        double [] rate = new double[1];
        rate[0] = 0.0;
        pi[0] = 0.0;
        pi[1] = 0.0;
        int totalSites = 0;

        int [] count = new int[3];
        count[0] = count[1] = count[2] = 0;
        for(Alignment alg : alignments) {
            for(String s : alg.getAlignment().values()) {
                for(int i = 0 ; i < s.length() ; i++) {
                    count[s.charAt(i) - '0']++;
                }
            }
            totalSites += alg.getSiteCount() * alg.getTaxonSize();
        }

        if(Utils._DIPLOID_SPECIES == null) {
            pi[1] = 1.0 * count[1] / totalSites;
        } else {
            pi[1] = 1.0 * (count[1] + 2 * count[2]) / (totalSites * 2);
        }
        pi[0] = 1 - pi[1];

        rate[0] = 1 / (2*pi[0]);

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(pi, rate);

        return BAGTRModel;
    }

    static public Map<RPattern, double[]> diploidSequenceToPatterns(Map<String, String> alleles2species, List<Alignment> alignments) {
        Map<RPattern, double[]> result = new HashMap<>();
        for(Alignment aln : alignments) {
            for(int i = 0 ; i < aln.getSiteCount() ; i++) {
                Map<String, Tuple<int[], int[]>> currentPattern = new HashMap<>(); //item1[0]: n, item2[0]: R0
                double weight = 1.0;
                for(String allele : aln.getTaxaNames()) {
                    char c = aln.getAlignment().get(allele).charAt(i);
                    String species;
                    if(alleles2species == null)
                        species = allele;
                    else
                        species = alleles2species.get(allele);
                    if(!currentPattern.containsKey(species))
                        currentPattern.put(species, new Tuple<>(new int[]{0}, new int[]{0}));
                    if(c == '-') continue;
                    else if(c == '0') currentPattern.get(species).Item2[0] += 2;
                    else if(c == '1') currentPattern.get(species).Item2[0] += 1;
                    currentPattern.get(species).Item1[0] += 2;

                    if(c == '1')
                        weight *= 2;
                }

                boolean notGood = false;
                Map<String, R> newPattern = new HashMap<>();
                for(String species : currentPattern.keySet()) {
                    if(currentPattern.get(species).Item1[0] == 0) notGood = true;   //ignore the case that there is a branch contains no data
                    weight /= ArithmeticUtils.binomialCoefficient(currentPattern.get(species).Item1[0], currentPattern.get(species).Item2[0]);
                    newPattern.put(species, new R(currentPattern.get(species).Item1[0], currentPattern.get(species).Item2));
                }

                if(notGood) continue;

                RPattern rpattern = new RPattern(newPattern);
                if(!result.containsKey(rpattern))
                    result.put(rpattern, new double[]{0.0, 0.0});
                result.get(rpattern)[0] += 1.0;
                result.get(rpattern)[1] += Math.log(weight);
            }
        }
        return result;
    }

    static public Map<RPattern, double[]> haploidSequenceToPatterns(Map<String, String> alleles2species, List<Alignment> alignments) {
        Map<RPattern, double[]> result = new HashMap<>();
        for(Alignment aln : alignments) {
            for(int i = 0 ; i < aln.getSiteCount() ; i++) {
                Map<String, Tuple<int[], int[]>> currentPattern = new HashMap<>();
                double weight = 1.0;
                for(String allele : aln.getTaxaNames()) {
                    char c = aln.getAlignment().get(allele).charAt(i);
                    String species;
                    if(alleles2species == null)
                        species = allele;
                    else
                        species = alleles2species.get(allele);
                    if(!currentPattern.containsKey(species))
                        currentPattern.put(species, new Tuple<>(new int[]{0}, new int[]{0}));
                    if(c == '-') continue;
                    else if(c == '0') currentPattern.get(species).Item2[0] += 1;
                    currentPattern.get(species).Item1[0] += 1;

                }

                boolean notGood = false;
                Map<String, R> newPattern = new HashMap<>();
                for(String species : currentPattern.keySet()) {
                    if(currentPattern.get(species).Item1[0] == 0) notGood = true;
                    weight /= ArithmeticUtils.binomialCoefficient(currentPattern.get(species).Item1[0], currentPattern.get(species).Item2[0]);
                    newPattern.put(species, new R(currentPattern.get(species).Item1[0], currentPattern.get(species).Item2));
                }

                if(notGood) continue;

                RPattern rpattern = new RPattern(newPattern);
                if(!result.containsKey(rpattern))
                    result.put(rpattern, new double[]{0.0, 0.0});
                result.get(rpattern)[0] += 1.0;
                result.get(rpattern)[1] += Math.log(weight);
            }
        }
        return result;
    }

    static public double computeSNAPPLikelihoodST(Network network, Map<RPattern, double[]> patterns, BiAllelicGTR BAGTRModel) {
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

        Double theta = null;
        if(Utils._CONST_POP_SIZE)
            theta = network.getRoot().getRootPopSize();

        String netstring = net.toString();

        double sum = 0.0;
        R.maxLineages = patterns.keySet().iterator().next().sumLineages();
        Network cloneNetwork = Networks.readNetwork(netstring);
        cloneNetwork.getRoot().setRootPopSize(network.getRoot().getRootPopSize());
        SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, BAGTRModel, theta);
        //long start = System.currentTimeMillis();
        for(RPattern pattern : patterns.keySet()) {
            double count = patterns.get(pattern)[0];
            double correction = patterns.get(pattern)[1];

            double likelihood = 0;
            //long start = System.currentTimeMillis();
            try {
                likelihood = run.getProbability(pattern);
                sum += Math.log(likelihood) * count + correction;
                //System.out.println(represent + " " + likelihood  + " " + count );
            } catch(Exception e) {
                e.printStackTrace();
                System.out.println("Exceptional network" + netstring);
            }
            //double time1 = (System.currentTimeMillis()-start)/1000.0;
            //System.out.println(time1);
        }
        //System.out.println((System.currentTimeMillis()-start)/1000.0);
        return sum;
    }

    static public double computeSNAPPLikelihoodMT(Network network, Map<RPattern, double[]> patterns, BiAllelicGTR BAGTRModel) {
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

        final Double theta = Utils._CONST_POP_SIZE ? network.getRoot().getRootPopSize(): null;

        String netstring = net.toString();

        ExecutorService executor = Executors.newFixedThreadPool(Utils._NUM_THREADS);

        double sum = 0.0;
        DoubleAdder adder = new DoubleAdder();

        for(RPattern pattern : patterns.keySet()) {
            executor.execute(new Runnable() {
                public void run() {
                    double count = patterns.get(pattern)[0];
                    double correction = patterns.get(pattern)[1];
                    Network cloneNetwork = Networks.readNetwork(netstring);
                    cloneNetwork.getRoot().setRootPopSize(network.getRoot().getRootPopSize());
                    R.maxLineages = pattern.sumLineages();
                    SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, BAGTRModel, theta);
                    double likelihood = 0;
                    try {
                        likelihood = run.getProbability(pattern);
                        adder.add(Math.log(likelihood) * count + correction);
                        //System.out.println(represent + " " + likelihood  + " " + count );
                    } catch (Exception e) {
                        e.printStackTrace();
                        System.out.println("Exceptional network" + netstring);
                    }
                }
            });
        }
        try {
            executor.shutdown();
            executor.awaitTermination(1000, TimeUnit.SECONDS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        sum = adder.sumThenReset();
        return sum;
    }

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

            if(ALGORITHM == -1) {
                for(int i = 0 ; i < alg.getSiteCount() ; i++) {
                    final int cursite = i;
                    //executor.execute(new Runnable() {
                    //    public void run() {
                            Map<String, Character> colorMap = new HashMap<String, Character>();
                            for(String taxon : alg.getAlignment().keySet()) {
                                colorMap.put(taxon, alg.getAlignment().get(taxon).charAt(cursite));
                            }
                            OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);
                            Network cloneNetwork = Networks.readNetwork(netstring);
                            cloneNetwork.getRoot().setRootPopSize(network.getRoot().getRootPopSize());
                            SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, alleles2species, BAGTRModel, (Utils._CONST_POP_SIZE ? network.getRoot().getRootPopSize() : null));
                            double likelihood = 0;
                            try {
                                likelihood = run.getProbability(converter);
                            } catch (Exception e) {
                                e.printStackTrace();
                                System.out.println("Exceptional network: " + netstring);
                            }
                            double logL = Math.log(likelihood);
                            adder.add(logL);
                    //    }
                    //});
                }
                if(useOnlyPolymorphic) {
                    for(int i = 0 ; i <= 1 ; i++) {
                        final char curNum = i == 0 ? '0' : '1';
                    //    executor.execute(new Runnable() {
                    //        public void run() {
                                Map<String, Character> colorMap = new HashMap<String, Character>();
                                for(String taxon : alg.getAlignment().keySet()) {
                                    colorMap.put(taxon, curNum);
                                }
                                OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);
                                Network cloneNetwork = Networks.readNetwork(netstring);
                                cloneNetwork.getRoot().setRootPopSize(network.getRoot().getRootPopSize());
                                SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, BAGTRModel, (Utils._CONST_POP_SIZE ? network.getRoot().getRootPopSize() : null));
                                double likelihood = 0;
                                try {
                                    likelihood = run.getProbability(converter);
                                    if(curNum == '0') P0.add(likelihood);
                                    else P1.add(likelihood);
                                } catch (Exception e) {
                                    e.printStackTrace();
                                    System.out.println("Exceptional network: " + netstring);
                                }
                            }
                    //    });
                    //}
                }

            } else if(ALGORITHM == 0) {
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

    private static void testMultipleAllelesBiallelic() {
        R.dims = 1;
        double sum =0;
        char[] nucleotides = {'0','1', '2'};
        int index = 0;
        HashSet<String> tried = new HashSet<String>();
        Map<String, String> colorMap = new HashMap<>();
        colorMap.put("a", "" );
        colorMap.put("b1", "" );
        colorMap.put("b2", "" );
        colorMap.put("b3", "" );
        colorMap.put("c", "" );
        for (char a : nucleotides)
            for (char b1: nucleotides)
                for (char b2: nucleotides)
                    for (char b3: nucleotides)
                        for (char c: nucleotides)
                        {
                            //if(index++<4)continue;

                            colorMap.put("a", colorMap.get("a") + a);
                            colorMap.put("b1", colorMap.get("b1") + b1);
                            colorMap.put("b2", colorMap.get("b2") + b2);
                            colorMap.put("b3", colorMap.get("b3") + b3);
                            colorMap.put("c", colorMap.get("c") + c);



                            System.out.println("Observation:  " + a +',' + b1 + ','+ b2 + "," +  b3 + "," + c);


                        }


        Network speciesNetworkTopology = Networks.readNetwork("(A:.5,(B:.2,C:.2)n1:.3)n0;");
        //Network speciesNetworkTopology = Networks.readNetwork("((((C:0.07,(B:0.06)I7#H1:0.01::0.2)I4:0.02,(I7#H1:0.02::0.8)I6#H2:0.01::0.2)I3:0.02,(I6#H2:0.02::0.8)I5#H3:0.01::0.2)I2:0.07,(I5#H3:0.06::0.8,A:0.16)I1:0.02)I0;");
        //Network speciesNetworkTopology = Networks.readNetwork("((((C:0.07)I4:0.02,((B:0.06)I7:0.02)I6#H2:0.01::0.2)I3:0.02,(I6#H2:0.02::0.8)I5#H3:0.01::0.2)I2:0.07,(I5#H3:0.06::0.8,A:0.16)I1:0.02)I0;");
        //Network speciesNetworkTopology = Networks.readNetwork("((A:.2,X#H1:0::0.2)n2:0.3,((B:.1)X#H1:0::0.8,C:.2)n1:.3)root;");
        speciesNetworkTopology.getRoot().setRootPopSize(0.001);

        Map<String, String> allele2species = new HashMap<String, String>();
        allele2species.put("a","A");
        allele2species.put("b1","B");
        allele2species.put("b2","B");
        allele2species.put("b3","B");
        allele2species.put("c","C");
        //allele2species = null;

        BiAllelicGTR gtrModel = new BiAllelicGTR(new double[] {0.5, 0.5}, new double[] {1.0});

        Alignment aln = new Alignment(colorMap);
        List<Alignment> alns = new ArrayList<>();
        alns.add(aln);

        sum += computeSNAPPLikelihoodST(speciesNetworkTopology, diploidSequenceToPatterns(allele2species, alns), gtrModel);

        System.out.println("Sum: " + sum);




    }

    private static void testMultipleAllelesBiallelic2() {
        R.dims = 1;
        double sum =0;
        char[] nucleotides = {'0','1', '2'};
        int index = 0;
        HashSet<String> tried = new HashSet<String>();

        for (char a : nucleotides)
            for (char b1: nucleotides)
                for (char b2: nucleotides)
                    for (char b3: nucleotides)
                        for (char c: nucleotides)
                        {
                            long start = System.currentTimeMillis();
                            //if(index++<4)continue;
                            Map<String, String> colorMap = new HashMap<>();
                            colorMap.put("a", "" );
                            colorMap.put("b1", "" );
                            colorMap.put("b2", "" );
                            colorMap.put("b3", "" );
                            colorMap.put("c", "" );

                            colorMap.put("a", colorMap.get("a") + a);
                            colorMap.put("b1", colorMap.get("b1") + b1);
                            colorMap.put("b2", colorMap.get("b2") + b2);
                            colorMap.put("b3", colorMap.get("b3") + b3);
                            colorMap.put("c", colorMap.get("c") + c);



                            System.out.println("Observation:  " + a +',' + b1 + ','+ b2 + "," +  b3 + "," + c);

                            //Network speciesNetworkTopology = Networks.readNetwork("(A:.5,(B:.2,C:.2)n1:.3)n0;");
                            Network speciesNetworkTopology = Networks.readNetwork("((((C:0.07,(B:0.06)I7#H1:0.01::0.2)I4:0.02,(I7#H1:0.02::0.8)I6#H2:0.01::0.2)I3:0.02,(I6#H2:0.02::0.8)I5#H3:0.01::0.2)I2:0.07,(I5#H3:0.06::0.8,A:0.16)I1:0.02)I0;");
                            //Network speciesNetworkTopology = Networks.readNetwork("((((C:0.07)I4:0.02,((B:0.06)I7:0.02)I6#H2:0.01::0.2)I3:0.02,(I6#H2:0.02::0.8)I5#H3:0.01::0.2)I2:0.07,(I5#H3:0.06::0.8,A:0.16)I1:0.02)I0;");
                            //Network speciesNetworkTopology = Networks.readNetwork("((A:.2,X#H1:0::0.2)n2:0.3,((B:.1)X#H1:0::0.8,C:.2)n1:.3)root;");
                            speciesNetworkTopology.getRoot().setRootPopSize(0.001);

                            Map<String, String> allele2species = new HashMap<String, String>();
                            allele2species.put("a","A");
                            allele2species.put("b1","B");
                            allele2species.put("b2","B");
                            allele2species.put("b3","B");
                            allele2species.put("c","C");
                            //allele2species = null;

                            BiAllelicGTR gtrModel = new BiAllelicGTR(new double[] {0.5, 0.5}, new double[] {1.0});

                            Alignment aln = new Alignment(colorMap);
                            List<Alignment> alns = new ArrayList<>();
                            alns.add(aln);

                            sum += Math.exp(computeSNAPPLikelihoodST(speciesNetworkTopology, diploidSequenceToPatterns(allele2species, alns), gtrModel));
                            System.out.println((System.currentTimeMillis()-start)/1000.0);
                        }




        System.out.println("Sum: " + sum);




    }

    public static void main(String[] args) {

        testMultipleAllelesBiallelic2();
    }
}
