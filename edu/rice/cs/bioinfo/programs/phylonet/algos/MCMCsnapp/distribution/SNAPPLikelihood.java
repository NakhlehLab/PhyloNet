package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.MarkerSeq;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import org.apache.commons.math3.util.ArithmeticUtils;

import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.DoubleAdder;
import java.util.concurrent.atomic.LongAdder;

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
    public static boolean usePseudoLikelihood = false;
    public static boolean useApproximateBayesian = false;
    public static boolean timeSavingMode = false;
    public static boolean debugMode = false ;
    public static LongAdder workloadCounter = null;
    public static SNAPPPseudoLikelihood pseudoLikelihood = null;
    public static ApproximateBayesian approximateBayesian = null;

    static public BiAllelicGTR getModel(List<MarkerSeq> markerSeqs) {
        double [] pi =  new double[2];
        double [] rate = new double[1];
        rate[0] = 0.0;
        pi[0] = 0.0;
        pi[1] = 0.0;
        int totalSites = 0;

        int [] count = new int[3];
        count[0] = count[1] = count[2] = 0;
        for(MarkerSeq alg : markerSeqs) {
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

    static public Map<String, String> randomlyPhasing(Map<String, String> sequence, Random random) {
        if(random == null)
            random = new Random();
        Map<String, StringBuilder> newbuilder = new HashMap<>();
        for(String taxon : sequence.keySet()) {
            newbuilder.put(taxon, new StringBuilder());
            for(int i = 0 ; i < sequence.get(taxon).length() ; i++) {
                char c = sequence.get(taxon).charAt(i);
                if(c == '1') {
                    if(random.nextDouble() < 0.5) c = '0';
                    else c = '2';
                }

                if(c == '2')
                    c = '1';
                newbuilder.get(taxon).append(c);
            }
        }

        Map<String, String> ret = new HashMap<>();
        for(String taxon : newbuilder.keySet()) {
            ret.put(taxon, newbuilder.get(taxon).toString());
        }
        return ret;
    }

    // TODO: not fully tested
    static public Map<RPattern, double[]> polyploidSequenceToPatterns(Map<String, String> alleles2species, List<MarkerSeq> markerSeqs, int totalPerSpecies) {
        Map<RPattern, double[]> result = new HashMap<>();
        Map<String, Integer> maxLineages = new HashMap<>();
        for(MarkerSeq aln : markerSeqs) {
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
                    if(c == '-' || c == '?') continue;
                    else currentPattern.get(species).Item2[0] += totalPerSpecies - (c - '0');
                    currentPattern.get(species).Item1[0] += totalPerSpecies;

                    weight *= ArithmeticUtils.binomialCoefficient(totalPerSpecies, c - '0');
                }

                boolean notGood = false;
                Map<String, R> newPattern = new HashMap<>();
                for(String species : currentPattern.keySet()) {
                    if(currentPattern.get(species).Item1[0] == 0) notGood = true;   //ignore the case that there is a branch contains no data
                    weight /= ArithmeticUtils.binomialCoefficient(currentPattern.get(species).Item1[0], currentPattern.get(species).Item2[0]);
                    newPattern.put(species, new R(currentPattern.get(species).Item1[0], currentPattern.get(species).Item2));
                }

                if(notGood) continue;

                for(String species : newPattern.keySet()) {
                    if(!maxLineages.containsKey(species))
                        maxLineages.put(species, 0);
                    maxLineages.put(species, Math.max(maxLineages.get(species), newPattern.get(species).getN()));
                }

                RPattern rpattern = new RPattern(newPattern);
                if(useOnlyPolymorphic && rpattern.isMonomorphic()) continue;

                if(!result.containsKey(rpattern))
                    result.put(rpattern, new double[]{0.0, 0.0});
                result.get(rpattern)[0] += 1.0;
                result.get(rpattern)[1] += Math.log(weight);
            }
        }

        if(useOnlyPolymorphic) {
            for(int i = 0 ; i <= R.dims ; i++) {
                Map<String, R> newPattern = new HashMap<>();

                for(String species : maxLineages.keySet()) {
                    int a[] = new int[R.dims];
                    if(i < R.dims)
                        a[i] = maxLineages.get(species);
                    newPattern.put(species, new R(maxLineages.get(species), a));
                }
                RPattern rpattern = new RPattern(newPattern);
                result.put(rpattern, new double[]{0.0, 0.0});
            }
        }

        return result;
    }

    static public Map<RPattern, double[]> diploidSequenceToPatterns(Map<String, String> alleles2species, List<MarkerSeq> markerSeqs) {
        Map<RPattern, double[]> result = new HashMap<>();
        Map<String, Integer> maxLineages = new HashMap<>();
        for(MarkerSeq aln : markerSeqs) {
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
                    if(c == '-' || c == '?') continue;
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

                for(String species : newPattern.keySet()) {
                    if(!maxLineages.containsKey(species))
                        maxLineages.put(species, 0);
                    maxLineages.put(species, Math.max(maxLineages.get(species), newPattern.get(species).getN()));
                }

                RPattern rpattern = new RPattern(newPattern);
                if(useOnlyPolymorphic && rpattern.isMonomorphic()) continue;

                if(!result.containsKey(rpattern))
                    result.put(rpattern, new double[]{0.0, 0.0});
                result.get(rpattern)[0] += 1.0;
                result.get(rpattern)[1] += Math.log(weight);
            }
        }

        if(useOnlyPolymorphic) {
            for(int i = 0 ; i <= R.dims ; i++) {
                Map<String, R> newPattern = new HashMap<>();

                for(String species : maxLineages.keySet()) {
                    int a[] = new int[R.dims];
                    if(i < R.dims)
                    a[i] = maxLineages.get(species);
                    newPattern.put(species, new R(maxLineages.get(species), a));
                }
                RPattern rpattern = new RPattern(newPattern);
                result.put(rpattern, new double[]{0.0, 0.0});
            }
        }

        return result;
    }

    static public Map<RPattern, double[]> haploidSequenceToPatterns(Map<String, String> alleles2species, List<MarkerSeq> markerSeqs) {
        Map<RPattern, double[]> result = new HashMap<>();
        Map<String, Integer> maxLineages = new HashMap<>();
        for(MarkerSeq aln : markerSeqs) {
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
                    if(c == '-' || c == '?') continue;
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

                for(String species : newPattern.keySet()) {
                    if(!maxLineages.containsKey(species))
                        maxLineages.put(species, 0);
                    maxLineages.put(species, Math.max(maxLineages.get(species), newPattern.get(species).getN()));
                }

                RPattern rpattern = new RPattern(newPattern);
                if(useOnlyPolymorphic && rpattern.isMonomorphic()) continue;

                if(!result.containsKey(rpattern))
                    result.put(rpattern, new double[]{0.0, 0.0});
                result.get(rpattern)[0] += 1.0;
                result.get(rpattern)[1] += Math.log(weight);
            }
        }

        if(useOnlyPolymorphic) {
            for(int i = 0 ; i <= R.dims ; i++) {
                Map<String, R> newPattern = new HashMap<>();

                for(String species : maxLineages.keySet()) {
                    int a[] = new int[R.dims];
                    if(i < R.dims)
                        a[i] = maxLineages.get(species);
                    newPattern.put(species, new R(maxLineages.get(species), a));
                }
                RPattern rpattern = new RPattern(newPattern);
                result.put(rpattern, new double[]{0.0, 0.0});
            }
        }

        return result;
    }

    static public double computeApproximateBayesian(Network network, Map<String, String> alleles2species, List<MarkerSeq> markerSeqs, BiAllelicGTR BAGTRModel, Map<String, Double> abcData) {
        if(approximateBayesian == null) {
            System.out.println("Initiating approximate bayesian");
            approximateBayesian = new ApproximateBayesian(alleles2species, markerSeqs, markerSeqs.get(0)._diploid, markerSeqs.get(0)._dominant != null, markerSeqs.get(0)._polyploid, useOnlyPolymorphic, BAGTRModel);
            System.out.println("Finished initiating approximate bayesian");
        }

        Network net = network.clone();

        int nameCount = 0;
        for(Object node : net.dfs()) {
            NetNode mynode = (NetNode) node;
            if(mynode.getName().equals("")) {
                mynode.setName("II" + nameCount);
                nameCount++;
            }
        }

        return approximateBayesian.computeApproximateBayesianMT(net, abcData);

    }

    static public double computeSNAPPPseudoLikelihood(Network network, Map<String, String> alleles2species, List<MarkerSeq> markerSeqs, BiAllelicGTR BAGTRModel) {
        if(pseudoLikelihood == null) {
            System.out.println("Initiating pseudo-likelihood");
            pseudoLikelihood = new SNAPPPseudoLikelihood(alleles2species, markerSeqs, markerSeqs.get(0)._diploid);
            System.out.println("Finished initiating pseudo-likelihood");
        }

        Network net = network.clone();

        int nameCount = 0;
        for(Object node : net.dfs()) {
            NetNode mynode = (NetNode) node;
            if(mynode.getName().equals("")) {
                mynode.setName("II" + nameCount);
                nameCount++;
            }
        }

        return pseudoLikelihood.computeSNAPPPseudoLogLikelihoodMT(net, alleles2species, BAGTRModel);

    }

    static public double computeSNAPPLikelihood(Network network, Map<RPattern, double[]> patterns, BiAllelicGTR BAGTRModel) {
        double prob;

        if(workloadCounter == null) {
            workloadCounter = new LongAdder();
            Runtime.getRuntime().addShutdownHook(new Thread(new Runnable() {
                public void run() {
                    System.out.println("Total number of network processed: " + workloadCounter.sum());
                    System.out.println("FMatrix cache hit rate: " + FMatrix.cacheHit.sum() / FMatrix.cacheAccess.sum());
                }
            }));
        }

        if(Utils._NUM_THREADS == 1)
            prob = SNAPPLikelihood.computeSNAPPLikelihoodST(network, patterns, BAGTRModel);
        else
            prob = SNAPPLikelihood.computeSNAPPLikelihoodMTC(network, patterns, BAGTRModel);
        if(!usePseudoLikelihood)
            workloadCounter.increment();
        return prob;
    }

    static public double computeSNAPPLikelihoodST(Network network, Map<RPattern, double[]> patterns, BiAllelicGTR BAGTRModel) {
        if(workloadCounter == null) {
            workloadCounter = new LongAdder();
            Runtime.getRuntime().addShutdownHook(new Thread(new Runnable() {
                public void run() {
                    System.out.println("Total number of network processed: " + workloadCounter.sum());
                    System.out.println("FMatrix cache hit rate: " + FMatrix.cacheHit.sum() / FMatrix.cacheAccess.sum());
                }
            }));
        }

        //System.out.println(network.toString());

        int nameCount = 0;
        Network net = Networks.readNetwork(network.toString());
        net.getRoot().setRootPopSize(network.getRoot().getRootPopSize());
        for(Object node : net.dfs()) {
            NetNode mynode = (NetNode) node;
            if(mynode.getName().equals("")) {
                mynode.setName("II" + nameCount);
                nameCount++;
            }
        }

        if(!Utils._ESTIMATE_POP_SIZE) {
            net.getRoot().setRootPopSize(Utils._POP_SIZE_MEAN);
        }

        Double theta = null;
        if(Utils._CONST_POP_SIZE) {
            // theta = network.getRoot().getRootPopSize();
            for(Object childObj : net.bfs()) {
                NetNode child = (NetNode) childObj;
                for(Object parentObj : child.getParents()) {
                    NetNode parent = (NetNode) parentObj;
                    child.setParentSupport(parent, net.getRoot().getRootPopSize());
                }
            }
        }

        String netstring = net.toString();

        double sum = 0.0;
        double sumMono = 0.0;
        double numsites = 0.0;
        int maxLineages = 0;

        for(RPattern pattern : patterns.keySet()) {
            maxLineages = Math.max(maxLineages, pattern.sumLineages());
        }

        /*if(Utils._TAXON_MAP != null) {
            maxLineages = 0;

            for(Object leafObj : network.getLeaves()) {
                NetNode leaf = (NetNode) leafObj;
                String species = leaf.getName();
                for (String allele : Utils._TAXON_MAP.get(species)) {
                    maxLineages ;
                }

            }
        }*/

        if(Algorithms.HAS_DOMINANT_MARKER) {
            maxLineages = maxLineages * 2;
        }

        Network cloneNetwork = Networks.readNetwork(netstring);
        cloneNetwork.getRoot().setRootPopSize(network.getRoot().getRootPopSize());
        SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, BAGTRModel, theta, maxLineages);
        long start = System.currentTimeMillis();
        for(RPattern pattern : patterns.keySet()) {
            double count = patterns.get(pattern)[0];
            double correction = patterns.get(pattern)[1];

            double likelihood = 0;
            //long start = System.currentTimeMillis();
            try {
                if(count > 0.0) {
                    likelihood = run.getProbability(pattern);
                    sum += Math.log(likelihood) * count + correction;
                } else {
                    sum += correction;
                }
                if(pattern.isMonomorphic() && pattern.sumLineages() == maxLineages)
                    sumMono += likelihood;
                numsites += count;
                //System.out.println(represent + " " + likelihood  + " " + count );
            } catch(Exception e) {
                e.printStackTrace();
                System.out.println("Exceptional network" + netstring);
            }
            //double time1 = (System.currentTimeMillis()-start)/1000.0;
            //System.out.println(time1);
        }

        if(useOnlyPolymorphic) {
            sum -= numsites * Math.log(1.0 - sumMono);
        }

        // System.out.println("Time to compute likelihood for trinet " + network.toString() + " " + (System.currentTimeMillis()-start)/1000.0);

        return sum;
    }

    static public Map<NetNode,Boolean> findNodes2Update(Network network, List<String> leaves2update) {
        Map<NetNode,Boolean> node2update = new HashMap<>();
        Set<NetNode> nodes = new HashSet<>();
        for(String name : leaves2update) {
            nodes.add(network.findNode(name));
        }

        node2update.clear();
        while(!nodes.isEmpty()) {
            NetNode node = nodes.iterator().next();
            node2update.put(node, true);
            nodes.remove(node);
            for(Object parentObject : node.getParents()) {
                NetNode parent = (NetNode) parentObject;
                nodes.add(parent);
            }
        }

        return node2update;
    }

    static public double computeSNAPPLikelihoodMTC(Network network, Map<RPattern, double[]> patterns, BiAllelicGTR BAGTRModel) {
        //System.out.println(network.toString());
        int nameCount = 0;
        Network net = Networks.readNetwork(network.toString());
        net.getRoot().setRootPopSize(network.getRoot().getRootPopSize());
        for(Object node : net.dfs()) {
            NetNode mynode = (NetNode) node;
            if(mynode.getName().equals("")) {
                mynode.setName("II" + nameCount);
                nameCount++;
            }
        }

        if(!Utils._ESTIMATE_POP_SIZE) {
            net.getRoot().setRootPopSize(Utils._POP_SIZE_MEAN);
        }

        Double theta = null;
        if(Utils._CONST_POP_SIZE) {
            // theta = network.getRoot().getRootPopSize();
            for(Object childObj : net.bfs()) {
                NetNode child = (NetNode) childObj;
                for(Object parentObj : child.getParents()) {
                    NetNode parent = (NetNode) parentObj;
                    child.setParentSupport(parent, net.getRoot().getRootPopSize());
                }
            }
        }

        String netstring = net.toString();

        ExecutorService executor = Executors.newFixedThreadPool(Utils._NUM_THREADS);

        double sum = 0.0;
        DoubleAdder adder = new DoubleAdder();
        DoubleAdder sitecounter = new DoubleAdder();
        DoubleAdder monoadder = new DoubleAdder();

        int maxLineages = 0;

        for(RPattern pattern : patterns.keySet()) {
            maxLineages = Math.max(maxLineages, pattern.sumLineages());
        }

        /*if(Utils._TAXON_MAP != null) {
            maxLineages = 0;
            for (String species : Utils._TAXON_MAP.keySet()) {
                for (String allele : Utils._TAXON_MAP.get(species)) {
                    maxLineages++;
                }
            }
        }*/

        if(Algorithms.HAS_DOMINANT_MARKER) {
            maxLineages = maxLineages * 2;
        }

        R.maxLineages = maxLineages;

        QParameters Q = new QParameters(BAGTRModel, R.maxLineages, theta);
        long start = System.currentTimeMillis();

        //find best leaf
        Map<NetNode,Boolean> node2update = new HashMap<>();
        List<String> names = patterns.keySet().iterator().next().getNames();
        int bestCost = 9999;
        String bestName = null;
        for(String name : names) {
            Set<NetNode> nodes = new HashSet<>();
            Network cloneNetwork = Networks.readNetwork(netstring);
            nodes.add(cloneNetwork.findNode(name));

            node2update.clear();
            while(!nodes.isEmpty()) {
                NetNode node = nodes.iterator().next();
                node2update.put(node, true);
                nodes.remove(node);
                for(Object parentObject : node.getParents()) {
                    NetNode parent = (NetNode) parentObject;
                    nodes.add(parent);
                }
            }

            int cost = 0;
            for(NetNode node : node2update.keySet()) {
                if(node.isNetworkNode()) cost += 1;
                else cost += 1;
            }

            if(bestCost > cost) {
                bestCost = cost;
                bestName = name;
            }
        }

        final List<List<Tuple<RPattern, double[]>>> jobLists = new ArrayList<>();
        Set<RPattern> patternSet = new HashSet<>( patterns.keySet());
        Map<RPattern, List<Tuple<RPattern, double[]>>> jobMapping = new HashMap<>();

        for(RPattern pattern : patterns.keySet()) {
            if(!patternSet.contains(pattern)) continue;
            Map<String, R> partialPattern = new HashMap<>();
            for(String name : names) {
                if(bestName.equals(name)) continue;
                partialPattern.put(name, pattern.getR(name));
            }
            RPattern partialPatternRP = new RPattern(partialPattern);
            if(!jobMapping.containsKey(partialPatternRP))
                jobMapping.put(partialPatternRP, new ArrayList<>());
            jobMapping.get(partialPatternRP).add(new Tuple<RPattern, double[]>(pattern, patterns.get(pattern)));
        }

        for(RPattern partialPatternRP : jobMapping.keySet()) {
            jobLists.add(jobMapping.get(partialPatternRP));
        }
        //System.out.println("Preprocessing " + (System.currentTimeMillis()-start)/1000.0);

        for(int i = 0 ; i < jobLists.size() ; i++) {
            final int jobIndex = i;
            final String bestName1 = bestName;
            final int maxLineages0 = maxLineages;
            executor.execute(new Runnable() {
                public void run() {

                    Network cloneNetwork = null;
                    SNAPPAlgorithm run = null;
                    List<Tuple<RPattern, double[]>> jobList = jobLists.get(jobIndex);

                    Tuple<RPattern, double[]> firstJob = jobList.get(0);

                    double count = firstJob.Item2[0];
                    double correction = firstJob.Item2[1];

                    cloneNetwork = Networks.readNetwork(netstring);
                    cloneNetwork.getRoot().setRootPopSize(network.getRoot().getRootPopSize());
                    run = new SNAPPAlgorithm(cloneNetwork, Q);

                    double likelihood = 0;
                    try {
                        likelihood = run.getProbability(firstJob.Item1);
                        adder.add(Math.log(likelihood) * count + correction);
                        if(useOnlyPolymorphic) {
                            sitecounter.add(count);
                            if (firstJob.Item1.isMonomorphic() && firstJob.Item1.sumLineages() == maxLineages0)
                                monoadder.add(likelihood);
                        }
                        //System.out.println(represent + " " + likelihood  + " " + count );
                    } catch (Exception e) {
                        e.printStackTrace();
                        System.out.println("Exceptional network" + netstring);
                    }

                    if (jobList.size() <= 1) return;

                    //Map<NetNode, Boolean> node2update = findNodes2Update(cloneNetwork, Collections.singletonList(bestName1));
                    Set<NetNode> nodes = new HashSet<>();
                    nodes.add(cloneNetwork.findNode(bestName1));

                    Map<NetNode,Boolean> node2update = new HashMap<>();
                    while(!nodes.isEmpty()) {
                        NetNode node = nodes.iterator().next();
                        node2update.put(node, true);
                        nodes.remove(node);
                        for(Object parentObject : node.getParents()) {
                            NetNode parent = (NetNode) parentObject;
                            nodes.add(parent);
                        }
                    }

                    for (int j = 1; j < jobList.size(); j++) {
                        Tuple<RPattern, double[]> nextJob = jobList.get(j);
                        likelihood = run.updateProbability(0, nextJob.Item1, null, node2update);
                        count = nextJob.Item2[0];
                        correction = nextJob.Item2[1];
                        adder.add(Math.log(likelihood) * count + correction);
                        if(useOnlyPolymorphic) {
                            sitecounter.add(count);
                            if (nextJob.Item1.isMonomorphic() && nextJob.Item1.sumLineages() == maxLineages0)
                                monoadder.add(likelihood);
                        }
                    }

                }
            });
        }
        try {
            executor.shutdown();
            if(timeSavingMode) {
                boolean finished = executor.awaitTermination(600, TimeUnit.SECONDS);
                if(!finished) {
                    executor.shutdownNow();
                    return Double.NEGATIVE_INFINITY;
                }
            } else {
                while (!executor.awaitTermination(1000, TimeUnit.SECONDS)) {
                    System.out.println("Super long: " + network.toString());
                }
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        sum = adder.sumThenReset();
        if(useOnlyPolymorphic) {
            sum -= sitecounter.sumThenReset() * Math.log(1.0 - monoadder.sumThenReset());
        }
        //System.out.println("Finished " + (System.currentTimeMillis()-start)/1000.0);
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

        MarkerSeq aln = new MarkerSeq(colorMap);
        List<MarkerSeq> alns = new ArrayList<>();
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

                            MarkerSeq aln = new MarkerSeq(colorMap);
                            List<MarkerSeq> alns = new ArrayList<>();
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
