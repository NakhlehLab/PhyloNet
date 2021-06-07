package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution.SNAPPLikelihood;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.MarkerSeq;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimSNPInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.RateModel;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 2019-04-21
 * Time: 18:48
 * To change this template use File | Settings | File Templates.
 */
public class AlgorithmMULTree {

    public static <T1, T2> Iterable<Map<T1, T2>> AllMappings(List<T1> from, List<T2> to) {
        return new Iterable<Map<T1, T2>>()
        {
            @Override
            public Iterator<Map<T1, T2>> iterator()
            {
                return new Iterator<Map<T1, T2>>()
                {

                    int index[] = null;
                    int size = from.size();

                    @Override
                    public void remove()
                    {
                        throw new UnsupportedOperationException();
                    }


                    @Override
                    public boolean hasNext()
                    {
                        if(index == null) return true;
                        for(int i = 0; i < index.length ; i++) {
                            if(index[i] != to.size() - 1) return true;
                        }
                        return false;
                    }

                    @Override
                    public Map<T1, T2> next()
                    {
                        if (index == null) {
                            index = new int[size];
                            for(int i = 0 ; i < size; i++)
                                index[i] = 0;
                        }
                        else
                        {
                            tryAdvance();
                        }

                        Map<T1, T2> temp = new HashMap<>();
                        for(int i = 0 ; i < size ; i++)
                            temp.put(from.get(i), to.get(index[i]));

                        return temp;
                    }

                    private void tryAdvance()
                    {
                        int j = size - 1;
                        index[j]++;
                        while(index[j] == to.size()) {
                            index[j] = 0;
                            j--;
                            index[j]++;
                        }

                    }
                };
            }
        };

    }

    // Not working!
    public static void compute(List<MarkerSeq> alns, Network theSpeciesNetwork, RateModel rModel, Double theta, int maxLineages, Map<String, List<String>> species2alleles, Map<RPattern, double[]> patterns) {
        long startTime = System.currentTimeMillis();

        Tuple<Network, Map<String, Double>> multreetuple = MULTreeUtils.GetMULTree(theSpeciesNetwork);
        Network multree = multreetuple.Item1;
        Map<String, Double> probs = multreetuple.Item2;
        System.out.println(Networks.getFullString(multree));


        List<String> bb = new ArrayList<>();

        for(Object leafObj : multree.getLeaves()) {
            NetNode leaf = (NetNode) leafObj;
            if(leaf.getName().startsWith("B")) {
                bb.add(leaf.getName());
            }
        }

        List<String> alleleOfB = species2alleles.get("B");
        double sum = 0;

        int index = 0;
        for(Map<String, String> mapping : AllMappings(alleleOfB, bb)) {
            Network cloneNetwork = multree.clone();
            SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, rModel, theta, maxLineages);


            Map<String, List<String>> newS2A = new HashMap<>();
            for(String species : species2alleles.keySet()) {
                newS2A.put(species, new ArrayList<>(species2alleles.get(species)));
            }
            newS2A.remove("B");


            double prob = 1.0;
            for(String allele : mapping.keySet()) {
                String taxon = mapping.get(allele);
                if(!newS2A.containsKey(taxon)) {
                    newS2A.put(taxon, new ArrayList<>());
                }
                newS2A.get(taxon).add(allele);
                prob *= probs.get(mapping.get(allele));
            }

            Map<String, String> a2s = new HashMap<>();
            for(String species : newS2A.keySet()) {
                for(String allele : newS2A.get(species)) {
                    a2s.put(allele, species);
                }
            }
            patterns = SNAPPLikelihood.haploidSequenceToPatterns(a2s, alns);

            for(RPattern pattern : patterns.keySet()) {
                double likelihood = run.getProbability(pattern);

                System.out.println(likelihood + " " + prob);
                sum += likelihood * prob;
            }
            index++;
        }

        System.out.println( "Done " + sum);

        System.out.println(String.format("Total elapsed time : %2.5f s\n",
                (double) (System.currentTimeMillis() - startTime) / 1000.0));

    }

    public static void main(String args[]) {
        double pi0 = 0.5;
        double pi1 = 1- pi0;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {pi0, pi1}, new double[] {1.0/ (2.0 * pi0)});


        double alpha1 = 0.05;
        double alpha2 = 0.04;
        double alpha3 = 0.06;
        double alpha4 = 0.08;
        double a = 0;
        double b = 0;
        double c = 0;
        double theta = 0.01;
        double gamma = 0.5;

        String netstring = String.format("[%f](((B:%f)I3#H1:%f::%f,A:%f)I1:%f,(C:%f,I3#H1:%f::%f)I2:%f)I0;", theta, alpha2, alpha1-alpha2, gamma, alpha1, alpha4-alpha1, alpha3, alpha3-alpha2, 1.0-gamma, alpha4-alpha3);
        //netstring = "[0.01]((((B:0.01:0.01)#H1:0.05:0.01:0.5)#H2:0.01:0.01:0.5,A:0.07:0.01):0.03:0.01,((C:0.02:0.01,#H1:0.01:0.01:0.5):0.06:0.01,#H2:0.02:0.01:0.5):0.02:0.01);";

        Network<NetNodeInfo> trueNetwork = Networks.readNetwork("I0;");
        trueNetwork = Networks.readNetworkWithRootPop(netstring);
        //Networks.autoLabelNodes(trueNetwork);
        for(Object nodeObj : trueNetwork.dfs()) {
            NetNode node = (NetNode) nodeObj;
            for(Object parentObj : node.getParents()) {
                NetNode parent = (NetNode) parentObj;
                node.setParentSupport(parent, trueNetwork.getRoot().getRootPopSize());
            }
        }

        System.out.println(Networks.getFullString(trueNetwork));
        Map<String, List<String>> species2alleles = new HashMap<>();
        species2alleles.put("A", new ArrayList<>());
        species2alleles.put("B", new ArrayList<>());
        species2alleles.put("C", new ArrayList<>());
        species2alleles.get("A").add("a");
        species2alleles.get("B").add("b0");
        species2alleles.get("B").add("b1");
//        species2alleles.get("B").add("b2");
//        species2alleles.get("B").add("b3");
//        species2alleles.get("B").add("b4");
//        species2alleles.get("B").add("b5");
//        species2alleles.get("B").add("b6");
//        species2alleles.get("B").add("b7");
//        species2alleles.get("B").add("b8");
//        species2alleles.get("B").add("b9");
//        species2alleles.get("B").add("b10");
//        species2alleles.get("B").add("b11");
        species2alleles.get("C").add("c");

        Map<String, String> allele2species = new HashMap<>();
        allele2species.put("a", "A");
        allele2species.put("b0", "B");
        allele2species.put("b1", "B");
//        allele2species.put("b2", "B");
//        allele2species.put("b3", "B");
//        allele2species.put("b4", "B");
//        allele2species.put("b5", "B");
//        allele2species.put("b6", "B");
//        allele2species.put("b7", "B");
//        allele2species.put("b8", "B");
//        allele2species.put("b9", "B");
//        allele2species.put("b10", "B");
//        allele2species.put("b11", "B");
        allele2species.put("c", "C");

        int numSites = 100;
        boolean useOnlyPolymorphic = false;
        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, 12345678L);
        Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, species2alleles, numSites, !useOnlyPolymorphic);
        onesnp.clear();
        onesnp.put("a", "0");
        onesnp.put("b0", "0");
        onesnp.put("b1", "0");
//        onesnp.put("b2", "0");
//        onesnp.put("b3", "0");
//        onesnp.put("b4", "0");
//        onesnp.put("b5", "0");
//        onesnp.put("b6", "0");
//        onesnp.put("b7", "0");
//        onesnp.put("b8", "0");
//        onesnp.put("b9", "0");
//        onesnp.put("b10", "0");
//        onesnp.put("b11", "0");
        onesnp.put("c", "0");


        List<MarkerSeq> alns = new ArrayList<>();
        MarkerSeq aln = new MarkerSeq(onesnp);
        alns.add(aln);
        aln._RPatterns = SNAPPLikelihood.haploidSequenceToPatterns(allele2species, alns);

        int maxLineages = 0;

        for(RPattern pattern : aln._RPatterns.keySet()) {
            maxLineages = Math.max(maxLineages, pattern.sumLineages());
        }
        R.maxLineages = maxLineages;


        long startTime = System.currentTimeMillis();

        Algorithms.targetSplittingIndices = new ArrayList<>();
        Algorithms.targetSplittingIndices.add(Arrays.asList(0));
        Algorithms.targetSplittingIndices.add(Arrays.asList(2));
        Algorithms.targetSplittingIndices.add(Arrays.asList(4));
        Algorithms.targetSplittingIndices.add(Arrays.asList(7));
        Algorithms.targetSplittingIndices.add(Arrays.asList(11));

        double sum = SNAPPLikelihood.computeSNAPPLikelihood(trueNetwork, aln._RPatterns, BAGTRModel);

        System.out.println(Math.exp(sum));
        System.out.println(String.format("Total elapsed time : %2.5f s\n",
                (double) (System.currentTimeMillis() - startTime) / 1000.0));
        System.out.println();



        compute(alns, trueNetwork, BAGTRModel, theta, maxLineages, species2alleles, aln._RPatterns);
    }

}
