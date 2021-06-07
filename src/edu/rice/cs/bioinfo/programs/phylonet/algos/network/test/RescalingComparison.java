package edu.rice.cs.bioinfo.programs.phylonet.algos.network.test;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityIntegrated;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.MDCOnNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.jblas.util.Random;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by hunter on 9/18/18.
 *
 * Intention: compare MDC and NCM species tree "likelihood" (or parsimony) distributions
 */
public class RescalingComparison {
    public static void main(String[] args) {
        String[] taxa = {"A", "B", "C", "D"};//, "D"};
        // This method generates trees that are somehow malformed. Can be corrected by reading from their Newick string.
        List<Tree> badTrees = Trees.generateAllBinaryTrees(taxa);

        List<Tree> treesOnTaxa = new ArrayList();
        List<Network> speciesNetworks = new ArrayList();
        int treeNumber = 0;
        for (Tree t : badTrees) {
            ((MutableTree) t).getRoot().setName(Integer.toString(treeNumber++));
            treesOnTaxa.add(Trees.readTree(t.toNewick()));
            speciesNetworks.add(Networks.readNetwork(t.toNewick()));
        }

//        for (Network n : speciesNetworks) {
//            System.out.println("\ngenerating virtual gene trees from network: " + n);
//            GeneTreeProbabilityIntegrated gtp = new GeneTreeProbabilityIntegrated();
//            gtp.setBranchLengthExponentialPrior(1, Double.POSITIVE_INFINITY);
//            List<Double> probabilities = gtp.calculateGTDistribution(n, treesOnTaxa, null, false);
//            System.out.println("probabilities: " + probabilities);
//            List<MutableTuple<Tree, Double>> dataDistribution = new ArrayList<>();
//            for (int i = 0; i < treesOnTaxa.size(); i++) {
//                dataDistribution.add(new MutableTuple(treesOnTaxa.get(i), probabilities.get(i)));
//            }
//
//            compare(speciesNetworks, dataDistribution);
//        }

        int randomSettings = 1;
        for (int i = 0; i < randomSettings; i++) {
            List<MutableTuple<Tree, Double>> dataDistribution = new ArrayList<>();

            List<Double> unnormalizedProbs = new ArrayList();
            Double normalization = 0.0;
            Double nextDouble;

            double[] possibleProbs = {1, 2, 4, 8, 16};

            int numAdds = 3;
//            double[] explicitProbs = new double[treesOnTaxa.size()];
//            for (int unused = 0; unused < numAdds; unused++) {
//                int explicitAdd = Random.nextInt(explicitProbs.length);
//                explicitProbs[explicitAdd] += 1;
//            }
            double[] explicitProbs = {0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};
            int explicitIndex = 0;

            for (Tree t: treesOnTaxa) {
//                nextDouble = possibleProbs[Random.nextInt(possibleProbs.length)]; //Random.nextDouble();
                nextDouble = explicitProbs[explicitIndex++];
                unnormalizedProbs.add(nextDouble);
                normalization += nextDouble;
            }
            for (int j = 0; j < treesOnTaxa.size(); j++) {
                if (unnormalizedProbs.get(j) > 0) {
                    dataDistribution.add(new MutableTuple(treesOnTaxa.get(j), unnormalizedProbs.get(j) / normalization));
                }
            }
//            System.out.println("\n" + dataDistribution);
            compare(speciesNetworks, dataDistribution);
        }
    }

    /**
     * Compare the likelihood and parsimony distributions over the species networks, given the data as observed
     * @param speciesNetworks
     * @param dataDistribution
     */
    private static void compare(List<Network> speciesNetworks, List<MutableTuple<Tree, Double>> dataDistribution) {
        List<MutableTuple<Network, Double>> likelihoods = getIntegratedLogLikelihoodDistribution(speciesNetworks, dataDistribution);
        List<MutableTuple<Network, Double>> parsimonies = getParsimonyDistribution(speciesNetworks, dataDistribution);
//        System.out.println("likelihoods: " + likelihoods);
//        System.out.println("parsimonies: " + parsimonies);

        Double maxLogLikelihood = Double.NEGATIVE_INFINITY;
        Network maxLogLikelihoodNetwork = null;
        for (MutableTuple<Network, Double> t: likelihoods) {
//                maxLogLikelihood = Math.max(maxLogLikelihood, t.Item2);
            if (t.Item2 > maxLogLikelihood) {
                maxLogLikelihood = t.Item2;
                maxLogLikelihoodNetwork = t.Item1;
            }
        }
        Double minParsimony = Double.POSITIVE_INFINITY;
        Network minParsimonyNetwork = null;
        for (MutableTuple<Network, Double> t: parsimonies) {
//                minParsimony = Math.min(minParsimony, t.Item2);
            if (t.Item2 < minParsimony) {
                minParsimony = t.Item2;
                minParsimonyNetwork = t.Item1;
            }
        }
//        System.out.println("best likelihood: " + maxLogLikelihoodNetwork);
//        System.out.println("best parsimony: " + minParsimonyNetwork);
        if (! Networks.hasTheSameTopology(maxLogLikelihoodNetwork, minParsimonyNetwork)) {
            System.out.println("\nmismatch!\n" + maxLogLikelihoodNetwork + " from likelihood, " + minParsimonyNetwork + " from parsimony");
            System.out.println("data: " + dataDistribution);
            System.out.println("likelihoods: " + likelihoods);
            System.out.println("parsimonies: " + parsimonies);
        }
    }

    private static Double getIntegratedLogLikelihood(Network network, List<MutableTuple<Tree, Double>> treeDistribution) {
        List<Tree> gts_stripped = new ArrayList<>();
        for (MutableTuple<Tree, Double> tup: treeDistribution) {
            gts_stripped.add(tup.Item1);
        }

        GeneTreeProbabilityIntegrated gtp = new GeneTreeProbabilityIntegrated();
        gtp.setBranchLengthExponentialPrior(1, Double.POSITIVE_INFINITY);
        double prob = 0;
        List<Double> geneTreeProbabilities = gtp.calculateGTDistribution(network, gts_stripped, null, true);
        for (int i = 0; i < geneTreeProbabilities.size(); i++) {
            prob += Math.log(geneTreeProbabilities.get(i)) * treeDistribution.get(i).Item2; // multiply log probability by weight factor
//                    System.out.println("adding log of " + geneTreeProbabilities.get(i) + " times " + gts.get(i).Item2);
        }
        return prob;
    }

    private static List<MutableTuple<Network, Double>> getIntegratedLogLikelihoodDistribution(List<Network> networks, List<MutableTuple<Tree, Double>> treeDistribution) {
        List distribution = new ArrayList();
        for (Network n : networks) {
            Double prob = getIntegratedLogLikelihood(n, treeDistribution);
            distribution.add(new MutableTuple(n, prob));
        }
        return distribution;
    }

    private static List<MutableTuple<Network, Double>> getParsimonyDistribution(List<Network> networks, List<MutableTuple<Tree, Double>> treeDistribution) {
        List distribution = new ArrayList();
        for (Network n : networks) {
            MDCOnNetwork mdc = new MDCOnNetwork();
            List<Integer> deepCoalescentCount = mdc.countExtraCoal(n, treeDistribution, null);
            MutableTuple<Network, Double> parsimony = new MutableTuple(n, 0.0);
            for (int i = 0; i < treeDistribution.size(); i++) {
                parsimony.Item2 += treeDistribution.get(i).Item2 * deepCoalescentCount.get(i).intValue();
            }
            distribution.add(parsimony);
        }
        return distribution;
    }
}
