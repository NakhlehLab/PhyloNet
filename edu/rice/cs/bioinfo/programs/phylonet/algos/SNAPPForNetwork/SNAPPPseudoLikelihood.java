package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;

import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution.SNAPPLikelihood;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.RateModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.SuperNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.DoubleAdder;

import static edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.Algorithms.HAS_DOMINANT_MARKER;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 11/15/17
 * Time: 10:24 AM
 * To change this template use File | Settings | File Templates.
 */
public class SNAPPPseudoLikelihood {

    private QParameters Q;
    private Network speciesNetwork;
    private Map<String, String> allele2species;
    private List<Network> subNetworks;
    private RateModel rateModel;
    private List<Tuple3<String, String, String>> triplets;
    private boolean speciesInTriplets = false;
    public Map<Tuple3<String, String, String>, Map<RPattern, double[]>> triplets2patterns;
    public Map<Tuple3<String, String, String>, Double> cache = new HashMap<>();

    // by alleles
    private void initialize1(Map<String, String> alleles2species, List<Alignment> alignments, boolean diploid) {
        List<String> taxa = alignments.get(0).getTaxaNames();
        triplets = new ArrayList<>();
        triplets2patterns = new HashMap<>();

        for(int i = 0 ; i < taxa.size() ; i++) {
            for(int j = i + 1 ; j < taxa.size() ; j++) {
                for(int k = j + 1 ; k < taxa.size() ; k++) {
                    Tuple3<String, String, String> triplet = new Tuple3<>(taxa.get(i), taxa.get(j), taxa.get(k));
                    triplets.add(triplet);

                    Map<String, String> sequence = new HashMap<>();
                    sequence.put(triplet.Item1, alignments.get(0).getAlignment().get(triplet.Item1));
                    sequence.put(triplet.Item2, alignments.get(0).getAlignment().get(triplet.Item2));
                    sequence.put(triplet.Item3, alignments.get(0).getAlignment().get(triplet.Item3));
                    List<Alignment> alnwarp = new ArrayList<>();
                    alnwarp.add(new Alignment(sequence));

                    Map<RPattern, double[]> patterns;
                    if(diploid && !HAS_DOMINANT_MARKER) {
                        patterns = SNAPPLikelihood.diploidSequenceToPatterns(alleles2species, alnwarp);
                    } else {
                        patterns = SNAPPLikelihood.haploidSequenceToPatterns(alleles2species, alnwarp);
                    }
                    triplets2patterns.put(triplet, patterns);

                    System.out.println(triplet.Item1 + " " + triplet.Item2 + " " + triplet.Item3);

                }
            }
        }

        speciesInTriplets = false;
    }

    // by species
    void initialize2(Map<String, String> alleles2species, List<Alignment> alignments, boolean diploid) {
        triplets = new ArrayList<>();
        triplets2patterns = new HashMap<>();

        List<String> speciesList = new ArrayList<>();
        for(String name : alleles2species.values()) {
            if(!speciesList.contains(name)) {
                speciesList.add(name);
            }
        }
        Collections.sort(speciesList);

        Map<String, List<String>> species2alleles = new HashMap<>();
        for(String allele : alleles2species.keySet()) {
            String species = alleles2species.get(allele);
            if(!species2alleles.containsKey(species)) {
                species2alleles.put(species, new ArrayList<>());
            }
            species2alleles.get(species).add(allele);
        }

        for(int i = 0 ; i < speciesList.size() ; i++) {
            for (int j = i + 1; j < speciesList.size(); j++) {
                for (int k = j + 1; k < speciesList.size(); k++) {
                    Tuple3<String, String, String> triplet = new Tuple3<>(speciesList.get(i), speciesList.get(j), speciesList.get(k));
                    triplets.add(triplet);
                    Map<String, String> sequence = new HashMap<>();
                    for(String allele : species2alleles.get(speciesList.get(i))) {
                        sequence.put(allele, alignments.get(0).getAlignment().get(allele));
                    }
                    for(String allele : species2alleles.get(speciesList.get(j))) {
                        sequence.put(allele, alignments.get(0).getAlignment().get(allele));
                    }
                    for(String allele : species2alleles.get(speciesList.get(k))) {
                        sequence.put(allele, alignments.get(0).getAlignment().get(allele));
                    }

                    List<Alignment> alnwarp = new ArrayList<>();
                    alnwarp.add(new Alignment(sequence));

                    Map<RPattern, double[]> patterns;
                    if(diploid && !HAS_DOMINANT_MARKER) {
                        patterns = SNAPPLikelihood.diploidSequenceToPatterns(alleles2species, alnwarp);
                    } else {
                        patterns = SNAPPLikelihood.haploidSequenceToPatterns(alleles2species, alnwarp);
                    }
                    triplets2patterns.put(triplet, patterns);

                    System.out.println(triplet.Item1 + " " + triplet.Item2 + " " + triplet.Item3);
                }
            }
        }

        speciesInTriplets = true;
    }

    public SNAPPPseudoLikelihood(Map<String, String> alleles2species, List<Alignment> alignments, boolean diploid) {
        initialize2(alleles2species, alignments, diploid);
    }

    public SNAPPPseudoLikelihood(Network theSpeciesNetwork, RateModel rModel, Double theta) {
        SuperNetwork superNetwork = new SuperNetwork(new ArrayList<>());
        rateModel = rModel;
        subNetworks = superNetwork.genAllSubNetworks(theSpeciesNetwork, 3);
    }

    public double getProbability(RPattern pattern) {
        double product = 1.0;

        for(Network subNetwork : subNetworks) {
            Double theta = 0.036;
            Network cloneNetwork = Networks.readNetwork(subNetwork.toString());
            cloneNetwork.getRoot().setRootPopSize(subNetwork.getRoot().getRootPopSize());
            SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, rateModel, theta);
            double likelihood = run.getProbability(pattern);
            product *= likelihood;
            //System.out.println(subNetwork);
            //System.out.println(likelihood);
        }

        return product;
    }

    private Map<Tuple3<String, String, String>, Network> buildSubNetworks(Network network, Map<String, String> alleles2species) {
        SuperNetwork superNetwork = new SuperNetwork(new ArrayList<>());
        Map<Tuple3<String, String, String>, Network> triplets2subnetworks = new HashMap<>();
        for(Tuple3<String, String, String> triplet : triplets) {
            List<String> leaves = new ArrayList<>();
            if(alleles2species != null && !speciesInTriplets) {
                leaves.add(alleles2species.get(triplet.Item1));
                if(!leaves.contains(alleles2species.get(triplet.Item2))) {
                    leaves.add(alleles2species.get(triplet.Item2));
                }
                if(!leaves.contains(alleles2species.get(triplet.Item3))) {
                    leaves.add(alleles2species.get(triplet.Item3));
                }
            } else {
                leaves.add(triplet.Item1);
                leaves.add(triplet.Item2);
                leaves.add(triplet.Item3);
            }
            Network subNetwork = superNetwork.getSubNetwork(network, leaves, false, false);
            Networks.autoLabelNodes(subNetwork);
            triplets2subnetworks.put(triplet, subNetwork);
        }
        return triplets2subnetworks;
    }

    public double computeSNAPPPseudoLogLikelihood(Network network, Map<String, String> alleles2species, BiAllelicGTR BAGTRModel) {
        Map<Tuple3<String, String, String>, Network> triplets2subnetworks = buildSubNetworks(network, alleles2species);

        double sum = 0.0;
        //Utils._NUM_THREADS = 1;
        //Collections.shuffle(triplets);
        for(Tuple3<String, String, String> triplet : triplets) {
            if(HAS_DOMINANT_MARKER) {
                R.maxLineages = triplets2patterns.get(triplet).keySet().iterator().next().sumLineages() * 2;
            } else {
                R.maxLineages = triplets2patterns.get(triplet).keySet().iterator().next().sumLineages();
            }

            double logLikelihood = SNAPPLikelihood.computeSNAPPLikelihood(triplets2subnetworks.get(triplet), triplets2patterns.get(triplet), BAGTRModel);
            //System.out.println("Triplet " + triplet.toString() + " likelihood " + logLikelihood + " network " + triplets2subnetworks.get(triplet));
            //System.out.println(logLikelihood);
            sum += logLikelihood;
        }

        return sum;
    }

    public double computeSNAPPPseudoLogLikelihoodMT(Network network, Map<String, String> alleles2species, BiAllelicGTR BAGTRModel) {
        Map<Tuple3<String, String, String>, Network> triplets2subnetworks = buildSubNetworks(network, alleles2species);

        ExecutorService executor = Executors.newFixedThreadPool(Utils._NUM_THREADS);
        DoubleAdder adder = new DoubleAdder();

        for(Tuple3<String, String, String> triplet : triplets) {
            /*if(cache.containsKey(triplet)) {
                adder.add(cache.get(triplet));
            } else {
                double logLikelihood = SNAPPLikelihood.computeSNAPPLikelihoodST(triplets2subnetworks.get(triplet), triplets2patterns.get(triplet), BAGTRModel);
                adder.add(logLikelihood);
                // cache.put(triplet, logLikelihood);


            }*/

            executor.execute(new Runnable() {
                public void run() {
                    double logLikelihood = SNAPPLikelihood.computeSNAPPLikelihoodST(triplets2subnetworks.get(triplet), triplets2patterns.get(triplet), BAGTRModel);
                    adder.add(logLikelihood);
                }
            });

        }

        try {
            executor.shutdown();
            executor.awaitTermination(1000, TimeUnit.SECONDS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        SNAPPLikelihood.workloadCounter.increment();

        double sum = adder.sumThenReset();
        return sum;
    }


}
