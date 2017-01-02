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

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
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
                mynode.setName("I" + nameCount);
                nameCount++;
            }
        }
        String netstring = net.toString();


        ExecutorService executor = Executors.newFixedThreadPool(Utils._NUM_THREADS);

        for(Alignment alg : alignments) {
            List<String> names = alg.getTaxaNames();

            for(Integer represent : alg.getCache().keySet()) {
                executor.execute(new Runnable() {
                    public void run() {
                        int cur = represent;
                        Integer count = alg.getCache().get(cur);
                        Map<String, Character> colorMap = new HashMap<String, Character>();
                        for(int i = names.size() - 1 ; i >= 0 ; i--) {
                            Character site = cur % 2 == 1 ? '1' : '0';
                            colorMap.put(names.get(i), site);
                            cur /= 2;
                        }
                        OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);
                        Network cloneNetwork = Networks.readNetwork(netstring);
                        cloneNetwork.getRoot().setRootPopSize(network.getRoot().getRootPopSize());
                        SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, BAGTRModel, 1);
                        double likelihood = 0;
                        try {
                            likelihood = run.getProbability(converter);
                        } catch(Exception e) {
                            e.printStackTrace();
                            System.out.println("Exceptional network: " + netstring);
                        }
                        double logL = Math.log(likelihood);
                        adder.add(logL * count);
                    }
                });
                //executor.shutdown();




            }

            try {
                executor.shutdown();
                executor.awaitTermination(1000, TimeUnit.SECONDS);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            sum += adder.sumThenReset();
        }

        return sum;

    }
}
