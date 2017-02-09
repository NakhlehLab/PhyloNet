package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.core;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution.SNAPPLikelihood;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.SNAPPAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.clustering.DataGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.InferMLNetworkFromSequences;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimSNPInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.OneNucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.*;
import java.lang.reflect.Array;
import java.util.*;
import java.util.concurrent.Exchanger;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import javax.rmi.CORBA.Util;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.DocumentBuilder;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import org.w3c.dom.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 12/12/16
 * Time: 2:16 PM
 * To change this template use File | Settings | File Templates.
 */
public class Test {
    public static void main(String[] args) {
        //checkSNAPPSimulationResults();
        //testWithSNAPPSimulation(args);
        //testMCMC1();
        //testMultiThread();
        //testBranchlength();
        //testSecondData(args);
        testMosquitoData(args);
        //testWithDingqiaoNetwork(args);
        //checkDingqiaoNetworkResults();
        //generateSNPdata();
        //generateSecondSNPdata();
        //processMosquitoData();
        //verifyLikelihood();

    }

    public static void testMultiThread() {
        ExecutorService executor = Executors.newFixedThreadPool(1);
        executor.execute(new Runnable() {
            public void run() {
                for(int i = 0 ; i < 10 ; i++) {
                    System.out.println("A");
                    try {
                        Thread.sleep(1000);
                    } catch (InterruptedException e) {

                    }
                }
            }
        });
        executor.execute(new Runnable() {
            public void run() {
                for(int i = 0 ; i < 10 ; i++) {
                    System.out.println("B");
                    try {
                        Thread.sleep(1000);
                    } catch (InterruptedException e) {

                    }
                }
            }
        });
        executor.execute(new Runnable() {
            public void run() {
                for(int i = 0 ; i < 10 ; i++) {
                    System.out.println("C");
                    try {
                        Thread.sleep(1000);
                    } catch (InterruptedException e) {

                    }
                }
            }
        });
        try {
            executor.shutdown();
            executor.awaitTermination(1000, TimeUnit.SECONDS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    public static void testBranchlength() {
        SNAPPLikelihood.useOnlyPolymorphic = true;
        Network trueNetwork = Networks.readNetwork("(((A:0.1,B:0.1)I1:0.1,C:0.2)I2:0.1,D:0.3)I0;");
        //Network trueNetwork = Networks.readNetwork("(((A:10.1,B:10.1)I1:10.1,C:20.2)I2:10.1,D:30.3)I0;");

        trueNetwork.getRoot().setRootPopSize(0.04);

        Map<String, String> sequence2 = new HashMap<>();
        sequence2.put("A", "1111111111");
        sequence2.put("B", "0000000000");
        sequence2.put("C", "0000000000");
        sequence2.put("D", "0000000000");

        List<Map<String, String>> snpdata = new ArrayList<>();
        snpdata.add(sequence2);

        List<Alignment> alns = new ArrayList<>();
        for(Map<String, String> input : snpdata) {

            System.out.println();
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

            aln.setCache(cache);
            alns.add(aln);
        }

        double [] pi =  new double[2];
        double [] rate = new double[1];


        pi[1] = 0.5;
        pi[0] = 1 - pi[1];

        rate[0] = 1 / (2*pi[0]);

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(pi, rate);

        Network cloneNetwork = Networks.readNetwork(trueNetwork.toString());
        cloneNetwork.getRoot().setRootPopSize(trueNetwork.getRoot().getRootPopSize());
        System.out.println("True Likelihood = " + SNAPPLikelihood.computeSNAPPLikelihood(cloneNetwork, null, alns, BAGTRModel));
    }

    public static void testMCMC() {
        Utils._NUM_THREADS = 2;
        Utils._CHAIN_LEN = 20000;
        Utils._BURNIN_LEN = 10000;
        Utils._SAMPLE_FREQUENCY = 500;
        //Utils._NET_MAX_RETI = 2;
        SNAPPLikelihood.useOnlyPolymorphic = true;


        Utils._MC3_CHAINS = new ArrayList<>();
        //Utils._MC3_CHAINS.add(1.0);
        //Utils._MC3_CHAINS.add(2.0);
        //Utils._MC3_CHAINS.add(4.0);
        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {0.5, 0.5}, new double[] {1.0});

        double gamma = 1 - 0.3;

        double y = 1.0;
        double x = 1000.0;

        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, 12345678L + 100000);
        DataGenerator dataGenerator = new DataGenerator();
        Network trueNetwork = Networks.readNetwork(dataGenerator.getNewNetwork1().toString());
        //System.out.println("True Network: " + trueNetwork.toString());

        //SimGTInNetworkByMS gtsimulator = new SimGTInNetworkByMS();
        //SimGTInNetwork gtsimulator = new SimGTInNetwork();
        //List<Tree> gts = gtsimulator.generateGTs(trueNetwork, null, 10);


        trueNetwork = Networks.readNetwork("(((((A:0.7)I6#H1:1.3::0.8,Q:2.0)I4:1.0,L:3.0)I3:1.0,R:4.0)I2:1.0,(G:2.0,(I6#H1:0.7::0.2,C:1.4)I5:0.6)I1:3.0)I0;");
        //trueNetwork = Networks.readNetwork("(((((A:0.7)I6#H1:1.3::1.0,Q:2.0)I4:1.0,L:3.0)I3:1.0,R:4.0)I2:1.0,(G:2.0,(I6#H1:0.7::0.0,C:1.4)I5:0.6)I1:3.0)I0;");
        //trueNetwork = Networks.readNetwork("((d:2.6,((c:1.0,b:1.0)I4:1.5)I3#H1:0.1::0.7)I2:0.5,(a:2.6,I3#H1:0.1::0.30000000000000004)I1:0.5)I0;");
        //trueNetwork = Networks.readNetwork("((((b:1.0,c:1.0)I4:" + y + ")I3#H1:0.1::" + (1 - gamma) + ",a:" + (1.1 + y) + ")I1:" + x + ",(I3#H1:0.1::" + gamma + ",d:" + (1.1 + y) + ")I2:"+ x +")I0;");
        trueNetwork = Networks.readNetwork("(((((Q:2.0,A:2.0)I4:1.0,L:3.0)I3:0.5)I8#H1:0.5::0.7,R:4.0)I2:1.0,(I8#H1:0.5::0.3,(G:2.0,C:2.0)I1:2.0)I7:1.0)I0;");
        //trueNetwork = Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,A:1.0)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,(G:1.0,C:1.0)I1:4.0)I0;");
        //trueNetwork = Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,(A:0.5)I6#H2:0.5::0.8)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,((I6#H2:0.5::0.2,C:1.0)I5:1.0,G:2.0)I1:3.0)I0;");

        double constTheta = 0.006;
        for(Object nodeObject : Networks.postTraversal(trueNetwork)) {
            NetNode node = (NetNode) nodeObject;
            for(Object parentObject :  node.getParents()) {
                NetNode parent = (NetNode) parentObject;
                node.setParentSupport(parent, constTheta);
                node.setParentDistance(parent, node.getParentDistance(parent) * constTheta / 2.0);
            }
        }
        trueNetwork.getRoot().setRootPopSize(constTheta);

        //trueNetwork = Networks.readNetwork("((((C:0.003,(B:0.002)I7#H1:0.001::0.2)I4:0.002,(I7#H1:0.002::0.8)I6#H2:0.001::0.2)I3:0.002,(I6#H2:0.002::0.8)I5#H3:0.001::0.2)I2:0.003,(I5#H3:0.002::0.8,A:0.008)I1:0.002)I0;\n");
        //trueNetwork = Networks.readNetwork("((((C:0.07,(B:0.06)I7#H1:0.01::0.2)I4:0.02,(I7#H1:0.02::0.8)I6#H2:0.01::0.2)I3:0.02,(I6#H2:0.02::0.8)I5#H3:0.01::0.2)I2:0.07,(I5#H3:0.06::0.8,A:0.16)I1:0.02)I0;");
        //trueNetwork.getRoot().setRootPopSize(0.04);

        Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, 100000, false);
        // List<Map<String, String>> snpdata = simulator.generateGTs(Networks.readNetwork("((B:0.5)X#H1:1.5::0.5,((X#H1:0.5::0.5,A:1)n1:0.5,C:1.5)n2:0.5)root;"), null, 100, "/scratch/jz55/Luay/msdir/ms");
        List<Map<String, String>> snpdata = new ArrayList<>();
        snpdata.add(onesnp);
        List<Alignment> alns = new ArrayList<>();
        for(Map<String, String> input : snpdata) {
            for(String allele : input.keySet()) {
                System.out.println(allele + " " + input.get(allele));
            }
            System.out.println();
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

            aln.setCache(cache);
            alns.add(aln);
        }

        double [] pi =  new double[2];
        double [] rate = new double[1];

        int totalSites = 0;
        for(Alignment alg : alns) {
            int count = 0;
            for(String s : alg.getAlignment().values()) {
                for(int i = 0 ; i < s.length() ; i++) {
                    count += s.charAt(i) == '1' ? 1 : 0;
                }
            }
            totalSites += alg.getSiteCount() * alg.getTaxonSize();
            pi[1] += 1.0 * count;
        }
        pi[1] /= 1.0 * totalSites;
        pi[0] = 1 - pi[1];

        rate[0] = 1 / (2*pi[0]);

        BiAllelicGTR curBAGTRModel = new BiAllelicGTR(pi, rate);

        Network cloneNetwork = Networks.readNetwork(trueNetwork.toString());
        cloneNetwork.getRoot().setRootPopSize(trueNetwork.getRoot().getRootPopSize());
        System.out.println("True Likelihood = " + SNAPPLikelihood.computeSNAPPLikelihood(cloneNetwork, null, alns, BAGTRModel));

        //Utils._START_NET = "[0.036]" + trueNetwork.toString();
        /*String[] gtTaxa = new String[onesnp.size()];
        char[][] sequences = new char[onesnp.size()][500000];
        int count = 0;
        for(String taxon : onesnp.keySet()) {
            gtTaxa[count] = taxon;
            for(int i = 0 ; i < onesnp.get(taxon).length() ; i++) {
                sequences[count][i] = onesnp.get(taxon).charAt(i);
            }
            count++;
        }
        List<char[][]> allLoci = new ArrayList<char[][]>();
        allLoci.add(sequences);
        InferMLNetworkFromSequences inference = new InferMLNetworkFromSequences();
        inference.setParallel(3);
        inference.setStartingNetwork(Networks.readNetwork("(((((A:0.7)I6#H1:1.3::1.0,Q:2.0)I4:1.0,L:3.0)I3:1.0,R:4.0)I2:1.0,(G:2.0,(I6#H1:0.7::0.0,C:1.4)I5:0.6)I1:3.0)I0;"));
        System.out.println(inference.inferNetwork(gtTaxa, allLoci, null, BAGTRModel, 0.036, 4, 1));
*/

        long startTime = System.currentTimeMillis();

        MC3Core run = new MC3Core(alns, curBAGTRModel);
        run.run();


        System.out.println(String.format("Total elapsed time : %2.5f s\n",
                (double) (System.currentTimeMillis() - startTime) / 1000.0));
    }

    public static void testMCMC1() {
        Utils._NUM_THREADS = 2;
        //Utils._CHAIN_LEN = 16000;
        //Utils._BURNIN_LEN = 4000;
        Utils._SAMPLE_FREQUENCY = 500;
        Utils._NET_MAX_RETI = 2;
        Utils._MC3_CHAINS = new ArrayList<>();
        //Utils._MC3_CHAINS.add(1.0);
        //Utils._MC3_CHAINS.add(2.0);
        //Utils._MC3_CHAINS.add(4.0);
        Utils._DIAMETER_PRIOR = true;
        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {0.5, 0.5}, new double[] {1.0});
        SNAPPLikelihood.useOnlyPolymorphic = false;
        SNAPPLikelihood.debugMode = true;
        double gamma = 1 - 0.3;

        double y = 1.0;
        double x = 1000.0;

        int numSites = 100000;

        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, 12345678L + numSites);
        DataGenerator dataGenerator = new DataGenerator();
        Network trueNetwork = Networks.readNetwork(dataGenerator.getNewNetwork1().toString());
        //System.out.println("True Network: " + trueNetwork.toString());

        //SimGTInNetworkByMS gtsimulator = new SimGTInNetworkByMS();
        //SimGTInNetwork gtsimulator = new SimGTInNetwork();
        //List<Tree> gts = gtsimulator.generateGTs(trueNetwork, null, 10);


        trueNetwork = Networks.readNetwork("(((((A:0.7)I6#H1:1.3::0.8,Q:2.0)I4:1.0,L:3.0)I3:1.0,R:4.0)I2:1.0,(G:2.0,(I6#H1:0.7::0.2,C:1.4)I5:0.6)I1:3.0)I0;");
        //trueNetwork = Networks.readNetwork("(((((A:0.7)I6#H1:1.3::1.0,Q:2.0)I4:1.0,L:3.0)I3:1.0,R:4.0)I2:1.0,(G:2.0,(I6#H1:0.7::0.0,C:1.4)I5:0.6)I1:3.0)I0;");
        //trueNetwork = Networks.readNetwork("((d:2.6,((c:1.0,b:1.0)I4:1.5)I3#H1:0.1::0.7)I2:0.5,(a:2.6,I3#H1:0.1::0.30000000000000004)I1:0.5)I0;");
        //trueNetwork = Networks.readNetwork("((((b:1.0,c:1.0)I4:" + y + ")I3#H1:0.1::" + (1 - gamma) + ",a:" + (1.1 + y) + ")I1:" + x + ",(I3#H1:0.1::" + gamma + ",d:" + (1.1 + y) + ")I2:"+ x +")I0;");
        //trueNetwork = Networks.readNetwork("(((((Q:2.0,A:2.0)I4:1.0,L:3.0)I3:0.5)I8#H1:0.5::0.7,R:4.0)I2:1.0,(I8#H1:0.5::0.3,(G:2.0,C:2.0)I1:2.0)I7:1.0)I0;");
        trueNetwork = Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,A:1.0)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,(G:1.0,C:1.0)I1:4.0)I0;");
        //trueNetwork = Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,(A:0.5)I6#H2:0.5::0.8)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,((I6#H2:0.5::0.2,C:1.0)I5:1.0,G:2.0)I1:3.0)I0;");

        double constTheta = 0.036;
        for(Object nodeObject : Networks.postTraversal(trueNetwork)) {
            NetNode node = (NetNode) nodeObject;
            for(Object parentObject :  node.getParents()) {
                NetNode parent = (NetNode) parentObject;
                node.setParentSupport(parent, constTheta);
                node.setParentDistance(parent, node.getParentDistance(parent) * constTheta / 2.0);
            }
        }
        trueNetwork.getRoot().setRootPopSize(constTheta);


        Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, numSites, !SNAPPLikelihood.useOnlyPolymorphic);
        // List<Map<String, String>> snpdata = simulator.generateGTs(Networks.readNetwork("((B:0.5)X#H1:1.5::0.5,((X#H1:0.5::0.5,A:1)n1:0.5,C:1.5)n2:0.5)root;"), null, 100, "/scratch/jz55/Luay/msdir/ms");
        //for(String allele : onesnp.keySet()) {
        //    System.out.println(allele + " " + onesnp.get(allele).substring(0, 10) + "..."+onesnp.get(allele).substring(onesnp.get(allele).length() - 10, onesnp.get(allele).length()));
        //}
        Map<String, String> sequence = new HashMap<>();

        try {
            File file = new File("../data/networkC_" + numSites + "_ac_seed12345678.snp");
            Scanner scanner = new Scanner(file);

            while(scanner.hasNext()) {
                String taxon = scanner.next();
                if(taxon.length() == 0) break;
                String s = scanner.next();
                sequence.put(taxon, s);
                System.out.println(s.length());
            }

            scanner.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        for(String taxon : sequence.keySet()) {
            System.out.println(taxon + " " + sequence.get(taxon).substring(0, 10) + "..." + sequence.get(taxon).substring(sequence.get(taxon).length() - 10, sequence.get(taxon).length()));
        }
        List<Map<String, String>> snpdata = new ArrayList<>();
        //snpdata.add(sequence);

        Map<String, String> sequence2 = new HashMap<>();
        sequence2.put("A", onesnp.get("A"));
        sequence2.put("C", onesnp.get("C"));
        sequence2.put("G", onesnp.get("G"));
        sequence2.put("L", onesnp.get("L"));
        sequence2.put("Q", onesnp.get("Q"));
        sequence2.put("R", onesnp.get("R"));
        snpdata.add(sequence2);

        //snpdata.add(onesnp);
        List<Alignment> alns = new ArrayList<>();
        for(Map<String, String> input : snpdata) {

            System.out.println();
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

            aln.setCache(cache);
            alns.add(aln);
        }

        double [] pi =  new double[2];
        double [] rate = new double[1];

        int totalSites = 0;
        for(Alignment alg : alns) {
            int count = 0;
            for(String s : alg.getAlignment().values()) {
                for(int i = 0 ; i < s.length() ; i++) {
                    count += s.charAt(i) == '1' ? 1 : 0;
                }
            }
            totalSites += alg.getSiteCount() * alg.getTaxonSize();
            pi[1] += 1.0 * count;
        }
        pi[1] /= 1.0 * totalSites;
        pi[0] = 1 - pi[1];

        rate[0] = 1 / (2*pi[0]);

        BiAllelicGTR curBAGTRModel = new BiAllelicGTR(pi, rate);

        Network cloneNetwork = Networks.readNetwork(trueNetwork.toString());
        cloneNetwork.getRoot().setRootPopSize(trueNetwork.getRoot().getRootPopSize());
        System.out.println("True Likelihood = " + SNAPPLikelihood.computeSNAPPLikelihood(cloneNetwork, null, alns, BAGTRModel));

        /*String[] gtTaxa = new String[onesnp.size()];
        char[][] sequences = new char[onesnp.size()][500000];
        int count = 0;
        for(String taxon : onesnp.keySet()) {
            gtTaxa[count] = taxon;
            for(int i = 0 ; i < onesnp.get(taxon).length() ; i++) {
                sequences[count][i] = onesnp.get(taxon).charAt(i);
            }
            count++;
        }
        List<char[][]> allLoci = new ArrayList<char[][]>();
        allLoci.add(sequences);
        InferMLNetworkFromSequences inference = new InferMLNetworkFromSequences();
        inference.setParallel(3);
        inference.setStartingNetwork(Networks.readNetwork("(((((A:0.7)I6#H1:1.3::1.0,Q:2.0)I4:1.0,L:3.0)I3:1.0,R:4.0)I2:1.0,(G:2.0,(I6#H1:0.7::0.0,C:1.4)I5:0.6)I1:3.0)I0;"));
        System.out.println(inference.inferNetwork(gtTaxa, allLoci, null, BAGTRModel, 0.036, 4, 1));
*/

        long startTime = System.currentTimeMillis();
        //Utils._START_NET = "[0.00817528575792145](((G:0.040360738566830895,C:0.040360738566830895):0.002333190642093877)#H1:5.486184649373872::0.04647906647548805,(((R:0.08110800117484995,((L:0.015242059440912224)#H4:0.06507842315816675::0.4169585812227872,((Q:0.01964782667669728)#H2:0.01030261251091685::0.49982865039323054,(#H2:0.010276565369475205::0.5001713496067695,(A:8.140857204844255E-5)#H3:0.029842983474124046::0.44717984052179527):2.6047141441644384E-5):0.05037004341146485):7.875185757709735E-4):1.3144339721918308E-4,#H4:0.06599738513115691::0.5830414187772128):0.020102587039841108,(#H1:0.037216724597387355::0.953520933524512,#H3:0.07982924523426368::0.5528201594782047):0.021431377805598117):5.427536546970886);";
        //Utils._START_NET = "[6.14445682191189E-4]((((C:0.051483563409399594,G:0.051483563409399594):6.708677194252957E-4)#H1:0.05592294494608705::0.83923391801978,(((L:0.02859994721641137)#H2:0.034276001750649716::0.6764465662392649,#H1:0.010721517838236196::0.16076608198022002):8.530354707823734E-4,((R:0.036841234549774934)#H3:0.026416637056929024::0.32685285771114914,((Q:0.0016061554721252332)#H4:0.046795953601801::0.8826884447627434,A:0.04840210907392623):0.014855762532777726):4.711128311395002E-4):0.04434839163706848):5.586966979772012E-4,((#H4:0.034203820575851446::0.11731155523725656,#H2:0.00721002883156531::0.3235534337607351):0.05651079359803648,#H3:0.05547953509623822::0.6731471422888509):0.016315303126875985);";
        //Utils._START_NET = "[1.0415328744458263E-4](((((A:0.017884252774813307)#H2:0.01960559813680477::0.13301357860234186,(L:0.019926847488378936)#H3:0.017563003423239142::0.37368009305891103):0.0014012702562922752)#H1:0.01999008354722731::0.2699211215749967,((C:0.05161526536924467,G:0.05161526536924467):1.1233265413338606E-4)#H4:0.007153606691759608::0.8237656306881619):0.0514882256089412,(((#H3:0.04484845516709575::0.626319906941089,(#H2:0.03238610541864937::0.8669864213976581,Q:0.05027035819346268):0.014504944462012004):0.006644263118072757,#H4:0.019691967750169383::0.1762343693118381):0.02463170998006932,(#H1:0.04383574296997041::0.7300788784250033,R:0.08272686413788076):0.013324411615735998):0.014318154570462108);";
        //Utils._START_NET = "[0.04827255749109394]((((A:0.033117465035632704,Q:0.033117465035632704):0.018796310245952803,(L:0.008191985661944505)#H1:0.043721789619641::0.871229455857914):0.017834572221000204,((#H1:0.012420419281326305::0.128770544142086)#H2:0.0014545150950438686::0.2658977866450809,R:0.022066920038314678):0.04768142746427104):0.009908050948504665,(#H2:0.04227411005335714::0.7341022133549191,(C:0.029764458003666394,G:0.029764458003666394):0.03312205699296156):0.016769883454462425);";
        //Utils._START_NET = "[0.008567344576519902]((((A:0.004292752975453313)#H1:0.02751264488645601::0.22074415864518515,C:0.03180539786190932):0.02164213927789757,G:0.05344753713980689):0.0549169145146883,((L:0.07338498131303631,(#H1:0.04711299821013786::0.7792558413548148,Q:0.051405751185591174):0.021979230127445133):0.017530370949661617,R:0.09091535226269792):0.01744909939179727);";
        //Utils._START_NET = "[0.02465629296292851]((((G:0.07386815517968381,C:0.07386815517968381):0.032749731416242905)#H1:0.00781188013045489::0.9232151324251652,(((Q:0.05245820588880809,A:0.05245820588880809):0.04020161471864697,(R:0.03596367662813288)#H2:0.05669614397932218::0.9127111699408714):0.01036784618498815,L:0.10302766679244321):0.011402099933938398):7.057332883518841,(#H1:2.676717630468579::0.07678486757483483,#H2:2.747371840436373::0.08728883005912857):4.388427133180716);";
        Utils._START_NET = "[3.9858926527049085E-10](((((Q:0.031762044099368496)#H1:0.004355658854667328::0.665255971102495,A:0.036117702954035824):0.018404860764948502,L:0.054522563718984327):0.03506935078073144,(#H1:0.012462567694216886::0.334744028897505,R:0.04422461179358538):0.04536730270613038):0.016286332748302504,(C:0.0348545521461004,G:0.0348545521461004):0.07102369510191786);";
        MC3Core run = new MC3Core(alns, BAGTRModel);
        run.run();


        System.out.println(String.format("Total elapsed time : %2.5f s\n",
                (double) (System.currentTimeMillis() - startTime) / 1000.0));
    }

    public static void testSecondData(String []args) {
        int curPart = Integer.parseInt(args[0]);
        int numSites = 500;
        Utils._NUM_THREADS = 2;
        Utils._CHAIN_LEN = 1500000;
        Utils._BURNIN_LEN = 500000;

        //Utils._CHAIN_LEN = 160000;
        //Utils._BURNIN_LEN = 40000;
        Utils._SAMPLE_FREQUENCY = 500;
        //Utils._NET_MAX_RETI = 0;
        Utils._MC3_CHAINS = new ArrayList<>();
        //Utils._MC3_CHAINS.add(1.0);
        //Utils._MC3_CHAINS.add(2.0);
        //Utils._MC3_CHAINS.add(4.0);
        SNAPPLikelihood.useOnlyPolymorphic = false;
        SNAPPLikelihood.ALGORITHM = 0;
        Utils._DIAMETER_PRIOR = true;
        Utils._TIMES_EXP_PRIOR = true;
        Utils._ESTIMATE_POP_SIZE = false;
        //Utils._CONST_POP_SIZE = false;
        Utils._POP_SIZE_MEAN = 0.04;
        //Utils.EXP_PARAM = 0.1;
        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {0.5, 0.5}, new double[] {1.0});

        Utils.printSettings();


        int number = 0;
        List<String> labels = new ArrayList<>();
        labels.add("net_0.1");
        labels.add("net_0.2");
        labels.add("net_0.4");
        labels.add("net_1.0");

        //curPart = paths.indexOf("../data/networkA_1000000_ac_seed12345678.snp");
        System.out.println("Total number " + number);
        System.out.println("Current number " + curPart);
        System.out.println(labels.get(curPart));
        String label = labels.get(curPart);

        Map<String, String> networks = new HashMap<>();
        networks.put("net_0.1", "((((C:0.003,(B:0.002)I7#H1:0.001::0.2)I4:0.002,(I7#H1:0.002::0.8)I6#H2:0.001::0.2)I3:0.002,(I6#H2:0.002::0.8)I5#H3:0.001::0.2)I2:0.003,(I5#H3:0.002::0.8,A:0.008)I1:0.002)I0;");
        networks.put("net_0.2", "((((C:0.006,(B:0.004)I7#H1:0.002::0.2)I4:0.004,(I7#H1:0.004::0.8)I6#H2:0.002::0.2)I3:0.004,(I6#H2:0.004::0.8)I5#H3:0.002::0.2)I2:0.006,(I5#H3:0.004::0.8,A:0.016)I1:0.004)I0;");
        networks.put("net_0.4", "((((C:0.012,(B:0.008)I7#H1:0.004::0.2)I4:0.008,(I7#H1:0.008::0.8)I6#H2:0.004::0.2)I3:0.008,(I6#H2:0.008::0.8)I5#H3:0.004::0.2)I2:0.012,(I5#H3:0.008::0.8,A:0.032)I1:0.008)I0;");
        networks.put("net_1.0", "((((C:0.07,(B:0.06)I7#H1:0.01::0.2)I4:0.02,(I7#H1:0.02::0.8)I6#H2:0.01::0.2)I3:0.02,(I6#H2:0.02::0.8)I5#H3:0.01::0.2)I2:0.07,(I5#H3:0.06::0.8,A:0.16)I1:0.02)I0;");

        Network trueNetwork = Networks.readNetwork(networks.get(label));
        trueNetwork.getRoot().setRootPopSize(0.04);


        Map<String, String> sequence = new HashMap<>();

        String filename = "../seconddata/" + label + "_" + numSites + "_ac_seed12345678.snp";
        System.out.println(filename);

        try {
            File file = new File(filename);
            Scanner scanner = new Scanner(file);

            while(scanner.hasNext()) {
                String taxon = scanner.next();
                if(taxon.length() == 0) break;
                String s = scanner.next();
                sequence.put(taxon, s);
                System.out.println(s.length());
            }

            scanner.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        for(String taxon : sequence.keySet()) {
            System.out.println(taxon + " " + sequence.get(taxon).substring(0, 10) + "..." + sequence.get(taxon).substring(sequence.get(taxon).length() - 10, sequence.get(taxon).length()));
        }
        List<Map<String, String>> snpdata = new ArrayList<>();
        snpdata.add(sequence);

        List<Alignment> alns = new ArrayList<>();
        for(Map<String, String> input : snpdata) {

            System.out.println();
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

            aln.setCache(cache);
            alns.add(aln);
        }

        double [] pi =  new double[2];
        double [] rate = new double[1];

        int totalSites = 0;
        for(Alignment alg : alns) {
            int count = 0;
            for(String s : alg.getAlignment().values()) {
                for(int i = 0 ; i < s.length() ; i++) {
                    count += s.charAt(i) == '1' ? 1 : 0;
                }
            }
            totalSites += alg.getSiteCount() * alg.getTaxonSize();
            pi[1] += 1.0 * count;
        }
        pi[1] /= 1.0 * totalSites;
        pi[0] = 1 - pi[1];

        rate[0] = 1 / (2*pi[0]);

        pi[1] = 0.5;
        pi[0] = 1 - pi[1];

        rate[0] = 1 / (2*pi[0]);

        BiAllelicGTR curBAGTRModel = new BiAllelicGTR(pi, rate);

        Network cloneNetwork = Networks.readNetwork(trueNetwork.toString());
        cloneNetwork.getRoot().setRootPopSize(trueNetwork.getRoot().getRootPopSize());
        System.out.println("True Likelihood = " + SNAPPLikelihood.computeSNAPPLikelihood(cloneNetwork, null, alns, BAGTRModel));

        long startTime = System.currentTimeMillis();

        //Utils._START_NET = "[9.057126525354496E-12]((((((((C:6.071885919301598E-4,B:6.071885919301598E-4):0.004084977900907154)#H2:1.9180763597067933E-4::0.9984510547505954,A:0.004883974128807993):0.0461067287811829,#H2:0.04629853641715358::0.0015489452494046319):12.119610091552316)#H1:1.178512708785492::0.12210260828789488)#H3:0.0696553828071309::0.9838799395493458,#H1:1.248168091592623::0.8778973917121051):92.1852811291236,#H3:92.25493651193074::0.016120060450654172);";
        //Utils._START_NET = "[0.04](C:920.8473410747589,(B:626.3380424743829,A:626.3380424743829):294.50929860037604);";
        //Utils._START_NET = "[0.04](((C:0.15496968657284896)#H1:0.09976506959530898::0.6848580260435815,B:0.25473475616815794):0.5224823341454845,(#H1:7.311128982365467E-4::0.3151419739564185,A:0.1557007994710855):0.6215162908425569);";
        MC3Core run = new MC3Core(alns, curBAGTRModel);
        run.run();


        System.out.println(String.format("Total elapsed time : %2.5f s\n",
                (double) (System.currentTimeMillis() - startTime) / 1000.0));
    }

    public static void testMosquitoData(String []args) {
        Utils._NUM_THREADS = 8;
        //Utils._CHAIN_LEN = 2000500;
        //Utils._BURNIN_LEN = 500;
        Utils._CHAIN_LEN = 2000000;
        Utils._BURNIN_LEN = 500000;
        Utils._SAMPLE_FREQUENCY = 500;
        //Utils._NET_MAX_RETI = 2;
        Utils._DIAMETER_PRIOR = true;
        Utils._TIMES_EXP_PRIOR = true;
        Utils._MC3_CHAINS = new ArrayList<>();
        //Utils._MC3_CHAINS.add(1.0);
        //Utils._MC3_CHAINS.add(2.0);
        //Utils._MC3_CHAINS.add(4.0);
        Utils._CONST_POP_SIZE = true;
        Utils._POP_SIZE_MEAN = 0.014;

        String filename = "../mosquito/sampled.snp";

        SNAPPLikelihood.ALGORITHM = 2;
        SNAPPLikelihood.useOnlyPolymorphic = false;
        //SNAPPLikelihood.debugMode = true;

        if(args.length > 0) {
            filename = args[0];
        }

        if(args.length > 1) {
            Utils._SEED = Integer.parseInt(args[1]);
        }

        if(args.length > 2) {
            Utils._CHAIN_LEN = Integer.parseInt(args[2]);
        }

        String breakpoint = null;//"breakpoint.txt";
        if(breakpoint != null) {
            System.out.println("Breakpoint: " + breakpoint);
            Utils._START_NET = getBestNetwork(breakpoint);
            System.out.println("Start network: " + Utils._START_NET);
        }

        Utils.printSettings();

        Map<String, String> sequence = new HashMap<>();

        try {
            File file = new File(filename);
            System.out.println(file.getName());
            Scanner scanner = new Scanner(file);

            while(scanner.hasNext()) {
                String taxon = scanner.next();
                if(taxon.length() == 0) break;
                String s = scanner.next();
                sequence.put(taxon, s);
                System.out.println(s.length());
            }

            scanner.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        for(String taxon : sequence.keySet()) {
            System.out.println(taxon + " " + sequence.get(taxon).substring(0, 10) + "..." + sequence.get(taxon).substring(sequence.get(taxon).length() - 10, sequence.get(taxon).length()));
        }
        List<Map<String, String>> snpdata = new ArrayList<>();
        snpdata.add(sequence);

        List<Alignment> alns = new ArrayList<>();
        for(Map<String, String> input : snpdata) {

            System.out.println();
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

            aln.setCache(cache);
            alns.add(aln);
        }

        double [] pi =  new double[2];
        double [] rate = new double[1];

       /* int totalSites = 0;
        for(Alignment alg : alns) {
            int count = 0;
            for(String s : alg.getAlignment().values()) {
                for(int i = 0 ; i < s.length() ; i++) {
                    count += s.charAt(i) == '1' ? 1 : 0;
                }
            }
            totalSites += alg.getSiteCount() * alg.getTaxonSize();
            pi[1] += 1.0 * count;
        }
        pi[1] /= 1.0 * totalSites;
        pi[0] = 1 - pi[1];

        rate[0] = 1 / (2*pi[0]);*/

        pi[1] = 0.5;
        pi[0] = 1 - pi[1];

        rate[0] = 1.0 / (2*pi[0]);

        BiAllelicGTR curBAGTRModel = new BiAllelicGTR(pi, rate);

        //Network cloneNetwork = Networks.readNetwork(trueNetwork.toString());
        //cloneNetwork.getRoot().setRootPopSize(trueNetwork.getRoot().getRootPopSize());
        //System.out.println("True Likelihood = " + SNAPPLikelihood.computeSNAPPLikelihood(cloneNetwork, null, alns, BAGTRModel));

        long startTime = System.currentTimeMillis();

        //Utils._START_NET = "[0.013941695932206571]((L:0.0027483538928538406)#H1:0.024591705100111893::0.3457545488771119,(((#H1:0.006070791744061946::0.6542454511228881,((C:0.001705161774251553,G:0.001705161774251553):0.0036209530403215947,A:0.0053261148145731475):0.0034930308223426387):7.363852463073764E-4,Q:0.009555530883223163):0.005224376373616809,R:0.014779907256839971):0.012560151736125763);";
        //Utils._START_NET = "[0.013183212466841288]((R:0.005393926178755438)#H1:103.9561992239639743::0.018548690495628684,(L:0.015519221548218956,(#H1:0.005032251955269181::0.9814513095043713,(((G:0.002179801517554337,C:0.002179801517554337):0.003448879780084071,A:0.005628681297638408):0.0039931010651085815,Q:0.00962178236274699):8.043957712776297E-4):0.005093043414194337):103.9460739285945103);";
        //Utils._START_NET = "[0.01432278740805789](R:0.014623597111565564,((Q:0.00911822640342767,(A:0.005012995498108841,(G:0.0013615204510037562,C:0.0013615204510037562):0.0036514750471050845):0.004105230905318829):8.83545194053861E-4,(L:0.004912143753730035):0.005089627843751496):0.004621825514084034);";
        //Utils._START_NET = "[0.004730554505572216]((L:0.0018026486850655757)#H1:0.006868170254180248::0.37434235096700086,((Q:0.003017192479968964,(((C:8.059411556637647E-4,G:8.059411556637647E-4):7.375228243481855E-4,A:0.0015434639800119502):0.0010181992335458793,#H1:7.590145284922538E-4::0.6256576490329991):4.5552926641113456E-4):0.0018414552390651788,R:0.004858647719034143):0.003812171220211681);";

        MC3Core run = new MC3Core(alns, curBAGTRModel);
        run.run();


        System.out.println(String.format("Total elapsed time : %2.5f s\n",
                (double) (System.currentTimeMillis() - startTime) / 1000.0));
    }

    public static void testWithSNAPPSimulation(String[] args) {
        int curPart = Integer.parseInt(args[0]);
        int number = 0;
        Utils._NET_MAX_RETI = 0;
        Utils._TIMES_EXP_PRIOR = true;
        SNAPPLikelihood.useOnlyPolymorphic = true;
        List<String> paths = new ArrayList<>();
        for (File file : new File("/Users/zhujiafan/snappPaperVerification/PaperSimulations/Simulation1Clean").listFiles()) {
            if (file.isFile() && file.getName().endsWith(".xml")) {
                paths.add(file.getPath());
                number++;
            }
        }
        Collections.sort(paths);
        System.out.println("Total number of xmls = " + number);
        System.out.println("Current number " + curPart);
        System.out.println(paths.get(curPart));

        long startTime = System.currentTimeMillis();

        Map<String, String> sequence = new HashMap<>();
        try {
            File fXmlFile = new File(paths.get(curPart));
            DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
            Document doc = dBuilder.parse(fXmlFile);
            doc.getDocumentElement().normalize();

            System.out.println("Root element :" + doc.getDocumentElement().getNodeName());
            NodeList nList = doc.getElementsByTagName("sequence");
            for(int i = 0 ; i < nList.getLength() ; i++) {
                Node nNode = nList.item(i);
                System.out.println("Current Element: " + nNode.getNodeName());
                Element eElement = (Element) nNode;
                System.out.println(eElement.getAttribute("taxon"));

                System.out.println(nNode.getTextContent());
                String seqToParse = nNode.getTextContent();
                String seq = "";
                for(int c = 0 ; c < seqToParse.length() ; c++) {
                    if(seqToParse.charAt(c) == '0' || seqToParse.charAt(c) == '1')
                        seq += seqToParse.charAt(c);
                }
                sequence.put(eElement.getAttribute("taxon"), seq);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        for(String taxon : sequence.keySet()) {
            System.out.println(taxon + " " + sequence.get(taxon));
        }
        List<Map<String, String>> snpdata = new ArrayList<>();
        snpdata.add(sequence);
        List<Alignment> alns = new ArrayList<>();

        for(Map<String, String> input : snpdata) {
            for(String allele : input.keySet()) {
                System.out.println(allele + " " + input.get(allele));
            }
            System.out.println();
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

            aln.setCache(cache);
            alns.add(aln);
        }

        double [] pi =  new double[2];
        double [] rate = new double[1];

        int totalSites = 0;
        for(Alignment alg : alns) {
            int count = 0;
            for(String s : alg.getAlignment().values()) {
                for(int i = 0 ; i < s.length() ; i++) {
                    count += s.charAt(i) == '1' ? 1 : 0;
                }
            }
            totalSites += alg.getSiteCount() * alg.getTaxonSize();
            pi[1] += 1.0 * count;
        }
        pi[1] /= 1.0 * totalSites;
        pi[0] = 1 - pi[1];

        rate[0] = 1 / (2*pi[1]);

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(pi, rate);

        MC3Core run = new MC3Core(alns, BAGTRModel);
        run.run();

        System.out.println(String.format("Total elapsed time : %2.5f s\n",
                (double) (System.currentTimeMillis() - startTime) / 1000.0));
    }

    public static String getBestNetwork(String filename) {
        String lastNetwork = null;
        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            String line;
            String inputfile = "";
            int seqLength;
            String allowConstantLabel;
            String currentChain = "";
            List<Network> trueBackbones = new ArrayList<>();
            boolean ready = false;
            boolean finished = false;
            double bestLogPosterior = -1e18;

            while ((line = br.readLine()) != null) {
                if(line.contains("/") && inputfile.equals("")) {
                    inputfile = line.substring(line.lastIndexOf('/') + 1);
                    System.out.println("Last Run Input: " + inputfile);
                }
                else if(line.startsWith("Temp")) {
                    currentChain = line;
                    ready = true;
                }
                else if(ready) {
                    Scanner lineScanner = new Scanner(line);
                    lineScanner.useDelimiter(";    ");
                    int numSamples = lineScanner.nextInt();
                    double logPosterior = lineScanner.nextDouble();
                    if(currentChain.contains("main") && bestLogPosterior < logPosterior) {
                        bestLogPosterior = logPosterior;
                        lastNetwork = br.readLine();
                    }
                    ready = false;
                }
                else if(line.startsWith("Rank = 0")){
                    System.out.println("Last Run Finished ");
                    String str = line.substring(line.indexOf(":") + 1);
                    lastNetwork = str.substring(0, str.indexOf(";") + 1);
                }
            }

            System.out.println("Best log posterior: " + bestLogPosterior);

        } catch (Exception e) {
            e.printStackTrace();
        }

        return lastNetwork;
    }

    public static void testWithDingqiaoNetwork(String[] args) {
        //String whichToRun = "seconddata_mono";
        //String whichToRun = "firstdata_mono_firstrun";
        String whichToRun = "firstdata_mono_secondrun";
        int curPart = Integer.parseInt(args[0]);

        String path = "";
        Utils._MC3_CHAINS = new ArrayList<>();
        Map<String, String> networks = new HashMap<>();
        String breakpoint = null;//"/scratch/jz55/useOnlyPolymorphic/run/slurm-3175687_";

        if(whichToRun.equals("firstdata_mono_firstrun")) {
            path = "/scratch/jz55/usePolyMono/data";
            Utils._NUM_THREADS = 8;
            Utils._NET_MAX_RETI = 2;
            if(curPart >= 12)
                Utils._NET_MAX_RETI = 3;

            Utils._CHAIN_LEN = 300000;
            Utils._BURNIN_LEN = 200000;
            Utils._DIAMETER_PRIOR = true;
            Utils._TIMES_EXP_PRIOR = true;

            Utils._MC3_CHAINS.add(2.0);
            Utils._MC3_CHAINS.add(4.0);

            SNAPPLikelihood.useOnlyPolymorphic = false;
            SNAPPLikelihood.ALGORITHM = 2;

            networks.put("networkA", "(((((A:0.7)I6#H1:1.3::0.8,Q:2.0)I4:1.0,L:3.0)I3:1.0,R:4.0)I2:1.0,(G:2.0,(I6#H1:0.7::0.2,C:1.4)I5:0.6)I1:3.0)I0;");
            networks.put("networkB", "(((((Q:2.0,A:2.0)I4:1.0,L:3.0)I3:0.5)I8#H1:0.5::0.7,R:4.0)I2:1.0,(I8#H1:0.5::0.3,(G:2.0,C:2.0)I1:2.0)I7:1.0)I0;");
            networks.put("networkC", "(((((Q:0.5)I8#H1:0.5::0.7,A:1.0)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,(G:1.0,C:1.0)I1:4.0)I0;");
            networks.put("networkD", "(((((Q:0.5)I8#H1:0.5::0.7,(A:0.5)I6#H2:0.5::0.8)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,((I6#H2:0.5::0.2,C:1.0)I5:1.0,G:2.0)I1:3.0)I0;");
        }
        else if(whichToRun.equals("firstdata_mono_secondrun")){
            path = "/scratch/jz55/usePolyMono/data";
            Utils._NUM_THREADS = 8;

            Utils._CHAIN_LEN = 1500000;
            Utils._BURNIN_LEN = 500000;
            Utils._DIAMETER_PRIOR = true;
            Utils._TIMES_EXP_PRIOR = true;
            //Utils._ESTIMATE_POP_PARAM = false;
            Utils._ESTIMATE_POP_SIZE = false;

            SNAPPLikelihood.useOnlyPolymorphic = false;
            SNAPPLikelihood.ALGORITHM = 2;

            breakpoint = "/scratch/jz55/usePolyMono/run/slurm-3199866_";

            networks.put("networkA", "(((((A:0.7)I6#H1:1.3::0.8,Q:2.0)I4:1.0,L:3.0)I3:1.0,R:4.0)I2:1.0,(G:2.0,(I6#H1:0.7::0.2,C:1.4)I5:0.6)I1:3.0)I0;");
            networks.put("networkB", "(((((Q:2.0,A:2.0)I4:1.0,L:3.0)I3:0.5)I8#H1:0.5::0.7,R:4.0)I2:1.0,(I8#H1:0.5::0.3,(G:2.0,C:2.0)I1:2.0)I7:1.0)I0;");
            networks.put("networkC", "(((((Q:0.5)I8#H1:0.5::0.7,A:1.0)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,(G:1.0,C:1.0)I1:4.0)I0;");
            networks.put("networkD", "(((((Q:0.5)I8#H1:0.5::0.7,(A:0.5)I6#H2:0.5::0.8)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,((I6#H2:0.5::0.2,C:1.0)I5:1.0,G:2.0)I1:3.0)I0;");

        }
        else if(whichToRun.equals("seconddata_mono")) {
            path = "/scratch/jz55/seconddata_mono";
            Utils._NUM_THREADS = 2;

            Utils._CHAIN_LEN = 1500000;
            Utils._BURNIN_LEN = 500000;
            SNAPPLikelihood.useOnlyPolymorphic = false;
            SNAPPLikelihood.ALGORITHM = 2;
            Utils._ESTIMATE_POP_SIZE = false;
            Utils._POP_SIZE_MEAN = 0.04;
            Utils._DIAMETER_PRIOR = true;
            Utils._TIMES_EXP_PRIOR = true;

            networks.put("net_0.1", "((((C:0.003,(B:0.002)I7#H1:0.001::0.2)I4:0.002,(I7#H1:0.002::0.8)I6#H2:0.001::0.2)I3:0.002,(I6#H2:0.002::0.8)I5#H3:0.001::0.2)I2:0.003,(I5#H3:0.002::0.8,A:0.008)I1:0.002)I0;");
            networks.put("net_0.2", "((((C:0.006,(B:0.004)I7#H1:0.002::0.2)I4:0.004,(I7#H1:0.004::0.8)I6#H2:0.002::0.2)I3:0.004,(I6#H2:0.004::0.8)I5#H3:0.002::0.2)I2:0.006,(I5#H3:0.004::0.8,A:0.016)I1:0.004)I0;");
            networks.put("net_0.4", "((((C:0.012,(B:0.008)I7#H1:0.004::0.2)I4:0.008,(I7#H1:0.008::0.8)I6#H2:0.004::0.2)I3:0.008,(I6#H2:0.008::0.8)I5#H3:0.004::0.2)I2:0.012,(I5#H3:0.008::0.8,A:0.032)I1:0.008)I0;");
            networks.put("net_1.0", "((((C:0.07,(B:0.06)I7#H1:0.01::0.2)I4:0.02,(I7#H1:0.02::0.8)I6#H2:0.01::0.2)I3:0.02,(I6#H2:0.02::0.8)I5#H3:0.01::0.2)I2:0.07,(I5#H3:0.06::0.8,A:0.16)I1:0.02)I0;");

        } else {
            System.out.println("Which to run?");
            return;
        }



        //String path = "../data";
        //Utils._NUM_THREADS = 1;

        //if(Runtime.getRuntime().availableProcessors() == 32) {
        //    Utils._NUM_THREADS = 16;
        //    path = "/scratch/jz55/data";
        //}

        //Utils._MC3_CHAINS.add(1.0);
        //Utils._MC3_CHAINS.add(2.0);
        //Utils._MC3_CHAINS.add(4.0);
        //Utils._NET_MAX_RETI = 2;


        //if(curPart >= 18)
        //Utils._NET_MAX_RETI = 3;


        Utils.printSettings();

        int number = 0;
        List<String> paths = new ArrayList<>();
        for (File file : new File(path).listFiles()) {
            if (file.isFile() && file.getName().endsWith(".snp")) {
                paths.add(file.getPath());
                number++;
            }
        }
        Collections.sort(paths);
        //curPart = paths.indexOf("../data/networkA_1000000_ac_seed12345678.snp");
        System.out.println("Total number " + number);
        System.out.println("Current number " + curPart);
        System.out.println(paths.get(curPart));

        Network trueNetwork = null;
        for(String s : networks.keySet()) {
            if(paths.get(curPart).contains(s)) {
                trueNetwork = Networks.readNetwork(networks.get(s));
                break;
            }
        }
        double constTheta = 0.036;
        double mu = 1.0;
        for(Object nodeObject : Networks.postTraversal(trueNetwork)) {
            NetNode node = (NetNode) nodeObject;
            for(Object parentObject :  node.getParents()) {
                NetNode parent = (NetNode) parentObject;
                node.setParentSupport(parent, constTheta);
                node.setParentDistance(parent, node.getParentDistance(parent) * constTheta / 2.0);
            }
        }
        trueNetwork.getRoot().setRootPopSize(constTheta);
        int nameCount = 0;
        for(Object node : trueNetwork.dfs()) {
            NetNode mynode = (NetNode) node;
            if(mynode.getName().equals("")) {
                mynode.setName("II" + nameCount);
                nameCount++;
            }
        }

        if(breakpoint != null) {
            breakpoint = breakpoint + curPart + ".out";
            System.out.println("Breakpoint: " + breakpoint);
            Utils._START_NET = getBestNetwork(breakpoint);
            System.out.println("Start network: " + Utils._START_NET);
        }

        long startTime = System.currentTimeMillis();

        Map<String, String> sequence = new HashMap<>();

        try {
            File file = new File(paths.get(curPart));
            Scanner scanner = new Scanner(file);

            while(scanner.hasNext()) {
                String taxon = scanner.next();
                if(taxon.length() == 0) break;
                String s = scanner.next();
                sequence.put(taxon, s);
            }

            scanner.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        for(String taxon : sequence.keySet()) {
            System.out.println(taxon + " " + sequence.get(taxon).substring(0, 10) + "..." + sequence.get(taxon).substring(sequence.get(taxon).length() - 10, sequence.get(taxon).length()));
        }
        List<Map<String, String>> snpdata = new ArrayList<>();
        snpdata.add(sequence);
        List<Alignment> alns = new ArrayList<>();

        for(Map<String, String> input : snpdata) {
            //for(String allele : input.keySet()) {
            //System.out.println(allele + " " + input.get(allele));
            //}
            System.out.println();
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

            aln.setCache(cache);
            alns.add(aln);
        }

        double [] pi =  new double[2];
        double [] rate = new double[1];

        int totalSites = 0;
        for(Alignment alg : alns) {
            int count = 0;
            for(String s : alg.getAlignment().values()) {
                for(int i = 0 ; i < s.length() ; i++) {
                    count += s.charAt(i) == '1' ? 1 : 0;
                }
            }
            totalSites += alg.getSiteCount() * alg.getTaxonSize();
            pi[1] += 1.0 * count;
        }
        pi[1] /= 1.0 * totalSites;
        pi[1] = 0.5;
        pi[0] = 1 - pi[1];

        rate[0] = 1 / (2*pi[0]);

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(pi, rate);

        System.out.println("True Network: " + trueNetwork.toString());
        Network cloneNetwork = Networks.readNetwork(trueNetwork.toString());
        cloneNetwork.getRoot().setRootPopSize(trueNetwork.getRoot().getRootPopSize());
        System.out.println("True Likelihood = " + SNAPPLikelihood.computeSNAPPLikelihood(cloneNetwork, null, alns, BAGTRModel));


        MC3Core run = new MC3Core(alns, BAGTRModel);
        run.run();

        System.out.println(String.format("Total elapsed time : %2.5f s\n",
                (double) (System.currentTimeMillis() - startTime) / 1000.0));
    }

    public static void checkSNAPPSimulationResults() {

        int number = 0;
        int correct = 0;
        List<String> paths = new ArrayList<>();
        for (File file : new File("/Users/zhujiafan/snappPaperVerification").listFiles()) {
            if (file.isFile() && file.getName().endsWith(".out")) {
                paths.add(file.getPath());
                number++;
            }
        }
        System.out.println("Total number = " + number);
        Collections.sort(paths);
        for(String path : paths) {
            try {
                BufferedReader br = new BufferedReader(new FileReader(path));
                String line;
                Network currentNetwork = null;
                Network trueNetwork = null;
                String inputfile = "";
                boolean correctness = false;
                while ((line = br.readLine()) != null) {
                    if(line.startsWith("/scratch/jz55/")) {
                        inputfile = line;
                        if(line.contains("easy8") || line.contains("hard8")) {
                            trueNetwork = Networks.readNetwork("((((A,B),C),D),((E,F),(G,H)));");
                        } else if(line.contains("easy4") || line.contains("hard4")) {
                            trueNetwork = Networks.readNetwork("(((A,B),C),D);");
                        }
                    }
                    else if(line.startsWith("Rank = ")){
                        String str = line.substring(line.indexOf(":") + 1);
                        currentNetwork = Networks.readNetwork(str);
                        if(trueNetwork != null) {
                            double [] diff = Networks.computeNormalizedTreeDistance(currentNetwork, trueNetwork);
                            if(diff[2] == 0) {
                                correctness = true;
                                break;
                            }
                        }
                    }
                }
                if(correctness) {
                    correct++;
                } else {
                    System.out.println("Incorrect: " + inputfile + " " + path);
                }

            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        System.out.println("Number of correct = " + correct);
    }

    public static void checkDingqiaoNetworkResults() {

        //Network n1 = Networks.readNetwork("(((((A:0.7)I6#H1:1.3::0.8,Q:2.0)I4:1.0,L:3.0)I3:1.0,R:4.0)I2:1.0,(G:2.0,(I6#H1:0.7::0.2,C:1.4)I5:0.6)I1:3.0)I0;");
        //Network n2 = Networks.readNetwork("(((((Q:2.0,A:2.0)I4:1.0,L:3.0)I3:0.5)I8#H1:0.5::0.7,R:4.0)I2:1.0,(I8#H1:0.5::0.3,(G:2.0,C:2.0)I1:2.0)I7:1.0)I0;");
        //boolean testing = Networks.hasTheSameTopology(n1, n2);

        //if(true) {
        //    return;
        //}

        int number = 0;
        List<String> paths = new ArrayList<>();
        for (File file : new File("/Users/zhujiafan/ABCD_firstRun").listFiles()) {
            if (file.isFile() && file.getName().endsWith(".out")) {
                paths.add(file.getPath());
                number++;
            }
        }
        Collections.sort(paths);
        System.out.println("Total number " + number);

        Map<String, Network> networks = new HashMap<>();
        networks.put("networkA", Networks.readNetwork("(((((A:0.7)I6#H1:1.3::0.8,Q:2.0)I4:1.0,L:3.0)I3:1.0,R:4.0)I2:1.0,(G:2.0,(I6#H1:0.7::0.2,C:1.4)I5:0.6)I1:3.0)I0;"));
        networks.put("networkB", Networks.readNetwork("(((((Q:2.0,A:2.0)I4:1.0,L:3.0)I3:0.5)I8#H1:0.5::0.7,R:4.0)I2:1.0,(I8#H1:0.5::0.3,(G:2.0,C:2.0)I1:2.0)I7:1.0)I0;"));
        networks.put("networkC", Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,A:1.0)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,(G:1.0,C:1.0)I1:4.0)I0;"));
        networks.put("networkD", Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,(A:0.5)I6#H2:0.5::0.8)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,((I6#H2:0.5::0.2,C:1.0)I5:1.0,G:2.0)I1:3.0)I0;"));

        int correct = 0;

        for(String path : paths) {
            try {
                BufferedReader br = new BufferedReader(new FileReader(path));
                String line;
                Network currentNetwork = null;
                Network trueNetwork = null;
                String inputfile = "";
                String networkLabel;
                int seqLength;
                String allowConstantLabel;
                List<Network> trueBackbones = new ArrayList<>();
                boolean correctNetwork = false;
                boolean correctBackbone = false;
                boolean containOther = false;
                boolean finished = false;
                boolean skip = false;

                while ((line = br.readLine()) != null) {
                    if(line.startsWith("/scratch/jz55/")) {
                        inputfile = line.substring(line.lastIndexOf('/') + 1);

                        Scanner lineScanner = new Scanner(inputfile);
                        lineScanner.useDelimiter("_");

                        networkLabel = lineScanner.next();

                        seqLength = Integer.parseInt(lineScanner.next());
                        allowConstantLabel = lineScanner.next();

                        if(!networkLabel.equals("networkB") || !allowConstantLabel.equals("nc"))
                            skip = true;

                        trueNetwork = networks.get(networkLabel).clone();
                        for (Object ntObject : Networks.getTrees(trueNetwork)) {
                            NetworkTree<Object> nt = (NetworkTree<Object>) ntObject;
                            trueBackbones.add(Networks.readNetwork(nt.makeTree().toNewick()));
                        }
                    }
                    else if(line.startsWith("Rank = ")){
                        String str = line.substring(line.indexOf(":") + 1);
                        str = str.substring(0, str.indexOf(";") + 1);
                        currentNetwork = Networks.readNetwork(str);
                        finished = true;

                        if(trueNetwork != null) {
                            if(Networks.hasTheSameTopology(currentNetwork, trueNetwork)) {
                                correctNetwork = true;
                            } else {
                                for(Network backbone : trueBackbones) {
                                    if(Networks.hasTheSameTopology(currentNetwork, backbone)){
                                        correctBackbone = true;
                                        break;
                                    }
                                }
                                if(!correctBackbone)
                                    containOther = true;
                            }
                        }
                    }
                }

                if(skip) continue;

                System.out.println(inputfile);
                System.out.println("Finished " + finished);
                System.out.println("Contain Correct Network " + correctNetwork);
                System.out.println("Contain Correct Backbone " + correctBackbone);
                System.out.println("Contain Other " + containOther);
                if(correctNetwork) {
                    correct++;
                }

            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    public static void generateSNPdata() {
        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {0.5, 0.5}, new double[] {1.0});

        Map<String, String> networks = new HashMap<>();
        networks.put("networkA", "(((((A:0.7)I6#H1:1.3::0.8,Q:2.0)I4:1.0,L:3.0)I3:1.0,R:4.0)I2:1.0,(G:2.0,(I6#H1:0.7::0.2,C:1.4)I5:0.6)I1:3.0)I0;");
        networks.put("networkB", "(((((Q:2.0,A:2.0)I4:1.0,L:3.0)I3:0.5)I8#H1:0.5::0.7,R:4.0)I2:1.0,(I8#H1:0.5::0.3,(G:2.0,C:2.0)I1:2.0)I7:1.0)I0;");
        networks.put("networkC", "(((((Q:0.5)I8#H1:0.5::0.7,A:1.0)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,(G:1.0,C:1.0)I1:4.0)I0;");
        networks.put("networkD", "(((((Q:0.5)I8#H1:0.5::0.7,(A:0.5)I6#H2:0.5::0.8)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,((I6#H2:0.5::0.2,C:1.0)I5:1.0,G:2.0)I1:3.0)I0;");


        int lengtharray[] = new int[]{50000};//

//        int lengtharray[] = new int[]{100, 1000, 10000, 100000, 250000, 500000, 1000000, 2000000, 4000000};//

        for(int i = 0 ; i < lengtharray.length ; i++) {
            System.out.println("Generating length=" + lengtharray[i]);
            for(String label : networks.keySet()) {
                //if(lengtharray[i] != 250000 || !label.equals("networkC")) continue;

                Network trueNetwork = Networks.readNetwork(networks.get(label));

                double constTheta = 0.036;
                for(Object nodeObject : Networks.postTraversal(trueNetwork)) {
                    NetNode node = (NetNode) nodeObject;
                    for(Object parentObject :  node.getParents()) {
                        NetNode parent = (NetNode) parentObject;
                        node.setParentSupport(parent, constTheta);
                        node.setParentDistance(parent, node.getParentDistance(parent) * constTheta / 2.0);
                    }
                }
                trueNetwork.getRoot().setRootPopSize(constTheta);


                String filename = "../data/" + label + "_" + lengtharray[i] + "_ac_seed12345678.snp";
                SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, Utils._SEED + lengtharray[i]);
                Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, lengtharray[i], true);
                List<String> taxa = new ArrayList<>(onesnp.keySet());
                Collections.sort(taxa);
                try {
                    PrintWriter out = new PrintWriter(filename);
                    for (String taxon : taxa) {
                        out.println(taxon + " " + onesnp.get(taxon));
                    }
                    out.close();
                } catch (IOException ioe) {
                    ioe.printStackTrace();
                }

            }
        }
    }

    public static void generateSecondSNPdata() {
        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {0.5, 0.5}, new double[] {1.0});

        Map<String, String> networks = new HashMap<>();
        networks.put("net_0.1", "((((C:0.003,(B:0.002)I7#H1:0.001::0.2)I4:0.002,(I7#H1:0.002::0.8)I6#H2:0.001::0.2)I3:0.002,(I6#H2:0.002::0.8)I5#H3:0.001::0.2)I2:0.003,(I5#H3:0.002::0.8,A:0.008)I1:0.002)I0;");
        networks.put("net_0.2", "((((C:0.006,(B:0.004)I7#H1:0.002::0.2)I4:0.004,(I7#H1:0.004::0.8)I6#H2:0.002::0.2)I3:0.004,(I6#H2:0.004::0.8)I5#H3:0.002::0.2)I2:0.006,(I5#H3:0.004::0.8,A:0.016)I1:0.004)I0;");
        networks.put("net_0.4", "((((C:0.012,(B:0.008)I7#H1:0.004::0.2)I4:0.008,(I7#H1:0.008::0.8)I6#H2:0.004::0.2)I3:0.008,(I6#H2:0.008::0.8)I5#H3:0.004::0.2)I2:0.012,(I5#H3:0.008::0.8,A:0.032)I1:0.008)I0;");
        networks.put("net_1.0", "((((C:0.07,(B:0.06)I7#H1:0.01::0.2)I4:0.02,(I7#H1:0.02::0.8)I6#H2:0.01::0.2)I3:0.02,(I6#H2:0.02::0.8)I5#H3:0.01::0.2)I2:0.07,(I5#H3:0.06::0.8,A:0.16)I1:0.02)I0;");


        int lengtharray[] = new int[]{100, 500, 1000, 50000};//

        for(int i = 0 ; i < lengtharray.length ; i++) {
            System.out.println("Generating length=" + lengtharray[i]);
            for(String label : networks.keySet()) {
                //if(lengtharray[i] != 250000 || !label.equals("networkC")) continue;

                Network trueNetwork = Networks.readNetwork(networks.get(label));

                double constTheta = 0.04;
                for(Object nodeObject : Networks.postTraversal(trueNetwork)) {
                    NetNode node = (NetNode) nodeObject;
                    for(Object parentObject :  node.getParents()) {
                        NetNode parent = (NetNode) parentObject;
                        node.setParentSupport(parent, constTheta);
                        //node.setParentDistance(parent, node.getParentDistance(parent) * constTheta / 2.0);
                    }
                }
                trueNetwork.getRoot().setRootPopSize(constTheta);
                System.out.println(trueNetwork.toString());


                String filename = "../seconddata/" + label + "_" + lengtharray[i] + "_ac_seed12345678.snp";
                SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, Utils._SEED + lengtharray[i]);
                Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, lengtharray[i], true);
                List<String> taxa = new ArrayList<>(onesnp.keySet());
                Collections.sort(taxa);
                try {
                    PrintWriter out = new PrintWriter(filename);
                    for (String taxon : taxa) {
                        out.println(taxon + " " + onesnp.get(taxon));
                    }
                    out.close();
                } catch (IOException ioe) {
                    ioe.printStackTrace();
                }

            }
        }
    }

    static public void processMosquitoData(){
        Random random = new Random(Utils._SEED);

        Map<String, StringBuilder> samples = new HashMap<>();
        samples.put("C", new StringBuilder());
        samples.put("G", new StringBuilder());
        samples.put("A", new StringBuilder());
        samples.put("R", new StringBuilder());
        samples.put("Q", new StringBuilder());
        samples.put("L", new StringBuilder());


        List<String> chromosomeLabels = new ArrayList<>();
        chromosomeLabels.add("2L");
        chromosomeLabels.add("2R");
        chromosomeLabels.add("3L");
        chromosomeLabels.add("3R");

        int nPolymorphic = 0;

        for(String chromosomeLabel : chromosomeLabels) {
            System.out.println("Chromosome: " + chromosomeLabel);
            int number = 0;
            Map<Integer, String> paths = new TreeMap<>();
            for (File file : new File("../mosquito/" + chromosomeLabel + "/").listFiles()) {
                if (file.isFile() && file.getName().endsWith(".fa")) {
                    String filename = file.getName();
                    Scanner scanner = new Scanner(filename);
                    scanner.useDelimiter("\\.");
                    scanner.next();
                    int position = Integer.parseInt(scanner.next());
                    int length = Integer.parseInt(scanner.next());
                    paths.put(position, file.getAbsolutePath());
                    number++;
                }
            }
            int nextPos = 0;
            int lastPos = 0;
            for(Integer position : paths.keySet()) {
                Map<String, String> sequence = new HashMap<>();
                int length = 0;

                try {
                    File file = new File(paths.get(position));
                    Scanner scanner = new Scanner(file);

                    while (scanner.hasNext()) {
                        String taxon = scanner.next();
                        if (taxon.length() == 0) break;
                        if(taxon.charAt(0) != '>') {
                            throw new IOException("There may be problem about taxon format!");
                        }
                        taxon = taxon.substring(1);
                        String s = scanner.next();
                        sequence.put(taxon, s);
                        if(length != 0 && length != s.length())
                            throw new IOException("Sequence length inconsistant!");
                        else if(length == 0)
                            length = s.length();
                    }

                    scanner.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }

                int start = (nextPos < position ? 0 : nextPos - position);
                for(int i = start ; i < length ; i++) {
                    Map<String, Character> site = new HashMap<>();
                    Map<Character, Integer> states = new TreeMap<>();
                    for(String taxon : samples.keySet()) {
                        Character c = sequence.get(taxon).charAt(i);
                        if(!states.containsKey(c))
                            states.put(c, 0);

                        states.put(c, states.get(c) + 1);
                        site.put(taxon, c);
                    }
                    System.out.println(states.size() + " " + states.keySet().iterator().next());
                    if(states.size() <= 2 && !states.containsKey('-') && !states.containsKey('N')) {
                        List<Character> statesList = new ArrayList<>(states.keySet());
                        Collections.shuffle(statesList, random);
                        int draw = random.nextInt(2);

                        Map<Character, Character> image = new HashMap<>();
                        image.put(statesList.get(0), (draw == 0 ? '0': '1'));
                        //if(states.get(statesList.get(0)) <= 1) continue;
                        if(statesList.size() >= 2) nPolymorphic++;
                        for(int j = 1 ; j < statesList.size() ; j++) {
                            image.put(statesList.get(j), (draw == 0 ? '1' : '0'));
                            //if(states.get(statesList.get(1)) <= 1) continue;
                        }

                        for(String taxon : site.keySet()) {
                            if(!site.get(taxon).equals('A') && !site.get(taxon).equals('T') && !site.get(taxon).equals('C') && !site.get(taxon).equals('G'))
                                throw new RuntimeException("error char " + site.get(taxon));
                            samples.get(taxon).append(image.get(site.get(taxon)));
                        }
                        System.out.println("Sampled at " + (position + i) + " (+" + (position + i - lastPos) + ")" + " " + statesList.size());
                        lastPos = position + i;
                        nextPos = position + i + 2500;
                        break;
                    }
                }
            }

        }

        System.out.println("Polymorphic sites: " + nPolymorphic);

        try {
            PrintWriter out = new PrintWriter("../mosquito/sampled_2k5.snp");
            List<String> taxa = new ArrayList<>(samples.keySet());
            //Collections.sort(taxa);
            for (String taxon : taxa) {
                out.println(taxon + " " + samples.get(taxon).toString());
                System.out.println(taxon + " " + samples.get(taxon).length());
            }
            out.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    static public void verifyLikelihood() {
        double pi0 = 0.5;
        double pi1 = 1- pi0;

        SNAPPLikelihood.useOnlyPolymorphic = true;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {pi0, pi1}, new double[] {1.0/ (2.0 * pi0)});

        double gamma = 1 - 0.3;

        double y = 1.5;
        double x = 0.5;

        Network<Object> trueNetwork;
        trueNetwork = Networks.readNetwork("((((b:1.0,c:1.0)I4:" + y + ")I3#H1:0.1::" + (1 - gamma) + ",a:" + (1.1 + y) + ")I1:" + x + ",(I3#H1:0.1::" + gamma + ",d:" + (1.1 + y) + ")I2:"+ x +")I0;");
        trueNetwork = Networks.readNetwork("(((((A:0.7)I6#H1:1.3::0.8,Q:2.0)I4:1.0,L:3.0)I3:1.0,R:4.0)I2:1.0,(G:2.0,(I6#H1:0.7::0.2,C:1.4)I5:0.6)I1:3.0)I0;");
        trueNetwork = Networks.readNetwork("(((((Q:2.0,A:2.0)I4:1.0,L:3.0)I3:0.5)I8#H1:0.5::0.7,R:4.0)I2:1.0,(I8#H1:0.5::0.3,(G:2.0,C:2.0)I1:2.0)I7:1.0)I0;");
        //trueNetwork = Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,A:1.0)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,(G:1.0,C:1.0)I1:4.0)I0;");
        //trueNetwork = Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,(A:0.5)I6#H2:0.5::0.8)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,((I6#H2:0.5::0.2,C:1.0)I5:1.0,G:2.0)I1:3.0)I0;");

        System.out.println(trueNetwork.toString());

        double constTheta = 0.036;
        double mu = 1.0;
        for(Object nodeObject : Networks.postTraversal(trueNetwork)) {
            NetNode node = (NetNode) nodeObject;
            for(Object parentObject :  node.getParents()) {
                NetNode parent = (NetNode) parentObject;
                node.setParentSupport(parent, constTheta);
                node.setParentDistance(parent, node.getParentDistance(parent) * constTheta / 2.0);
            }
        }
        trueNetwork.getRoot().setRootPopSize(constTheta);


        trueNetwork = Networks.readNetwork("((((A:0.09791508158706645,Q:0.09791508158706645):0.04737303547719145,L:0.1452881170642579):0.023386428388076835,R:0.16867454545233473):0.016009187542328024,(G:0.10020765998425706,C:0.10020765998425706):0.0844760730104057);");
        trueNetwork.getRoot().setRootPopSize(0.036);
        int nameCount = 0;
        for(Object nodeObject : Networks.postTraversal(trueNetwork)) {
            NetNode node = (NetNode) nodeObject;
            for(Object parentObject :  node.getParents()) {
                NetNode parent = (NetNode) parentObject;
                node.setParentSupport(parent, 0.036);
            }
        }
        for(Object node : trueNetwork.dfs()) {
            NetNode mynode = (NetNode) node;
            if(mynode.getName().equals("")) {
                mynode.setName("I" + nameCount);
                nameCount++;
            }
        }


        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, null);
        Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, 100000, true);

        List<Map<String, String>> snpdata = new ArrayList<>();
        snpdata.add(onesnp);
        List<Alignment> alns = new ArrayList<>();
        for(Map<String, String> input : snpdata) {
            for(String allele : input.keySet()) {
                System.out.println(allele + " " + input.get(allele));
            }
            System.out.println();
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

            aln.setCache(cache);
            alns.add(aln);

        }

        SNAPPLikelihood.ALGORITHM = 0;
        System.out.println(SNAPPLikelihood.computeSNAPPLikelihood(trueNetwork, null, alns, BAGTRModel));
        //SNAPPLikelihood.ALGORITHM = 1;
        //System.out.println(SNAPPLikelihood.computeSNAPPLikelihood(trueNetwork, null, alns, BAGTRModel));
        SNAPPLikelihood.ALGORITHM = 2;
        System.out.println(SNAPPLikelihood.computeSNAPPLikelihood(trueNetwork, null, alns, BAGTRModel));
    }
}
