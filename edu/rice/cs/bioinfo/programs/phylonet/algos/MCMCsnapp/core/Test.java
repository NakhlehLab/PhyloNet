package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.core;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution.SNAPPLikelihood;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.clustering.DataGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.InferMLNetworkFromSequences;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimSNPInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.*;
import java.lang.reflect.Array;
import java.util.*;
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
        //checkResults();
        //testWithPopSize(args);
        //testMCMC();
        testWithDingqiaoNetwork(args);
        //checkDingqiaoNetworkResults();
        //generateSNPdata();
    }

    public static void testMCMC() {
        /*Utils._MC3_CHAINS = new ArrayList<>();
        Utils._MC3_CHAINS.add(1.0);
        Utils._MC3_CHAINS.add(2.0);
        Utils._MC3_CHAINS.add(4.0);*/
        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {0.5, 0.5}, new double[] {1.0});

        double gamma = 1 - 0.3;

        double y = 1.0;
        double x = 1000.0;

        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, null);
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


        Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, 500000, true);
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
        MC3Core run = new MC3Core(alns, BAGTRModel);
        run.run();
    }

    public static void testWithSNAPPSimulation(String[] args) {
        int curPart = Integer.parseInt(args[0]);
        int number = 0;
        List<String> paths = new ArrayList<>();
        for (File file : new File("/scratch/jz55/PaperSimulations/Simulation1Clean").listFiles()) {
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

    public static void testWithDingqiaoNetwork(String[] args) {
        int curPart = Integer.parseInt(args[0]);
        int number = 0;
        List<String> paths = new ArrayList<>();
        for (File file : new File("/scratch/jz55/data").listFiles()) {
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

        rate[0] = 1 / (2*pi[0]);

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(pi, rate);

        MC3Core run = new MC3Core(alns, BAGTRModel);
        run.run();

        System.out.println(String.format("Total elapsed time : %2.5f s\n",
                (double) (System.currentTimeMillis() - startTime) / 1000.0));
    }

    public static void checkSNAPPSimulationResults() {

        int number = 0;
        int correct = 0;
        List<String> paths = new ArrayList<>();
        for (File file : new File("/scratch/jz55").listFiles()) {
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
                    System.out.println("Incorrect: " + inputfile);
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
        for (File file : new File("/Users/zhujiafan/").listFiles()) {
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

                        if(!networkLabel.equals("networkC") || !allowConstantLabel.equals("ac"))
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

        Map<String, Network> networks = new HashMap<>();
        networks.put("networkA", Networks.readNetwork("(((((A:0.7)I6#H1:1.3::0.8,Q:2.0)I4:1.0,L:3.0)I3:1.0,R:4.0)I2:1.0,(G:2.0,(I6#H1:0.7::0.2,C:1.4)I5:0.6)I1:3.0)I0;"));
        networks.put("networkB", Networks.readNetwork("(((((Q:2.0,A:2.0)I4:1.0,L:3.0)I3:0.5)I8#H1:0.5::0.7,R:4.0)I2:1.0,(I8#H1:0.5::0.3,(G:2.0,C:2.0)I1:2.0)I7:1.0)I0;"));
        networks.put("networkC", Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,A:1.0)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,(G:1.0,C:1.0)I1:4.0)I0;"));
        networks.put("networkD", Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,(A:0.5)I6#H2:0.5::0.8)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,((I6#H2:0.5::0.2,C:1.0)I5:1.0,G:2.0)I1:3.0)I0;"));


        int lengtharray[] = new int[]{100, 1000, 10000, 100000, 250000, 500000, 1000000, 2000000, 4000000};//

        for(int i = 0 ; i < lengtharray.length ; i++) {
            System.out.println("Generating length=" + lengtharray[i]);
            for(String label : networks.keySet()) {
                Network trueNetwork = networks.get(label).clone();

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

                String filename = "../data/" + label + "_" + lengtharray[i] + "_nc_seed12345678.snp";
                SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, Utils._SEED + lengtharray[i]);
                Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, lengtharray[i], false);
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
}
