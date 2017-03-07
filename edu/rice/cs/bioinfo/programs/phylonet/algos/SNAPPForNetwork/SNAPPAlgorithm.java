package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm.NucleotideProbabilityAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.GTRModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.RateModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.OneNucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import org.apache.commons.math3.util.ArithmeticUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.io.File;
import java.io.StringReader;
import java.util.*;

public class SNAPPAlgorithm extends NucleotideProbabilityAlgorithm {

    //private MatrixQ Q;
    private QParameters Q;
    private Network speciesNetwork;
    private Map<String, String> allele2species;

    public SNAPPAlgorithm(Network theSpeciesNetwork, RateModel rModel, Double theta){

        speciesNetwork = theSpeciesNetwork;
        Q = new QParameters(rModel, R.maxLineages, theta);
        if(speciesNetwork.getRoot().getData()==null) {
            Networks.removeBinaryNodes(speciesNetwork);
            for (Object node : Networks.postTraversal(speciesNetwork)) {
                ((NetNode<SNAPPData[]>) node).setData(new SNAPPData[1]);
            }
        }
        allele2species = null;
    }

    public SNAPPAlgorithm(Network theSpeciesNetwork, Map<String, String> alleleMapping, RateModel rModel){

        new SNAPPAlgorithm(theSpeciesNetwork, alleleMapping, rModel, 0.02);
    }

    public SNAPPAlgorithm(Network theSpeciesNetwork, Map<String, String> alleleMapping, RateModel rModel, Double theta){
        speciesNetwork = theSpeciesNetwork;
        if(speciesNetwork.getRoot().getData()==null) {
            Networks.removeBinaryNodes(speciesNetwork);
            for (Object node : Networks.postTraversal(speciesNetwork)) {
                ((NetNode<SNAPPData[]>) node).setData(new SNAPPData[1]);
            }
        }
        allele2species = alleleMapping;
        if(allele2species == null) {
            Q = new QParameters(rModel, speciesNetwork.getLeafCount(), theta);
        }else{
            Q = new QParameters(rModel, alleleMapping.size(), theta);
        }
    }


    /*
    public SNAPPAlgorithm(Network theSpeciesNetwork, Map<String, String> alleleMapping, RateModel rModel, double theta){
        allele2species = alleleMapping;
        speciesNetwork = theSpeciesNetwork;
        if(allele2species == null) {
            Q = new MatrixQ(rModel, speciesNetwork.getLeafCount(), theta);
        }else{
            Q = new MatrixQ(rModel, alleleMapping.size(), theta);
        }
    }
    */

    public double getProbability(NucleotideObservation dna) {
        return getProbability(dna, 0);
    }

    public double getProbability(RPattern pattern) {
        return getProbability(pattern, 0);
    }

    public double getProbability(NucleotideObservation dna, int siteID) {
        double resultVal= Algorithms.getProbabilityObservationGivenNetwork(speciesNetwork, allele2species, dna, Q, siteID);

        if (Double.isNaN(resultVal)||Double.isInfinite(resultVal))
            throw new RuntimeException("Result is not finite");
        return resultVal;
    }

    public double getProbability(RPattern pattern, int siteID) {
        Q._M = pattern.sumLineages();
        double resultVal= Algorithms.getProbabilityObservationGivenNetwork(speciesNetwork, pattern, Q, siteID);

        if (Double.isNaN(resultVal)||Double.isInfinite(resultVal))
            throw new RuntimeException("Result is not finite");
        return resultVal;
    }



    public double updateProbability(int siteID, NucleotideObservation obs, Map<String, String> allele2species, Map<NetNode,Boolean> node2update) {

        double resultVal= Algorithms.updateProbabilityObservationGivenNetwork(speciesNetwork, Q, siteID, obs, allele2species, node2update);

        if (Double.isNaN(resultVal)||Double.isInfinite(resultVal))
            throw new RuntimeException("Result is not finite");
        return resultVal;
    }



    public static void main(String[] args)
    {

/*        String tree = "(((1:0.0)n2:7.425226489143906,(2:0.0)n3:7.425226489143906)n1:3.5972158884074084,(3:0.0)n6:11.022442377551315)n0";

        Tree theTree = Trees.toTree(tree);

        GTRModel model = new GTRModel(new double[]{.1, .2, .3, .4},new double[]{1,1.5,2,2.5,3,3.5});
        SNAPPAlgorithm algo = new SNAPPAlgorithm(theTree,model);

        Map<String,Character> map = new HashMap<>();

        map.put("1",'A');
        map.put("2",'A');
        map.put("3",'A');

        System.out.println(algo.getProbability(new OneNucleotideObservation(map)));



        R.dims = 1;
        testBiallelic();
*/
        //R.dims = 3;

        //testSpeed();
        //testFourAlleles();
        //testFourAllelesTwoUnderReticulation();
        //testBiallelicLarge();
        //testBiallelic();
        //checkWithOriginalPaperData();
        //testMultipleAllelesBiallelic2();
        //testMultipleAlleles();
        testMultipleAllelesBiallelic();
        //R.dims = 3;

        //testFourAllelesMedium();


        //R.dims = 3;
        //testFourAllelesLarge();

    }

    private static void testFourAlleles() {
        R.dims = 3;
        double sum =0;
        char[] nucleotides = {'A','C','T','G'};
        for (char a : nucleotides)
            for (char b: nucleotides)
                for (char c: nucleotides)
                {
                    //Network speciesNetworkTopology = Networks.readNetwork("((A:0.04,X#H1:0.02::0.7)n3:0.02,((B:0.02)X#H1:0.02::0.3,C:0.04)n1:0.02)root;");
                    //Network speciesNetworkTopology = Networks.readNetwork("((B:0.09719404370089872)#H1:0.46845495567048107::0.6928009147744886,((#H1:0.30331736190632846::0.3071990852255114,A:0.4005114056072272):0.018286593135532436,C:0.41879799874275964):0.14685100062862017);");
                    Network speciesNetworkTopology = Networks.readNetwork("((((C:0.07,(B:0.06)I7#H1:0.01::0.2)I4:0.02,(I7#H1:0.02::0.8)I6#H2:0.01::0.2)I3:0.02,(I6#H2:0.02::0.8)I5#H3:0.01::0.2)I2:0.07,(I5#H3:0.06::0.8,A:0.16)I1:0.02)I0;");

                    //Network speciesNetworkTopology = Networks.string2network("(A:2,(B:1,C:1)n1:1)root;");
                    //Network speciesNetworkTopology = Networks.string2network("(A:0.04,(B:0.02,C:0.02)n1:0.02)root;");
                    //Network speciesNetworkTopology = Networks.string2network("(A:.5,(B:.2,C:.2)n1:.3)n0;");
                    //Network speciesNetworkTopology = Networks.string2network("((A:.2,X#H1:0::0.2)n2:0.3,((B:.1)X#H1:0::0.8,C:.2)n1:.3)root;");
                    Map<String, Character> colorMap = new HashMap<String, Character>();
                    colorMap.put("A", a);
                    colorMap.put("B", b);
                    colorMap.put("C", c);

                    OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);
                    GTRModel gtrModel = new GTRModel(new double[]{.1, .2, .3, .4},new double[]{1,1.5,2,2.5,3,3.5});

                    SNAPPAlgorithm run = new SNAPPAlgorithm(speciesNetworkTopology, gtrModel, 0.02);
                    //System.out.println("Observation:  " + a +',' + b + ','+ c);
                    System.out.println("Probability:  " + run.getProbability(converter));
                    sum += run.getProbability(converter);
                }


        System.out.println("Sum: " + sum);

    }

    private static void testMultipleAlleles() {
        R.dims = 3;
        double sum =0;
        char[] nucleotides = {'A','C','T','G'};
        int index = 0;
        HashSet<String> tried = new HashSet<String>();
        for (char a : nucleotides)
            for (char b1: nucleotides)
                for (char b2: nucleotides)
                    for (char b3: nucleotides)
                        for (char c: nucleotides)
                        {
                            //if(index++<4)continue;
                            //Network speciesNetworkTopology = Networks.string2network("(A:.5,(B:.2,C:.2)n1:.3)n0;");
                            //Network speciesNetworkTopology = Networks.string2network("((A:0.3,X#H1:0::0.5)n2:0.1,(((B:0.1,C:0.1)n3:0.1)X#H1:0::0.5,D:0.3)n1:0.1)root;");
                            Network speciesNetworkTopology = Networks.readNetwork("((A:.2,X#H1:0::0.2)n2:0.3,((B:.1)X#H1:0::0.8,C:.2)n1:.3)root;");
                            //Network speciesNetworkTopology = Networks.readNetwork("((((C:0.07,(B:0.06)I7#H1:0.01::0.2)I4:0.02,(I7#H1:0.02::0.8)I6#H2:0.01::0.2)I3:0.02,(I6#H2:0.02::0.8)I5#H3:0.01::0.2)I2:0.07,(I5#H3:0.06::0.8,A:0.16)I1:0.02)I0;");
                            //Network speciesNetworkTopology = Networks.string2network("((a:.2,X#H1:0::0.2)n2:0.3,(((b1:0,b2:0)n3:.1)X#H1:0::0.8,c:.2)n1:.3)root;");

                    /*
                    for (Object node : Networks.postTraversal(speciesNetworkTopology)) {
                        for (Object child : ((NetNode)node).getChildren()) {
                            ((NetNode)child).setParentProbability((NetNode)node, NetNode.NO_PROBABILITY);
                            ((NetNode)child).setParentDistance((NetNode)node, NetNode.NO_DISTANCE);
                            ((NetNode)child).setParentSupport((NetNode)node, NetNode.NO_PROBABILITY);
                        }
                    }
                    System.out.println(Networks.network2string(speciesNetworkTopology));
                    */
                        /*
                        String seq = String.valueOf(a)+String.valueOf(b2)+String.valueOf(b1)+String.valueOf(c);
                        //System.out.println(seq);
                        if(tried.contains(seq) && b1!=b2){
                            System.out.println(seq + " vs. " + a+b1+b2+c+"");
                            continue;
                        }
                        else{
                            tried.add(String.valueOf(a)+String.valueOf(b1)+String.valueOf(b2)+String.valueOf(c));
                        }
                        */

                            Map<String, Character> colorMap = new HashMap<String, Character>();
                            colorMap.put("a", a);
                            colorMap.put("b1", b1);
                            colorMap.put("b2", b2);
                            colorMap.put("b3", b3);
                            colorMap.put("c", c);

                            Map<String, String> allele2species = new HashMap<String, String>();
                            allele2species.put("a","A");
                            allele2species.put("b1","B");
                            allele2species.put("b2","B");
                            allele2species.put("b3","B");
                            allele2species.put("c","C");
                            //allele2species = null;

                            OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);


                            GTRModel gtrModel = new GTRModel(new double[]{.1, .2, .3, .4},new double[]{1,1.5,2,2.5,3,3.5});



                            SNAPPAlgorithm run = new SNAPPAlgorithm(speciesNetworkTopology, allele2species, gtrModel, 1.0);
                            System.out.println("Observation:  " + a +',' + b1 + ','+ b2 + "," +  b3 + "," + c);
                            //System.out.println("Probability:  " + run.getProbability(converter) + "\n");
                            sum += run.getProbability(converter);

                        }


        System.out.println("Sum: " + sum);




    }

    private static void testMultipleAllelesBiallelic() {
        R.dims = 1;
        double sum =0;
        char[] nucleotides = {'0','1'};
        int index = 0;
        HashSet<String> tried = new HashSet<String>();
        for (char a : nucleotides)
            for (char b1: nucleotides)
                for (char b2: nucleotides)
                    for (char b3: nucleotides)
                        for (char c: nucleotides)
                        {
                            //if(index++<4)continue;
                            //Network speciesNetworkTopology = Networks.readNetwork("(A:.5,(B:.2,C:.2)n1:.3)n0;");
                            Network speciesNetworkTopology = Networks.readNetwork("((((C:0.07,(B:0.06)I7#H1:0.01::0.2)I4:0.02,(I7#H1:0.02::0.8)I6#H2:0.01::0.2)I3:0.02,(I6#H2:0.02::0.8)I5#H3:0.01::0.2)I2:0.07,(I5#H3:0.06::0.8,A:0.16)I1:0.02)I0;");
                            //Network speciesNetworkTopology = Networks.string2network("((A:0.3,X#H1:0::0.5)n2:0.1,(((B:0.1,C:0.1)n3:0.1)X#H1:0::0.5,D:0.3)n1:0.1)root;");
                            //Network speciesNetworkTopology = Networks.readNetwork("((A:.2,X#H1:0::0.2)n2:0.3,((B:.1)X#H1:0::0.8,C:.2)n1:.3)root;");
                            //Network speciesNetworkTopology = Networks.string2network("((a:.2,X#H1:0::0.2)n2:0.3,(((b1:0,b2:0)n3:.1)X#H1:0::0.8,c:.2)n1:.3)root;");

                    /*
                    for (Object node : Networks.postTraversal(speciesNetworkTopology)) {
                        for (Object child : ((NetNode)node).getChildren()) {
                            ((NetNode)child).setParentProbability((NetNode)node, NetNode.NO_PROBABILITY);
                            ((NetNode)child).setParentDistance((NetNode)node, NetNode.NO_DISTANCE);
                            ((NetNode)child).setParentSupport((NetNode)node, NetNode.NO_PROBABILITY);
                        }
                    }
                    System.out.println(Networks.network2string(speciesNetworkTopology));
                    */
                        /*
                        String seq = String.valueOf(a)+String.valueOf(b2)+String.valueOf(b1)+String.valueOf(c);
                        //System.out.println(seq);
                        if(tried.contains(seq) && b1!=b2){
                            System.out.println(seq + " vs. " + a+b1+b2+c+"");
                            continue;
                        }
                        else{
                            tried.add(String.valueOf(a)+String.valueOf(b1)+String.valueOf(b2)+String.valueOf(c));
                        }
                        */

                            Map<String, Character> colorMap = new HashMap<String, Character>();
                            colorMap.put("a", a);
                            colorMap.put("b1", b1);
                            colorMap.put("b2", b2);
                            colorMap.put("b3", b3);
                            colorMap.put("c", c);

                            Map<String, String> allele2species = new HashMap<String, String>();
                            allele2species.put("a","A");
                            allele2species.put("b1","B");
                            allele2species.put("b2","B");
                            allele2species.put("b3","B");
                            allele2species.put("c","C");
                            //allele2species = null;

                            OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);


                            //GTRModel gtrModel = new GTRModel(new double[]{.1, .2, .3, .4},new double[]{1,1.5,2,2.5,3,3.5});
                            BiAllelicGTR gtrModel = new BiAllelicGTR(new double[] {0.5, 0.5}, new double[] {1.0});



                            SNAPPAlgorithm run = new SNAPPAlgorithm(speciesNetworkTopology, allele2species, gtrModel, 1.0);
                            System.out.println("Observation:  " + a +',' + b1 + ','+ b2 + "," +  b3 + "," + c);
                            //System.out.println("Probability:  " + run.getProbability(converter) + "\n");
                            sum += run.getProbability(converter);

                        }


        System.out.println("Sum: " + sum);




    }

    private static void testMultipleAllelesBiallelic2() {
        R.dims = 1;
        double sum =0;
        char[] nucleotides = {'0','1'};
        int index = 0;
        HashSet<String> tried = new HashSet<String>();
        for (char a : nucleotides)
            for (char b1: nucleotides)
                for (char b2: nucleotides)
                    for (char b3: nucleotides)
                        for (char c1: nucleotides)
                            for(char c2 : nucleotides)
                            {
                                //if(index++<4)continue;
                                Network speciesNetworkTopology = Networks.readNetwork("(A:.5,(B:.2,C:.2)n1:.3)n0;");
                                //Network speciesNetworkTopology = Networks.string2network("((A:0.3,X#H1:0::0.5)n2:0.1,(((B:0.1,C:0.1)n3:0.1)X#H1:0::0.5,D:0.3)n1:0.1)root;");
                                //Network speciesNetworkTopology = Networks.readNetwork("((A:.2,X#H1:0::0.2)n2:0.3,((B:.1)X#H1:0::0.8,C:.2)n1:.3)root;");
                                //Network speciesNetworkTopology = Networks.string2network("((a:.2,X#H1:0::0.2)n2:0.3,(((b1:0,b2:0)n3:.1)X#H1:0::0.8,c:.2)n1:.3)root;");

                        /*
                        for (Object node : Networks.postTraversal(speciesNetworkTopology)) {
                            for (Object child : ((NetNode)node).getChildren()) {
                                ((NetNode)child).setParentProbability((NetNode)node, NetNode.NO_PROBABILITY);
                                ((NetNode)child).setParentDistance((NetNode)node, NetNode.NO_DISTANCE);
                                ((NetNode)child).setParentSupport((NetNode)node, NetNode.NO_PROBABILITY);
                            }
                        }
                        System.out.println(Networks.network2string(speciesNetworkTopology));
                        */
                            /*
                            String seq = String.valueOf(a)+String.valueOf(b2)+String.valueOf(b1)+String.valueOf(c);
                            //System.out.println(seq);
                            if(tried.contains(seq) && b1!=b2){
                                System.out.println(seq + " vs. " + a+b1+b2+c+"");
                                continue;
                            }
                            else{
                                tried.add(String.valueOf(a)+String.valueOf(b1)+String.valueOf(b2)+String.valueOf(c));
                            }
                            */

                                Map<String, Character> colorMap = new HashMap<String, Character>();
                                colorMap.put("a", a);
                                colorMap.put("b1", b1);
                                colorMap.put("b2", b2);
                                colorMap.put("b3", b3);
                                colorMap.put("c1", c1);
                                colorMap.put("c2", c2);

                                Map<String, String> allele2species = new HashMap<String, String>();
                                allele2species.put("a","A");
                                allele2species.put("b1","B");
                                allele2species.put("b2","B");
                                allele2species.put("b3","B");
                                allele2species.put("c1","C");
                                allele2species.put("c2","C");
                                //allele2species = null;

                                double modifier = 1.0;
                                int countB = (b1 == '1' ? 1 : 0) + (b2 == '1' ? 1 : 0) + (b3 == '1' ? 1 : 0);
                                int countC = (c1 == '1' ? 1 : 0) + (c2 == '1' ? 1 : 0);
                                //modifier /= ArithmeticUtils.binomialCoefficient(3, countB) * ArithmeticUtils.binomialCoefficient(2, countC);

                                OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);


                                //GTRModel gtrModel = new GTRModel(new double[]{.1, .2, .3, .4},new double[]{1,1.5,2,2.5,3,3.5});
                                BiAllelicGTR gtrModel = new BiAllelicGTR(new double[] {0.5, 0.5}, new double[] {1.0});



                                SNAPPAlgorithm run = new SNAPPAlgorithm(speciesNetworkTopology, allele2species, gtrModel, 1.0);
                                System.out.println("Observation:  " + a +',' + b1 + ','+ b2 + "," +  b3 + "," + c1 + "," + c2);
                                double prob = run.getProbability(converter) * modifier;
                                System.out.println("Probability:  " + prob + "\n");
                                sum += prob;

                            }


        System.out.println("Sum: " + sum);




    }

    private static void testFourAllelesTwoUnderReticulation() {
        R.dims = 3;
        double sum =0;
        char[] nucleotides = {'A','C','T','G'};
        for (char a : nucleotides)
            for (char b: nucleotides)
                for (char c: nucleotides)
                    for (char d: nucleotides)
                    {
                        //Network speciesNetworkTopology = Networks.string2network("(A:.5,(B:.2,C:.2)n1:.3)n0;");
                        //Network speciesNetworkTopology = Networks.string2network("((A:0.3,X#H1:0::0.5)n2:0.1,(((B:0.1,C:0.1)n3:0.1)X#H1:0::0.5,D:0.3)n1:0.1)root;");
                        //Network speciesNetworkTopology = Networks.readNetwork("((A:0.02,X#H1:0.02::0.3)n2:0.02,(((B:0.02,C:0.02)n3:0.02)X#H1:0.02::0.7,D:0.02)n1:0.02)root;");
                        Network speciesNetworkTopology = Networks.readNetwork("((A:1,X#H1:1::0.3)n2:1,(((B:1,C:1)n3:1)X#H1:1::0.7,D:1)n1:1)root;");
                        double constTheta = 0.036;
                        for(Object nodeObject : Networks.postTraversal(speciesNetworkTopology)) {
                            NetNode node = (NetNode) nodeObject;
                            for(Object parentObject :  node.getParents()) {
                                NetNode parent = (NetNode) parentObject;
                                node.setParentSupport(parent, constTheta);
                                node.setParentDistance(parent, node.getParentDistance(parent) * constTheta / 2.0);
                            }
                        }
                        speciesNetworkTopology.getRoot().setRootPopSize(constTheta);
                    /*
                    for (Object node : Networks.postTraversal(speciesNetworkTopology)) {
                        for (Object child : ((NetNode)node).getChildren()) {
                            ((NetNode)child).setParentProbability((NetNode)node, NetNode.NO_PROBABILITY);
                            ((NetNode)child).setParentDistance((NetNode)node, NetNode.NO_DISTANCE);
                            ((NetNode)child).setParentSupport((NetNode)node, NetNode.NO_PROBABILITY);
                        }
                    }
                    System.out.println(Networks.network2string(speciesNetworkTopology));
                    */

                        Map<String, Character> colorMap = new HashMap<String, Character>();
                        colorMap.put("A", a);
                        colorMap.put("B", b);
                        colorMap.put("C", c);
                        colorMap.put("D", d);


                        OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);


                        GTRModel gtrModel = new GTRModel(new double[]{.1, .2, .3, .4},new double[]{1,1.5,2,2.5,3,3.5});
                        //GTRModel gtrModel = new GTRModel(new double[]{1, 0, 0, 0},new double[]{0,0,0,0,0,0});



                        SNAPPAlgorithm run = new SNAPPAlgorithm(speciesNetworkTopology, gtrModel, 0.02);
                        System.out.println("Observation:  " + a +',' + b + ','+ c + "," + d);
                        System.out.println("Probability:  " + run.getProbability(converter) + "\n");
                        sum += run.getProbability(converter);

                    }


        System.out.println("Sum: " + sum);




    }




    private static void testFourAllelesMedium() {
        R.dims = 3;
        double sum =0;
        long start = System.currentTimeMillis();
        char[] nucleotides = {'A','C','T','G'};
        for (char a : nucleotides)
            for (char b: nucleotides)
                for (char c: nucleotides)
                    for(char d: nucleotides)
                        for(char e: nucleotides)
                        {
                            //Network speciesNetworkTopology = Networks.string2network("((((B:1)X#H1:0::0.2,A:1)i1:1,X#H1:1::0.8)i2:1,(((D:1)Y#H2:0::0.4,E:1)i3:1,(C:1,Y#H2:0::0.6)i4:1)i6:1)i5;");
                            Network speciesNetworkTopology = Networks.readNetwork("((((C:1.0)Y#H1:1.0::0.5,(B:1.0)X#H2:1.0::0.5)i4:1.0,(A:1.0,X#H2:1.0::0.5)i1:1.0)i5:1.0,((D:1.0,E:1.0)i2:1.0,Y#H1:1.0::0.5)i3:1.0)i6;");

                            //Network speciesNetworkTopology = Networks.string2network("((((B:1)X#H1:1::0,A:1)i1:1,X#H1:1::1)i2:1,((D:1)Y#H2:1::0,E:1)i3:1,(C:1,Y#H2:1::1)i4:1)i5;");
                            //Network speciesNetworkTopology = Networks.string2network("((A:.2, B:.2)n1:.3,((C:.1,D:.1)n0:.3, E:.4)n2:.1)n3");
                            Map<String, Character> colorMap = new HashMap<String, Character>();
                            colorMap.put("A", a);
                            colorMap.put("B", b);
                            colorMap.put("C", c);
                            colorMap.put("D", d);
                            colorMap.put("E", e);


                            OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);


                            GTRModel gtrModel = new GTRModel(new double[]{.1, .2, .3, .4},new double[]{1,1.5,2,2.5,3,3.5});



                            SNAPPAlgorithm run = new SNAPPAlgorithm(speciesNetworkTopology, gtrModel, 1.0);
                            System.out.println("Observation:  " + a +',' + b + ','+ c +','+ d + ','+ e);
                            //System.out.println("Probability:  " + run.getProbability(converter) + "\n");
                            sum += run.getProbability(converter);
                        }


        System.out.println((System.currentTimeMillis()-start)/1000.0);
        System.out.println("Sum: " + sum);




    }

    private static void testFourAllelesLarge() {
        int i = 0;
        double sum =0;
        Network speciesNetworkTopology = Networks.readNetwork("(((A:.2, B:.2)n1:.3,((C:.1,D:.1)n0:.3, E:.4)n2:.1)n3:.1,(F:.2, G:.2)n4:.4)");
        GTRModel gtrModel = new GTRModel(new double[]{.1, .2, .3, .4},new double[]{1,1.5,2,2.5,3,3.5});

        SNAPPAlgorithm run = new SNAPPAlgorithm(speciesNetworkTopology, gtrModel, 1.0);

        long last = System.currentTimeMillis();
        char[] nucleotides = {'A','C','T','G'};
        for (char a : nucleotides)
            for (char b: nucleotides)
                for (char c: nucleotides)
                    for(char d: nucleotides)
                        for(char e: nucleotides)
                            for (char f : nucleotides)
                                for (char g: nucleotides)
                                {
                                    Map<String, Character> colorMap = new HashMap<String, Character>();
                                    colorMap.put("A", a);
                                    colorMap.put("B", b);
                                    colorMap.put("C", c);
                                    colorMap.put("D", d);
                                    colorMap.put("E", e);
                                    colorMap.put("F", f);
                                    colorMap.put("G", g);

                                    OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);

                                    //System.out.println("Observation:  " + a +',' + b + ','+ c +','+ d + ','+ e+ ','+ f+ ','+ g);
                                    double result = run.getProbability(converter);
                                    //System.out.println("Probability:  " + result + "\n");
                                    sum += result;
                                    i++;

                                    if (i%2000 == 0)
                                    {
                                        long cur = System.currentTimeMillis();
                                        System.out.println(cur - last);
                                        last = cur;

                                    }

                                }
        System.out.println("Sum: " + sum);
    }

    private static void testBiallelic() {
        double sum =0;
        for (int i = 0; i<8; i++)
        {
            char a = ((i&4) != 0) ? '1' :'0';
            char b = ((i&2) != 0) ? '1' :'0';
            char c = ((i&1) != 0) ? '1' :'0';

            //Network speciesNetworkTopology = Networks.readNetwork("(A:.5,(B:.2,C:.2)n1:.3)n0;");
            //Network speciesNetworkTopology = Networks.readNetwork("((A:0.04,X#H1:0.02::0.7)n3:0.02,((B:0.02)X#H1:0.02::0.3,C:0.04)n1:0.02)root;");
            //Network speciesNetworkTopology = Networks.readNetwork("((B:0.5)X#H1:1.5::0.5,((X#H1:0.5::0.5,A:1)n1:0.5,C:1.5)n2:0.5)root;");
            //Network speciesNetworkTopology = Networks.readNetwork("(((C:0.09,(B:0.07)I6#H2:0.01::0.2)I3:0.02,(I6#H2:0.02::0.8)I5#H3:0.01::0.2)I2:0.07,(I5#H3:0.06::0.8,A:0.16)I1:0.02)I0;");
            //Network speciesNetworkTopology = Networks.readNetwork("(((C:0.09)I3:0.02,((B:0.07)I6:0.01)I5#H3:0.01::0.2)I2:0.07,(I5#H3:0.06::0.8,A:0.16)I1:0.02)I0;");
            Network speciesNetworkTopology = Networks.readNetwork("((A:1.5,((B:0.5,C:0.5)I3:0.5)I2#H1:0.5::0.5)I1:0.5,I2#H1:1.0::0.5)I0;");
            //Network speciesNetworkTopology = Networks.readNetwork("((A:1.5,C:1.5):0.5,B:2.0);");

            double constTheta = 0.036;
            for(Object nodeObject : Networks.postTraversal(speciesNetworkTopology)) {
                NetNode node = (NetNode) nodeObject;
                for(Object parentObject :  node.getParents()) {
                    NetNode parent = (NetNode) parentObject;
                    node.setParentSupport(parent, constTheta);
                    //node.setParentDistance(parent, node.getParentDistance(parent) * constTheta / 2.0);
                }
            }
            speciesNetworkTopology.getRoot().setRootPopSize(constTheta);


            Map<String, Character> colorMap = new HashMap<String, Character>();
            colorMap.put("A", a);
            colorMap.put("B", b);
            colorMap.put("C", c);

            OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);
            BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {0.5, 0.5}, new double[] {1.0});

            //System.out.println(BAGTRModel.getRateMatrix());


            SNAPPAlgorithm run = new SNAPPAlgorithm(speciesNetworkTopology, BAGTRModel, 1.0);
            System.out.println("Observation:  " + a + "," + b + "," + c);
            System.out.println("Probability:  " + run.getProbability(converter) + "\n");
            sum += run.getProbability(converter);
        }
        System.out.println("Sum: " + sum);
    }

    private static void testBiallelicLarge() {
        double sum =0;


        //Network speciesNetworkTopology = Networks.readNetwork("(A:.5,(B:.2,C:.2)n1:.3)n0;");
        //Network speciesNetworkTopology = Networks.readNetwork("((A:0.04,X#H1:0.02::0.7)n3:0.02,((B:0.02)X#H1:0.02::0.3,C:0.04)n1:0.02)root;");
        Network speciesNetworkTopology = Networks.readNetwork("((((((((G:0.22950212771454837)#H3:0.08670753342120574::0.00282793363498135,A:0.3162096611357541):0.017035797590358892,(Q:0.19012318627246375)#H1:0.14312227245364925::0.5228892817568864):0.010328596967308479,C:0.3435740556934215):0.0012203081657338188,(R:0.1195921448597946)#H2:0.22520221899936071::0.026849304992340173):0.00859180946928989,L:0.3533861733284452):1.008223109245107,(#H1:0.9660558616585218::0.4771107182431136,#H2:1.036586903071191::0.9731506950076598):0.2054302346425667):0.06680556134344817,#H3:1.198912716202452::0.9971720663650187);");
        //Network speciesNetworkTopology = Networks.readNetwork("((A:1.5,C:1.5):0.5,B:2.0);");

        speciesNetworkTopology = Networks.readNetwork("(((((A:0.7)I6#H1:1.3::0.8,Q:2.0)I4:1.0,L:3.0)I3:1.0,R:4.0)I2:1.0,(G:2.0,(I6#H1:0.7::0.2,C:1.4)I5:0.6)I1:3.0)I0;");
        speciesNetworkTopology = Networks.readNetwork("(((((Q:2.0,A:2.0)I4:1.0,L:3.0)I3:0.5)I8#H1:0.5::0.7,R:4.0)I2:1.0,(I8#H1:0.5::0.3,(G:2.0,C:2.0)I1:2.0)I7:1.0)I0;");
        //speciesNetworkTopology = Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,A:1.0)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,(G:1.0,C:1.0)I1:4.0)I0;");
        //speciesNetworkTopology = Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,(A:0.5)I6#H2:0.5::0.8)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,((I6#H2:0.5::0.2,C:1.0)I5:1.0,G:2.0)I1:3.0)I0;");
        speciesNetworkTopology = Networks.readNetwork("(((((R:0.10472155953678641,(L:0.08956088683777966,(Q:0.061388364670773195,A:0.061388364670773195)II8:0.028172522167006463)II7:0.015160672699006755)II6:0.015145851252580833,((C:0.07315556526977926,G:0.07315556526977926)II3:0.03223357104021027)II2#H1:0.014478274479377717::0.8669225099912875)II5:2.152380577149497)II4#H2:0.1822392555742689::0.22368190104831598,II2#H1:2.3490981072031434::0.1330774900087125)II1:6.446533712494073,II4#H2:6.6287729680683425::0.776318098951684)II0;");
        //speciesNetworkTopology = Networks.readNetwork("((((L:0.009304825035069735)#H1:14.370818397517256::0.020382388550110098)#H2:25.83351345488807::0.28805683793703096,((#H1:0.0013174924892173801::0.9796176114498899,(Q:0.009166552336431677,(A:0.005296210274212645,(G:0.0014429534014067207,C:0.0014429534014067207):0.003853256872805924):0.0038703420622190326):0.0014557651878554373):0.0047109417058845745,R:0.015333259230171689):40.19830341821022):75.28824647371039,#H2:101.12175992859846::0.711943162062969);");
        double constTheta = 0.036;
        for(Object nodeObject : Networks.postTraversal(speciesNetworkTopology)) {
            NetNode node = (NetNode) nodeObject;
            for(Object parentObject :  node.getParents()) {
                NetNode parent = (NetNode) parentObject;
                node.setParentSupport(parent, constTheta);
                //node.setParentDistance(parent, node.getParentDistance(parent) * constTheta / 2.0);
            }
        }
        speciesNetworkTopology.getRoot().setRootPopSize(constTheta);

        int nameCount = 0;
        for(Object node : speciesNetworkTopology.dfs()) {
            NetNode mynode = (NetNode) node;
            if(mynode.getName().equals("")) {
                mynode.setName("I" + nameCount);
                nameCount++;
            }
        }

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {0.5, 0.5}, new double[] {1.0});

        //System.out.println(BAGTRModel.getRateMatrix());


        SNAPPAlgorithm run = new SNAPPAlgorithm(speciesNetworkTopology, BAGTRModel, 1.0);
        long last = System.currentTimeMillis();

        for(int repeat = 0 ; repeat < 20 ; repeat++) {
            sum = 0;
            for (int i = 0; i < 64; i++) {

                char a = ((i & 1) != 0) ? '1' : '0';
                char b = ((i & 2) != 0) ? '1' : '0';
                char c = ((i & 4) != 0) ? '1' : '0';
                char d = ((i & 8) != 0) ? '1' : '0';
                char e = ((i & 16) != 0) ? '1' : '0';
                char f = ((i & 32) != 0) ? '1' : '0';


                Map<String, Character> colorMap = new HashMap<String, Character>();
                colorMap.put("A", a);
                colorMap.put("C", b);
                colorMap.put("G", c);
                colorMap.put("L", d);
                colorMap.put("Q", e);
                colorMap.put("R", f);


                OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);
                //System.out.println("Observation:  " + a + "," + b + "," + c+","+d+","+e+"," +f);
                //System.out.println("Probability:  " + run.getProbability(converter) + "\n");
                sum += run.getProbability(converter);
                //break;
            }
            System.out.println("Sum: " + sum);
        }


        System.out.println(String.format("Total elapsed time : %2.5f s\n",
                (double) (System.currentTimeMillis() - last) / 1000.0));
    }

    private static void checkWithOriginalPaperData() {
        String filename = "/Users/zhujiafan/Downloads/PaperSimulations/Simulation1Clean/tree_easy4_100_aa_seed200/tree_easy4_100_aa.xml";

        long startTime = System.currentTimeMillis();

        Map<String, String> sequence = new HashMap<>();
        try {
            File fXmlFile = new File(filename);
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
            cache.put(0, 0);
            cache.put((1 << 4) - 1, 0);

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

        pi[1] = 0.5;
        pi[0] = 1 - pi[1];

        rate[0] = 1 / (2*pi[0]);

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(pi, rate);
        Network trueNetwork;
        trueNetwork = Networks.readNetwork("(((A:0.0057:0.006,B:0.0057:0.006)I1:0.0045:0.005,C:0.0102:0.006)I2:0.0138:0.005,D:0.024:0.006)I3;");
        trueNetwork.getRoot().setRootPopSize(0.006);

        trueNetwork = Networks.readNetwork("(((A:0.013922554907530933:0.0016612632913445985,B:0.01392255490753094:0.004606359528304711)I1:0.0102260158920536:0.006014080275576303,C:0.02414857079958455:0.0029854913730910735)I2:0.015961937881117386:0.007008380469480965,D:0.04011050868070192:0.007741741979758894)I0;");
        trueNetwork.getRoot().setRootPopSize(0.009213561041854954);

        trueNetwork = Networks.readNetwork("(((A:0.003408359879116683:4.4775539481432705E-4,B:0.0034083598791166897:0.006364282857113693)I1:0.005552517402502629:0.0036728018937483005,C:0.00896087728161933:0.0018937185668051557)I2:0.005843886090002517:0.005523809655807946,D:0.014804763371621833:0.003751804012158342)I0;");
        trueNetwork.getRoot().setRootPopSize(0.006176353468824599);


        for(Alignment alg : alns) {
            List<String> names = alg.getTaxaNames();
            double sum = 0;

            double P0 = 0;
            double P1 = 0;
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
                    if(represent == ((1 << 4) - 1)) P1 = likelihood;
                    sum += Math.log(likelihood) * count;
                    System.out.println(represent + " " + likelihood + " " + count);
                } catch(Exception e) {
                    e.printStackTrace();
                    System.out.println("Exceptional network");
                }


            }
            sum -= 100 * Math.log(1.0 - P0 - P1);
            System.out.println("sum " + sum);


        }

        System.out.println(String.format("Total elapsed time : %2.5f s\n",
                (double) (System.currentTimeMillis() - startTime) / 1000.0));
    }


    private static void testSpeed(){
        long startTime = System.currentTimeMillis();
        GTRModel gtrModel = new GTRModel(new double[]{.1, .2, .3, .4},new double[]{1,1.5,2,2.5,3,3.5});
        Network speciesTreeTopology = Networks.readNetwork("((A:0.006,B:0.006)i4:0.006,(C:0.006,(D:0.006,((E:0.006,G:0.006)i0:0.5,F:0.006)i1:0.006)i2:0.006)i3:0.006)i5;");
        Map<String, Character> colorMap = new HashMap<String, Character>();
        colorMap.put("A", 'A');
        colorMap.put("B", 'A');
        colorMap.put("C", 'A');
        colorMap.put("D", 'A');
        colorMap.put("E", 'A');
        colorMap.put("F", 'A');
        colorMap.put("G", 'A');
        for(int i=0; i<87; i++) {
            SNAPPAlgorithm run = new SNAPPAlgorithm(speciesTreeTopology, gtrModel, 1.0);
            OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);
            run.getProbability(converter);
        }
        System.out.println(System.currentTimeMillis()-startTime);
    }

    private static Tree toTree(String exp){
        Tree tree = null;
        try{
            NewickReader nr = new NewickReader(new StringReader(exp));
            tree = nr.readTree();
        }catch (Exception e){

        }
        return tree;
    }


}
