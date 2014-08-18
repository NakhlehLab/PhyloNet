package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;

import edu.rice.cs.bioinfo.library.programming.Tuple;
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
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.StringReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

public class SNAPPAlgorithm extends NucleotideProbabilityAlgorithm {

    private MatrixQ Q;
    private Network speciesNetwork;
    private Map<String, String> allele2species;

    public SNAPPAlgorithm(Network theSpeciesNetwork, RateModel rModel, double theta){

        speciesNetwork = theSpeciesNetwork;
        Q = new MatrixQ(rModel, speciesNetwork.getLeafCount(), theta);
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

    public SNAPPAlgorithm(Network theSpeciesNetwork, Map<String, String> alleleMapping, RateModel rModel, double theta){
        speciesNetwork = theSpeciesNetwork;
        if(speciesNetwork.getRoot().getData()==null) {
            Networks.removeBinaryNodes(speciesNetwork);
            for (Object node : Networks.postTraversal(speciesNetwork)) {
                ((NetNode<SNAPPData[]>) node).setData(new SNAPPData[1]);
            }
        }
        allele2species = alleleMapping;
        if(allele2species == null) {
            Q = new MatrixQ(rModel, speciesNetwork.getLeafCount(), theta);
        }else{
            Q = new MatrixQ(rModel, alleleMapping.size(), theta);
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



    public double getProbability(NucleotideObservation dna, int siteID) {
        double resultVal= Algorithms.getProbabilityObservationGivenNetwork(speciesNetwork, allele2species, dna, Q, siteID);

        if (Double.isNaN(resultVal)||Double.isInfinite(resultVal))
            throw new RuntimeException("Result is not finite");
        return resultVal;
    }



    public double updateProbability(int siteID, Map<NetNode,Boolean> node2update) {

        double resultVal= Algorithms.updateProbabilityObservationGivenNetwork(speciesNetwork, Q, siteID, node2update);

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

        R.dims = 3;
*/
        //testSpeed();
        //testFourAlleles();
        testFourAllelesTwoUnderReticulation();
/*


        R.dims = 3;

        testFourAllelesMedium();


        R.dims = 3;
        testFourAllelesLarge();
*/
    }

    private static void testFourAlleles() {
        double sum =0;
        char[] nucleotides = {'A','C','T','G'};
        for (char a : nucleotides)
            for (char b: nucleotides)
                for (char c: nucleotides)
                {
                    Network speciesNetworkTopology = Networks.readNetwork("((A:0.04,X#H1:0.02::0.7)n3:0.02,((B:0.02)X#H1:0.02::0.3,C:0.04)n1:0.02)root;");

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



                        SNAPPAlgorithm run = new SNAPPAlgorithm(speciesNetworkTopology, allele2species, gtrModel, 1);
                        System.out.println("Observation:  " + a +',' + b1 + ','+ b2 + "," +  b3 + "," + c);
                        //System.out.println("Probability:  " + run.getProbability(converter) + "\n");
                        sum += run.getProbability(converter);

                    }


        System.out.println("Sum: " + sum);




    }



    private static void testFourAllelesTwoUnderReticulation() {
        double sum =0;
        char[] nucleotides = {'A','C','T','G'};
        for (char a : nucleotides)
            for (char b: nucleotides)
                for (char c: nucleotides)
                    for (char d: nucleotides)
                {
                    //Network speciesNetworkTopology = Networks.string2network("(A:.5,(B:.2,C:.2)n1:.3)n0;");
                    //Network speciesNetworkTopology = Networks.string2network("((A:0.3,X#H1:0::0.5)n2:0.1,(((B:0.1,C:0.1)n3:0.1)X#H1:0::0.5,D:0.3)n1:0.1)root;");
                    Network speciesNetworkTopology = Networks.readNetwork("((A:0.02,X#H1:0.02::0.3)n2:0.02,(((B:0.02,C:0.02)n3:0.02)X#H1:0.02::0.7,D:0.02)n1:0.02)root;");
                    //Network speciesNetworkTopology = Networks.string2network("((A:1,X#H1:1::0.3)n2:1,(((B:1,C:1)n3:1)X#H1:1::0.7,D:1)n1:1)root;");

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


                    //GTRModel gtrModel = new GTRModel(new double[]{.1, .2, .3, .4},new double[]{1,1.5,2,2.5,3,3.5});
                    GTRModel gtrModel = new GTRModel(new double[]{1, 0, 0, 0},new double[]{0,0,0,0,0,0});



                    SNAPPAlgorithm run = new SNAPPAlgorithm(speciesNetworkTopology, gtrModel, 0.02);
                    System.out.println("Observation:  " + a +',' + b + ','+ c + "," + d);
                    System.out.println("Probability:  " + run.getProbability(converter) + "\n");
                    sum += run.getProbability(converter);

                }


        System.out.println("Sum: " + sum);




    }


    private static void testFourAllelesMedium() {
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



                            SNAPPAlgorithm run = new SNAPPAlgorithm(speciesNetworkTopology, gtrModel, 1);
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

        SNAPPAlgorithm run = new SNAPPAlgorithm(speciesNetworkTopology, gtrModel, 1);

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

            Network speciesNetworkTopology = Networks.readNetwork("(A:.5,(B:.2,C:.2)n1:.3)n0");
            Map<String, Character> colorMap = new HashMap<String, Character>();
            colorMap.put("A", a);
            colorMap.put("B", b);
            colorMap.put("C", c);

            OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);
            BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {5/6.0, 1/6.0}, new double[] {0.1});

            //System.out.println(BAGTRModel.getRateMatrix());


            SNAPPAlgorithm run = new SNAPPAlgorithm(speciesNetworkTopology, BAGTRModel, 1);
            System.out.println("Observation:  " + a + "," + b + "," + c);
            System.out.println("Probability:  " + run.getProbability(converter) + "\n");
            sum += run.getProbability(converter);
        }
        System.out.println("Sum: " + sum);
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
            SNAPPAlgorithm run = new SNAPPAlgorithm(speciesTreeTopology, gtrModel, 1);
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
