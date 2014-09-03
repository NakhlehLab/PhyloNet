package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.prototypes;

import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.SNAPPAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.JCModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.RateModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.OneNucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.HmmNetworkUtils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SimpleTest
{

    public static void main(String[] args)
    {

        bestTest();
        simpleTestOne();
        simpleTestTwo();
        simpleTestThree();
        simpleTestFour();

    }

    private static void simpleTestOne()
    {
        Network<String> net = HmmNetworkUtils.fromENewickString("((A:1,B:1):1,C:1);");
        Map<String,String> alleleToSpecies = new HashMap<String,String>();
        alleleToSpecies.put("A","A");
        alleleToSpecies.put("B","B");
        alleleToSpecies.put("C","C");


        RateModel model = new JCModel(.2);
        SNAPPAlgorithm algo2 = new SNAPPAlgorithm(net,model,0.2);


        double sum = 0;

        char[] alleles = {'A','C','T','G'};

        for (char a : alleles)
            for (char b: alleles)
                for (char c: alleles)
                {
                    Map<String,Character> obsMap = new HashMap<String,Character>();
                    obsMap.put("A",a);
                    obsMap.put("B",b);
                    obsMap.put("C",c);

                    sum += algo2.getProbability(new OneNucleotideObservation(obsMap));
                }

        System.out.println(sum);
    }

    private static void simpleTestTwo()
    {
        Network<String> net = Networks.readNetwork("((A:3,ANC#H1:2.25::.2):2,((B:.75)ANC#H1:0::.8,C:.75):4.25);");
        Map<String,String> alleleToSpecies = new HashMap<String,String>();
        alleleToSpecies.put("A","A");
        alleleToSpecies.put("B","B");
        alleleToSpecies.put("C","C");


        RateModel model = new JCModel(.2);
        SNAPPAlgorithm algo2 = new SNAPPAlgorithm(net,alleleToSpecies,model);


        double sum = 0;

        char[] alleles = {'A','C','T','G'};

        for (char a : alleles)
            for (char b: alleles)
                for (char c: alleles)
                {
                    Map<String,Character> obsMap = new HashMap<String,Character>();
                    obsMap.put("A",a);
                    obsMap.put("B",b);
                    obsMap.put("C",c);

                    sum += algo2.getProbability(new OneNucleotideObservation(obsMap));
                }

        System.out.println(sum);
    }

    private static void simpleTestThree()
    {
        Network<String> net = Networks.readNetwork("((A:3,ANC#H1:2.25::0.2):2,((B:.75)ANC#H1:0::0.8,C:.75):4.25);");
        Map<String,String> alleleToSpecies = new HashMap<String,String>();
        alleleToSpecies.put("A","A");
        alleleToSpecies.put("B1","B");
        alleleToSpecies.put("B2","B");
        alleleToSpecies.put("C","C");


        RateModel model = new JCModel(.2);
        SNAPPAlgorithm algo2 = new SNAPPAlgorithm(net,alleleToSpecies,model);


        double sum = 0;

        char[] alleles = {'A','C','T','G'};

        for (char a : alleles)
            for (char b1: alleles)
                for (char b2: alleles)
                for (char c: alleles)
                {
                    Map<String,Character> obsMap = new HashMap<String,Character>();
                    obsMap.put("A",a);
                    obsMap.put("B1",b1);
                    obsMap.put("B2",b2);
                    obsMap.put("C",c);

                    sum += algo2.getProbability(new OneNucleotideObservation(obsMap));
                }

        System.out.println(sum);
    }

    private static void simpleTestFour()
    {
        Network<String> net = Networks.readNetwork("((A:3,ANC#H1:2.25::0.2):2,((B:.75)ANC#H1:0::0.8,C:.75):4.25);");
        Map<String,List<String>> alleleToSpecies = new HashMap<String,List<String>>();
        alleleToSpecies.put("A",Arrays.asList("A"));
        alleleToSpecies.put("B",Arrays.asList("B1","B2"));
        alleleToSpecies.put("C",Arrays.asList("C"));


        List<STITree<String>> trees = HmmNetworkUtils.generateParentalTrees(net, alleleToSpecies);

        RateModel model = new JCModel(.2);
        edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPP.SNAPPAlgorithm algo2 = new edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPP.SNAPPAlgorithm(trees.get(0),model);


        double sum = 0;

        char[] alleles = {'A','C','T','G'};

        for (char a : alleles)
            for (char b1: alleles)
                for (char b2: alleles)
                    for (char c: alleles)
                    {
                        Map<String,Character> obsMap = new HashMap<String,Character>();
                        obsMap.put("A",a);
                        obsMap.put("B1",b1);
                        obsMap.put("B2",b2);
                        obsMap.put("C",c);

                        sum += algo2.getProbability(new OneNucleotideObservation(obsMap));
                    }

        System.out.println(sum);
    }

    private static void bestTest()
    {
        Network<String> net = Networks.readNetwork("((A:0,ANC#H1:0::0.5):500,((B:0)ANC#H1:0::0.5,C:0):500);");
        Map<String,String> alleleToSpecies = new HashMap<String,String>();
        alleleToSpecies.put("A","A");
        alleleToSpecies.put("B1","B");
        alleleToSpecies.put("B2","B");
        alleleToSpecies.put("C","C");


        RateModel model = new JCModel(.0001);
        SNAPPAlgorithm algo2 = new SNAPPAlgorithm(net,alleleToSpecies,model);



        char[] alleles = {'A','C'};

        char a='A';
        char c='C';
            for (char b1: alleles)
                for (char b2: alleles)
                    {
                        Map<String,Character> obsMap = new HashMap<String,Character>();
                        obsMap.put("A",a);
                        obsMap.put("B1",b1);
                        obsMap.put("B2",b2);
                        obsMap.put("C",c);

                        System.out.println(obsMap+ "," + algo2.getProbability(new OneNucleotideObservation(obsMap)));
                    }

    }







}
