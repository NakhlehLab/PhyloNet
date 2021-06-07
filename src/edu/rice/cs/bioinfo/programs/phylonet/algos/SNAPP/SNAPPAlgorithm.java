package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPP;

import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.OneNucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.GTRModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm.NucleotideProbabilityAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.RateModel;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.HmmTreeUtils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.StringReader;
import java.util.HashMap;
import java.util.Map;

public class SNAPPAlgorithm extends NucleotideProbabilityAlgorithm {

    private MatrixQ Q;
    private Tree speciesTree;

    public SNAPPAlgorithm(Tree theSpeciesTree, RateModel rModel, double theta){
        speciesTree = theSpeciesTree;

        Q = new MatrixQ(rModel, speciesTree.getLeafCount(), theta);
    }


    public SNAPPAlgorithm(Tree theSpeciesTree, RateModel rModel){
        speciesTree = theSpeciesTree;

        Q = new MatrixQ(rModel, speciesTree.getLeafCount(), 0.02);
    }


    @Override
    public double getProbability(NucleotideObservation dna) {

        Tree result = HmmTreeUtils.flattenTree(speciesTree);

        double resultVal= Algorithms.getProbabilityObservationGivenTree(result, dna, Q);

        if (Double.isNaN(resultVal) || Double.isInfinite(resultVal))
            throw new RuntimeException("SNAPP result is not finite");
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
        testFourAlleles();
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
                    Tree speciesTreeTopology = toTree("(A:1,(B:.4,C:.4)n1:.6)n0");
                    Map<String, Character> colorMap = new HashMap<String, Character>();
                    colorMap.put("A", a);
                    colorMap.put("B", b);
                    colorMap.put("C", c);


                    OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);


                    GTRModel gtrModel = new GTRModel(new double[]{.1, .2, .3, .4},new double[]{1,1.5,2,2.5,3,3.5});

                    SNAPPAlgorithm run = new SNAPPAlgorithm(speciesTreeTopology, gtrModel, 2);
                    System.out.println("Observation:  " + a +',' + b + ','+ c);
                    System.out.println("Probability:  " + run.getProbability(converter) + "\n");
                    sum += run.getProbability(converter);

                }


        System.out.println("Sum: " + sum);




    }


    private static void testSpeed1(){
        long startTime = System.currentTimeMillis();
        GTRModel gtrModel = new GTRModel(new double[]{.1, .2, .3, .4},new double[]{1,1.5,2,2.5,3,3.5});
        //Tree speciesTreeTopology = toTree("((A:0.006,B:0.006)i4:0.006,(C:0.006,(D:0.006,((E:0.006,G:0.006)i0:0.5,F:0.006)i1:0.006)i2:0.006)i3:0.006)i5;");
        //Tree speciesTreeTopology = toTree("((A:0.006,B:0.006)i4:0.006,(C:0.006,(D:0.006,E:0.006):0.006):0.006)i5;");
        Tree speciesTreeTopology = toTree("(((1:5.436891196340337,((2:1.9514374480674375,3:1.9514374480674375):3.0851202667288735,(4:1.944547302080008,5:1.944547302080008):3.092010412716303):0.400333481544026):7.327154442501835,(6:8.25561963243531,(7:1.0546860644694407,8:1.0546860644694407):7.20093356796587):4.508426006406862):3.235954361157828,(((9:10.202434506154612,(10:9.917387872449376,(11:8.77311617198174,(12:7.74802492571023,(13:5.593016807909792,14:5.593016807909792):2.1550081178004383):1.0250912462715092):1.1442717004676373):0.28504663370523564):1.8592676286438592,((15:2.6943829262972905,16:2.6943829262972905):2.2917498878556506,17:4.986132814152941):7.07556932064553):1.6998791684653605,(18:11.889096083654344,(19:2.70030386259598,20:2.70030386259598):9.188792150606252):1.872485219609489):2.2384186967361677);");

        Map<String, Character> colorMap = new HashMap<String, Character>();
        for(int i=1; i<=20; i++){
            colorMap.put(""+i, 'G');
        }
        //colorMap.put("E", 'G');
        //colorMap.put("F", 'A');
        //colorMap.put("G", 'A');
        SNAPPAlgorithm run = new SNAPPAlgorithm(speciesTreeTopology, gtrModel, 1);
        OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);
        run.getProbability(converter);
        System.out.println((System.currentTimeMillis() - startTime) / 1000.0);
    }


    private static void testSpeed2(){
        long startTime = System.currentTimeMillis();
        GTRModel gtrModel = new GTRModel(new double[]{.1, .2, .3, .4},new double[]{1,1.5,2,2.5,3,3.5});
        //Tree speciesTreeTopology = toTree("((A:0.006,B:0.006)i4:0.006,(C:0.006,(D:0.006,((E:0.006,G:0.006)i0:0.5,F:0.006)i1:0.006)i2:0.006)i3:0.006)i5;");
        //Tree speciesTreeTopology = toTree("((A:0.006,B:0.006)i4:0.006,(C:0.006,(D:0.006,E:0.006):0.006):0.006)i5;");
        Tree speciesTreeTopology = toTree("((A:0.006,B:0.006)i4:0.006,(C:0.006,D:0.006):0.006)i5;");

        Map<String, Character> colorMap = new HashMap<String, Character>();
        colorMap.put("A", 'A');
        colorMap.put("B", 'A');
        colorMap.put("C", 'A');
        colorMap.put("D", 'T');
        //colorMap.put("E", 'G');
        //colorMap.put("F", 'A');
        //colorMap.put("G", 'A');
        for(int j=0; j<100; j++) {
            startTime = System.currentTimeMillis();
            for (int i = 0; i < 80; i++) {
                SNAPPAlgorithm run = new SNAPPAlgorithm(speciesTreeTopology, gtrModel, 1);
                OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);
                run.getProbability(converter);
            }
            System.out.println((System.currentTimeMillis() - startTime) / 1000.0);

        }

    }

    private static void testFourAllelesMedium() {
        double sum =0;
        char[] nucleotides = {'A','C','T','G'};
        for (char a : nucleotides)
            for (char b: nucleotides)
                for (char c: nucleotides)
                    for(char d: nucleotides)
                        for(char e: nucleotides)
                        {
                            Tree speciesTreeTopology = toTree("((A:.2, B:.2)n1:.3,((C:.1,D:.1)n0:.3, E:.4)n2:.1)n3");
                            Map<String, Character> colorMap = new HashMap<String, Character>();
                            colorMap.put("A", a);
                            colorMap.put("B", b);
                            colorMap.put("C", c);
                            colorMap.put("D", d);
                            colorMap.put("E", e);


                            OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);


                            GTRModel gtrModel = new GTRModel(new double[]{.1, .2, .3, .4},new double[]{1,1.5,2,2.5,3,3.5});



                            SNAPPAlgorithm run = new SNAPPAlgorithm(speciesTreeTopology, gtrModel, 1);
                            System.out.println("Observation:  " + a +',' + b + ','+ c +','+ d + ','+ e);
                            System.out.println("Probability:  " + run.getProbability(converter) + "\n");
                            sum += run.getProbability(converter);
                        }



        System.out.println("Sum: " + sum);




    }

    private static void testFourAllelesLarge() {
        int i = 0;
        double sum =0;
        Tree speciesTreeTopology = toTree("(((A:.2, B:.2)n1:.3,((C:.1,D:.1)n0:.3, E:.4)n2:.1)n3:.1,(F:.2, G:.2)n4:.4)");
        GTRModel gtrModel = new GTRModel(new double[]{.1, .2, .3, .4},new double[]{1,1.5,2,2.5,3,3.5});

        SNAPPAlgorithm run = new SNAPPAlgorithm(speciesTreeTopology, gtrModel, 1);

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

            Tree speciesTreeTopology = toTree("(A:.5,(B:.2,C:.2)n1:.3)n0");
            Map<String, Character> colorMap = new HashMap<String, Character>();
            colorMap.put("A", a);
            colorMap.put("B", b);
            colorMap.put("C", c);

            OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);
            BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {5/6.0, 1/6.0}, new double[] {0.1});

            //System.out.println(BAGTRModel.getRateMatrix());


            SNAPPAlgorithm run = new SNAPPAlgorithm(speciesTreeTopology, BAGTRModel, 1);
            System.out.println("Observation:  " + a + "," + b + "," + c);
            System.out.println("Probability:  " + run.getProbability(converter) + "\n");
            sum += run.getProbability(converter);
        }
        System.out.println("Sum: " + sum);
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
