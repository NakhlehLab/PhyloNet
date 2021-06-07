package edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.integration.FullLikelihoodFromSequence;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.SubstitutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.ExNewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.io.IOException;
import java.io.StringReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class IntegrationAlgorithm extends NucleotideProbabilityAlgorithm {

    Network sourceTree;
    Tree geneTreeTopology;
    SubstitutionModel model;

    double theta;

    public IntegrationAlgorithm(Tree sourceTree,Tree geneTreeTopology, SubstitutionModel model, double theta)
    {
        this.sourceTree = toNetwork(sourceTree);
        this.geneTreeTopology = geneTreeTopology;
        this.model = model;
        this.theta = theta;
    }

    public IntegrationAlgorithm(Tree sourceTree,Tree geneTreeTopology, SubstitutionModel model)
    {
        this(sourceTree,geneTreeTopology,model,1);
    }

    /**
     * Creates a network from a tree.
     *
     * @param speciesTree The tree to convert.
     * @return The tree in network form. Note the parametrized type (so it can work with the GeneTreeProbability class).
     */
    private static <T> Network<T> toNetwork(Tree speciesTree) {
        String newickString = speciesTree.toNewickWD();

        ExNewickReader<T> reader = new ExNewickReader<T>(new StringReader(newickString));

        try {
            return reader.readNetwork();
        } catch (IOException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }




    /**
     * Run the Felstein algorithm and obtain the likelihood of the whole tree.
     * @param dna A mapping of leaf node names from the gene tree to 'A','C','T', or 'G'.
     * @return The likelihood.
     */
    @Override
    public double getProbability(NucleotideObservation dna) {

        String[] leaves = geneTreeTopology.getLeaves();
        Map<String,List<String>> mapp = new HashMap<String,List<String>>();


        char[] obs = new char[leaves.length];

        for (int i = 0; i < leaves.length; i++)
        {
            obs[i] = dna.getObservationForAllele(leaves[i]);
            mapp.put(leaves[i],Arrays.asList(leaves[i]));
        }

        //FullLikelihoodFromSequence l = new FullLikelihoodFromSequence(Arrays.asList(this.geneTreeTopology),geneTreeTopology.getLeaves(),Arrays.asList(new Tuple<>(obs,1)),model,theta, 0.1, 20);

        //l.computeGTLikelhoodWithIntegral(100000,2);
        //return Math.exp(l.computeNetworkLikelhoodWithIntegral(sourceTree,mapp));
        return 0;
    }

    private static Tree toTree(String s)
    {
        NewickReader a = new NewickReader(new StringReader(s));
        try {
            return a.readTree();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    /*public static void main(String[] args)
    {
        Map<String,String> dna = new HashMap<>();
        dna.put("1","AAA");
        dna.put("2","ACA");
        dna.put("3","CAA");

        Tree topology = toTree("((1:1.0001438964487415,2:2.999144327357899)n1:0.9998732541373248,3:2.0000160579853397)n0");
        List<String> geneTreesStrings = Arrays.asList("(1:5.364,(2:1.176,3:1.176):4.188)","(1:6.138,(2:1.176,3:1.176):4.962)","(3:5.985,(1:3.891,2:3.891):2.094)");
        List<Tree> geneTrees = geneTreesStrings.stream().map(Trees::toTree).collect(Collectors.toList());

        NucleotideObservation obs = NucleotideObservation.getObservations(dna).get(0);

        long seed = 1010123012;
        for (Tree geneTree : geneTrees)
        {
            GeneTreeProbabilityYF f = new GeneTreeProbabilityYF();
            double[] result = new double[1];

            f.calculateGTDistribution(toNetwork(topology),Arrays.asList(geneTree),null,result);


            System.out.println(geneTree);
            double foo = new FelsensteinAlgorithm(geneTree,new JCModel(1)).getProbability(obs);
            //System.out.println(foo);
            //System.out.println(result[0]);
            //System.out.println(foo *result[0]);

            IntegrationAlgorithm algo = new IntegrationAlgorithm(topology,geneTree,new JCModel(1),seed);
            System.out.println(algo.getProbability(obs));
            //System.out.println(algo.getProbability(obs));
            //System.out.println(algo.getProbability(obs));

        }

    }*/
}
