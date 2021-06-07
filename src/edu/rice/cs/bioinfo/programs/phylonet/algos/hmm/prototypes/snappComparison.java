package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.prototypes;

import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm.IntegrationAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.JCModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.ExNewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.HmmTreeUtils;

import java.io.IOException;
import java.io.StringReader;
import java.util.*;

public class snappComparison {

    public long seed;
    public  Tree topology;
    private static Map<String,String> baDATA = new HashMap<String,String>();


    public snappComparison(Tree specTreeTopology, long numSeed){
        this.seed = numSeed;
        this.topology = specTreeTopology;
    }


    public static void main(String[] args){
        baDATA.put("A", "AATCGGGT");
        baDATA.put("B", "ACTAATGT");
        baDATA.put("C", "ATGCCTGT");

        double allSum = 0;
        for(int i =0; i<8; i++) {
            NucleotideObservation obs = NucleotideObservation.getObservations(baDATA).get(i);
            System.out.println("\nOutput for LOCI #" + i);
            System.out.println("THE OBSERVATION for A,B,C = " + obs.getObservationForAllele("A") + obs.getObservationForAllele("B") + obs.getObservationForAllele("C"));

            List<String> geneTreesStrings = Arrays.asList("(A,(B,C)n1)n0", "(B,(A,C)n1)n0", "(C,(B,A)n1)n0");

            List<Tree> geneTrees = new ArrayList<Tree>();
            for (String str: geneTreesStrings)
                geneTrees.add(HmmTreeUtils.toTree(str));

            Tree speciesTreeTopology = toTree("(A:.5,(B:.2,C:.2)n1:.3)n0");

            double sum = 0;
            for (Tree geneTree : geneTrees) {
                snappComparison intAlg = new snappComparison(speciesTreeTopology, 34563456345123490L);
                double result = intAlg.getProbability(obs, geneTree);
 //               System.out.println(result);
 //               System.out.println(geneTree);
                sum+=result;
            }
            System.out.println(sum);
            allSum += sum;
        }
        System.out.println(allSum);

    }

    //GTRModel gtrModel = new GTRModel(new double[]{.1, .2, .3, .4},new double[]{1,1.5,2,2.5,3,3.5});
    JCModel gtrModel = new JCModel(.2);

    private double getProbability(NucleotideObservation dna, Tree geneTreeTopology) {
        return new IntegrationAlgorithm(topology,geneTreeTopology,gtrModel,2).getProbability(dna);
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






}
