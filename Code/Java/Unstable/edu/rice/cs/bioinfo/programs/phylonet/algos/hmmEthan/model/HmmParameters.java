package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model;


import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.Configuration;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.ConfigurationBuilder;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.HmmNetworkUtils;

import java.util.*;

/**
 * This is a class which facilitates converting the hidden markov parameters to a double[] and back again.
 */
public class HmmParameters {


    int numberOfTaxa;
    Network<String> net;
    private int numberOfSpeciesTrees;
    private int numberOfGeneTrees;
    Configuration config;

    /**
     * Creates a new HmmParameters for a given number of taxa.
     * The number of taxa is important for determining the number of gene tree lengths.
     *
     * @param numberOfTaxa The number of taxa.
     * @param config
     */
    public HmmParameters(int numberOfTaxa,Network<String> net,Map<String,List<String>> speciesOptions, int numGeneTrees, Configuration config) {
        this.numberOfTaxa = numberOfTaxa;
        this.net = net;
        this.config = config;
        numberOfSpeciesTrees = HmmNetworkUtils.generateParentalTrees(net, speciesOptions).size();
        numberOfGeneTrees = numGeneTrees;
    }

    public int getNumberOfSpeciesTrees()
    {
        return numberOfSpeciesTrees;
    }

    /**
     * Simply scales the equilibrium frequencies so that they sum to 1.
     *
     * @param preScaledEquilibriumFrequencies The previous unscaled frequencies.
     * @return A new array with the new scaled frequencies.
     */
    double[] scaleEquilibriumFrequencies(double[] preScaledEquilibriumFrequencies) {

        double[] temp = preScaledEquilibriumFrequencies.clone();
        double sum = 0;
        for (double v : preScaledEquilibriumFrequencies) {
            sum += v;
        }

        for (int i = 0; i < preScaledEquilibriumFrequencies.length; i++) {
            temp[i] /= sum;
        }
        return temp;
    }

    /**
     * Decodes an input array of doubles into an easier to use Data class.
     *
     * @param input The input array.
     * @return A Data object with the correct elements.
     */
    public Data decode(double[] input) {
        Data result = new Data();

        if (input.length != numberOfParameters())
            throw new RuntimeException("Wrong length of input");


        ArrayReader foo = new ArrayReader(input);
        result.scale = foo.getNext();
        result.equilibriumFrequencies = scaleEquilibriumFrequencies(foo.getNext(4));
        result.transitionFrequencies = scaleEquilibriumFrequencies(foo.getNext(6));
        result.speciesStayProbs = foo.getNext(numberOfSpeciesStayProbs());
        result.geneStayProb = foo.getNext();
        result.speciesNetworkBranchLengths = foo.getNext(numberOfNetworkLengths());

        if (config.ALGORITHM == Configuration.AlgorithmType.NORMAL)
            result.allGeneTreeLengths = foo.getNext(numberOfGeneTreeLengths());

        return result;
    }

    /**
     * Encodes a Data object into an array of doubles for use with apache commons.
     *
     * @param input A Data object with the correct information.
     * @return A double array with the information from the Data object.
     */
    public double[] encode(Data input) {
        double[] result = new double[numberOfParameters()];


        ArrayWriter foo = new ArrayWriter(result);
        foo.putNext(input.scale);
        foo.putNext(input.equilibriumFrequencies);
        foo.putNext(input.transitionFrequencies);
        foo.putNext(input.speciesStayProbs);
        foo.putNext(input.geneStayProb);
        foo.putNext(input.speciesNetworkBranchLengths);

        if (config.ALGORITHM == Configuration.AlgorithmType.NORMAL)
            foo.putNext(input.allGeneTreeLengths);

        return result;
    }

    /**
     * Calculates the number of gene tree lengths in this model.
     *
     * @return The number of gene tree lengths.
     */
    public int numberOfGeneTreeLengths() {
        return numberOfTreeLengthsPerGeneTree() * numberOfGeneTrees();
    }


    public int numberOfSpeciesStayProbs()
    {
        if (config.USEMULTIPLESTAYPROBABILITIES)
            return numberOfSpeciesTrees;
        else
            return 1;
    }

    /**
     * Calculates the number of parameters (IE, the number of elements in the double arrays).
     *
     * @return The number of parameters.
     */
    public int numberOfParameters() {


        switch(config.ALGORITHM)
        {
            case INTEGRATION:
            case SNAPP:
                return 12 + numberOfNetworkLengths() + numberOfSpeciesStayProbs();
            case NORMAL:
                return 12 + numberOfGeneTreeLengths() + numberOfSpeciesStayProbs() + numberOfNetworkLengths();
        }
        throw new RuntimeException("Invalid algorithm type");
    }

    /**
     * Calculates the number of tree lengths per each gene tree.
     *
     * @return The number of tree lengths per each gene tree.
     */
    public int numberOfTreeLengthsPerGeneTree() {
        return numberOfTaxa * 2 - 2;
    }

    /**
     * Calculates the number of gene trees for the given model.
     *
     * @return The number of gene trees.
     */
    public int numberOfGeneTrees() {
       return numberOfGeneTrees;
    }

    <T> int sizeOfIterable(Iterable<T> iter)
    {
        return ((Collection<T>) iter).size();
    }

    public int numberOfNetworkLengths()
    {
        return sizeOfIterable(net.getTreeNodes()) + sizeOfIterable(net.getNetworkNodes())- sizeOfIterable(net.getLeaves());
    }

    /**
     * Gets the state id for given gene and species indices.
     *
     * @param geneTreeIndex    The index of the gene tree.
     * @param speciesTreeIndex The index of the species tree.
     * @return The index as an integer.
     */
    public int getStateIndex(int geneTreeIndex, int speciesTreeIndex) {
        return geneTreeIndex * numberOfSpeciesTrees + speciesTreeIndex;
    }

    /**
     * A simple data class to store all of the parameters.
     */
    public static class Data {
        public double scale;
        public double[] equilibriumFrequencies;
        public double[] transitionFrequencies;
        public double[] speciesStayProbs;
        public double geneStayProb;
        public double[] allGeneTreeLengths;
        public double[] speciesNetworkBranchLengths;

        @Override
        public String toString() {
            return "Data{" +
                    "scale=" + scale +
                    "equilibriumFrequencies=" + Arrays.toString(equilibriumFrequencies) +
                    ", transitionFrequencies=" + Arrays.toString(transitionFrequencies) +
                    ", speciesStayProbs=" + Arrays.toString(speciesStayProbs) +
                    ", geneStayProb=" + geneStayProb +
                    ", allGeneTreeLengths=" + Arrays.toString(allGeneTreeLengths) +
                    ", speciesNetworkBranchLengths=" + Arrays.toString(speciesNetworkBranchLengths) +
                    '}';
        }
    }


    public static void main(String[] args)
    {
        Network<String> foo = HmmNetworkUtils.fromENewickString("(A:1,B:1);");
        Map<String,List<String>> alleleMapping = new HashMap<String,List<String>>();
        alleleMapping.put("A",Arrays.asList("A"));
        alleleMapping.put("B",Arrays.asList("B"));
        HmmParameters p = new HmmParameters(2,foo,alleleMapping,1, ConfigurationBuilder.getSNAPP().build(""));

        Random rand = new Random();
        double[] arr = new double[p.numberOfParameters()];
        for (int i = 0 ;i < arr.length;i++)
            arr[i] = rand.nextDouble();

        System.out.println(Arrays.toString(arr));

        Data fooData = p.decode(arr);
        double[] arrTwo = p.encode(fooData);
        System.out.println(Arrays.toString(arrTwo));

        double[] arrThree = p.encode(p.decode(arrTwo));
        System.out.println(Arrays.toString(arrThree));


    }

}

class ArrayWriter
{
    int index = 0;
    final double[] arr;

    public ArrayWriter(double[] arr)
    {
        this.arr = arr;
    }

    void putNext(double num)
    {
        arr[index++] = num;
    }

    void putNext(double[] nums)
    {
        System.arraycopy(nums,0,arr,index,nums.length);
        index+=nums.length;
    }


}

class ArrayReader
{
    int index = 0;
    final double[] arr;

    public ArrayReader(double[] arr)
    {
        this.arr = arr;
    }

    double getNext()
    {
        return arr[index++];
    }

    double[] getNext(int num)
    {
        double[] result =  Arrays.copyOfRange(arr,index,index+num);
        index+=num;
        return result;
    }
}