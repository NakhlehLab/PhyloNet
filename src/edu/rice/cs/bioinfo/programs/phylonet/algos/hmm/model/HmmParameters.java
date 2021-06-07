package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.model;


import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.Configuration;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.ConfigurationBuilder;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.GTRModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.JCModel;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.HmmNetworkUtils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * This is a class which facilitates converting the hidden markov parameters to a double[] and back again.
 */
public class HmmParameters {


    int numberOfTaxa;
    Network<String> net;
    Map<String, List<String>> species2alleles;
    private int numberOfAlleleMappings;
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
        numberOfAlleleMappings = HmmNetworkUtils.generateParentalTrees(net, speciesOptions).size();
        numberOfGeneTrees = numGeneTrees;
        species2alleles = speciesOptions;
    }

    public int getNumberOfAlleleMappings()
    {
        return numberOfAlleleMappings;
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
        switch (config.MODEL)
        {
            case GTR:
                result.equilibriumFrequencies = scaleEquilibriumFrequencies(foo.getNext(4));
                result.transitionFrequencies = scaleEquilibriumFrequencies(foo.getNext(6));
        }

        result.speciesNetworkParameters = foo.getNext(numberOfNetworkParameters());
        result.geneStayProb = foo.getNext();
        result.speciesStayProb = foo.getNext();
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

        if(input.getSize() != numberOfParameters()){
            throw new RuntimeException("Size not correct!");
        }

        ArrayWriter foo = new ArrayWriter(result);
        foo.putNext(input.scale);
        switch (config.MODEL)
        {
            case GTR:
                foo.putNext(input.equilibriumFrequencies);
                foo.putNext(input.transitionFrequencies);
        }
        foo.putNext(input.speciesNetworkParameters);
        foo.putNext(input.geneStayProb);
        foo.putNext(input.speciesStayProb);
        return result;
    }



    /**
     * Calculates the number of parameters (IE, the number of elements in the double arrays).
     *
     * @return The number of parameters.
     */
    public int numberOfParameters() {

        int numParametersForSubstitutionModel = 0;
        switch (config.MODEL)
        {
            case GTR:
                numParametersForSubstitutionModel = 11;
                break;
            case JC:
                numParametersForSubstitutionModel = 1;
        }

        switch(config.ALGORITHM)
        {
            case INTEGRATION:
            case SNAPP:
                return numParametersForSubstitutionModel + 1 + numberOfNetworkParameters() + (numberOfAlleleMappings>1?1:0);
            case NORMAL:
                return numParametersForSubstitutionModel + 1 + numberOfNetworkParameters() + (numberOfAlleleMappings>1?1:0);
        }
        throw new RuntimeException("Invalid algorithm type");
    }


    /**
     * Calculates the number of gene trees for the given model.
     *
     * @return The number of gene trees.
     */
    public int numberOfGeneTrees() {
        return numberOfGeneTrees;
    }

/*
    public int numberOfReticulations() {
        return sizeOfIterable(net.getNetworkNodes());
    }
*/

    <T> int sizeOfIterable(Iterable<T> iter)
    {
        return ((Collection<T>) iter).size();
    }

    public int numberOfNetworkParameters()
    {
        int numParameters = 0;
        Map<NetNode, Set<String>> node2leaves = new HashMap<>();
        for(NetNode node: Networks.postTraversal(net)){
            if(node.isRoot())break;
            Set<String> leaves = new HashSet<>();
            node2leaves.put(node, leaves);
            if(node.isLeaf()){
                leaves.addAll(species2alleles.get(node.getName()));
            }
            else{
                for(Object childO: node.getChildren()){
                    leaves.addAll(node2leaves.get(childO));
                }
            }
            /*
            if(node.isNetworkNode()){
                numParameters++;
            }
            */
            if(leaves.size()>1){
                numParameters += sizeOfIterable(node.getParents());
            }
        }
        return numParameters;
    }



    /**
     * Gets the state id for given gene and species indices.
     *
     * @param geneTreeIndex    The index of the gene tree.
     * @param speciesTreeIndex The index of the species tree.
     * @return The index as an integer.
     */
    public int getStateIndex(int geneTreeIndex, int speciesTreeIndex) {
        return geneTreeIndex * numberOfAlleleMappings + speciesTreeIndex;
    }

    /**
     * A simple data class to store all of the parameters.
     */
    public static class Data {
        public double scale;
        public double[] equilibriumFrequencies;
        public double[] transitionFrequencies;
        public double[] speciesNetworkParameters;
        public double geneStayProb;
        public double speciesStayProb;


        public int getSize(){
            int size = 3 + speciesNetworkParameters.length;
            if(equilibriumFrequencies!=null){
                size += equilibriumFrequencies.length;
                size += transitionFrequencies.length;
            }
            return size;
        }

        @Override
        public String toString() {
            return "Data{" +
                    "scale=" + scale +
                    ", equilibriumFrequencies=" + Arrays.toString(equilibriumFrequencies) +
                    ", transitionFrequencies=" + Arrays.toString(transitionFrequencies) +
                    ", speciesNetworkParameters=" + Arrays.toString(speciesNetworkParameters) +
                    ", geneStayProb=" + geneStayProb +
                    ", speciesStayProb=" + speciesStayProb +
                    //", speciesNetworkInheritanceProbabilities=" + Arrays.toString(speciesNetworkInheritanceProbabilities) +
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