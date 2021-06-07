package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.model;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPP.SNAPPAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.fasttrees.FastTrees;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.Configuration;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYFForMulTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.GTRModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.JCModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.SubstitutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.util.LimitedRunMultivariateFunction;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.util.NetworkLengthApplier;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeWithBranchLengthProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.HmmNetworkUtils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.HmmTreeUtils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.apache.commons.math3.analysis.MultivariateFunction;

import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * This class does the task of creating the hmm model and observations, and then finding the total likelihood of the whole model.
 */
public class HmmOptimizationFunction extends LimitedRunMultivariateFunction implements MultivariateFunction, AutoCloseable
{

    List<NucleotideObservation> observations;
    List<JahmmNucleotideObservation> jahmmObservations;

    List<Map<String,List<String>>> allAlleleMappings;
    List<Tree> allGeneTrees;

    List<double[]> geneTreeProbabilities;

    HmmParameters params;

    Network<String> net;
    Map<String, List<String>> speciesToAlleles;

    ExecutorService threadPool;
    Configuration config;

    /**
     * @param net              A species network which represents the genetic history of a species including introgression.
     * @param speciesToAlleles A mapping of a species name to a list of alleles for that species.
     * @param observations     The nucleotide observations to work with.
     */
    public HmmOptimizationFunction(Network<String> net, Map<String, List<String>> speciesToAlleles, List<NucleotideObservation> observations, Configuration config)
    {
        super(config.ITERATIONS);
        this.config = config;

        this.net = net;
        this.speciesToAlleles = speciesToAlleles;
        this.observations = observations;

        createGeneTrees();

        int numberOfAlleles = observations.get(0).getAlleles().size();

        params = new HmmParameters(numberOfAlleles, net, speciesToAlleles, allGeneTrees.size(),config);
        jahmmObservations = JahmmNucleotideObservation.wrapObservations(observations);
        threadPool = Executors.newFixedThreadPool(config.threads);
    }

    /**
     * Closes and cleans up the thread pool.
     */
    private void shutdownThreadPool()
    {
        threadPool.shutdown();
        try
        {
            threadPool.awaitTermination(2, TimeUnit.SECONDS);
        } catch (InterruptedException e)
        {
            e.printStackTrace();
        }
    }

    /**
     * Creates all distinct gene trees to use in the calculations.
     */
    private void createGeneTrees()
    {
        if (config.USEFASTTREES)
        {
            FastTrees f = new FastTrees(observations);
            System.out.println(f.calculateTrees());

            allGeneTrees = f.calculateTrees();
        }
        else if (config.ALGORITHM == Configuration.AlgorithmType.SNAPP)
            allGeneTrees = Arrays.asList((Tree)new STITree<>()); //Create a dummy gene tree as SNAPP does not use gene trees.
        else
        {
            Collection<String> alleles = observations.get(0).getAlleles();
            allGeneTrees = Trees.generateAllBinaryTrees(alleles.toArray(new String[alleles.size()]));
            /*
            for(Tree gt: allGeneTrees){
                for(TNode node: gt.getNodes()){
                    node.setParentDistance(1.0);
                }
            }
            */
        }
    }




    private void computeGeneTreeProbabilities(Tree mulTree, Map<NetNode, List<TNode>> netNode2treeNodes)
    {
        geneTreeProbabilities = new ArrayList<>();

        for (Map<String,List<String>> alleleMapping: allAlleleMappings)
        {
            geneTreeProbabilities.add(computeGeneTreeProbForOneAlleleMapping(mulTree, alleleMapping, netNode2treeNodes));

        }
    }

    private double[] computeGeneTreeProbForOneAlleleMapping(Tree mulTree, Map<String,List<String>> alleleMapping, Map<NetNode, List<TNode>> netNode2treeNodes)
    {
        if (config.USEGENETREELENGTHS)
        {
            return computeGeneTreeProbWithLengths((STITree)mulTree);
        } else
        {
            return computeGeneTreeProbWithTopologies(mulTree, alleleMapping, netNode2treeNodes);
        }
    }

    //TODO: Not correct
    private double[] computeGeneTreeProbWithLengths(STITree<String> speciesTree)
    {
        Network<Integer> test = HmmNetworkUtils.toNetwork(speciesTree);

        double[] result = new double[allGeneTrees.size()];

       GeneTreeWithBranchLengthProbabilityYF g = new GeneTreeWithBranchLengthProbabilityYF(test,allGeneTrees,null);

        g.calculateGTDistribution(result);
        return result;
    }



    private double[] computeGeneTreeProbWithTopologies(Tree mulTree, Map<String,List<String>> alleleMapping, Map<NetNode, List<TNode>> netNode2treeNodes)
    {
        double[] result = new double[allGeneTrees.size()];

        int index = 0;
        for(Tree gt: allGeneTrees) {
            GeneTreeProbabilityYFForMulTree calculator = new GeneTreeProbabilityYFForMulTree();
            result[index++] = calculator.calculateProbability((Network)net, mulTree, netNode2treeNodes, alleleMapping, gt);
        }
        return result;
    }

    /**
     * Returns the list of observations.
     */
    public List<NucleotideObservation> getObservations()
    {
        return observations;
    }

    /**
     * Returns the list of jahmm observations.
     */
    public List<JahmmNucleotideObservation> getJahmmObservations()
    {
        return jahmmObservations;
    }

    /**
     * Computes the number of states for the HMM.
     *
     * @return The number of states.
     */
    private int numberOfStates()
    {
        return allGeneTrees.size() * allAlleleMappings.size();
    }

    /**
     * Creates a hidden markov model from the given input.
     *
     * @param input A double array of parameters, see {@link HmmParameters} for format.
     * @return A hidden markov model.
     */
    private Hmm<JahmmNucleotideObservation> getHmm(double[] input)
    {

        HmmParameters.Data data = params.decode(input);

        return createHmm(
                getModel(data.scale, data.equilibriumFrequencies, data.transitionFrequencies),
                data.speciesNetworkParameters,
                data.geneStayProb,
                data.speciesStayProb);

    }

    private SubstitutionModel getModel(double scale, double[] equilibriumFrequencies, double[] transitionFrequencies)
    {


        switch (config.MODEL)
        {
            case GTR:
                return new GTRModel(scale, equilibriumFrequencies, transitionFrequencies);
            case JC:
                return new JCModel(scale);
        }

        throw new RuntimeException("This should never be used. GTR and JC Model selection is incorrect. See CreateOptimumHMM.java" +
                "or HmmOptimizationFunction.java");
    }

    public Hmm<JahmmNucleotideObservation> getBestHmm()
    {
        return getHmm(getBestInput());
    }


    /**
     * Creates a hidden markov model given the following parameters.
     *
     * @param model           A {@link edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.SubstitutionModel} which represents the nucleotide changes.
     * @return A hidden markov model for those parameters.
     */
    private Hmm<JahmmNucleotideObservation> createHmm(SubstitutionModel model, double[] networkParameters, double geneStayProb, double speciesStayProb)
    {

        new NetworkLengthApplier(net, networkParameters, speciesToAlleles).apply();

        Map<NetNode, List<TNode>> netNode2treeNodes = new HashMap<>();
        allAlleleMappings = new ArrayList<>();
        Tree mulTree = HmmNetworkUtils.generateMulTree(net, speciesToAlleles, netNode2treeNodes, allAlleleMappings);


        if (config.ALGORITHM == Configuration.AlgorithmType.NORMAL)
        {
            computeGeneTreeProbabilities(mulTree, netNode2treeNodes);
        }

        //TODO Comment me a lot more.
        double[] initialStateProbabilities = new double[numberOfStates()];
        double[][] transitionProbabilities = new double[numberOfStates()][numberOfStates()];

        List<HmmState> states = new ArrayList<>();
        List<NucleotideProbabilityAlgorithm> algorithms = new ArrayList<>();

        //Create all HMM states.
        for (int geneTreeIndex = 0; geneTreeIndex < allGeneTrees.size(); geneTreeIndex++)
        {

            Tree geneTree = allGeneTrees.get(geneTreeIndex);

            NucleotideProbabilityAlgorithm algo = null;

            if (config.ALGORITHM == Configuration.AlgorithmType.NORMAL)
            {
                algo = new FelsensteinAlgorithmExtended(geneTree, model);
                //algo = new FelsensteinAlgorithm(geneTree, model);
                algorithms.add(algo);
            }

            for (int alleleMappingIndex = 0; alleleMappingIndex < allAlleleMappings.size(); alleleMappingIndex++)
            {
                Map<String,List<String>> alleleMapping = allAlleleMappings.get(alleleMappingIndex);
                //Compute an initial state probability.
                if (config.ALGORITHM == Configuration.AlgorithmType.NORMAL)
                    initialStateProbabilities[params.getStateIndex(geneTreeIndex, alleleMappingIndex)] = geneTreeProbabilities.get(alleleMappingIndex)[geneTreeIndex];
                else
                    initialStateProbabilities[params.getStateIndex(geneTreeIndex, alleleMappingIndex)] = 1;

                double geneTreeSum = 0;
                for (int nextStateAlleleMappingIndex = 0; nextStateAlleleMappingIndex < allAlleleMappings.size(); nextStateAlleleMappingIndex++)
                {
                    geneTreeSum = 0;
                    for (int nextStateGeneTreeIndex = 0; nextStateGeneTreeIndex < allGeneTrees.size(); nextStateGeneTreeIndex++)
                    {

                        double prob = 1;

                        if (config.ALGORITHM == Configuration.AlgorithmType.NORMAL)
                            prob = geneTreeProbabilities.get(nextStateAlleleMappingIndex)[nextStateGeneTreeIndex];

                        if (config.MARKOV && nextStateAlleleMappingIndex == alleleMappingIndex)
                        {
                            //Have a special probability when switching to another gene tree
                            double geneSwitchFactor = (nextStateGeneTreeIndex == geneTreeIndex) ? geneStayProb : (1 - geneStayProb) / (allGeneTrees.size() - 1);
                            prob *= geneSwitchFactor;
                        }

                        geneTreeSum += prob;
                        transitionProbabilities[params.getStateIndex(geneTreeIndex, alleleMappingIndex)][params.getStateIndex(nextStateGeneTreeIndex, nextStateAlleleMappingIndex)] = prob;
                    }

                    /*
                    if(geneTreeSum != 0) {
                        for (int nextStateGeneTreeIndex = 0; nextStateGeneTreeIndex < allGeneTrees.size(); nextStateGeneTreeIndex++) {
                            transitionProbabilities[params.getStateIndex(geneTreeIndex, alleleMappingIndex)][params.getStateIndex(nextStateGeneTreeIndex, nextStateAlleleMappingIndex)] /= geneTreeSum;
                        }
                    }
                    */


                    double speciesSwitchFactor = (nextStateAlleleMappingIndex == alleleMappingIndex) ? speciesStayProb : (1 - speciesStayProb) / (allAlleleMappings.size() - 1);
                    for (int otherGeneTreeIndex = 0; otherGeneTreeIndex < allGeneTrees.size(); otherGeneTreeIndex++)
                    {
                        if(geneTreeSum != 0)
                            transitionProbabilities[params.getStateIndex(geneTreeIndex, alleleMappingIndex)][params.getStateIndex(otherGeneTreeIndex, nextStateAlleleMappingIndex)] /= geneTreeSum;

                        // Add species tree dependency
                        transitionProbabilities[params.getStateIndex(geneTreeIndex, alleleMappingIndex)][params.getStateIndex(otherGeneTreeIndex, nextStateAlleleMappingIndex)] *= speciesSwitchFactor;
                        //geneTreeSum += transitionProbabilities[params.getStateIndex(geneTreeIndex, alleleMappingIndex)][params.getStateIndex(otherGeneTreeIndex, nextStateAlleleMappingIndex)];
                    }
                }


/*
                if(geneTreeSum != 0) {
                    for (int nextStateAlleleMappingIndex = 0; nextStateAlleleMappingIndex < allAlleleMappings.size(); nextStateAlleleMappingIndex++) {
                        for (int nextStateGeneTreeIndex = 0; nextStateGeneTreeIndex < allGeneTrees.size(); nextStateGeneTreeIndex++) {

                            transitionProbabilities[params.getStateIndex(geneTreeIndex, alleleMappingIndex)][params.getStateIndex(nextStateGeneTreeIndex, nextStateAlleleMappingIndex)] /= geneTreeSum;
                        }
                    }
                }
*/
/*
                if (config.ALGORITHM == Configuration.AlgorithmType.INTEGRATION)
                {
                    //Noting the seed to that it is different from the CreateOptimium hmm seed.
                    algo = new IntegrationAlgorithm(speciesTree, geneTree, model);
                    algorithms.add(algo);
                }

                if (config.ALGORITHM == Configuration.AlgorithmType.SNAPP)
                {
                    algo = new SNAPPAlgorithm(speciesTree, model, 0.018);
                    algorithms.add(algo);
                }
*/
                states.add(new HmmState(alleleMapping, geneTree, model, algo));

            }

        }


        MultithreadedNucleotideProbabilityAlgorithmPrecompute.precomputeUniqueObservations(algorithms, observations, threadPool);

        return new Hmm<>(initialStateProbabilities, transitionProbabilities, states);
    }




    @Override
    protected double calculateValue(double[] input)
    {
        Hmm<JahmmNucleotideObservation> hmm = getHmm(input);

        return hmm.lnProbability(jahmmObservations);
    }

    public HmmParameters getParams()
    {
        return params;
    }

    public List<Tree> getTrees()
    {
        return allGeneTrees;
    }

    public Tree getAlleleMappings(List<Map<String,List<String>>> allAlleleMappings)
    {
        Tree mulTree = HmmNetworkUtils.generateMulTree(net, speciesToAlleles, new HashMap<NetNode, List<TNode>>(), allAlleleMappings);
        for(TNode node: mulTree.getNodes()) {
            if (!node.isLeaf()) {
                ((STINode) node).setName(TNode.NO_NAME);
            }
        }
        return mulTree;
        /*

        }
        List<STITree<String>> allSpeciesTrees = new ArrayList<>();
        for(Map<String,List<String>> alleleMapping: allAlleleMappings) {
            Tree st = Trees.readTree(mulTree.toString());
            ((STITree)st).constrainByLeaves(alleleMapping.keySet());
            Trees.removeBinaryNodes((MutableTree)st);
            for(TNode node: st.getNodes()){
                if(node.isLeaf()){
                    String[] fields = ;
                    ((STINode)node).setName(fields[1]);
                }
            }
            allSpeciesTrees.add((STITree)st);
        }
        return allSpeciesTrees;
        */
    }

    @Override
    public void close()
    {
        shutdownThreadPool();
    }
}
