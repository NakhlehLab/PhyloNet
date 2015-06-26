package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPP.SNAPPAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.fasttrees.FastTrees;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.Configuration;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.GTRModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.JCModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.SubstitutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.util.LimitedRunMultivariateFunction;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.util.NetworkLengthApplier;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeWithBranchLengthProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.HmmNetworkUtils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
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

    List<STITree<String>> allSpeciesTrees;
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
     * Creates all of the gene trees to use in the calculations.
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
        }


    }

    private void computeGeneTreeProbabilities()
    {
        geneTreeProbabilities = new ArrayList<>();

        for (STITree<String> speciesTree: allSpeciesTrees)
        {
            geneTreeProbabilities.add(computeGeneTreeProbForSpeciesTree(speciesTree));

        }
    }

    private double[] computeGeneTreeProbForSpeciesTree(STITree<String> speciesTree)
    {
        if (config.USEGENETREELENGTHS)
        {
            return computeGeneTreeProbWithLengths(speciesTree);
        } else
        {
            return computeGeneTreeProbWithTopologies(speciesTree);
        }
    }

    private double[] computeGeneTreeProbWithLengths(STITree<String> speciesTree)
    {
        Network<Integer> test = HmmNetworkUtils.toNetwork(speciesTree);

        double[] result = new double[allGeneTrees.size()];

       GeneTreeWithBranchLengthProbabilityYF g = new GeneTreeWithBranchLengthProbabilityYF(test,allGeneTrees,null);

        g.calculateGTDistribution(result);
        return result;
    }

    private double[] computeGeneTreeProbWithTopologies(STITree<String> speciesTree)
    {
        GeneTreeProbabilityYF g = new GeneTreeProbabilityYF();

        Network test = HmmNetworkUtils.toNetwork(speciesTree);

        double[] result = new double[allGeneTrees.size()];

        g.calculateGTDistribution(test, allGeneTrees, null, result);
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
        return allGeneTrees.size() * allSpeciesTrees.size();
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
                data.speciesStayProbs,
                data.geneStayProb,
                data.allGeneTreeLengths,
                data.speciesNetworkBranchLengths);

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

    public void applyGeneTreeLengths(double[] allGeneTreeLengths)
    {

        for (int i = 0; i < allGeneTrees.size(); i++)
        {
            HmmTreeUtils.setParentDistances(allGeneTrees.get(i), getGeneTreeLengthsForIndex(allGeneTreeLengths, i));
        }
    }

    /**
     * Creates a hidden markov model given the following parameters.
     *
     * @param model           A {@link edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.SubstitutionModel} which represents the nucleotide changes.
     * @param speciesStayProbs The probability that a state will stay in its current species tree.
     * @param geneStayProb    The probability that a state will stay in its current gene tree (only if also staying in the same species tree).
     * @return A hidden markov model for those parameters.
     */
    private Hmm<JahmmNucleotideObservation> createHmm(SubstitutionModel model, double[] speciesStayProbs, double geneStayProb, double[] allGeneTreeLengths, double[] networkLengths)
    {

        new NetworkLengthApplier(net,networkLengths).apply();

        allSpeciesTrees = HmmNetworkUtils.generateParentalTrees(net, speciesToAlleles);

        if (config.ALGORITHM == Configuration.AlgorithmType.NORMAL)
        {
            applyGeneTreeLengths(allGeneTreeLengths);
            computeGeneTreeProbabilities();
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
                algorithms.add(algo);
            }

            for (int speciesTreeIndex = 0; speciesTreeIndex < allSpeciesTrees.size(); speciesTreeIndex++)
            {
                STITree<String> speciesTree = allSpeciesTrees.get(speciesTreeIndex);
                //Compute an initial state probability.
                if (config.ALGORITHM == Configuration.AlgorithmType.NORMAL)
                    initialStateProbabilities[params.getStateIndex(geneTreeIndex, speciesTreeIndex)] = geneTreeProbabilities.get(speciesTreeIndex)[geneTreeIndex];
                else
                    initialStateProbabilities[params.getStateIndex(geneTreeIndex, speciesTreeIndex)] = 1;


                for (int nextStateSpeciesTreeIndex = 0; nextStateSpeciesTreeIndex < allSpeciesTrees.size(); nextStateSpeciesTreeIndex++)
                {
                    double geneTreeSum = 0;
                    for (int nextStateGeneTreeIndex = 0; nextStateGeneTreeIndex < allGeneTrees.size(); nextStateGeneTreeIndex++)
                    {

                        double prob = 1;

                        if (config.ALGORITHM == Configuration.AlgorithmType.NORMAL)
                            prob = geneTreeProbabilities.get(nextStateSpeciesTreeIndex)[nextStateGeneTreeIndex];

                        // Add gene tree dependency
                        if (config.MARKOV && nextStateSpeciesTreeIndex == speciesTreeIndex)
                        {
                            //Have a special probability when switching to another gene tree
                            double geneSwitchFactor = (nextStateGeneTreeIndex == geneTreeIndex) ? geneStayProb : (1 - geneStayProb) / (allGeneTrees.size() - 1);
                            prob *= geneSwitchFactor;
                        }
                        geneTreeSum += prob;
                        transitionProbabilities[params.getStateIndex(geneTreeIndex, speciesTreeIndex)][params.getStateIndex(nextStateGeneTreeIndex, nextStateSpeciesTreeIndex)] = prob;
                    }
                    for (int otherGeneTreeIndex = 0; otherGeneTreeIndex < allGeneTrees.size(); otherGeneTreeIndex++)
                    {
                        transitionProbabilities[params.getStateIndex(geneTreeIndex, speciesTreeIndex)][params.getStateIndex(otherGeneTreeIndex, nextStateSpeciesTreeIndex)] /= geneTreeSum;

                        // Add species tree dependency
                        double speciesSwitchFactor = getSpeciesSwitchProbability(speciesStayProbs, speciesTreeIndex, nextStateSpeciesTreeIndex);
                        transitionProbabilities[params.getStateIndex(geneTreeIndex, speciesTreeIndex)][params.getStateIndex(otherGeneTreeIndex, nextStateSpeciesTreeIndex)] *= speciesSwitchFactor;
                    }
                }


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

                states.add(new HmmState(speciesTree, geneTree, model, algo));

            }

        }

        double initialStateProbSum = 0;
        for (double initialStateProbability : initialStateProbabilities)
        {
            initialStateProbSum += initialStateProbability;
        }

        for (int i = 0; i < initialStateProbabilities.length; i++)
        {
            initialStateProbabilities[i] /= initialStateProbSum;
        }


        MultithreadedNucleotideProbabilityAlgorithmPrecompute.precomputeUniqueObservations(algorithms, observations, threadPool);

        return new Hmm<>(initialStateProbabilities, transitionProbabilities, states);
    }

    private double getSpeciesStayProb(double speciesStayProbs[],int speciesTreeIndex, Configuration config)
    {
        if (config.USEMULTIPLESTAYPROBABILITIES)
            return speciesStayProbs[speciesTreeIndex];
        else
            return speciesStayProbs[0];
    }

    private double getSpeciesSwitchProbability(double speciesStayProbs[], int speciesTreeIndex, int otherSpeciesTreeIndex)
    {
        if (config.MARKOV)
        {
            double speciesStayProb = getSpeciesStayProb(speciesStayProbs,speciesTreeIndex,config);
            return (otherSpeciesTreeIndex == speciesTreeIndex) ? speciesStayProb : (1 - speciesStayProb) / (allSpeciesTrees.size() - 1);
        }
        else
            return 1.0 / (allSpeciesTrees.size());
    }

    private double[] getGeneTreeLengthsForIndex(double[] allGeneTreeLengths, int geneTreeIndex)
    {
        return Arrays.copyOfRange(
                allGeneTreeLengths,
                geneTreeIndex * params.numberOfTreeLengthsPerGeneTree(),
                (geneTreeIndex + 1) * params.numberOfTreeLengthsPerGeneTree());
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

    public List<STITree<String>> getSpeciesTrees()
    {
        return allSpeciesTrees;
    }

    @Override
    public void close()
    {
        shutdownThreadPool();
    }
}
