package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan;


import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model.HmmOptimizationFunction;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model.HmmParameters;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optim.*;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.random.MersenneTwister;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * The main class for the problem we are working on. Alex Yozzo & Ethan Steinberg
 * Creates a HMM to find the probability of each state at each position of a genome.
 */
public class CreateOptimumHMM {

    Configuration myConfiguration;


    HmmParameters.Data lower = new HmmParameters.Data();
    HmmParameters.Data upper = new HmmParameters.Data();

    HmmParameters params;
    HmmOptimizationFunction f;
    MersenneTwister rand;



    public CreateOptimumHMM(Network<String> net, Map<String, List<String>> speciesToAlleles,
                            List<NucleotideObservation> observations, Configuration config) {

        myConfiguration = config;
        //myConfiguration.seed = Long.parseLong("8699692860394415879");
        myConfiguration.seed = ConfigurationBuilder.createSeed();
        System.out.println("\nThe seed was " + myConfiguration.seed);

        rand = new MersenneTwister(myConfiguration.seed);
        f = new HmmOptimizationFunction(net, speciesToAlleles, observations, config);
        params = f.getParams();

        System.out.println("\nStates: " + params.numberOfGeneTrees() * params.getNumberOfSpeciesTrees() + " Params: " + params.numberOfParameters());
    }

    public long getSeed()
    {
        return myConfiguration.seed;
    }

    double[] initialGuess()
    {
        HmmParameters.Data initial = new HmmParameters.Data();


        initial.scale = 1;
        initial.equilibriumFrequencies = new double[] {.25,.25,.25,.25};
        initial.transitionFrequencies = new double[6];

        //Initialize Random Values for transition frequencies and gene tree lengths.
        for (int i = 0; i < initial.transitionFrequencies.length; i++) {
            initial.transitionFrequencies[i] = 0.1 + .9 * rand.nextDouble();
        }

        initial.speciesStayProbs = new double[params.numberOfSpeciesStayProbs()];
        Arrays.fill(initial.speciesStayProbs,.99);

        initial.geneStayProb = .5;
        initial.allGeneTreeLengths = new double[params.numberOfGeneTreeLengths()];

        for (int i = 0; i < initial.allGeneTreeLengths.length; i++) {
            initial.allGeneTreeLengths[i] = .1 + .9 * rand.nextDouble();
        }

        initial.speciesNetworkBranchLengths = new double[params.numberOfNetworkLengths()];
        for (int i = 0; i < initial.speciesNetworkBranchLengths.length; i++) {
            initial.speciesNetworkBranchLengths[i] = .001 + rand.nextDouble() * 5;
        }

        return params.encode(initial);
    }

    SimpleBounds bounds()
    {
        lower.scale = 0.001;
        upper.scale = Double.POSITIVE_INFINITY;

        lower.equilibriumFrequencies =new double[] {0.1,0.1,0.1,0.1};
        upper.equilibriumFrequencies = new double[] {1,1,1,1};

        lower.transitionFrequencies = new double[] {.1,.1,.1,.1,.1,.1};
        upper.transitionFrequencies = new double[] {1,  1, 1, 1, 1, 1};

        lower.speciesStayProbs = new double[params.numberOfSpeciesStayProbs()];
        upper.speciesStayProbs = new double[params.numberOfSpeciesStayProbs()];

        Arrays.fill(lower.speciesStayProbs,.95);
        Arrays.fill(upper.speciesStayProbs,1);


        lower.geneStayProb =  0.0001; // Lower bound cannot be zero as results in a divide by zero and NaNs
        upper.geneStayProb =  1;

        lower.allGeneTreeLengths = new double[params.numberOfGeneTreeLengths()];
        upper.allGeneTreeLengths = new double[params.numberOfGeneTreeLengths()];
        Arrays.fill(lower.allGeneTreeLengths,.00001);
        Arrays.fill(upper.allGeneTreeLengths,Double.POSITIVE_INFINITY);

        lower.speciesNetworkBranchLengths = new double[params.numberOfNetworkLengths()];
        upper.speciesNetworkBranchLengths = new double[params.numberOfNetworkLengths()];
        Arrays.fill(lower.speciesNetworkBranchLengths,.00001);
        Arrays.fill(upper.speciesNetworkBranchLengths,Double.POSITIVE_INFINITY);

        return new SimpleBounds(params.encode(lower),params.encode(upper));
    }

    public HmmOptimizationFunction solve()
    {
        MultivariateOptimizer optimizer = new BOBYQAOptimizer(params.numberOfParameters()+5);

        try {
            PointValuePair result = optimizer.optimize(
                    new ObjectiveFunction(f),
                    new InitialGuess(initialGuess()),
                    GoalType.MAXIMIZE,
                    bounds(),
                    new MaxEval(myConfiguration.ITERATIONS),
                    new MaxIter(1));
            System.out.println( "\n Solving " + Arrays.toString(result.getPoint()));
        }
        catch (TooManyEvaluationsException e){
            System.out.println("Optimization finished.");
        }

        System.out.println("Best Input = " + Arrays.toString(f.getBestInput()));
        System.out.println("The seed for pseudo-random numbers was: " + myConfiguration.seed);

        return f;
    }
}


