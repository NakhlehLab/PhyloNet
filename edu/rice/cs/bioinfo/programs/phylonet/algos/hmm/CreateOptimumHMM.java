package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm;


import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.model.HmmOptimizationFunction;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.model.HmmParameters;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.GTRModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.JCModel;
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
import java.util.Random;

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
        //Random random = new Random();
        myConfiguration.seed = ConfigurationBuilder.createSeed();
        System.out.println("\nThe seed was " + myConfiguration.seed);

        rand = new MersenneTwister(myConfiguration.seed);
        f = new HmmOptimizationFunction(net, speciesToAlleles, observations, config);
        params = f.getParams();

        System.out.println("\nStates: " + params.numberOfGeneTrees() * params.getNumberOfAlleleMappings() + " Params: " + params.numberOfParameters());
    }

    public long getSeed()
    {
        return myConfiguration.seed;
    }

    double[] initialGuess()
    {
        HmmParameters.Data initial = new HmmParameters.Data();

        initial.scale = 1;
        switch (myConfiguration.MODEL)
        {
            case GTR:
                initial.equilibriumFrequencies = new double[] {.25,.25,.25,.25};
                initial.transitionFrequencies = new double[6];

                //Initialize Random Values for transition frequencies and gene tree lengths.
                for (int i = 0; i < initial.transitionFrequencies.length; i++) {
                    initial.transitionFrequencies[i] = 0.1 + .9 * rand.nextDouble();
                    //initial.transitionFrequencies[i] = 1/6.0;
                }
        }

        initial.speciesNetworkParameters = new double[params.numberOfNetworkParameters()];
        //int numReticulations = params.numberOfReticulations();
        for (int i = 0; i < initial.speciesNetworkParameters.length; i++) {
            initial.speciesNetworkParameters[i] = .001 + rand.nextDouble() * 5;
            /*
            if(i < numReticulations){
                initial.speciesNetworkParameters[i] = rand.nextDouble();
                //initial.speciesNetworkParameters[i] = 0.5;
            }
            else {
                initial.speciesNetworkParameters[i] = .001 + rand.nextDouble() * 5;
                //initial.speciesNetworkParameters[i] = 3;
            }
            */
        }

        initial.geneStayProb = 0.5;
        initial.speciesStayProb = 0.5;

        return params.encode(initial);
    }

    SimpleBounds bounds()
    {
        lower.scale = 0.001;
        upper.scale = Double.POSITIVE_INFINITY;
        //lower.scale = 0.99;
        //upper.scale = 1.01;
        switch (myConfiguration.MODEL)
        {
            case GTR:
                lower.equilibriumFrequencies =new double[] {0.1,0.1,0.1,0.1};
                upper.equilibriumFrequencies = new double[] {1,1,1,1};
                //lower.equilibriumFrequencies =new double[] {0.24,0.24,0.24,0.24};
                //upper.equilibriumFrequencies = new double[] {0.26,0.26,0.26,0.26};

                lower.transitionFrequencies = new double[] {.1,.1,.1,.1,.1,.1};
                upper.transitionFrequencies = new double[] {1, 1, 1, 1, 1, 1};
        }

        //lower.transitionFrequencies = new double[] {1/6.0-0.01,1/6.0-0.01,1/6.0-0.01,1/6.0-0.01,1/6.0-0.01,1/6.0-0.01};
        //upper.transitionFrequencies = new double[] {1/6.0+0.01,1/6.0+0.01,1/6.0+0.01,1/6.0+0.01,1/6.0+0.01,1/6.0+0.01};

        int numNetworkParameters = params.numberOfNetworkParameters();
        lower.speciesNetworkParameters = new double[numNetworkParameters];
        upper.speciesNetworkParameters = new double[numNetworkParameters];
        //int numReticulations = params.numberOfReticulations();
        for(int i = 0; i < numNetworkParameters; i++){
            lower.speciesNetworkParameters[i] = .00001;
            upper.speciesNetworkParameters[i] = Double.POSITIVE_INFINITY;
            /*
            if(i < numReticulations){
                lower.speciesNetworkParameters[i] = 0;
                upper.speciesNetworkParameters[i] = 1;
                //lower.speciesNetworkParameters[i] = 0.49;
                //upper.speciesNetworkParameters[i] = 0.51;
            }
            else{
                lower.speciesNetworkParameters[i] = .00001;
                upper.speciesNetworkParameters[i] = Double.POSITIVE_INFINITY;
                //lower.speciesNetworkParameters[i] = 2.9;
                //upper.speciesNetworkParameters[i] = 3.1;
            }
            */
        }

        lower.geneStayProb =  0.0001; // Lower bound cannot be zero as results in a divide by zero and NaNs
        upper.geneStayProb =  1;

        lower.speciesStayProb =  0.0001; // Lower bound cannot be zero as results in a divide by zero and NaNs
        upper.speciesStayProb =  1;

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


