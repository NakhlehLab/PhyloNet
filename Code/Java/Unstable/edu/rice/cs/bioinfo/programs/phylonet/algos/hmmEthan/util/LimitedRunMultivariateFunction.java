package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.util;

import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model.HmmOptimizationFunction;
import org.apache.commons.math3.analysis.MultivariateFunction;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

public abstract class LimitedRunMultivariateFunction implements MultivariateFunction
{
    double bestSoFar = Double.NEGATIVE_INFINITY;
    double[] bestInput;
    List<Double> probScoreData;
    int run = 0;
    final int maxIterations;
    Logger logg = Logger.getLogger(HmmOptimizationFunction.class.toString());

    public LimitedRunMultivariateFunction(int maxIterations)
    {
        probScoreData = new ArrayList<Double>();
        this.maxIterations = maxIterations;
    }

    /**
     * Returns the input that had the highest likelihood when {@link #value(double[])} was called.
     *
     * @return The best input.
     */
    public double[] getBestInput()
    {
        checkRanOnce();
        return bestInput;
    }

    /**
     * Returns the score that the best input recieved.
     * Is equivalent to value(getBestInput()).
     *
     * @return The best score.
     */
    public double getBestInputScore()
    {
        checkRanOnce();
        return bestSoFar;
    }

    private void checkRanOnce()
    {
        if (run ==0)
            throw new RuntimeException("You are trying to get the best result on a function that didn't even run once.");
    }

    /**
     * Calculates the likelihood of a model specified by the input parameters.
     *
     * @param input A double array of parameters, see {@link edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.model.HmmParameters} for format.
     * @return The likelihood of the model specified by the input parameters.
     */
    @Override
    public double value(double[] input)
    {

        double result = calculateValue(input);

        if (result > 0)
        {
            throw new RuntimeException("The result is greater than 0, thus p>1.");
        }

        if (Double.isNaN(result) || Double.isInfinite(result))
        {
            throw new RuntimeException("Non finite result " + result + " The input was " + Arrays.toString(input));
        }

        probScoreData.add(result);             //Store the result, to be plotted.
        if (result > bestSoFar)
        {
            bestSoFar = result;
            bestInput = input.clone();
        }


        //Print to the console.
        int printMod = maxIterations/100;
        if (printMod == 0)
            printMod = 1;

        if (run % printMod == 0)
            logg.log(Level.INFO, "Iteration : {0}, Prob : {1}", new Object[]{run, result});

        run++;

        return result;
    }

    protected abstract double calculateValue(double[] input);

    public List<Double> getProbScoreData()
    {
        return probScoreData;
    }
}
