package edu.rice.cs.bioinfo.programs.phylonet.algos.integration;

import edu.rice.cs.bioinfo.library.programming.Func1;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by yunyu on 6/2/14.
 */
public class MultivariateMonteCarloIntegral{

    // bins is the number of divisions that each axis it split into (so the full number of bins is this to the
    // power of the dimension of the function
    // sampleSize is the number of samples PER BIN

    public MultivariateMonteCarloIntegral(int sampleSize, int bins) {
        this.sampleSize = sampleSize;
        this.bins = bins;
    }

    public MultivariateMonteCarloIntegral(int sampleSize) {
        this(sampleSize, 1);
    }

    /**
     * @return the approximate integral of the given function
     * within the given range using simple monte carlo integration.
     * @param f the function whose integral is of interest
     * @param mins the minimum value of the function
     * @param maxes the  upper limit of the function
     */
    public double integrate(Func1<double[],Double> f, int dim, double[] mins, double[] maxes) {

        int totalBins = (int)Math.pow(bins, dim);
        double[] steps = new double[dim];
        double totalArea=1;

        double[] currentCorner = new double[dim];
        for(int i=0; i<dim; i++){
            steps[i] = (maxes[i]-mins[i])/bins;
            totalArea *= maxes[i]-mins[i];
            currentCorner[i] = mins[i];
        }

        List<double[]> binCorners = new ArrayList<double[]>();

        /*
        int[] binStepCounter = new int[dim];
        for(int i=0; i<totalBins; i++){
            double[] corner = new double[dim];
            for(int j=0; j<dim; j++){
                corner[j] = mins[j] + steps[j] * binStepCounter[j];
            }
            binCorners.add(corner);
            increaseBinStepCounter(binStepCounter, bins);
        }
        */

        for(int i=0; i<totalBins; i++){
            binCorners.add(currentCorner.clone());
            for(int j=dim-1; j>=0; j--){
                currentCorner[j] += steps[j];
                if(Math.abs(currentCorner[j]-maxes[j])<0.00000000001){
                    currentCorner[j] = mins[j];
                }
                else{
                    break;
                }
            }
        }


        double integral = 0.0;
        for(double[] sampleBase: binCorners){
            for (int j=1; j <= sampleSize; j++) {
                double[] sample = new double[dim];
                for(int k=0; k<sample.length; k++){
                    sample[k] = sampleBase[k] + Math.random()*(steps[k]);
                }
                integral += f.execute(sample);
            }
        }

        integral *= totalArea/((double)sampleSize*totalBins);
        return integral;
    }

    private void increaseBinStepCounter(int[] binStepCounter, int maxBin){
        for(int i=binStepCounter.length-1; i>=0; i--){
            binStepCounter[i]++;
            if(binStepCounter[i] == maxBin){
                binStepCounter[i] = 0;
            }
            else{
                break;
            }
        }
    }

    protected int getSampleSize(){
        return sampleSize;
    }
    protected int getBins(){
        return bins;
    }

    private int sampleSize;
    private int bins;

}
