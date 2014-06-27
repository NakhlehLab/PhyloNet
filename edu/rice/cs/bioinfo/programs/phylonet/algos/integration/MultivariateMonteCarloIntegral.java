package edu.rice.cs.bioinfo.programs.phylonet.algos.integration;

import edu.rice.cs.bioinfo.library.programming.Func1;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Created by yunyu on 6/2/14.
 */
public class MultivariateMonteCarloIntegral{
    private boolean _printDetails = false;



    // bins is the number of divisions that each axis it split into (so the full number of bins is this to the
    // power of the dimension of the function
    // sampleSize is the number of samples PER BIN

    public MultivariateMonteCarloIntegral(int sampleSize, int bins) {
        this._sampleSize = sampleSize;
        this._bins = bins;
    }

    public void setPrintDetails(boolean printDetails){
        _printDetails = printDetails;
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
        Random random;
        if(_seed == null){
            random = new Random();
        }
        else{
            random = new Random(_seed);
        }
        int totalBins = (int)Math.pow(_bins, dim);
        double[] steps = new double[dim];
        double totalArea=1;

        double[] currentCorner = new double[dim];
        for(int i=0; i<dim; i++){
            steps[i] = (maxes[i]-mins[i])/_bins;
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
            for (int j=1; j <= _sampleSize; j++) {
                double[] sample = new double[dim];
                for(int k=0; k<sample.length; k++){
                    double randomNumber = random.nextDouble();
                    //System.out.println(randomNumber);
                    sample[k] = sampleBase[k] + randomNumber * (steps[k]);
                }
                integral += f.execute(sample);
            }
        }
        //System.out.println(integral +"*"+ totalArea +"/"+ _sampleSize*totalBins);
        integral *= totalArea/((double)_sampleSize*totalBins);
        return integral;
    }


    public List<double[]> getSamplePoints(int dim, double min, double max) {
        Random random;
        if(_seed == null){
            random = new Random();
        }
        else{
            random = new Random(_seed);
        }

        int totalBins = (int)Math.pow(_bins, dim);

        double[] steps = new double[dim];
        Arrays.fill(steps, (max - min)/_bins);

        double[] currentCorner = new double[dim];
        Arrays.fill(currentCorner, min);

        List<double[]> binCorners = new ArrayList<double[]>();

         for(int i=0; i<totalBins; i++){
            binCorners.add(currentCorner.clone());
            for(int j=dim-1; j>=0; j--){
                currentCorner[j] += steps[j];
                if(Math.abs(currentCorner[j]-max)<0.00000000001){
                    currentCorner[j] = min;
                }
                else{
                    break;
                }
            }
        }


        List<double[]> sampleList = new ArrayList<double[]>();
        for(double[] sampleBase: binCorners){
            for (int j=1; j <= _sampleSize; j++) {
                double[] sample = new double[dim];
                //do {
                    for (int k = 0; k < sample.length; k++) {
                        double randomNumber = random.nextDouble();
                        //System.out.println(randomNumber);
                        sample[k] = sampleBase[k] + randomNumber * (steps[k]);
                    }
                //}while(!validator.execute(sample));
                sampleList.add(sample);
                //printArray(sample);
            }
        }

        return sampleList;
    }

    private void printArray(double array[]){
        for(double value: array){
            System.out.print(value + " ");
        }
        System.out.println();
    }

    public static double integrate(int dim, double min, double max, double[] evaluatedValue){
        double totalArea = Math.pow(max - min, dim);
        double totalValue = 0;
        for(double value: evaluatedValue){
            totalValue += value;
        }
        //System.out.println(totalValue +"*"+ totalArea +"/"+ evaluatedValue.length);
        return totalValue * totalArea / evaluatedValue.length;
    }

    /*
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
    */

    public int getSampleSize(){
        return _sampleSize;
    }


    public int getBins(){
        return _bins;
    }

    public void setRandomSeed(Long seed){
        _seed = seed;
    }

    private int _sampleSize;
    private int _bins;
    private Long _seed = null;
}
