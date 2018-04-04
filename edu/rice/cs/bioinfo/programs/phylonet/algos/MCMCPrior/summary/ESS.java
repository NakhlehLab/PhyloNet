package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.summary;

import java.util.ArrayList;
import java.util.List;

/**
 * Compute ESS effective sample size at runtime
 * The method is adapted from BEAST2
 *
 * Created by wendingqiao on 11/1/14.
 */
public class ESS {

    // list of values
    private List<Double> values;

    // dividend of pho_k: auto correlation for lag k
    // sum::x_i * x_i+k
    private List<Double> sumOfLagged;

    private double sum = 0.0;

    /**
     * Calculate effective sample size based on Auto Correlation(AC)
     */
    public ESS() {
        values = new ArrayList<>();
        sumOfLagged = new ArrayList<>();
    }

    final static int MAX_LAG = 2000;

    /**
     * Add new value to the list
     * @param val   new value
     * @return  ess value : the number of post burn-in GTT samples M / (1 + 2 * sum_pho_k)
     */
    public double add(double val) {

        values.add(val);
        sumOfLagged.add(0.0); // lag = M-1
        sum += val;
        int M = values.size();
        // calculate auto correlation for selected lag times.
        int nMax = Math.min( MAX_LAG, M );

        double ave = sum / (double)(M);
        double[] fAC = new double[M];
        double fSum1 = sum;
        double fSum2 = sum;
        for(int lag = 0; lag < nMax; lag++) {
            sumOfLagged.set(lag, sumOfLagged.get(lag) + values.get(M-1) * values.get(M-1-lag));
            fAC[lag] = sumOfLagged.get(lag) - (fSum1+fSum2) * ave;
            fAC[lag] /= (double)(M-lag);
            fAC[lag] += ave * ave;
            fSum1 -= values.get(M-1-lag);
            fSum2 -= values.get(lag);
        }
        // calculate auto correlation time fACT
        double integralACT = fAC[0];
        for(int lag = 1; lag < nMax; lag++) {
            if(lag % 2 == 0) {
                // fancy stopping criterion
                double tmp = fAC[lag - 1] + fAC[lag];
                if(tmp < 0) break;
                integralACT += 2.0 * tmp;
            }
//            if(lag > 5 && fAC[lag] / fAC[0] < 0.1) break;
        }
        double fACT = integralACT / fAC[0];
        // return ESS
        return M == 1 ? 1.0 : (double)(M) / fACT;
    }

}
