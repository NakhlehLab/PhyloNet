package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.summary;


import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Compute ESS effective sample size at runtime
 * The method is adapted from BEAST2 dr.inference.trace/TraceCorrelation
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

//    public double add(double val) {
//
//        values.add(val);
//        sumOfLagged.add(0.0); // lag = M-1
//        sum += val;
//        int M = values.size();
//        // calculate auto correlation for selected lag times.
//        int nMax = Math.min( MAX_LAG, M );
//
//        double ave = sum / (double)(M);
//        double[] fAC = new double[M];
//        double fSum1 = sum;
//        double fSum2 = sum;
//        for(int lag = 0; lag < nMax; lag++) {
////            sumOfLagged.set(lag, sumOfLagged.get(lag) + values.get(M-1) * values.get(M-1-lag));
////            fAC[lag] = sumOfLagged.get(lag) - (fSum1+fSum2) * ave;
////            fAC[lag] /= (double)(M-lag);
////            fAC[lag] += ave * ave;
////            fSum1 -= values.get(M-1-lag);
////            fSum2 -= values.get(lag);
//            for (int j = 0; j < M - lag; j++) {
//                final double del1 = values.get(j) - ave;
//                final double del2 = values.get(j + lag) - ave;
//                fAC[lag] += (del1 * del2);
//                //varGammaStat[lag] += (del1*del1*del2*del2);
//            }
//
//            fAC[lag] /= ((double) (M - lag));
//        }
//        // calculate auto correlation time fACT
//        double integralACT = fAC[0];
//        for(int lag = 1; lag < nMax; lag++) {
//            if(lag % 2 == 0) {
//                // fancy stopping criterion
//                double tmp = fAC[lag - 1] + fAC[lag];
//                if(tmp < 0) break;
//                integralACT += 2.0 * tmp;
//            }
////            if(lag > 5 && fAC[lag] / fAC[0] < 0.1) break;
//        }
//        double fACT = integralACT / fAC[0];
//        // return ESS
////        System.out.println(M+"/"+fACT);
//        return M == 1 ? 1.0 : (double)(M) / fACT;
//    }


    public static void main(String[] args) {
        ESS ess = new ESS();

//        double[] arr = {-47906.93633, -47901.64519, -47887.74341, -47887.27738, -47912.89766, -47927.58793, -47920.97556, -47936.8555, -47911.79779, -47911.18417};
//        for (double val: arr){
//            System.out.println(ess.add(val));
//        }


        try {
            String filepath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/butterfly/data/combined.log";
//            String filepath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/butterfly/data/trinets13/gamma/slurm-5510908_0_posterior.log";
            BufferedReader br = new BufferedReader(new FileReader(filepath));
            String s;
            String[] ss;
//            List<Double> posterior_list = new ArrayList<>();
            while((s = br.readLine()) != null) {
//                System.out.println(s);
                if (s.startsWith("Sample") || s.startsWith("sample")) continue;
                ss = s.trim().split("\\s+");
                double posterior = Double.parseDouble(ss[1]);
//                posterior_list.add(posterior);
                System.out.println(posterior+","+ess.add(posterior));
            }
        }catch (IOException e){

        }


    }
}
