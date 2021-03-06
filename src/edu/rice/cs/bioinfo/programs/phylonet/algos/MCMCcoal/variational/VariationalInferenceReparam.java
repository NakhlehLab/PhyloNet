package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational;

import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm.HmmBuilder;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.hmm.HmmCore;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.Prior;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.VariationalVariable;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Variational inference procedure on the reparametrized model
 *
 * Created by Xinhao Liu on 1/24/21.
 */
public class VariationalInferenceReparam {
    private VariationalModelReparam variationalPosterior;
    private Prior prior;

    public VariationalInferenceReparam(ModelTree model, Prior prior) {
        variationalPosterior = new VariationalModelReparam(model);
        this.prior = prior;
    }

    public void run() {
        System.out.println("");
        System.out.println("Inference starts...");
        for (int i = 0; i < Utils.nIterations; i++) {
            if (i == Utils.nIterations - 1) {
                Utils.lastiter = true;
            }

            System.out.println("=====================================ITERATION " + i +"=====================================");
            List<Tuple3<VariationalVariable, Double, Double>> BLGradientList = new ArrayList<>();
            List<Tuple3<VariationalVariable, Double, Double>> NGradientList = new ArrayList<>();

            List<Map<VariationalVariable, Double>> samples = new ArrayList<>();
            for (int s = 0; s < Utils.nSamples; s++) {
                samples.add(variationalPosterior.sample());
            }

            // evaluate coalhmm likelihood for each sample
            List<Double> logJoints = new ArrayList<>();
            for (int s = 0; s < samples.size(); s++) {
                Map<VariationalVariable, Double> sample = samples.get(s);
                variationalPosterior.setTreeBySample(sample);
                ModelTree sampledVariationalPosterior = variationalPosterior.getModel();
                // code below is to circumvent a strange bug causing likelihood to return positive value.
                // not sure why. maybe beagle screwed up.
                double likelihood = 1;
                while (likelihood >= 0) {
                    HmmBuilder builder = new HmmBuilder(sampledVariationalPosterior.getTree(), sampledVariationalPosterior.getRecombRate());
                    // log time used by builder.build(), add up to time used by building HMM by simulation
                    long buildingStartTime = System.currentTimeMillis();
                    HmmCore hmm = builder.build();
                    Utils.buildingTime += System.currentTimeMillis() - buildingStartTime;
                    // log time used by hmm.logLikelihood(), add up to time used by forward algorithm
                    long likelihoodStartTime = System.currentTimeMillis();
                    likelihood = hmm.logLikelihood();
                    Utils.likelihoodTime += System.currentTimeMillis() - likelihoodStartTime;
                }
                double prior = this.prior.logPrior(sampledVariationalPosterior);
                logJoints.add(likelihood + prior);
            }

            // loop over branchlengthvariables and popsizevariables, compute the gradient for mean and stddev of each separately
            // Branch length variables
            for (VariationalVariable var:variationalPosterior.getBranchLengthVariableList()) {
                // mean
                List<Double> fmean = new ArrayList<>();
                List<Double> hmean = new ArrayList<>();
                for (int s = 0; s < samples.size(); s++) {
                    double sample = samples.get(s).get(var);
                    double fmean_s = var.scoreMean(sample) * (logJoints.get(s) - var.logDensity(sample));
                    double hmean_s = var.scoreMean(sample);
                    fmean.add(fmean_s);
                    hmean.add(hmean_s);
                }
                double amean = estimateCovariance(fmean, hmean) / estimateVariance(hmean);

                double sumMean = 0.0;
                for (int s = 0; s < samples.size(); s++) {
                    sumMean += fmean.get(s) - amean * hmean.get(s);
                }
                double meanGradient = sumMean / samples.size();

                // standard deviation
                List<Double> fStdDev = new ArrayList<>();
                List<Double> hStdDev = new ArrayList<>();
                for (int s = 0; s < samples.size(); s++) {
                    double sample = samples.get(s).get(var);
                    double fStdDev_s = var.scoreStandardDeviation(sample) * (logJoints.get(s) - var.logDensity(sample));
                    double hStdDev_s = var.scoreStandardDeviation(sample);
                    fStdDev.add(fStdDev_s);
                    hStdDev.add(hStdDev_s);
                }
                double aStdDev = estimateCovariance(fStdDev, hStdDev) / estimateVariance(hStdDev);

                double sumStdDev = 0.0;
                for (int s = 0; s < samples.size(); s++) {
                    sumStdDev += fStdDev.get(s) - aStdDev * hStdDev.get(s);
                }
                double stdDevGradient = sumStdDev / samples.size();

                BLGradientList.add(new Tuple3<>(var, meanGradient, stdDevGradient));
            }

            // Population size variables
            for (VariationalVariable var:variationalPosterior.getPopSizeVariableList()) {
                // mean
                List<Double> fmean = new ArrayList<>();
                List<Double> hmean = new ArrayList<>();
                for (int s = 0; s < samples.size(); s++) {
                    double sample = samples.get(s).get(var);
                    double fmean_s = var.scoreMean(sample) * (logJoints.get(s) - var.logDensity(sample));
                    double hmean_s = var.scoreMean(sample);
                    fmean.add(fmean_s);
                    hmean.add(hmean_s);
                }
                double amean = estimateCovariance(fmean, hmean) / estimateVariance(hmean);

                double sumMean = 0.0;
                for (int s = 0; s < samples.size(); s++) {
                    sumMean += fmean.get(s) - amean * hmean.get(s);
                }
                double meanGradient = sumMean / samples.size();

                // standard deviation
                List<Double> fStdDev = new ArrayList<>();
                List<Double> hStdDev = new ArrayList<>();
                for (int s = 0; s < samples.size(); s++) {
                    double sample = samples.get(s).get(var);
                    double fStdDev_s = var.scoreStandardDeviation(sample) * (logJoints.get(s) - var.logDensity(sample));
                    double hStdDev_s = var.scoreStandardDeviation(sample);
                    fStdDev.add(fStdDev_s);
                    hStdDev.add(hStdDev_s);
                }
                double aStdDev = estimateCovariance(fStdDev, hStdDev) / estimateVariance(hStdDev);

                double sumStdDev = 0.0;
                for (int s = 0; s < samples.size(); s++) {
                    sumStdDev += fStdDev.get(s) - aStdDev * hStdDev.get(s);
                }
                double stdDevGradient = sumStdDev / samples.size();

                NGradientList.add(new Tuple3<>(var, meanGradient, stdDevGradient));
            }

            double OLD_BRANCH_LENGTH_MEAN_LEARNING_RATE = 0.0;
            double OLD_POP_SIZE_MEAN_LEARNING_RATE = 0.0;
            double OLD_BRANCH_LENGTH_STDDEV_LEARNING_RATE = 0.0;
            double OLD_POP_SIZE_STDDEV_LEARNING_RATE = 0.0;
            // Also decrease learning rate when i == 0 because the gradient will always be +/- 1 in the first iteration of AdaGrad
            if (Utils.ILLEGAL_SAMPLE_GENERATED || i == 0) {
                if (i > 0) {
                    System.out.println("!!!");
                    System.out.println("INFO: Illegal value sampled this iteration. Automatically decreasing learning rate...");
                    System.out.println("If you see too many iterations with this warning, please restart with lower learning rates and smaller standard deviations.");
                    System.out.println("!!!");
                }
                OLD_BRANCH_LENGTH_MEAN_LEARNING_RATE = Utils.BRANCH_LENGTH_MEAN_LEARNING_RATE;
                Utils.BRANCH_LENGTH_MEAN_LEARNING_RATE /= 1000;

                OLD_POP_SIZE_MEAN_LEARNING_RATE = Utils.POP_SIZE_MEAN_LEARNING_RATE;
                Utils.POP_SIZE_MEAN_LEARNING_RATE /= 1000;

                OLD_BRANCH_LENGTH_STDDEV_LEARNING_RATE = Utils.BRANCH_LENGTH_STDDEV_LEARNING_RATE;
                Utils.BRANCH_LENGTH_STDDEV_LEARNING_RATE /= 1000;

                OLD_POP_SIZE_STDDEV_LEARNING_RATE = Utils.POP_SIZE_STDDEV_LEARNING_RATE;
                Utils.POP_SIZE_STDDEV_LEARNING_RATE /= 1000;
            }

            // gradient update
            System.out.println("********************* Branch lengths in postorder traversal *********************");
            for (Tuple3<VariationalVariable, Double, Double> gradientTuple:BLGradientList) {
                VariationalVariable var = gradientTuple.Item1;
                System.out.println("Current mean: " + var.getMean() + ", Current standard deviation: " + var.getStandardDeviation());
                double meanGradient = gradientTuple.Item2;
                double stdDevGradient = gradientTuple.Item3;
                var.meanGradientUpdate(meanGradient);
                var.standardDeviationGradientUpdate(stdDevGradient);
                System.out.println("-----");
            }
            System.out.println("********************* Population sizes in postorder traversal *********************");
            for (Tuple3<VariationalVariable, Double, Double> gradientTuple:NGradientList) {
                VariationalVariable var = gradientTuple.Item1;
                System.out.println("Current mean: " + var.getMean() + ", Current standard deviation: " + var.getStandardDeviation());
                double meanGradient = gradientTuple.Item2;
                double stdDevGradient = gradientTuple.Item3;
                var.meanGradientUpdate(meanGradient);
                var.standardDeviationGradientUpdate(stdDevGradient);
                System.out.println("-----");
            }

            if (Utils.ILLEGAL_SAMPLE_GENERATED || i == 0) {
                Utils.ILLEGAL_SAMPLE_GENERATED = false;
                Utils.BRANCH_LENGTH_MEAN_LEARNING_RATE = OLD_BRANCH_LENGTH_MEAN_LEARNING_RATE;
                Utils.POP_SIZE_MEAN_LEARNING_RATE = OLD_POP_SIZE_MEAN_LEARNING_RATE;
                Utils.BRANCH_LENGTH_STDDEV_LEARNING_RATE = OLD_BRANCH_LENGTH_STDDEV_LEARNING_RATE;
                Utils.POP_SIZE_STDDEV_LEARNING_RATE = OLD_POP_SIZE_STDDEV_LEARNING_RATE;
            }
        }
    }

    private static double estimateCovariance(List<Double> x, List<Double> y) {
        double meanX = x.stream().mapToDouble(a -> a).average().getAsDouble();
        double meanY = y.stream().mapToDouble(a -> a).average().getAsDouble();

        double sum = 0.0;
        for (int i = 0; i < x.size(); i++) {
            sum += (x.get(i) - meanX) * (y.get(i) - meanY);
        }

        return sum / (x.size() - 1);
    }

    private static double estimateVariance(List<Double> x) {
        double meanX = x.stream().mapToDouble(a -> a).average().getAsDouble();

        double sum = 0.0;
        for (double num:x) {
            sum += (num - meanX) * (num - meanX);
        }

        return sum / (x.size() - 1);
    }

    public VariationalModelReparam getVariationalPosterior() {
        return variationalPosterior;
    }
}
