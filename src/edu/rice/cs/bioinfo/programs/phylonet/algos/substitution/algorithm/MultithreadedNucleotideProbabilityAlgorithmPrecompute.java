package edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm;

import com.google.common.collect.Sets;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

public class MultithreadedNucleotideProbabilityAlgorithmPrecompute
{
    public static void precomputeUniqueObservations(List<NucleotideProbabilityAlgorithm> algorithms, List<NucleotideObservation> observations, ExecutorService threadPool) {
        Set<NucleotideObservation> uniqueObservations = Sets.newHashSet(observations);

        PrecomputeTask[] tasks = new PrecomputeTask[uniqueObservations.size() * algorithms.size()];
        int i = 0;
        for (NucleotideObservation obs: uniqueObservations)
        {
            for (NucleotideProbabilityAlgorithm algo : algorithms)
            {
                PrecomputeTask task = new PrecomputeTask(obs,algo);
                tasks[i++] = task;
            }
        }

        try {
            List<Future<PrecomputeResult>> results = threadPool.invokeAll(Arrays.asList(tasks));

            for (Future<PrecomputeResult>  result : results)
            {
                result.get().seedResult();
            }
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    private static class PrecomputeTask implements Callable<PrecomputeResult>
    {
        NucleotideObservation obs;
        NucleotideProbabilityAlgorithm algo;


        public PrecomputeTask(NucleotideObservation obs, NucleotideProbabilityAlgorithm algo)
        {
            this.obs = obs;
            this.algo = algo;
        }


        @Override
        public PrecomputeResult call() throws Exception {
            return new PrecomputeResult(algo.getProbability(obs),obs,algo);
        }
    }

    private static class PrecomputeResult
    {
        double value;
        NucleotideObservation obs;
        NucleotideProbabilityAlgorithm algo;

        private PrecomputeResult(double value, NucleotideObservation obs, NucleotideProbabilityAlgorithm algo) {
            this.value = value;
            this.obs = obs;
            this.algo = algo;
        }

        public void seedResult()
        {
            if (Double.isNaN(value) || Double.isInfinite(value))
                throw new RuntimeException("Bad result");
            algo.seedProbability(obs,value);
        }
    }

}
