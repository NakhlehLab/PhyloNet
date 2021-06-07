package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm;

public class Configuration
{
    public final int ITERATIONS;                 //The number of iterations to run the optimizer.
    public final boolean MARKOV;                //Choose between Markov and non-Markov models.
    public final boolean USEGENETREELENGTHS;
    public final boolean USEFASTTREES;
    public final SubstitutionModelType MODEL;
    public final AlgorithmType ALGORITHM;
    public int numberOfRuns;
    public long seed;
    public final int threads;
    public final boolean PPATraining;
    public final int numberOfFolds;
    public double threshold = 0.0;

    public final String inputFile;

    public Configuration(int ITERATIONS, boolean MARKOV, boolean USEGENETREELENGTHS, boolean USEFASTTREES,
                         SubstitutionModelType MODEL, AlgorithmType ALGORITHM, int numberOfRuns, long seed, int threads,
                         boolean PPATraining, int numberOfFolds, String inputFile, double threshold)
    {
        this.ITERATIONS = ITERATIONS;
        this.MARKOV = MARKOV;
        this.USEGENETREELENGTHS = USEGENETREELENGTHS;
        this.USEFASTTREES = USEFASTTREES;
        this.MODEL = MODEL;
        this.ALGORITHM = ALGORITHM;
        this.numberOfRuns = numberOfRuns;
        this.seed = seed;
        this.threads = threads;
        this.PPATraining = PPATraining;
        this.numberOfFolds = numberOfFolds;
        this.inputFile = inputFile;
        this.threshold = threshold;
    }

    public enum AlgorithmType {
        INTEGRATION,
        SNAPP,
        NORMAL
    }

    public enum SubstitutionModelType {
        GTR,
        JC
    }
}
