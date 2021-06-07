package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm;

import java.security.SecureRandom;

public class ConfigurationBuilder
{
    private int ITERATIONS;                 //The number of iterations to run the optimizer.
    private boolean MARKOV;                //Choose between Markov and non-Markov models.
    private boolean USEGENETREELENGTHS;
    private boolean USEFASTTREES;
    private Configuration.SubstitutionModelType MODEL;
    private Configuration.AlgorithmType ALGORITHM;
    private int numberOfRuns;
    private long seed;
    private int threads;

    private int numberOfFolds;
    private boolean PPATraining;
    private double threshold = 0.0;

    private ConfigurationBuilder()
    {}

    public Configuration build(String inputFile)
    {
        return new Configuration(ITERATIONS,MARKOV,USEGENETREELENGTHS,USEFASTTREES,MODEL,ALGORITHM,
                numberOfRuns, seed, threads, PPATraining, numberOfFolds, inputFile, threshold);
    }

    public ConfigurationBuilder withThreshold(double thresh)
    {
        this.threshold = thresh;
        return this;
    }

    public ConfigurationBuilder withPPATraining(int numberOfFolds)
    {
        this.numberOfFolds = numberOfFolds;
        this.PPATraining = true;
        return this;
    }

    public ConfigurationBuilder withThreads(int threads)
    {
        this.threads = threads;
        return this;
    }

    public ConfigurationBuilder withSeed(int seed)
    {
        this.seed = seed;
        return this;
    }

    public ConfigurationBuilder withIterations(int iterations)
    {
        this.ITERATIONS = iterations;
        return this;
    }

    public ConfigurationBuilder withNumberOfRuns(int num)
    {
        this.numberOfRuns = num;
        return this;
    }

    public ConfigurationBuilder withGTR()
    {
        this.MODEL = Configuration.SubstitutionModelType.GTR;
        return this;
    }

    public static ConfigurationBuilder getSNAPP()
    {
        return getDefault().withSNAPP();
    }

    public static ConfigurationBuilder getNormal()
    {
        return getDefault().withNormal();
    }



    private ConfigurationBuilder withNormal()
    {
        this.ALGORITHM = Configuration.AlgorithmType.NORMAL;
        return this;
    }

    private ConfigurationBuilder withSNAPP()
    {
        this.ALGORITHM = Configuration.AlgorithmType.SNAPP;
        return this;
    }

    private static ConfigurationBuilder getDefault()
    {
        ConfigurationBuilder config = new ConfigurationBuilder();

        config.ITERATIONS = 300;
        config.MARKOV = true;

        config.USEGENETREELENGTHS = false;
        config.USEFASTTREES = false;

        config.MODEL = Configuration.SubstitutionModelType.JC;
        config.numberOfRuns = 10;
        config.seed = createSeed();
        config.threads = 1;
        config.PPATraining = false;
        config.threshold = 0.0;

        return config;
    }

    public static long createSeed()
    {
        return new SecureRandom().nextLong();
    }

}
