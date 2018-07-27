package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util;

import java.util.Random;

/**
 * Randomizer for the whole program. Default seed is 12345678.
 * Created by wendingqiao on 2/15/16.
 */
public class Randomizer {

    final private static Random _rand = new Random(Utils._SEED);

    private Randomizer() {}

    public static int getRandomInt(int upperBound) {
        return _rand.nextInt(upperBound);
    }

    public static double getRandomDouble() {
        return _rand.nextDouble();
    }

    public static void setSeed(long seed) {
        _rand.setSeed(seed);
    }

    public static Random getRandom() { return _rand; }
}
