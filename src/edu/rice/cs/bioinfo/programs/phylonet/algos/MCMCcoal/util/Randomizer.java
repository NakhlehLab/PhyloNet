package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util;

import java.util.Random;

/**
 * Randomizer for the whole program. Default seed is 12345678.
 * Created by Xinhao Liu on 1/1/2020.
 * Modified from edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Randomizer
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

    /**
     * Returns a pseudo-random number between min and max, inclusive.
     * The difference between min and max can be at most Integer.MAX_VALUE - 1
     */
    public static int getIntRange(int min, int max) {
        return _rand.nextInt((max - min) + 1) + min;
    }

    public static void setSeed(long seed) {
        _rand.setSeed(seed);
    }

    public static Random getRandom() { return _rand; }


    public static void main(String[] args) {
        //System.out.println(Randomizer.getRandomInt(100));
        System.out.println(Randomizer.getIntRange(500, 50000));
        System.out.println(Randomizer.getIntRange(500, 50000));
        System.out.println(Randomizer.getIntRange(500, 50000));

    }
}
