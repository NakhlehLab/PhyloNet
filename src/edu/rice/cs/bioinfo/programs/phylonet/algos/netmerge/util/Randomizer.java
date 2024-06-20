package edu.rice.cs.bioinfo.programs.phylonet.algos.netmerge.util;
/*
 * @ClassName:   Randomizer
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        9/14/23 2:32 PM
 */

import java.util.Random;

public class Randomizer {

    /* Constructor */
    final private static Random _rand = new Random(Utils._SEED);

    private Randomizer() {}

    public static int getRandomInt(int upperBound) {
        return _rand.nextInt(upperBound);
    }

    public static double getRandomDouble() {
//        System.out.println(_rand.);
        return _rand.nextDouble();
    }


    public static boolean getRandomBoolean() {
        return _rand.nextBoolean();
    }

    public static void setSeed(long seed) {
        _rand.setSeed(seed);
    }

    public static Random getRandom() { return _rand; }

}
