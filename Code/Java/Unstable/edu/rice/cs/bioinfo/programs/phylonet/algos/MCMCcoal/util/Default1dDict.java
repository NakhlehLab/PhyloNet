package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util;

import java.util.HashMap;

/**
 * An 1-dimensional defaultdict mapping int to int.
 * Created by Xinhao Liu on 1/4/20.
 */
public class Default1dDict extends HashMap<Integer, Integer> {
    private int defaultValue;

    public Default1dDict(int value) {
        this.defaultValue = value;
    }

    public Integer get(int key) {
//        V returnValue = super.get(key);
//        if (returnValue == null) {
//            returnValue = defaultValue;
//            this.put((K) key, returnValue);
//        }
//        return returnValue;
        return super.getOrDefault(key, defaultValue);
    }

    // test
    public static void main(String[] args) {
        Default1dDict test = new Default1dDict(0);
        System.out.println(test.get(90));
        test.put(1, test.get(1) + 10);
        System.out.println(test.get(1));
        test.put(1, test.get(1) + 20);
        System.out.println(test.get(1));
        test.put(100, 22);
        System.out.println(test.get(100));

        System.out.println(test.size());

    }
}

