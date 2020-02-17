package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util;

import java.util.HashMap;

/**
 * A 2-dimensional defaultdict mapping (int, int) pair to int.
 * Created by Xinhao Liu on 1/4/20.
 */
public class Default2dDict extends HashMap<Integer, Default1dDict> {

    private Default1dDict get(int key) {
        super.putIfAbsent(key, new Default1dDict(0));
        return super.get(key);
    }

    public int get(int key1, int key2) {
//        Default1dDict returnValue = super.get(key);
//        if (returnValue == null) {
//            returnValue = new Default1dDict(0);
//            this.put((int) key, returnValue);
//        }
        //super.putIfAbsent(key1, new Default1dDict(0));
        return this.get(key1).get(key2);
    }

    public void put(int key1, int key2, int value) {
        this.get(key1).put(key2, value);
    }

    // test
    public static void main(String[] args) {
        Default2dDict test = new Default2dDict();

        System.out.println(test.get(4, 4));

        test.put(4, 5, test.get(4, 5) + 1);
        System.out.println(test.get(4, 5));

        System.out.println(test.get(5, 5));
        System.out.println(test.get(5, 5));
        test.put(5, 5, 22);
        System.out.println(test.get(5, 5));
    }
}
