package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 2/19/17
 * Time: 3:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class SplittingIndexVector {
    private int[] _arr;

    public SplittingIndexVector(int [] arr) {
        _arr = Arrays.copyOf(arr, arr.length);
    }

    @Override
    public int hashCode() {
        int primes[] = {2, 3, 5, 7, 11};
        int ret = 0;
        for(int i = 0 ; i < _arr.length ; i++) {
            ret += _arr[i] * primes[i];
        }
        return ret;
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof SplittingIndexVector))
            return false;
        if (obj == this)
            return true;

        SplittingIndexVector rhs = (SplittingIndexVector) obj;

        if(this._arr.length != rhs._arr.length)
            return false;

        for(int i = 0 ; i < _arr.length ; i++)
            if(this._arr[i] != rhs._arr[i])
                return false;

        return true;

    }

    private void getAllCompatibleHelper(int [] arr, List<SplittingIndexVector> ret, int i) {
        if(i == arr.length) {
            ret.add(new SplittingIndexVector(arr));
            return;
        }

        if(arr[i] == 0)
            getAllCompatibleHelper(arr, ret, i + 1);
        else {
            int tmp = arr[i];
            arr[i] = 0;
            getAllCompatibleHelper(arr, ret, i + 1);
            arr[i] = tmp;
            getAllCompatibleHelper(arr, ret, i + 1);
        }
    }

    public List<SplittingIndexVector> getAllCompatible() {
        List<SplittingIndexVector> ret = new ArrayList<>();
        int[] arr = Arrays.copyOf(_arr, _arr.length);
        getAllCompatibleHelper(arr, ret, 0);
        return ret;
    }
}
