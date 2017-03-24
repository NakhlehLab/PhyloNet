package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 2/15/17
 * Time: 10:59 AM
 * To change this template use File | Settings | File Templates.
 */
public class RPattern {
    private Map<String, R> pattern;
    private List<String> sortedKeys;

    public RPattern(Map<String, R> pattern) {
        this.pattern = pattern;
        this.sortedKeys = new ArrayList<>(pattern.keySet());
        Collections.sort(sortedKeys);
    }

    @Override
    public int hashCode() {
        int ret = 0;
        for(String species : pattern.keySet())
            ret += species.hashCode() * pattern.get(species).hashCode();
        return ret;
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof RPattern))
            return false;
        if (obj == this)
            return true;

        RPattern rhs = (RPattern) obj;

        for(String species : pattern.keySet()) {
            if(!rhs.getPattern().containsKey(species))
                return false;
            if(!rhs.getR(species).equals(pattern.get(species)))
                return false;
        }

        return true;

    }

    public Map<String, R> getPattern() {
        return pattern;
    }

    public List<String> getNames() { return sortedKeys; }

    public R getR(String species) {
        return pattern.get(species);
    }

    public int sumLineages() {
        int sum = 0;
        for(String species : pattern.keySet()) {
            sum += pattern.get(species).n;
        }
        return sum;
    }

    public int diff(RPattern p2) {
        int count = 0;
        for(String species : pattern.keySet()) {
            if(!p2.getR(species).equals(pattern.get(species)))
                count++;
        }
        return count;
    }

    public List<String> leaves2update(RPattern p2) {
        List<String> ret = new ArrayList<>();
        for(String species : pattern.keySet()) {
            if(!p2.getR(species).equals(pattern.get(species)))
                ret.add(species);
        }
        return ret;
    }

    public int compareTo(RPattern p2) {
        for(String species : sortedKeys) {
            int k = this.pattern.get(species).compareTo(p2.pattern.get(species));
            if(k != 0)
                return k;
        }
        return 0;
    }

    public boolean isMonomorphic() {
        int count[] = new int[R.dims];
        int n = sumLineages();
        for(String species : pattern.keySet()) {
            R r = pattern.get(species);
            for(int i = 0 ; i < R.dims ; i++) {
                count[i] += r.values[i];
            }
        }

        boolean flag = true;
        for(int i = 0 ; i < R.dims ; i++) {
            if(count[i] == n)
                return true;
            if(count[i] != 0)
                flag = false;
        }
        return flag;
    }
}
