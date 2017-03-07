package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;

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

    public RPattern(Map<String, R> pattern) {
        this.pattern = pattern;
    }

    @Override
    public int hashCode() {
        int ret = 0;
        for(String species : pattern.keySet())
            ret += pattern.get(species).hashCode();
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
}
