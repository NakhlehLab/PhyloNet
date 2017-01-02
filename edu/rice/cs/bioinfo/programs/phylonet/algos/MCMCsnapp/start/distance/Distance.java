package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.start.distance;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Given sequence alignment, calculate the distance between each pair of taxa.
 * Created by wendingqiao on 1/27/16.
 */
public abstract class Distance {

    protected Map<String, Map<String, Double>> distanceMap;
    protected Map<String, String> alignment;
    protected Set<String> taxa;

    public Distance(Map<String, String> sequences) {
        distanceMap = new HashMap<>();
        alignment = sequences;
        taxa = sequences.keySet();
    }

    protected void buildDistanceMap() {
        for(String taxon : taxa) {
            distanceMap.put(taxon, new HashMap<String, Double>());
        }
        for(String taxon1 : taxa) {
            for(String taxon2 : taxa) {
                if(distanceMap.get(taxon1).containsKey(taxon2)) continue;
                if(taxon1.compareTo(taxon2) == 0) {
                    distanceMap.get(taxon1).put(taxon2, 0.0);
                } else {
                    double dist = getPairwiseDistance(taxon1, taxon2);
                    distanceMap.get(taxon1).put(taxon2, dist);
                    distanceMap.get(taxon2).put(taxon1, dist);
                }
            }
        }
    }

    public double getDistance(String taxon1, String taxon2) {
        return distanceMap.get(taxon1).get(taxon2);
    }

    protected abstract double getPairwiseDistance(String taxon1, String taxon2);

    public Set<String> getTaxa() {
        return taxa;
    }

    public Map<String, String> getAlignment() {
        return alignment;
    }
}
