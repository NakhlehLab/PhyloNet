package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.start.distance;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.start.distance.Distance;

import java.util.*;

/**
 * P-distance
 * Created by wendingqiao on 1/26/16.
 */
public class PDistance extends Distance {

    private static Set<Character> missingMark = new HashSet<>(Arrays.asList(new Character[]{'-', '?', 'N'}));

    public PDistance(Map<String, String> sequences) {
        super(sequences);
        buildDistanceMap();
    }

    protected double getPairwiseDistance(String taxon1, String taxon2) {
        if(!alignment.containsKey(taxon1) || !alignment.containsKey(taxon2)) {
            System.err.println("cannot find sequence of taxa " + taxon1 + ", " + taxon2 + " from map");
            System.exit(1);
        }
        String seq1 = alignment.get(taxon1);
        String seq2 = alignment.get(taxon2);
        char c1, c2;

        double sumDistance = 0.0;
        for (int i = 0; i < seq1.length(); i++) {
            c1 = seq1.charAt(i);
            c2 = seq2.charAt(i);
            if(!missingMark.contains(c1) && !missingMark.contains(c2) && c1 != c2) {
                sumDistance += 1.0;
            }
        }
        return sumDistance / seq1.length();
    }
}
