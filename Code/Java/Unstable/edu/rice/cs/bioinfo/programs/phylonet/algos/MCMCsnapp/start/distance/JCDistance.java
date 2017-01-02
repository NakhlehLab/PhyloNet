package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.start.distance;

import java.util.Map;

/**
 * Jukes-Cantor distance
 * Created by wendingqiao on 1/26/16.
 */
public class JCDistance extends Distance {

    private static double MIN_DISTANCE = 0.01;
    private static double MAX_DISTANCE = 0.75;
    private PDistance _pDistance;

    public JCDistance(Map<String, String> sequences) {
        super(sequences);
        _pDistance = new PDistance(sequences);
        buildDistanceMap();
    }

    @Override
    protected double getPairwiseDistance(String taxon1, String taxon2) {
        double pDist = _pDistance.getDistance(taxon1, taxon2);
        if (pDist <= MIN_DISTANCE) return MIN_DISTANCE;
        if (pDist >= MAX_DISTANCE) return MAX_DISTANCE;

        double expDist = -MAX_DISTANCE * Math.log(1.0 - pDist / MAX_DISTANCE);
        return expDist < MAX_DISTANCE ? expDist : MAX_DISTANCE;
    }
}
