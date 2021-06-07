package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.model;

import java.util.ArrayList;
import java.util.List;

public class HmmStateUtilities {

    public static List<Integer> getGeneTrees(int[] hmmStateSequence, HmmParameters params) {
        //Decomposing HMM States to geneTrees.
        List<Integer> geneTrees = new ArrayList<Integer>();
        for (int i : hmmStateSequence) {
            geneTrees.add(i / params.getNumberOfAlleleMappings());
        }
        return geneTrees;
    }

    public static List<Integer> getAlleleMappings(int[] hmmStateSequence, HmmParameters params) {
        //Decomposing HMM States to speciesTrees.
        List<Integer> alleleMappings = new ArrayList<Integer>();
        for (int i : hmmStateSequence) {
            alleleMappings.add(i % params.getNumberOfAlleleMappings());
        }
        return alleleMappings;
    }


}
