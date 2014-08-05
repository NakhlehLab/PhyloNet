package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.model;

import java.util.ArrayList;
import java.util.List;

public class HmmStateUtilities {

    public static List<Integer> getGeneTrees(int[] hmmStateSequence, HmmParameters params) {
        //Decomposing HMM States to geneTrees.
        List<Integer> geneTrees = new ArrayList<>();
        for (int i : hmmStateSequence) {
            geneTrees.add(i / params.getNumberOfSpeciesTrees());
        }
        return geneTrees;
    }

    public static List<Integer> getSpeciesTrees(int[] hmmStateSequence, HmmParameters params) {
        //Decomposing HMM States to speciesTrees.
        List<Integer> speciesTrees = new ArrayList<>();
        for (int i : hmmStateSequence) {
            speciesTrees.add(i % params.getNumberOfSpeciesTrees());
        }
        return speciesTrees;
    }


}
