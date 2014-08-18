/**
 * This file is part of PhyloNet-HMM.
 *
 * Copyright © 2013-2014 Kevin Liu, Jingxuan Dai, Kathy Truong, 
 * Ying Song, Michael H. Kohn, and Luay Nakhleh. <http://bioinfo.cs.rice.edu/>
 * 
 * PhyloNet-HMM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * PhyloNet-HMM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* jahmm package - v0.6.1 */

/*
  *  Copyright (c) 2004-2006, Jean-Marc Francois.
 *
 *  This file is part of Jahmm.
 *  Jahmm is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Jahmm is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Jahmm; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

 */

package edu.rice.cs.bioinfo.programs.phmm.src.be.ac.ulg.montefiore.run.jahmm;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

// kliu - these appear to be new
import edu.rice.cs.bioinfo.programs.phmm.src.runHmm.MemoryReport;
import edu.rice.cs.bioinfo.programs.phmm.src.containers.IArrayFactory;
import edu.rice.cs.bioinfo.programs.phmm.src.containers.IntArrayFactory;
import edu.rice.cs.bioinfo.programs.phmm.src.containers.ListOfArrays;


/**
 * This class can be used to compute the most probable state sequence matching
 * a given observation sequence (given an HMM).
 */
public class ViterbiCalculator
{
    /*
     * The psy and delta values, as described in Rabiner and Juand classical
     * papers.
     */
    private double[][] delta;
    // kliu - not sure why switched from original int indexing to
    // Integer indexing?
    private List<ListOfArrays<Integer>> Psy;
    private List<Integer> StateSequence;
    private double lnProbability;


    /**
     * Computes the most likely state sequence matching an observation
     * sequence given an HMM.
     *
     * @param hmm A Hidden Markov Model;
     * @param oseq An observations sequence.
     */
    public <O extends Observation>
    ViterbiCalculator(List<? extends O> oseq, Hmm<O> hmm)
    {
        if (oseq.isEmpty())
            throw new IllegalArgumentException("Invalid empty sequence");

        final int growthSize = 10000;

        IArrayFactory<Integer> arrayFactory = new IntArrayFactory();
        delta = new double[2][hmm.nbStates()];
        Psy = new ArrayList<ListOfArrays<Integer>>();
        StateSequence = new ListOfArrays<Integer>(arrayFactory, growthSize);

        //building Psy array
        System.out.println("Building Psy Array");
        for (int i = 0; i < oseq.size(); i++) {
            // kliu - added Integer type argument to fix compilation error
            Psy.add(new ListOfArrays<Integer>(arrayFactory, hmm.nbStates()));
            if ((i%10000) == 0) {
                System.out.print("\r"+ (((double)i/(double)oseq.size()) * 100.0) + "%");
            }
        }

        System.out.print("\r100% Done!                             ");
        System.out.println("");

        MemoryReport.report();

        for (int i = 0; i < hmm.nbStates(); i++) {
            delta[0][i] = -Math.log(hmm.getPi(i)) -
                    Math.log(hmm.getOpdf(i).probability(oseq.get(0)));
            Psy.get(0).add(i, 0);
        }

        System.out.println("Begin computing Viterbi's algorithm : FORWARD PART ");
        Iterator<? extends O> oseqIterator = oseq.iterator();
        if (oseqIterator.hasNext())
            oseqIterator.next();

        int t = 1;
        while (oseqIterator.hasNext()) {
            O observation = oseqIterator.next();

            for (int i = 0; i < hmm.nbStates(); i++) {

                computeStep(hmm, observation, t, i);

            }
            for (int k = 0; k < delta[0].length; k++) {
                delta[0][k] = delta[1][k];
            }

            if ((t%10000) == 0) {
                System.out.print("\r"+ (((double)t/(double)oseq.size()) * 100.0) + "%");
            }

            t++;

        }

        System.out.print("\r 100% DONE!                                \n");

        lnProbability = Double.MAX_VALUE;

        int index = -1;
        for (int i = 0; i < hmm.nbStates(); i++) {
            double thisProbability = delta[0][i];

            if (lnProbability > thisProbability) {
                lnProbability = thisProbability;
                index = i;
            }
        }
        StateSequence.add(index);

        lnProbability = -lnProbability;


        System.out.println("Begin computing Viterbi's algorithm : BACKWARD PART ");
        int ss = 0;
        for (int t2 = oseq.size() - 2; t2 >= 0; t2--) {
            StateSequence.add(Psy.get(t2+1).get(StateSequence.get(ss)));

            if ((ss%10000) == 0) {
                System.out.print("\r"+ (((double)ss/(double)oseq.size()) * 100.0) + "%");
            }

            ss++;
        }

        System.out.print("\r 100% DONE!                                   \n");


    }


    /*
     * Computes delta and psy[t][j] (t > 0)
     */
    private <O extends Observation> void
    computeStep(Hmm<O> hmm, O o, int t, int j)
    {
        double minDelta = Double.MAX_VALUE;
        int min_psy = 0;

        for (int i = 0; i < hmm.nbStates(); i++) {
            double thisDelta = delta[0][i] - Math.log(hmm.getAij(i, j));


            if (minDelta > thisDelta) {
                minDelta = thisDelta;
                min_psy = i;
            }
        }


        delta[1][j] = minDelta - Math.log(hmm.getOpdf(j).probability(o));
        Psy.get(t).add(min_psy);

    }


    /**
     * Returns the neperian logarithm of the probability of the given
     * observation sequence on the most likely state sequence of the given
     * HMM.
     *
     * @return <code>ln(P[O,S|H])</code> where <code>O</code> is the given
     *         observation sequence, <code>H</code> the given HMM and
     *         <code>S</code> the most likely state sequence of this observation
     *         sequence given this HMM.
     */
    public double lnProbability()
    {
        return lnProbability;
    }


    /**
     * Returns a (clone of) the array containing the computed most likely
     * state sequence.
     *
     * @return The state sequence; the i-th value of the array is the index
     *         of the i-th state of the state sequence.
     */
    public int[] stateSequence()
    {
        int[] statesseq = new int[StateSequence.size()];
        int j = 0;
        for (int i = StateSequence.size() - 1; i >= 0; i--) {
            statesseq[j] = StateSequence.get(i);
            j++;
        }
        return statesseq;
    }

    /**
     * Print state sequence out to a file
     *
     */
    public void printStateSequence(String filename) {
        System.out.println("Begin Printing State Sequence to File.");
        try {
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(filename, true)));
            for (int j = StateSequence.size() -1; j >= 0; j--) {
                out.print(StateSequence.get(j));

                if (j % 1024 == 0) {
                    out.flush();
                }

                if ((j %10000) == 0) {
                    System.out.print("\r"+ (((double)j/(double)StateSequence.size()) * 100.0) + "%");
                }
            }

            System.out.print("\r 100% DONE!                                \n");

            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
