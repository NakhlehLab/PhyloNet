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

package be.ac.ulg.montefiore.run.jahmm.toolbox;

import java.util.ArrayList;
import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.*;


/**
 * Generates sequences of markovian observations given a HMM.
 */
public class MarkovGenerator<O extends Observation>
{
    private final Hmm<O> hmm;
    private int stateNb;


    /**
     * Initializes a Markovian generator.
     *
     * @param hmm An Hidden Markov Model that perfectly models the sequences
     *            generated by this object.
     */
    public MarkovGenerator(Hmm<O> hmm)
    {
        if (hmm == null)
            throw new IllegalArgumentException("Invalid null HMM");

        this.hmm = hmm;
        newSequence();
    }


    /**
     * Generates a new (pseudo) random observation.
     *
     * @return The generated observation.
     */
    public O observation()
    {
        O o = hmm.getOpdf(stateNb).generate();
        double rand = Math.random();

        for (int j = 0; j < hmm.nbStates()-1; j++)
            if ((rand -= hmm.getAij(stateNb, j)) < 0) {
                stateNb = j;
                return o;
            }

        stateNb = hmm.nbStates() - 1;
        return o;
    }


    /**
     * Generates a new (pseudo) random observation sequence and start
     * a new one.
     *
     * @param length The length of the sequence.
     * @return An observation sequence.
     */
    public List<O> observationSequence(int length)
    {
        if (length <= 0)
            throw new IllegalArgumentException("Positive length required");

        ArrayList<O> sequence = new ArrayList<O>();
        while (length-- > 0)
            sequence.add(observation());
        newSequence();

        return sequence;
    }


    /**
     * Finds a new state according to the initial (pi) probabilities of each
     * state.
     */
    public void newSequence()
    {
        double rand = Math.random(), current = 0.;

        for (int i = 0; i < hmm.nbStates() - 1; i++) {
            current += hmm.getPi(i);

            if (current > rand) {
                stateNb = i;
                return;
            }
        }

        stateNb = hmm.nbStates() - 1;
    }


    /**
     * Returns the state number of the current state.
     *
     * @return A state number.
     */
    public int stateNb()
    {
        return stateNb;
    }
}
