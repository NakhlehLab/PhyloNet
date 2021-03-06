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

import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Collection;

import edu.rice.cs.bioinfo.programs.phmm.src.be.ac.ulg.montefiore.run.distributions.MultiGaussianDistribution;



/**
 * This class represents a multivariate gaussian distribution function.
 */
public class OpdfMultiGaussian
implements Opdf<ObservationVector>
{
    private MultiGaussianDistribution distribution;


    /**
     * Builds a new gaussian probability distribution with zero mean and
     * identity covariance matrix.
     *
     * @param dimension The dimension of the vectors.
     */
    public OpdfMultiGaussian(int dimension)
    {
        distribution = new MultiGaussianDistribution(dimension);
    }


    /**
     * Builds a new gaussian probability distribution with a given mean and
     * covariance matrix.
     *
     * @param mean The distribution's mean.
     * @param covariance The distribution's covariance matrix.
     */
    public OpdfMultiGaussian(double[] mean, double[][] covariance)
    {
        if (covariance.length == 0 || mean.length != covariance.length ||
                covariance.length != covariance[0].length)
            throw new IllegalArgumentException();

        distribution = new MultiGaussianDistribution(mean, covariance);
    }


    /**
     * Returns (a copy of) this distribution's mean vector.
     *
     * @return The mean vector.
     */
    public double[] mean()
    {
        return distribution.mean();
    }


    /**
     * Returns (a copy of) this distribution's covariance matrix.
     *
     * @return The covariance matrix.
     */
    public double[][] covariance()
    {
        return distribution.covariance();
    }


    /**
     * Returns the dimension of the vectors handled by this distribution.
     *
     * @return The dimension of the vectors handled by this distribution.
     */
    public int dimension()
    {
        return distribution.dimension();
    }


    public double probability(ObservationVector o)
    {
        if (o.dimension() != distribution.dimension())
            throw new IllegalArgumentException("Vector has a wrong " +
            "dimension");

        return distribution.probability(o.value);
    }


    public ObservationVector generate()
    {
        return new ObservationVector(distribution.generate());
    }


    public void fit(ObservationVector... oa)
    {
        fit(Arrays.asList(oa));
    }


    public void fit(Collection<? extends ObservationVector> co)
    {
        if (co.isEmpty())
            throw new IllegalArgumentException("Empty observation set");

        double[] weights = new double[co.size()];
        Arrays.fill(weights, 1. / co.size());

        fit(co, weights);
    }


    public void fit(ObservationVector[] o, double[] weights)
    {
        fit(Arrays.asList(o), weights);
    }


    public void fit(Collection<? extends ObservationVector> co,
            double[] weights)
    {
        if (co.isEmpty() || co.size() != weights.length)
            throw new IllegalArgumentException();

        // Compute mean
        double[] mean = new double[dimension()];
        for (int r = 0; r < dimension(); r++) {
            int i = 0;

            for (ObservationVector o : co)
                mean[r] += o.value[r] * weights[i++];
        }

        // Compute covariance
        double[][] covariance = new double[dimension()][dimension()];
        int i = 0;
        for (ObservationVector o : co) {
            double[] obs = o.value;
            double[] omm = new double[obs.length];

            for (int j = 0; j < obs.length; j++)
                omm[j] = obs[j] - mean[j];

            for (int r = 0; r < dimension(); r++)
                for (int c = 0; c < dimension(); c++)
                    covariance[r][c] += omm[r] * omm[c] * weights[i];

            i++;
        }

        distribution = new MultiGaussianDistribution(mean, covariance);
    }


    public OpdfMultiGaussian clone()
    {
        try {
            return (OpdfMultiGaussian) super.clone();
        } catch(CloneNotSupportedException e) {
            throw new AssertionError(e);
        }
    }


    public String toString()
    {
        return toString(NumberFormat.getInstance());
    }


    public String toString(NumberFormat numberFormat)
    {
        String s = "Multi-variate Gaussian distribution --- Mean: [ ";
        double[] mean = distribution.mean();

        for (int i = 0; i < mean.length; i++)
            s += numberFormat.format(mean[i]) + " ";

        return s + "]";
    }


    private static final long serialVersionUID = 1L;
}
