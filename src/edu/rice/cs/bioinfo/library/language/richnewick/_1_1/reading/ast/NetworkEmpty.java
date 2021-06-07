/*
 * Copyright (c) 2013 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/21/13
 * Time: 5:34 PM
 * To change this template use File | Settings | File Templates.
 */
public final class NetworkEmpty implements Network {

    public static final NetworkEmpty Singleton = new NetworkEmpty();

    private NetworkEmpty()
    {

    }

    public <R, E extends Exception> R execute(NetworkAlgo<R, E> algo) throws E {
        return algo.forNetworkEmpty(this);
    }

    public <R, T, E extends Exception> R execute(edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.NetworkAlgo<R, T, E> algo, T input) throws E {
        return algo.forNetworkEmpty(edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.NetworkEmpty.Singleton, input);
    }
}
