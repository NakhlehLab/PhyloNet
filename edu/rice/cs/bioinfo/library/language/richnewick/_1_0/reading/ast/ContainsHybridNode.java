/*
 * Copyright (c) 2012 Rice University.
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

package edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/17/11
 * Time: 6:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class ContainsHybridNode implements NetworkAlgo<Boolean, Object, RuntimeException>
{
    public static final ContainsHybridNode Singleton = new ContainsHybridNode();

    private ContainsHybridNode()
    {

    }

    public Boolean forNetworkEmpty(NetworkEmpty network, Object input) throws RuntimeException {
       return false;
    }

    public Boolean forNetworkNonEmpty(NetworkNonEmpty network, Object input) throws RuntimeException {

        if(containsHybridNode(network.PrincipleInfo))
        {
            return true;
        }
        else
        {
            return containsHybridNode(network.PrincipleDescendants.Subtrees);
        }


    }

    private Boolean containsHybridNode(NetworkInfo principleInfo) {

        return principleInfo.HybridNodeQualifier.execute(new HybridNodeQualifierAlgo<Boolean, Object, RuntimeException>() {
            public Boolean forHybridNodeQualifierEmpty(HybridNodeQualifierEmpty qualifier, Object input) throws RuntimeException {

                return false;
            }

            public Boolean forHybridNodeQualifierNonEmpty(HybridNodeQualifierNonEmpty qualifier, Object input) throws RuntimeException {

                return true;
            }

            public Boolean forHybridNodeQualifierWithType(HybridNodeQualifierWithType qualifier, Object input) throws RuntimeException {

                return true;
            }
        }, null);

    }

    private Boolean containsHybridNode(Iterable<Subtree> subtrees) {

        for(Subtree st : subtrees)
        {
            if(containsHybridNode(st.NetworkInfo))
            {
                return true;
            }

            if(containsHybridNode(st.Descendants.Subtrees))
            {
                return true;
            }
        }

        return false;
    }
}
