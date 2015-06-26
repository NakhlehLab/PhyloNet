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

package edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;

import java.util.HashMap;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 3/15/12
 * Time: 9:16 AM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkFactoryFromRNNetwork
{
    // kliu - same as below, except using interface implemented by
    // Matt's network representation class 
    public <T> BniNetwork<T> makeNetwork(Network astNetwork)
    {
        return astNetwork.execute(new NetworkAlgo<BniNetwork<T>, Object, RuntimeException>() {

            public BniNetwork<T> forNetworkEmpty(NetworkEmpty networkEmpty, Object o) throws RuntimeException {

                return new BniNetwork<T>();
            }

            public BniNetwork<T> forNetworkNonEmpty(NetworkNonEmpty networkNonEmpty, Object o) throws RuntimeException {

                return makeNetwork(networkNonEmpty);
            }
        }, null);
    }

    // kliu - the code to change Matt's network representation
    // into Cuong's network representation
    public <T> BniNetwork<T> makeNetwork(NetworkNonEmpty astNetwork)
    {
        BniNetwork<T> result = new BniNetwork<T>();
        Map<Integer, BniNetNode<T>> hybridNodeIndexToNetNode = new HashMap<Integer, BniNetNode<T>>();

        BniNetNode root = new BniNetNode();
        populateChild(root, astNetwork.PrincipleInfo, null);
        makeNetworkHelp(root, astNetwork.PrincipleDescendants, hybridNodeIndexToNetNode);

        return new BniNetwork<T>(root);



    }

    private <T> void makeNetworkHelp(BniNetNode<T> parent, DescendantList principleDescendants, Map<Integer, BniNetNode<T>> hybridNodeIndexToNetNode)
    {
        for(Subtree st : principleDescendants.Subtrees)
        {
            BniNetNode node = createOrGetNetNode(st.NetworkInfo, hybridNodeIndexToNetNode, parent);
            makeNetworkHelp(node, st.Descendants, hybridNodeIndexToNetNode);
        }
    }

    // kliu - for all the cruft, I think the two network
    // representations are the same??
    // in-memory node representation with list of references to child only??
    // not interested enough to look at code to confirm
    private <T> BniNetNode createOrGetNetNode(final NetworkInfo principleInfo, final Map<Integer, BniNetNode<T>> hybridNodeIndexToNetNode, final BniNetNode parent)
    {
        BniNetNode<T> result = principleInfo.HybridNodeQualifier.execute(new HybridNodeQualifierAlgo<BniNetNode, Object, RuntimeException>() {
            public BniNetNode forHybridNodeQualifierEmpty(HybridNodeQualifierEmpty hybridNodeQualifierEmpty, Object o) throws RuntimeException {

                BniNetNode<T> newNode = new BniNetNode<T>();
                return newNode;

            }

            public BniNetNode forHybridNodeQualifierNonEmpty(HybridNodeQualifierNonEmpty hybridNodeQualifierNonEmpty, Object o) throws RuntimeException {

                return forHybridNode(principleInfo, new Integer(hybridNodeQualifierNonEmpty.HybridNodeIndex.Content));
            }

            public BniNetNode forHybridNodeQualifierWithType(HybridNodeQualifierWithType hybridNodeQualifierWithType, Object o) throws RuntimeException {

                return forHybridNode(principleInfo, new Integer(hybridNodeQualifierWithType.HybridNodeIndex.Content));
            }

            private BniNetNode<T> forHybridNode(NetworkInfo principleInfo, Integer hybidNodeIndex)
            {
                if(hybridNodeIndexToNetNode.containsKey(hybidNodeIndex))
                {
                    return hybridNodeIndexToNetNode.get(hybidNodeIndex);
                }
                else
                {
                    BniNetNode<T> hybridNode = new BniNetNode<T>();
                    hybridNodeIndexToNetNode.put(hybidNodeIndex, hybridNode);
                    return hybridNode;
                }
            }

        }, null);

        parent.adoptChild(result, NetNode.NO_DISTANCE);
        populateChild(result, principleInfo, parent);
        return result;
    }

    private <T> void populateChild(final BniNetNode<T> child, NetworkInfo childInfo, final BniNetNode parent)
    {
        childInfo.NodeLabel.execute(new NodeLabelAlgo<Object, Object, RuntimeException>() {
            public Object forNodeLabelNonEmpty(NodeLabelNonEmpty nodeLabelNonEmpty, Object o) throws RuntimeException {
                child.setName(nodeLabelNonEmpty.Label.Content);
                return null;
            }

            public Object forNodeLabelEmpty(NodeLabelEmpty nodeLabelEmpty, Object o) throws RuntimeException {
                return null;
            }
        }, null);

        if(parent != null)
        {
            childInfo.BranchLength.execute(new BranchLengthAlgo<Object, Object, RuntimeException>() {

                public Object forBranchLengthEmpty(BranchLengthEmpty branchLengthEmpty, Object o) throws RuntimeException {
                    return null;  //To change body of implemented methods use File | Settings | File Templates.
                }

                public Object forBranchLengthNonEmpty(BranchLengthNonEmpty branchLengthNonEmpty, Object o) throws RuntimeException {
                    child.setParentDistance(parent, Double.parseDouble(branchLengthNonEmpty.Length.Content));
                    return null;  //To change body of implemented methods use File | Settings | File Templates.
                }
            }, null);

            childInfo.Support.execute(new SupportAlgo<Object, Object, RuntimeException>() {

                public Object forSupportNonEmpty(SupportNonEmpty supportNonEmpty, Object o) throws RuntimeException {
                    child.setParentSupport(parent, Double.parseDouble(supportNonEmpty.SupportValue.Content));
                    return null;  //To change body of implemented methods use File | Settings | File Templates.
                }

                public Object forSupportEmpty(SupportEmpty supportEmpty, Object o) throws RuntimeException {
                    return null;  //To change body of implemented methods use File | Settings | File Templates.
                }

            }, null);

            childInfo.Probability.execute(new ProbabilityAlgo<Object, Object, RuntimeException>() {
                public Object forProbabilityEmpty(ProbabilityEmpty probabilityEmpty, Object o) throws RuntimeException {
                    return null;
                }

                public Object forProbabilityNonEmpty(ProbabilityNonEmpty probabilityNonEmpty, Object o) throws RuntimeException {
                    child.setParentProbability(parent, Double.parseDouble(probabilityNonEmpty.ProbabilityValue.Content) );
                    return null;
                }
            }, null);
        }



    }
}
