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

package edu.rice.cs.bioinfo.programs.phylonet;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.*;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.BlockContents;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.RichNewickAssignment;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.ContainsHybridNode;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RootageQualifierAlgo;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RootageQualifierEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RootageQualifierNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.csa.ASTContextAnalyser;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.csa.CSAError;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.math.BigDecimal;
import java.util.HashSet;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 5:07 PM
 * To change this template use File | Settings | File Templates.
 */
class ContextSensitiveAnalyser {

    public static void analyseNetworks(Map<String,NetworkNonEmpty> sourceIdentToNetwork, BlockContents blockContents, BigDecimal hybridSumTolerance, final Proc3<String,Integer,Integer> errorDetected)
    {
       for(String rNewickStringIdent : sourceIdentToNetwork.keySet())
       {
           final NetworkNonEmpty network = sourceIdentToNetwork.get(rNewickStringIdent);

           for(CSAError error : new ASTContextAnalyser().analyse(network, hybridSumTolerance))
           {
               errorDetected.execute(error.Message, -1, -1);
           }

           final RichNewickAssignment assigment = blockContents.getRichNewickAssigment(rNewickStringIdent);

           if(assigment.isDefinedByNetworksBlock())
           {
               network.RootageQualifier.execute(new RootageQualifierAlgo<Object, Object, RuntimeException>() {
                   public Object forEmptyQualifier(RootageQualifierEmpty rootageQualifierEmpty, Object o) throws RuntimeException {
                       return null;  //To change body of implemented methods use File | Settings | File Templates.
                   }

                   public Object forNonEmptyQualifier(RootageQualifierNonEmpty rootageQualifierNonEmpty, Object o) throws RuntimeException {

                       if(!rootageQualifierNonEmpty.isRooted())
                       {
                           errorDetected.execute("Networks may not be unrooted.", assigment.getRichNewickStringLine(), assigment.getRichNewickStringColumn());
                       }
                       return null;
                   }
               }, null);
           }


       }
    }

    public static boolean checkforDuplicateAssignmentIdentifiers(Blocks blocks, final Proc3<String, Integer, Integer> errorDetected)
    {
         final HashSet<String> seenIdents = new HashSet<String>();

        Boolean observedDup = false;
        for(Block b: blocks.Contents)
        {
            observedDup = observedDup || b.execute(new BlockAlgo<Boolean, Object, RuntimeException>() {

                public Boolean forTreesBlock(TreesBlockBody treesBlockBody, Object o) throws RuntimeException {

                    boolean seenDup = false;
                    for(TreeAssignment ta : treesBlockBody.getAssignments())
                    {
                        String ident = ta.Assignment.LHSIdentifier.Content;
                        if(seenIdents.contains(ident))
                        {
                            seenDup = true;
                            errorDetected.execute("Duplicate identifier '" + ident + "'.", ta.Assignment.LHSIdentifier.Line, ta.Assignment.LHSIdentifier.Col);
                        }
                        else
                        {
                            seenIdents.add(ident);
                        }
                    }

                    return seenDup;  //To change body of implemented methods use File | Settings | File Templates.
                }

                public Boolean forNetworksBlock(NetworksBlockBody networksBlockBody, Object o) throws RuntimeException {

                    boolean seenDup = false;
                    for(NetworkAssignment na : networksBlockBody.getAssignments())
                    {
                        String ident = na.Assignment.LHSIdentifier.Content;
                        if(seenIdents.contains(ident))
                        {
                             seenDup = true;
                             errorDetected.execute("Duplicate identifier '" + ident + "'.", na.Assignment.LHSIdentifier.Line, na.Assignment.LHSIdentifier.Col);
                        }
                        else
                        {
                            seenIdents.add(ident);
                        }
                    }

                    return seenDup;  //To change body of implemented methods use File | Settings | File Templates.
                }

                public Boolean forPhylonetBlockBody(PhyloNetBlockBody phyloNetBlockBody, Object o) throws RuntimeException {
                    return false;  //To change body of implemented methods use File | Settings | File Templates.
                }

                public Boolean forDataBlock(DataBlockBody dataBlock, Object o) throws RuntimeException {
                  return false;
                }

            }, null);
        }

        return observedDup;
    }

    public static void checkforHybridNodesInTrees(Map<String,NetworkNonEmpty> sourceIdentToNetwork, BlockContents blockContents, Proc3<String, Integer, Integer> errorDetected) {

        for(String assignmentIdent : blockContents.getRickNewickAssignmentIdentifiers())
        {
            RichNewickAssignment assigment = blockContents.getRichNewickAssigment(assignmentIdent);
            NetworkNonEmpty network = sourceIdentToNetwork.get(assignmentIdent);

            if(assigment.isDefinedByTreesBlock() && network.execute(ContainsHybridNode.Singleton, null))
            {
                errorDetected.execute("Trees may not contain hybrid nodes.", assigment.getRichNewickStringLine(), assigment.getRichNewickStringColumn());
            }
        }

    }
}
