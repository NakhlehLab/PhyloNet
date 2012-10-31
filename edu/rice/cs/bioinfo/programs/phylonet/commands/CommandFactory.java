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

package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.util.ArrayList;
import java.util.Map;
import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 5:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class CommandFactory {



    public static Command make(SyntaxCommand directive, Map<String,NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected, Random rand)
    {
        String lowerCommandName = directive.getName().toLowerCase();

        ArrayList<Parameter> params = new ArrayList<Parameter>();
        for(Parameter p : directive.getParameters())
        {
            params.add(p);
        }

        if(lowerCommandName.equals("symmetricdifference") || lowerCommandName.equals("rf"))
        {
            return new SymmetricDifference(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("lca"))
        {
            return new LCA(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("mast"))
        {
            return new MAST(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("riatahgt"))
        {
            return new RIATAHGT(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("deepcoalcount_tree"))
        {
            return new DeepCoalCount(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("deepcoalcount_network"))
        {
            return new CountXLInNetwork(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("processgt"))
        {
            return new ProcessGT(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("charnet"))
        {
            return new CharNet(directive, params, sourceIdentToNetwork, errorDetected);
        }
         else if(lowerCommandName.equals("cmpnets"))
        {
            return new CmpNets(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("simgtinnetwork"))
        {
            return new SimGTinNetwork(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("countcoal"))
        {
           return new CountCoal(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("genst"))
        {
            return new GenST(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("gencplex"))
        {
            return new GenCPLEX(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("infer_st_mdc"))
        {
            return new InferST_MDC(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("infer_st_mdc_time"))
        {
            return new InferST_MDC_Time(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("infer_st_mdc_ur"))
        {
            return new InferST_MDC_UR(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("infer_st_glass"))
        {
            return new InfterST_MDC_GLASS(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("infer_st_dv"))
        {
            return new InferST_DV(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("infer_st_mc"))
        {
            return new InferST_MC(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("infer_st_bootstrap"))
        {
            return new InferST_Bootstrap(directive, params, sourceIdentToNetwork, errorDetected, rand);
        }
        else if(lowerCommandName.equals("infer_st_mdc_ilp"))
        {
            return new InferST_MDC_ILP(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("nninterchange"))
        {
            return new NearestNeighborInterchange(directive, params, sourceIdentToNetwork, errorDetected, rand);
        }
        else if(lowerCommandName.equals("spregraft"))
        {
            return new SubtreePruneAndRegraft(directive, params, sourceIdentToNetwork, errorDetected, rand);
        }
        else if(lowerCommandName.equals("tbreconnection"))
        {
            return new TreeBisectionAndReconnection(directive, params, sourceIdentToNetwork, errorDetected, rand);
        }
        else if(lowerCommandName.equals("nexus_out"))
        {
            return new NexusOut(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("calgtprob"))
        {
            return new CalGTProbInNetwork(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("optimizecontinuousnetworkmodelparameters"))
        {
            return new OptimizeContinuousNetworkModelParameters(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else if(lowerCommandName.equals("assignbranchlengthsmaxgtprob"))
        {
            return new AssignBranchLengthsMaxGTProb(directive, params, sourceIdentToNetwork, errorDetected);
        }
        else
        {
             throw new IllegalArgumentException(String.format("Unknown command name '%s'.", directive.getName()));
        }

    }
}
