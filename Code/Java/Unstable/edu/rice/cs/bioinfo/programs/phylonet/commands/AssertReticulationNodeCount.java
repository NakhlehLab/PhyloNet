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
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdent;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentList;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc1;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/29/12
 * Time: 2:25 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("AssertReticulationNodeCount")
public class AssertReticulationNodeCount extends CommandBase
{

    private LinkedList<NetworkNonEmpty> _nets;

    private int _expectedReticulationNodeCount = -1;

    public AssertReticulationNodeCount(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                                       Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams() {
        return 2;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 2;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand() {

         boolean noError = true;

        ParameterIdentList netList = this.assertParameterIdentList(0);
        noError = noError && netList != null;
        _nets = this.assertNetworksExist(netList);
        noError = noError && _nets != null;

        ParameterIdent expectedRetNodeCount = this.assertParameterIdent(1);
        noError = noError && expectedRetNodeCount != null;

        if(expectedRetNodeCount != null)
        {
            try
            {
                _expectedReticulationNodeCount = Integer.parseInt(expectedRetNodeCount.Content);
            }
            catch(NumberFormatException e)
            {
                this.errorDetected.execute("Expected number, found " + expectedRetNodeCount.Content, expectedRetNodeCount.getLine(), expectedRetNodeCount.getColumn());
                noError = false;
            }
        }

        return noError;


    }

    @Override
    protected void executeCommandHelp(Proc1<String> displayResult) throws IOException {

        for(NetworkNonEmpty net : _nets)
        {

            Network<Object> pnet =  new NetworkFactoryFromRNNetwork().makeNetwork(net);

            int foundReticulationNodeCount = 0;
            for(NetNode child : pnet.bfs())
            {
                if(child.isRoot())
                    continue;

                if(child.getParentCount() > 1)
                {
                    foundReticulationNodeCount++;
                }
            }

            if(foundReticulationNodeCount != _expectedReticulationNodeCount)
                throw new RuntimeException("Expected " + _expectedReticulationNodeCount + " reticulation nodes but found " + foundReticulationNodeCount);
        }
    }
}
