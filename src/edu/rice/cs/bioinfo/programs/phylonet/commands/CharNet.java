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

import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTripartition;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.ExNewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/27/11
 * Time: 4:04 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("charnet")
public class CharNet extends CommandBaseFileOut {

    enum Method
    {
        Tree,Tri,Cluster;
    }

    private Method _method;

    private NetworkNonEmpty _inputNetwork;

    public CharNet(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                   Proc3<String, Integer, Integer> errorDetected, RichNewickReader<edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    public int getMinNumParams() {
        return 3;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public int getMaxNumParams() {
        return 4;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand()
    {
        boolean noError = true;

        _inputNetwork = this.assertAndGetNetwork(0);
        noError = noError && _inputNetwork != null;

        ParamExtractor mSwitch = new ParamExtractor("m", this.params, this.errorDetected);
        boolean sawDashM = mSwitch.ContainsSwitch;

        if(sawDashM)
        {
             if(mSwitch.PostSwitchValue != null)
             {
                 String simplePlusOneValueLower = mSwitch.PostSwitchValue.toLowerCase();

                 if(simplePlusOneValueLower.equals("tree"))
                 {
                     _method = Method.Tree;
                 }
                 else if(simplePlusOneValueLower.equals("tri"))
                 {
                     _method = Method.Tri;
                 }
                 else if(simplePlusOneValueLower.equals("cluster"))
                 {
                     _method = Method.Cluster;
                 }
                 else
                 {
                     this.errorDetected.execute("Unknown method '" + mSwitch.PostSwitchValue + "'", mSwitch.PostSwitchParam.getLine(), mSwitch.PostSwitchParam.getColumn());
                     noError = false;
                 }
             }
            else
             {
                 this.errorDetected.execute("Expected subsequent method value 'tree', 'tri' or 'cluster'.", mSwitch.SwitchParam.getLine(), mSwitch.SwitchParam.getColumn());
                    noError = false;

             }
        }
        else
        {
             this.errorDetected.execute("Expected mandatory parameter '-m [tree|tri|cluster]'", this.getDefiningSyntaxCommand().getLine(),
                                                                                                   this.getDefiningSyntaxCommand().getColumn());
             noError = false;
        }

        this.checkAndSetOutFile(mSwitch);

          return noError;
    }

     @Override
    protected String produceResult()
    {
        if(_method == null)
        {
            throw new IllegalStateException("_method should be non-null.");
        }

        if(_inputNetwork == null)
        {
            throw new IllegalStateException("_inputNetwork should be non-null.");
        }

        String eNewickInputNetwork = NetworkTransformer.toENewick(_inputNetwork);

       // Read the input network.
		ExNewickReader<String> enr = new ExNewickReader<String>(new StringReader(eNewickInputNetwork));

        Network<String> net;
        try
        {
		    net = enr.readNetwork();
        }
        catch(Exception e)
        {
            throw new RuntimeException(e);
        }

        StringBuilder result = new StringBuilder();
        if(_method == Method.Tree)
        {
            boolean first = true;
            for (NetworkTree<String> nt : Networks.getTrees(net)) {
                String tree = nt.makeTree().toString();

                try
                {
                    String rTree = StringTransformer.toRNewickTree(tree);
				    result.append("\n" + rTree.toString());
                    this.richNewickGenerated(rTree);
                }
                catch(CoordinateParseErrorsException e)
                {
                    throw new RuntimeException(e);
                }
                first = false;
			}
        }
        else if(_method == Method.Tri)
        {

            for (NetworkTripartition<String> ntp : Networks.getTripartitions(net)) {
				for (int i = 0; i < ntp.getTripartitionNode().getIndeg(); i++) {
					result.append("\n" + ntp.toString());
				}
			}
        }
        else if(_method == Method.Cluster)
        {
	        for (NetworkCluster<String> nc : Networks.getSoftwiredClusters(net))
            {
				result.append("\n" + nc.toString());
			}
		}

        return result.toString();

    }
}
