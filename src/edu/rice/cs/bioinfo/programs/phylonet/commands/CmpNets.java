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
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetworkMetricNakhleh;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.characterization.NetworkTripartition;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.ExNewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/27/11
 * Time: 4:04 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("cmpnets")
public class CmpNets extends CommandBaseFileOut {

    enum Method
    {
        Tree,Tri,Cluster,Luay;
    }

    private Method _method;

    private NetworkNonEmpty _inputNetwork1;

    private NetworkNonEmpty _inputNetwork2;


    public CmpNets(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                   Proc3<String, Integer, Integer> errorDetected, RichNewickReader<edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams() {
        return 4;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 5;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand()
    {
        boolean noError = true;

        _inputNetwork1 = this.assertAndGetNetwork(0);
        noError = noError && _inputNetwork1 != null;

        _inputNetwork2 = this.assertAndGetNetwork(1);
        noError = noError && _inputNetwork2 != null;

        ParamExtractor mSwitch = new ParamExtractor("m", this.params, errorDetected);
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
                 else if(simplePlusOneValueLower.equals("luay"))
                 {
                     _method = Method.Luay;
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

        if(_inputNetwork1 == null)
        {
            throw new IllegalStateException("_inputNetwork1 should be non-null.");
        }

        if(_inputNetwork2 == null)
        {
            throw new IllegalStateException("_inputNetwork2 should be non-null.");
        }

        String eNewickInputNetwork1 = NetworkTransformer.toENewick(_inputNetwork1);

       // Read the input network.
		ExNewickReader<Object> enr = new ExNewickReader<Object>(new StringReader(eNewickInputNetwork1));

        Network<Object> net1;
        try
        {
		    net1 = enr.readNetwork();
        }
        catch(Exception e)
        {
            throw new RuntimeException(e);
        }


        String eNewickInputNetwork2 = NetworkTransformer.toENewick(_inputNetwork2);

       // Read the input network.
		enr = new ExNewickReader<Object>(new StringReader(eNewickInputNetwork2));

        Network<Object> net2;
        try
        {
		    net2 = enr.readNetwork();
        }
        catch(Exception e)
        {
            throw new RuntimeException(e);
        }

        StringBuilder result = new StringBuilder();
        List<NetworkCluster<Object>> clusters1 = null, clusters2 = null;

       if(_method == Method.Cluster)
        {
            if (clusters1 == null) {
				clusters1 = new LinkedList<NetworkCluster<Object>>();
				for (NetworkCluster<Object> nc : Networks.getSoftwiredClusters(net1)) {
					clusters1.add(nc);
				}
			}
			if (clusters2 == null) {
				clusters2 = new LinkedList<NetworkCluster<Object>>();
				for (NetworkCluster<Object> nc : Networks.getSoftwiredClusters(net2)) {
					clusters2.add(nc);
				}
			}

			double dist[] = Networks.computeSoftwiredClusterDistance(clusters1, clusters2);
			result.append("\n" + "The cluster-based distance between two networks: " + dist[0] + " " + dist[1] + " " + dist[2]);
		}
        else if(_method == Method.Tri)
        {
            List<NetworkTripartition<Object>> partitions1 = new LinkedList<NetworkTripartition<Object>>();
			List<NetworkTripartition<Object>> partitions2 = new LinkedList<NetworkTripartition<Object>>();

			for (NetworkTripartition<Object> ntp : Networks.getTripartitions(net1)) {
				partitions1.add(ntp);
			}
			for (NetworkTripartition<Object> ntp : Networks.getTripartitions(net2)) {
				partitions2.add(ntp);
			}

			double dist[] = Networks.computeTripartitionDistance(partitions1, partitions2);
			result.append("\n" + "The tripartition-based distance between two networks: " + dist[0] + " " + dist[1] + " " + dist[2]);
        }
        else if(_method == Method.Tree)
       {
          LinkedList<Tree> trees1 =new LinkedList<Tree>();
          LinkedList<Tree> trees2 =new LinkedList<Tree>();


           for (NetworkTree<Object> nt : Networks.getTrees(net1)) {
               trees1.add(nt.makeTree());
           }

           for (NetworkTree<Object> nt : Networks.getTrees(net2)) {
               trees2.add(nt.makeTree());
           }

           double dist[] = Networks.computeTreeDistance(trees1, trees2);
		   result.append("\nThe tree-based distance between two networks: " + dist[0] + " " + dist[1] + " " + dist[2]);
       }
       else if(_method == Method.Luay)
       {
           NetworkMetricNakhleh metric = new NetworkMetricNakhleh();
           double distance = metric.computeDistanceBetweenTwoNetworks(net1, net2);
           result.append("\nThe Luay's distance between two networks: " + distance);
       }

        return result.toString();

    }
}
