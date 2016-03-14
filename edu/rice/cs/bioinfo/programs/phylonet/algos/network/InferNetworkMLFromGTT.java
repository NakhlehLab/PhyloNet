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

package edu.rice.cs.bioinfo.programs.phylonet.algos.network;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_DP;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkRandomParameterNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NonUltrametricNetworkRandomParameterNeighbourGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class InferNetworkMLFromGTT extends InferNetworkMLFromGT {

    protected void findSingleAlleleSpeciesSet(Network speciesNetwork, Map<String,List<String>> species2alleles, Set<String> singleAlleleSpecies){
        for(Object node: speciesNetwork.getLeaves()){
            String species = ((NetNode)node).getName();
            if(species2alleles == null || species2alleles.get(species).size() == 1)
                singleAlleleSpecies.add(species);
        }
    }


    protected NetworkRandomParameterNeighbourGenerator getNetworkRandomParameterNeighbourGenerator(List dataForNetworkInference, Map<String,String> allele2species, Set<String> singleAlleleSpecies){
        return new NonUltrametricNetworkRandomParameterNeighbourGenerator(singleAlleleSpecies);
    }

}
