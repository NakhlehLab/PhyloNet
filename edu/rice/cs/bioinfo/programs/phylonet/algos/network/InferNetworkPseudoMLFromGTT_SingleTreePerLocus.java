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

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.util.List;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 2/11/13
 * Time: 11:40 AM
 * To change this template use File | Settings | File Templates.
 */
public class InferNetworkPseudoMLFromGTT_SingleTreePerLocus extends InferNetworkPseudoMLFromGTT {
    public InferNetworkPseudoMLFromGTT_SingleTreePerLocus(){
        _fullLikelihoodCalculator = new NetworkLikelihoodFromGTT_SingleTreePerLocus();
        _likelihoodCalculator = new NetworkPseudoLikelihoodFromGTT_SingleTreePerLocus();
        if(_batchSize!=0){
            ((NetworkPseudoLikelihoodFromGTT_SingleTreePerLocus)_likelihoodCalculator).setBatchSize(_batchSize);
        }
    }

    protected void summarizeData(List originalData, Map<String,String> allele2species, List dataForStartingNetwork, List dataForInferNetwork, List dataCorrespondences){
        NetworkPseudoLikelihoodFromGTT_SingleTreePerLocus likelihoodComputer = new NetworkPseudoLikelihoodFromGTT_SingleTreePerLocus();
        likelihoodComputer.summarizeData(originalData, allele2species, dataForStartingNetwork, dataForInferNetwork, dataCorrespondences);
    }


}
