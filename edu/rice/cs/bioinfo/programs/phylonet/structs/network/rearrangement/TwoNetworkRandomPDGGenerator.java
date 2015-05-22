package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by yunyu on 11/7/14.
 */
public class TwoNetworkRandomPDGGenerator extends NetworkNeighbourhoodGenerator {
    private Random _random;

    public TwoNetworkRandomPDGGenerator(Long seed){
        if (seed != null) {
            _random = new Random(seed);
        } else {
            _random = new Random();
        }
    }

    public void mutateNetwork(Network network){};

    public void undo(){};

    public boolean mutateNetwork(Network network1, Network network2){
        if(_printDetails){
            System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
            System.out.println("Select operation: "+ this.getClass().getSimpleName());
            System.out.println("Before operation: ");
            System.out.println("Parent #1: " + network1);
            System.out.println("Parent #2: " + network2);
        }
        List<NetNode> candidateNodes = new ArrayList<>();
        for(Object nodeO: Networks.getAllArticulationNodes(network1)){
            NetNode node = (NetNode)nodeO;
            if(!node.isLeaf() && !node.isRoot()){
                candidateNodes.add(node);
            }
        }
        if(candidateNodes.size()==0){
            if(_printDetails){
                System.out.println("Cannot do recombination!");
                System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
            }
            return false;
        }
        Set<String> allLeaves = new HashSet<>();
        for(Object leaf: network2.getLeaves()){
            allLeaves.add(((NetNode)leaf).getName());
        }
        NetNode targetNode = candidateNodes.get(_random.nextInt(candidateNodes.size()));
        Tuple<NetNode,NetNode> targetEdge = new Tuple<>((NetNode)targetNode.getParents().iterator().next(), targetNode);
        if(_printDetails){
            System.out.println("Select edge: "+printEdge(targetEdge));
        }
        PruneDeleteGraft pdg = new PruneDeleteGraft();
        pdg.setParameters(network2, targetEdge);
        pdg.performOperation();
        if(!isNetworkValid(network2, allLeaves)){
            throw new RuntimeException(network2.toString());
        }
        if(_printDetails){
            System.out.println("After operation: "+network2);
            System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
        }
        return true;
    }
}
