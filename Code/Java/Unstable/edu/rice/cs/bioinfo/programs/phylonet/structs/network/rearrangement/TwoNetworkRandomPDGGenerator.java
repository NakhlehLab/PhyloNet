package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by Yun Yu
 *
 * This class is a subclass of NetworkNeighbourhoodGenerator.
 * It generates a random network by combining two networks
 */
public class TwoNetworkRandomPDGGenerator extends NetworkNeighbourhoodGenerator {
    private Random _random;

    /**
     * Constructor of this class
     */
    public TwoNetworkRandomPDGGenerator(Long seed){
        if (seed != null) {
            _random = new Random(seed);
        } else {
            _random = new Random();
        }
    }

    /**
     * This function does nothing
     */
    public void mutateNetwork(Network network){};


    /**
     * This function does nothing
     */
    public void undo(){};


    /**
     * This function is for combining two networks <code>network1</code> and <code>network2</code>
     * Note that the resulting network is <code>network2</code>
     */
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
            System.out.println("Select edge: " + printEdge(targetEdge));
        }
        PruneDeleteGraft pdg = new PruneDeleteGraft();
        pdg.setParameters(network2, targetEdge);
        pdg.performOperation();
        if(!isNetworkValid(network2, allLeaves)){
            throw new RuntimeException(network2.toString());
        }
        if(_printDetails){
            System.out.println("After operation: " + network2);
            System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
        }
        return true;
    }



    public void resetList(){}

    public int getOperationID(){
        return -1;
    }
}
