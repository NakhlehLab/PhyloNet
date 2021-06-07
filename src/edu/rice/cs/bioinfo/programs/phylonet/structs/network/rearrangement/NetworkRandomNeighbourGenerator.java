package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;


import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.util.Random;

/**
 * Created by Yun Yu
 *
 * This class is a subclass of NetworkNeighbourhoodGenerator.
 * It generates a random neighbor of a given network (topology-wise or parameter-wise)
 */
public class NetworkRandomNeighbourGenerator extends NetworkNeighbourhoodGenerator{
    private NetworkNeighbourhoodGenerator _topologyMutate;
    private NetworkNeighbourhoodGenerator _parameterMutate;
    private double _operationProbability;
    private Random _random;
    private int _previousOperation;


    /**
     * This function is to undo the last rearrangement move
     */
    public void undo(){
        if(_previousOperation == 0){
            _topologyMutate.undo();
        }
        else{
            _parameterMutate.undo();
        }
    }


    /**
     * This function is to set the window size for changing parameters
     */
    public void setParameterWindowSize(double window){
        ((NetworkRandomParameterNeighbourGenerator)_parameterMutate).setWindowSize(window);
    }



    /**
     * Constructor of this class
     *
     * @param topologyMutate    the class for mutating topology of the network
     * @param weight1           the weight of choosing topology change
     * @param parameterMutate   the class for mutating parameters of the network
     * @param weight2           the weight of choosing parameter change
     * @param seed              the seed for randomness
     */
    public NetworkRandomNeighbourGenerator(NetworkNeighbourhoodGenerator topologyMutate, double weight1, NetworkNeighbourhoodGenerator parameterMutate, double weight2, Long seed)
    {
        _topologyMutate = topologyMutate;
        _parameterMutate = parameterMutate;
        _operationProbability = weight1/(weight1+weight2);
        if(seed!=null) {
            _random = new Random(seed);
        }else{
            _random = new Random();
        }
    }


    /**
     * This function is to randomly mutate the current network (topology or parameters)
     */
    public void mutateNetwork(Network network){
        if(_random.nextDouble()<_operationProbability){
            _topologyMutate.mutateNetwork(network);
            _previousOperation = 0;
        }
        else{
            _parameterMutate.mutateNetwork(network);
            _previousOperation = 1;
        }
    }


    /**
     * This function is to randomly mutate the parameters of the current network
     */
    public void mutateNetworkParameter(Network network){
        _parameterMutate.mutateNetwork(network);
        _previousOperation = 1;
    }


    /**
     * This function is to get the type of last mutation (topology or parameters)
     */
    public int getOperationType(){
        return _previousOperation;
    }


    /**
     * This function is to set the window size for changing parameters
     */
    public void resetList(){
        _topologyMutate.resetList();

    }

    public int getOperationID(){
        if(_previousOperation == 0){
            return _topologyMutate.getOperationID();

        }
        else{
            return 6;//length-1
        }
    }

}