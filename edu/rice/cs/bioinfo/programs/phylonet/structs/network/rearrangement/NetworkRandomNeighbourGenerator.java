package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;


import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 6/7/12
 * Time: 9:26 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkRandomNeighbourGenerator extends NetworkNeighbourhoodGenerator{
    private NetworkNeighbourhoodGenerator _topologyMutate;
    private NetworkNeighbourhoodGenerator _parameterMutate;
    private double _operationProbability;
    private Random _random;
    private int _previousOperation;


    public void undo(){
        if(_previousOperation == 0){
            _topologyMutate.undo();
        }
        else{
            _parameterMutate.undo();
        }
    }

    public void setParameterWindowSize(double window){
        ((NetworkRandomParameterNeighbourGenerator)_parameterMutate).setWindowSize(window);
    }

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

    public void mutateNetworkParameter(Network network){
        _parameterMutate.mutateNetwork(network);
        _previousOperation = 1;
    }

    public int getOperationType(){
        return _previousOperation;
    }

    public void setOperationType(int type){
        _previousOperation = type;
    }

}