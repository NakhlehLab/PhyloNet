package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Ref;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by Yun Yu
 *
 * This class is a subclass of NetworkNeighbourhoodGenerator.
 * It generates a random neighbor of a given network by changing its parameters
 */
public abstract class NetworkRandomParameterNeighbourGenerator extends NetworkNeighbourhoodGenerator {
    protected EdgeParameterChange _lengthChanger;
    protected EdgeParameterChange _inheriProbChanger = new EdgeInheritanceProbabilityChange();
    protected Func1<Integer, Integer> _numParametersToChange;
    protected List<Tuple<Tuple<Tuple<NetNode, NetNode>, Boolean>, Double>> _edgeChanged = new ArrayList<>();
    private Random _random;


    /**
     * Constructor of this class
     *
     * @param numParametersToChange      the function that calculates the number of parameters to change in one move
     * @param seed      the seed for controlling randomness
     */
    public NetworkRandomParameterNeighbourGenerator(Func1<Integer, Integer> numParametersToChange, Long seed) {
        _numParametersToChange = numParametersToChange;
        if (seed != null) {
            _random = new Random(seed);
        } else {
            _random = new Random();
        }
    }


    /**
     * Constructor of this class
     *
     * @param numParametersToChange      the function that calculates the number of parameters to change in one move
     */
    public NetworkRandomParameterNeighbourGenerator(Func1<Integer, Integer> numParametersToChange) {
        this(numParametersToChange, null);
    }


    /**
     * Constructor of this class
     * Note that here the default is set to change 2 parameters at a time
     */
    public NetworkRandomParameterNeighbourGenerator() {
        this(new Func1<Integer, Integer>() {
            @Override
            public Integer execute(Integer input) {
                return 2;
            }
        });
    }


    /**
     * This function is to set the function that calculates the number of parameters to change in one move
     */
    public void setParametersToChange(Func1<Integer, Integer> numParametersToChange){
        _numParametersToChange = numParametersToChange;
    }


    /**
     * This function is to set the size of the window when changing the parameters
     */
    public void setWindowSize(double window){
        _lengthChanger.setWindowSize(window);
        _inheriProbChanger.setWindowSize(window);
    }


    /**
     * This function is to mutate the network by changing its parameters
     */
    public void mutateNetwork(Network network) {
        if(_printDetails){
            System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
            System.out.println("Before rearrangement: "+ network.toString());
            System.out.println("Select operation: "+ this.getClass().getSimpleName());
        }
        if(_lengthChanger instanceof NodeHeightChange){
            ((NodeHeightChange)_lengthChanger).setInfo(network);
        }

        ArrayList<Tuple<Tuple<NetNode, NetNode>, Boolean>> allEdges = new ArrayList<>();
        getAllParameters(network, allEdges);
        _edgeChanged.clear();
        int numEdgeToChange = _numParametersToChange.execute(allEdges.size());

        for (int i = 0; i < numEdgeToChange; i++) {
            Tuple<Tuple<NetNode, NetNode>, Boolean> selectEdgeInfo = allEdges.remove(_random.nextInt(allEdges.size()));
            Tuple<NetNode, NetNode> selectEdge = selectEdgeInfo.Item1;
            if (selectEdgeInfo.Item2) {
                _edgeChanged.add(new Tuple<Tuple<Tuple<NetNode, NetNode>, Boolean>, Double>(selectEdgeInfo, selectEdge.Item2.getParentDistance(selectEdge.Item1)));
                _lengthChanger.setParameters(null, selectEdge, -1, -1);
                _lengthChanger.performOperation();
                if (_printDetails) {
                    System.out.println("Change length of: " + printEdge(selectEdgeInfo.Item1));
                }
            } else {
                _edgeChanged.add(new Tuple<Tuple<Tuple<NetNode, NetNode>, Boolean>, Double>(selectEdgeInfo, selectEdge.Item2.getParentProbability(selectEdge.Item1)));
                _inheriProbChanger.setParameters(null, selectEdge, -1, -1);
                _inheriProbChanger.performOperation();
                if (_printDetails) {
                    System.out.println("Change probability of: " + printEdge(selectEdgeInfo.Item1));
                }
            }
        }
        if (_printDetails) {
            System.out.println("After rearrangement: " + network.toString());
            System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
        }
    }


    /**
     * This function is to undo the last move
     */
    public void undo() {
        for(Tuple<Tuple<Tuple<NetNode, NetNode>, Boolean>, Double> selectEdgeInfo: _edgeChanged){
            if(selectEdgeInfo.Item1.Item2){
                _lengthChanger.setParameters(null, selectEdgeInfo.Item1.Item1, selectEdgeInfo.Item2, -1);
                _lengthChanger.performOperation();
            }
            else{
                _inheriProbChanger.setParameters(null, selectEdgeInfo.Item1.Item1, -1, selectEdgeInfo.Item2);
                _inheriProbChanger.performOperation();
            }
        }
    }


    /**
     * This function is to get all parameters of the network
     *
     * @param network   the species network
     * @param allEdges  all edges in the species network
     */
    abstract protected void getAllParameters(Network network, ArrayList<Tuple<Tuple<NetNode, NetNode>, Boolean>> allEdges);

}