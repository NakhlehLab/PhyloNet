package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move.network.param;

import edu.rice.cs.bioinfo.library.programming.Container;
import edu.rice.cs.bioinfo.library.programming.Proc;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.core.State;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.move.network.NetworkOperator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 12/10/17
 * Time: 6:19 PM
 * To change this template use File | Settings | File Templates.
 */

/**
 * This is NOT a MCMC move! Only work with MLE!
 */

public class OptimizeAll extends NetworkOperator {

    private int _maxRounds = 10;
    protected int _maxTryPerBranch = 10;
    protected double _improvementThreshold = 0.001;
    protected double _maxBranchLength = 6;
    protected double _Brent1 = 0.01;
    protected double _Brent2 = 0.001;
    private UltrametricNetwork _network;
    private State _state;
    private final Boolean _debugMode = false;

    private Map<NetNode, Double> _orginalHeight;
    private Map<NetNode, Map<NetNode, Double>> _originalInheriProb;
    private Map<NetNode, Map<NetNode, Double>> _originalPopSize;

    public OptimizeAll(UltrametricNetwork network) {
        super(network);
        _network = network;
    }

    public OptimizeAll(State state) {
        super(state.getUltrametricNetworkObject());
        _network = state.getUltrametricNetworkObject();
        _state = state;
    }

    @Override
    // returns new score
    public double propose() {
        _orginalHeight = new HashMap<>();
        _originalInheriProb = new HashMap<>();
        _originalPopSize = new HashMap<>();

        boolean continueRounds = true;
        double initialProb = _network.logDensity();

        System.out.println("Optimizer initial likelihood: " + initialProb);
        //if(true) {
        //    return initialProb;
        //}
        final Container<Double> lnGtProbOfSpeciesNetwork = new Container<Double>(initialProb);

        if(_debugMode && !_network.isUltrametric()) {
            System.out.println("Not ultrametric!!!!!!");
        }

        for(Object nodeObj : Networks.postTraversal(_network.getNetwork())) {
            NetNode<NetNodeInfo> node = (NetNode) nodeObj;

            if(node.isNetworkNode()) {
                Iterator<NetNode<NetNodeInfo>> hybridParents = node.getParents().iterator();
                final NetNode hybridParent1 = hybridParents.next();
                final NetNode hybridParent2 = hybridParents.next();

                _originalInheriProb.put(node, new HashMap<>());
                _originalInheriProb.get(node).put(hybridParent1, node.getParentProbability(hybridParent1));
                _originalInheriProb.get(node).put(hybridParent2, node.getParentProbability(hybridParent2));

                _originalPopSize.put(node, new HashMap<>());
                _originalPopSize.get(node).put(hybridParent1, node.getParentSupport(hybridParent1));
                _originalPopSize.get(node).put(hybridParent2, node.getParentSupport(hybridParent2));
            } else if(node.isRoot()) {
                _originalPopSize.put(node, new HashMap<>());
                _originalPopSize.get(node).put(node, node.getRootPopSize());
            } else {
                final NetNode parent = node.getParents().iterator().next();
                _originalPopSize.put(node, new HashMap<>());
                _originalPopSize.get(node).put(parent, node.getParentSupport(parent));
            }

            if (!node.isLeaf()) {
                _orginalHeight.put(node, node.getData().getHeight());
            }

        }

        if(_debugMode && !_network.isUltrametric()) {
            System.out.println("Not ultrametric!!!!!!");
        }

        int roundIndex = 0;
        for(; roundIndex <_maxRounds && continueRounds; roundIndex++)
        {
            System.out.println("Optimization Round: " + roundIndex + " from " + lnGtProbOfSpeciesNetwork.getContents());
            double lnGtProbLastRound = lnGtProbOfSpeciesNetwork.getContents();
            List<Proc> assigmentActions = new ArrayList<Proc>(); // store adjustment commands here.  Will execute them one by one later.

            // for inheritance probablitiy
            for(Object childObj : _network.getNetwork().getNetworkNodes()) // find every hybrid node
            {
                NetNode child = (NetNode) childObj;
                //System.out.println("Inheritance prob: " + child.getName());

                Iterator<NetNode<Object>> hybridParents = child.getParents().iterator();
                final NetNode hybridParent1 = hybridParents.next();
                final NetNode hybridParent2 = hybridParents.next();

                assigmentActions.add(new Proc()
                {
                    public void execute()
                    {
                        UnivariateFunction functionToOptimize = new UnivariateFunction() {
                            public double value(double suggestedProb) {
                                double incumbentHybridProbParent1 = child.getParentProbability(hybridParent1);
                                child.setParentProbability(hybridParent1, suggestedProb);
                                child.setParentProbability(hybridParent2, 1.0 - suggestedProb);

                                double lnProb =  _network.logDensity();
                                if(lnProb > lnGtProbOfSpeciesNetwork.getContents()) // change improved GTProb, keep it
                                {

                                    lnGtProbOfSpeciesNetwork.setContents(lnProb);
                                }
                                else // change did not improve, roll back
                                {
                                    child.setParentProbability(hybridParent1, incumbentHybridProbParent1);
                                    child.setParentProbability(hybridParent2, 1.0 - incumbentHybridProbParent1);
                                }
                                return lnProb;
                            }
                        };
                        BrentOptimizer optimizer = new BrentOptimizer(_Brent1, _Brent2); // very small numbers so we control when brent stops, not brent.

                        try
                        {
                            optimizer.optimize(_maxTryPerBranch, functionToOptimize, GoalType.MAXIMIZE, 0, 1.0);
                        }
                        catch(TooManyEvaluationsException e)  // _maxAssigmentAttemptsPerBranchParam exceeded
                        {
                        }
                    }
                });
            }

            // for population mutation rate
            for(Object childObj : _network.getNetwork().bfs()) // find every node
            {
                NetNode child = (NetNode) childObj;

                Iterator<NetNode> hybridParents = child.getParents().iterator();
                List<NetNode> parentList = new ArrayList<>();

                parentList.add(child.isRoot() ? child : hybridParents.next());
                if(child.isNetworkNode()) {
                    parentList.add(hybridParents.next());
                }

                if(Utils._CONST_POP_SIZE && !child.isRoot()) {
                    continue;
                }

                for(NetNode parent : parentList) {

                    assigmentActions.add(new Proc() {
                        public void execute() {
                            //System.out.println("Pop mutation rate: " + child.getName());

                            final double currentPopSize = parent == child ? child.getRootPopSize() : child.getParentSupport(parent);

                            final Container<Double> minPopSize = new Container<Double>(0.0);
                            final Container<Double> maxPopSize = new Container<Double>(Double.MAX_VALUE);

                            minPopSize.setContents(Math.max(currentPopSize - 0.5 * Utils._POP_SIZE_WINDOW_SIZE, 0.00001));
                            maxPopSize.setContents(currentPopSize + 0.5 * Utils._POP_SIZE_WINDOW_SIZE);

                            UnivariateFunction functionToOptimize = new UnivariateFunction() {
                                public double value(double suggestedPopSize) {
                                    double incumbentPopSize = parent == child ? child.getRootPopSize() : child.getParentSupport(parent);

                                    if(child == parent) {
                                        child.setRootPopSize(suggestedPopSize);
                                    } else {
                                        child.setParentSupport(parent, suggestedPopSize);
                                    }

                                    double lnProb = _network.logDensity();
                                    if (lnProb > lnGtProbOfSpeciesNetwork.getContents()) // change improved GTProb, keep it
                                    {

                                        lnGtProbOfSpeciesNetwork.setContents(lnProb);
                                    } else // change did not improve, roll back
                                    {
                                        if(child == parent) {
                                            child.setRootPopSize(incumbentPopSize);
                                        } else {
                                            child.setParentSupport(parent, incumbentPopSize);
                                        }
                                    }
                                    return lnProb;
                                }
                            };
                            BrentOptimizer optimizer = new BrentOptimizer(_Brent1, _Brent2); // very small numbers so we control when brent stops, not brent.

                            try {
                                optimizer.optimize(_maxTryPerBranch, functionToOptimize, GoalType.MAXIMIZE, minPopSize.getContents(), maxPopSize.getContents());

                            } catch (TooManyEvaluationsException e)  // _maxAssigmentAttemptsPerBranchParam exceeded
                            {
                            }
                        }
                    });
                }
            }

            if(_debugMode && !_network.isUltrametric()) {
                System.out.println("Not ultrametric!!!!!!");
            }

            for(Object nodeObj : Networks.postTraversal(_network.getNetwork())) {
                NetNode<NetNodeInfo> node = (NetNode) nodeObj;
                if (node.isLeaf()) {
                    continue;
                }

                if(_debugMode && !_network.isUltrametric()) {
                    System.out.println("Not ultrametric!!!!!!");
                }

                assigmentActions.add(new Proc() {
                    public void execute() {
                        //System.out.println("Height: " + node.getName());

                        ChangeTime changeTime = new ChangeTime(_network);
                        final Container<Double> minHeight = new Container<Double>(0.0);
                        final Container<Double> maxHeight = new Container<Double>(Double.MAX_VALUE);

                        double[] bounds = _network.getLowerAndUpperBound(node);
                        minHeight.setContents(bounds[0] + 0.00001);
                        maxHeight.setContents(Math.min(bounds[1], _maxBranchLength) - 0.00001);
                        if(maxHeight.getContents() < minHeight.getContents()) {
                            System.out.println("Invalid bound!!!!");
                        }

                        UnivariateFunction functionToOptimize = new UnivariateFunction() {

                            public double value(double suggestedHeight) {  // brent suggests a new branch length
                                double incumbentHeight = node.getData().getHeight();

                                if(_debugMode && !_network.isUltrametric()) {
                                    System.out.println("Not ultrametric!!!!!!");
                                }

                                changeTime.setNodeHeight(node, suggestedHeight);

                                if(_debugMode && !_network.isUltrametric()) {
                                    System.out.println("Not ultrametric!!!!!!");
                                }

                                double lnProb = _network.logDensity();
                                if(lnProb > lnGtProbOfSpeciesNetwork.getContents()) // did improve, keep change
                                {
                                    lnGtProbOfSpeciesNetwork.setContents(lnProb);
                                } else {
                                    changeTime.setNodeHeight(node, incumbentHeight);
                                }

                                if(_debugMode && !_network.isUltrametric()) {
                                    System.out.println("Not ultrametric!!!!!!");
                                }

                                return lnProb;
                            }
                        };

                        BrentOptimizer optimizer = new BrentOptimizer(_Brent1, _Brent2); // very small numbers so we control when brent stops, not brent.

                        try
                        {
                            optimizer.optimize(_maxTryPerBranch, functionToOptimize, GoalType.MAXIMIZE, minHeight.getContents(), maxHeight.getContents());
                        }
                        catch(TooManyEvaluationsException e) // _maxAssigmentAttemptsPerBranchParam exceeded
                        {
                        }
                    }
                });
            }

            // TODO: add random change when stuck in local optimal

            Collections.shuffle(assigmentActions, Randomizer.getRandom());

            for(Proc assigment : assigmentActions)   // for each change attempt, perform attempt
            {
                assigment.execute();
                //System.out.println(lnGtProbOfSpeciesNetwork.getContents());
            }
            if( ((double)lnGtProbOfSpeciesNetwork.getContents()) == lnGtProbLastRound)  // if no improvement was made wrt to last around, stop trying to find a better assignment
            {
                continueRounds = false;
            }
            else if (lnGtProbOfSpeciesNetwork.getContents() > lnGtProbLastRound) // improvement was made, ensure it is large enough wrt to improvement threshold to continue searching
            {

                double improvementPercentage = Math.pow(Math.E, (lnGtProbOfSpeciesNetwork.getContents() - lnGtProbLastRound)) - 1.0;  // how much did we improve over last round
                //System.out.println(improvementPercentage + " vs. " + _improvementThreshold);
                if(improvementPercentage < _improvementThreshold  )  // improved, but not enough to keep searching
                {
                    continueRounds = false;
                }
            }
            else
            {
                throw new IllegalStateException("Should never have decreased prob.");
            }
        }
        return lnGtProbOfSpeciesNetwork.getContents();
    }

    @Override
    public void undo() {

        if(_debugMode && _orginalHeight == null) {
            System.out.println("null pointer");
        }

        if(_debugMode && !_network.isUltrametric()) {
            System.out.println("Not ultrametric!!!!!!");
        }

        for(Object nodeObj : _network.getNetwork().bfs()) {
            NetNode<NetNodeInfo> node = (NetNode) nodeObj;
            if(_orginalHeight.containsKey(node)) {
                node.getData().setHeight(_orginalHeight.get(node));
            }
        }

        for(Object nodeObj : _network.getNetwork().bfs()) {
            NetNode<NetNodeInfo> node = (NetNode) nodeObj;
            for(Object parentObj : node.getParents()) {
                NetNode<NetNodeInfo> parent = (NetNode) parentObj;
                node.setParentDistance(parent, parent.getData().getHeight() - node.getData().getHeight());
                if(node.getParentDistance(parent) < 0) {
                    System.out.println("distance abnormal");
                }
            }
        }

        if(_debugMode && !_network.isUltrametric()) {
            System.out.println("Not ultrametric!!!!!!");
        }

        for(NetNode<NetNodeInfo> child : _originalInheriProb.keySet()) {
            for(NetNode<NetNodeInfo> parent : _originalInheriProb.get(child).keySet()) {
                child.setParentProbability(parent, _originalInheriProb.get(child).get(parent));
            }
        }

        for(NetNode<NetNodeInfo> child : _originalPopSize.keySet()) {
            for(NetNode<NetNodeInfo> parent : _originalPopSize.get(child).keySet()) {
                if(child != parent) {
                    child.setParentSupport(parent, _originalPopSize.get(child).get(parent));
                }
            }
            if(child.isRoot()) {
                child.setRootPopSize(_originalPopSize.get(child).get(child));
            }
        }

        if(_debugMode && !_network.isUltrametric()) {
            System.out.println("Not ultrametric!!!!!!");
        }

        _orginalHeight = null;
        _originalInheriProb = null;
    }

    @Override
    public String getName() {
        return "Optimize-All-NotMCMC";
    }
}
