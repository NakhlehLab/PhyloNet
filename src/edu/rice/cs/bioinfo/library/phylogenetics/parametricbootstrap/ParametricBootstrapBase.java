package edu.rice.cs.bioinfo.library.phylogenetics.parametricbootstrap;

import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.NetworkValidator;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Func3;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 3/25/13
 * Time: 2:13 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ParametricBootstrapBase<T1, T2, T3, S> implements ParametricBootstrap<T1, T2, T3, S>{
    protected Func2<T1, S, T2> _simulator;
    protected Func2<T1, List<T2>, T3> _estimator;

    public ParametricBootstrapBase(){};

    public ParametricBootstrapBase(Func2<T1, S, T2> simulator, Func2<T1, List<T2>, T3> estimator){
        _simulator = simulator;
        _estimator = estimator;
    }

    public T3 performBootstrapping(T1 model, S sampleSize, int numRound){
        List<T2> samplesDrawn = new ArrayList<T2>();
        for(int i=0; i<numRound; i++){
            T2 sample = _simulator.execute(model, sampleSize);
            samplesDrawn.add(sample);
        }
        return _estimator.execute(model, samplesDrawn);
    }
}
