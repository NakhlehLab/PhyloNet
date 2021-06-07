package edu.rice.cs.bioinfo.library.phylogenetics.parametricbootstrap.network;

import edu.rice.cs.bioinfo.library.phylogenetics.parametricbootstrap.ParametricBootstrapBase;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Func3;

import java.util.ArrayList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 3/25/13
 * Time: 3:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkParametricBootstrap<T1,T2,T3,S> extends ParametricBootstrapBase<T1,T2,T3,S> {
    Func1<T2, T1> _inference;
    Func2<T1, List<T1>, T3> _estimator;

    public NetworkParametricBootstrap(Func2<T1, S, T2> simulator, Func1<T2,T1> inference, Func2<T1, List<T1>, T3> estimator){
        _simulator = simulator;
        _estimator = estimator;
        _inference = inference;
    }

    public T3 performBootstrapping(T1 model, S sampleSize, int numRound){
        List<T1> inferredModels = new ArrayList<T1>();
        for(int i=0; i<numRound; i++){
            T2 sample = _simulator.execute(model, sampleSize);
            inferredModels.add(_inference.execute(sample));
        }
        return _estimator.execute(model, inferredModels);
    }
}
