package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Func4;
import org.joda.time.LocalDate;

import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/18/13
 * Time: 11:10 AM
 * To change this template use File | Settings | File Templates.
 */
public class SnitkinTransMapInferrer<P,S,D extends Comparable<D>> extends SnitkinTransMapInferrerTemplate<P,S,D>
{
    private final Func2<S,S,Integer> _getGeneticDistance;

    private final Func2<D,D,D> _add;

    private final Func2<D,D,D> _subtract;

    private final Func1<Integer,D> _makeDistance;

    private final Func4<P,P,Integer,Integer,SnitkinEdge<P,D>> _makeEdge;

    private final Func<D> _getMaxDistance;

    public SnitkinTransMapInferrer(Map<P, Map<LocalDate, Object>> patientTraces, Map<P, LocalDate> patientToFirstPositiveDate, Map<S, P> genomeToPatient,
                                   Func2<S, S, Integer> getGeneticDistance, Func2<D, D, D> add, Func2<D, D, D> subtract, Func1<Integer,D> makeDistance, Func<D> getMaxDistance,
                                   Func4<P, P, Integer, Integer, SnitkinEdge<P, D>> makeEdge) {
        super(patientTraces, patientToFirstPositiveDate, genomeToPatient);
        _getGeneticDistance = getGeneticDistance;
        _add = add;
        _subtract = subtract;
        _makeDistance = makeDistance;
        _getMaxDistance = getMaxDistance;
        _makeEdge = makeEdge;
    }


    @Override
    protected D getMaxDistance() {
        return _getMaxDistance.execute();
    }

    @Override
    protected int getGeneticDistance(S sequence1, S sequence2) {
        return _getGeneticDistance.execute(sequence1, sequence2);  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected D add(D a, D b) {
        return _add.execute(a, b);  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected D subtract(D a, D b) {
        return _subtract.execute(a, b);  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected D makeDistance(int intDistance) {
        return _makeDistance.execute(intDistance);  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected SnitkinEdge<P, D> makeEdge(P source, P destination, int geneticDistance, int epidemiologicalDistance) {
        return _makeEdge.execute(source,destination,geneticDistance,epidemiologicalDistance);  //To change body of implemented methods use File | Settings | File Templates.
    }
}
