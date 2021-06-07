package edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing;
import edu.rice.cs.bioinfo.library.phylogenetics.search.ObservableGenerationalScoringSearcher;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/8/12
 * Time: 1:34 PM
 * To change this template use File | Settings | File Templates.
 */
public interface HillClimberObservable<T,S> extends HillClimber<T,S>, ObservableGenerationalScoringSearcher<T,S>
{

}
