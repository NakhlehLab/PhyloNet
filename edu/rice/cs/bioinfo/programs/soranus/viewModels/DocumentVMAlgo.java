package edu.rice.cs.bioinfo.programs.soranus.viewModels;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/5/13
 * Time: 6:31 PM
 * To change this template use File | Settings | File Templates.
 */
public interface DocumentVMAlgo<R, E extends Exception>
{
    public <N,Ed> R forTransMapVM(TransMapVM<N,Ed> vm) throws E;

    public <N,Ed> R forNeighborJoiningVM(NeighborJoiningVM<N,Ed> vm) throws E;

    public R forXMLDataVM(XMLDataVM xmlDataVM);

    public R forVAALOutDataVM(VAALOutDataVM vaalOutDataVM);
}
