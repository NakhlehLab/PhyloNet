package edu.rice.cs.bioinfo.programs.soranus.views.swing;

import edu.rice.cs.bioinfo.programs.soranus.viewModels.NeighborJoiningVM;

import javax.swing.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/7/13
 * Time: 1:04 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class NeighborJoiningViewJFilesBase<N,E> extends JPanel
{
    private NeighborJoiningVM<N,E> _vm;

    public NeighborJoiningVM<N,E> getViewModel()
    {
        return _vm;
    }

    public NeighborJoiningViewJFilesBase(NeighborJoiningVM<N, E> vm) {
        _vm = vm;
    }
}
