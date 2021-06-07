package edu.rice.cs.bioinfo.programs.soranus;

import edu.rice.cs.bioinfo.programs.soranus.viewModels.WorkspaceVM;
import edu.rice.cs.bioinfo.programs.soranus.views.WorkspaceView;
import edu.rice.cs.bioinfo.programs.soranus.views.swing.WorkspaceViewSwing;

import javax.swing.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/1/13
 * Time: 3:03 PM
 * To change this template use File | Settings | File Templates.
 */
public class SwingApplication extends ApplicationBase
{
    @Override
    public void start()
    {
        SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                SwingApplication.super.start();
            }
        });
    }

    @Override
    protected WorkspaceView makeWorkspaceView(WorkspaceVM vm)
    {
        return new WorkspaceViewSwing(vm);
    }
}
