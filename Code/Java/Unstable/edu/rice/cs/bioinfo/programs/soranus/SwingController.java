package edu.rice.cs.bioinfo.programs.soranus;


import edu.rice.cs.bioinfo.programs.soranus.viewModels.WorkspaceVM;
import edu.rice.cs.bioinfo.programs.soranus.views.WorkspaceView;
import edu.rice.cs.bioinfo.programs.soranus.views.swing.WorkspaceViewSwing;

import javax.swing.*;

class SwingController extends ControllerBase
{
    @Override
    public void startView()
    {
        SwingUtilities.invokeLater(new Runnable() {
            public void run() {

                SwingController.super.startView();

            }
        });
    }

    @Override
    protected WorkspaceView makeWorkspaceView(WorkspaceVM vm)
    {
        return new WorkspaceViewSwing(vm);
    }

}
