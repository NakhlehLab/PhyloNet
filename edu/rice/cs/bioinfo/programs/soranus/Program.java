package edu.rice.cs.bioinfo.programs.soranus;

import edu.rice.cs.bioinfo.programs.soranus.viewModels.WorkspaceVM;
import edu.rice.cs.bioinfo.programs.soranus.views.swing.Workspace;

import javax.swing.*;

public class Program
{
    public static void main(String[] args)
    {
        new DefaultController().startView();

    }
}
